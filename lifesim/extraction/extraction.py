import numpy as np
from lifesim.extraction.extraction_auxiliary import *
from astropy import units as u
import matplotlib.pyplot as plt
import scipy as sp
from lifesim.util.radiation import black_body
from lifesim.core.modules import ExtractionModule
from tqdm import tqdm
import multiprocessing as multiprocessing
import mpmath as mp
from matplotlib.ticker import ScalarFormatter
import random as ran
import warnings

### Code Author: Philipp Binkert

class ML_Extraction(ExtractionModule):
    '''
    The ML extraction class is a signal extraction module based on the maximum likelihood method described by Dannert et
    al. (2022) in LIFE II:  Signal simulation, signal extraction and fundamental exoplanet parameters from single epoch
    observations. The code is based on the code for an older Lifesim version by Maurice Ottiger, as described in his
    Master Thesis "Exoplanet Detection Yield Simulations and Signal Analysis Studies for the LIFE Mission" (2020)
    '''

    def __init__(self, name: str):
        """
        Parameters
        ----------
        name : str
            Name of the module.
        """
        super().__init__(name=name)


    def create_dips(self):
        '''
        This function creates dips in the planet blackbody spectrum based on the self.atmospheric_scenario attributes

        :return isdip: np.ndarray of size n_molecules; denotes for each molecule whether a corresponding dip was
                        induced []
        :return SNR_ps_new: float; new SNR_ps value for the planet after a potential dip was induced []
        '''

        is_dip = np.zeros((self.n_molecules))

        perfect_blackbody_curve = self.single_data_row['planet_flux_use'][0]

        for i in range(is_dip.size):

            randomizer = ran.random()

            # if (randomizer <= self.dip_matrix[i][4]):
            if (randomizer <= self.atmospheric_scenario.dip_likelihoods[i]):
                # generate a Gaussian-shaped dip
                dip_center = self.atmospheric_scenario.dip_centers[i]
                dip_width = self.atmospheric_scenario.get_width(self.atmospheric_scenario.dip_molecules[i])
                dip_depth = self.atmospheric_scenario.get_depth(self.atmospheric_scenario.dip_molecules[i])

                gaussian_dip = sp.stats.norm.pdf(self.wl_bins, loc=dip_center, scale=dip_width)
                gaussian_dip_normalized = gaussian_dip / np.max(gaussian_dip) * \
                                          self.single_data_row['planet_flux_use'][0][np.abs(self.wl_bins -
                                                                                            dip_center).argmin()]

                dipped_curve = self.single_data_row['planet_flux_use'][0] - dip_depth * gaussian_dip_normalized

                is_dip[i] = True

                self.single_data_row['planet_flux_use'][0] = dipped_curve

        # calculate the new SNR_ps
        # this is required so that when you set the index to 'None' calling the next function it takes the right angsep;
        #   also, you have to reset the baseline
        self.data.single['angsep'] = self.single_data_row['angsep']
        self.data.inst['bl'] = self.single_data_row['baseline']

        transm_eff, transm_noise = self.run_socket(s_name='transmission',
                                                   method='transmission_efficiency',
                                                   index=None)

        noise_planet = 2 * self.single_data_row['planet_flux_use'][0] * transm_noise
        noise_tot = noise_planet + self.single_data_row['noise_astro'][0]

        flux_planet = self.single_data_row['planet_flux_use'][0] * transm_eff

        SNR_1h = np.sqrt((flux_planet ** 2 / noise_tot).sum())
        SNR_ps_new = SNR_1h * np.sqrt(self.single_data_row['int_time'] / 3600)

        if (self.plot == True):

            for i in range(is_dip.size):
                if (is_dip[i] == True):
                    print(self.atmospheric_scenario.dip_molecules[i], 'dip induced')

            plt.plot(self.wl_bins * 10 ** 6, perfect_blackbody_curve / self.wl_bin_widths,
                     label=r'perfect blackbody, $\mathrm{SNR_{ps}} =$' + str(
                         np.round(self.single_data_row['snr_current'], 1)), color='black')
            plt.plot(self.wl_bins * 10 ** 6, self.single_data_row['planet_flux_use'][0] / self.wl_bin_widths,
                     label=r'blackbody with dips, $\mathrm{SNR_{ps}} =$' + str(np.round(SNR_ps_new, 1)), color='red')
            plt.legend(loc='best')
            plt.grid()
            plt.title('Input blackbody curves')
            plt.xlabel(r"$\lambda$ [$\mu$m]", fontsize=12)
            plt.ylabel(r"Planet flux $F_\lambda$ [ph $\mathrm{s}^{-1}$m$^{-2}\mu \mathrm{m}^{-1}$]", fontsize=12)
            plt.show()

        return is_dip, SNR_ps_new


    def get_B_C(self):
        '''
        This functions calculates two matrices that aid in the calculation of the cost function J and the estimated
        planet flux based on the received signal and the transmission functions. See Appendix B of LIFE II for an
        explanation of the exact calculations

        Attributes
        --------------
        - self.B: np.ndarray of shape (L,radial_ang_px,n_steps); matrix used for the calculation of the estimated flux
                as described in Appendix B of LIFE II []
        - self.C: np.ndarray of shape (L,radial_ang_px,n_steps); matrix used for the calculation of the estimated flux
                and the cost function as described in Appendix B of LIFE II []
        '''

        # get dimensions of transmission map
        (n_l, n_r, n_p) = self.tm_chop.shape

        # get variance of received signals
        var = np.var(self.signals, axis=1, ddof=1)

        # calculate B; in this step, the time series is included by the np.repeat function
        B_vec = (self.tm_chop ** 2).sum(axis=-1) / var[:, np.newaxis]
        self.B = np.repeat(B_vec[:, :, np.newaxis], n_p, axis=2)

        # take T twice back to back
        T_exp = np.tile(self.tm_chop, 2)

        # calculate C; also here, the time series is included
        C = np.empty((n_l, n_r, n_p // 2))

        for i in range(n_p // 2):
            T_i = T_exp[:, :, n_p - i: 2 * n_p - i]
            C[:, :, i] = np.einsum("ij,ikj->ik", self.signals, T_i)

        # use antisymmetry of C to speed up calculation
        C = np.concatenate((C, -C), axis=2)
        self.C = (C.T / var).T

        return


    def get_B_C_whitened(self):
        '''
        This function calculates the matices B&C analogously to the previous functions, but now the attribute self.noise
        has been defined and quantities taking only the variance of the noise are calculated. These quantities are used
        in cost_func_MAP to calculate the various SNRs

        Attributes
        --------------
        - self.B_noise: np.ndarray of shape (L,radial_ang_px,n_steps); like B, but only for the noise (no planet signal)
        - self.C_noise_only: np.ndarray of shape (L,radial_ang_px,n_steps); like C, but only for the noise
                                        (no planet signal)
        - self.C_noise_var: np.ndarray of shape (L,radial_ang_px,n_steps); like C, but takes the variance only of
                                        the noise
        '''

        (n_l, n_r, n_p) = self.tm_chop.shape

        # get variance of received noise (without signal)
        var_noise = np.var(self.noise, axis=1, ddof=1)

        # calculate B analogously to before for the new variance
        B_vec_noise = (self.tm_chop ** 2).sum(axis=-1) / var_noise[:, np.newaxis]
        self.B_noise = np.repeat(B_vec_noise[:, :, np.newaxis], n_p, axis=2)

        # take T twice back to back
        T_exp = np.tile(self.tm_chop, 2)

        C = np.empty((n_l, n_r, n_p // 2))
        C_noise = np.empty((n_l, n_r, n_p // 2))

        for i in range(n_p // 2):
            T_i = T_exp[:, :, n_p - i: 2 * n_p - i]
            C[:, :, i] = np.einsum("ij,ikj->ik", self.signals, T_i)
            C_noise[:, :, i] = np.einsum("ij,ikj->ik", self.noise, T_i)

        # use antisymmetry of C(_noise) to speed up calculation
        C = np.concatenate((C, -C), axis=2)
        C_noise = np.concatenate((C_noise, -C_noise), axis=2)

        # calculate C_noise_only and C_noise_var; the latter used the normal signal (including planet) while the first
        #   is a purely noisy signal
        self.C_noise_only = (C_noise.T / var_noise).T
        self.C_noise_var = (C.T / var_noise).T

        return


    def get_transm_map(self):
        '''
        This function first calculates the transmission function by defining a polar coordinate system and then calling
        the corresponding function from the transmission module

        Attributes
        ----------------
        - self.tm_chop: np.ndarray of shape (L, self.radial_ang_pix, self.n_steps); contains the transmission function
                            values for each point in the transmission map
        '''

        # define the azimuthal coordinates
        phi_lin = np.linspace(0, 2 * np.pi, self.n_steps, endpoint=False)  # 1D array with azimuthal coordinates
        phi_mat = np.tile(phi_lin, (self.radial_ang_px, 1))  # 2D map with azimuthal coordinates

        # define the radial coordinates
        theta_lin = np.linspace(0, 1, self.radial_ang_px,
                                endpoint=False)  # 1D array with radial separation coord [radians]
        theta_lin += 1 * 0.5 / self.radial_ang_px  # shift the coordinates to the "center" of the bins
        theta_mat = np.tile(theta_lin, (self.n_steps, 1)).T  # 2D array with radial separation coordinates [radians]

        # define matrices which include the different wl parameters; the radial coordinate must be scaled by the hfov
        phi_mat_wl = np.zeros((self.L, self.radial_ang_px, self.n_steps))
        theta_mat_wl = np.zeros((self.L, self.radial_ang_px, self.n_steps))

        for i in range(self.L):
            phi_mat_wl[i, :, :] = phi_mat
            theta_mat_wl[i, :, :] = theta_mat * self.field_of_search

        # define the input parameters for the transmission function
        d_alpha = theta_mat_wl * np.cos(phi_mat_wl)
        d_beta = theta_mat_wl * np.sin(phi_mat_wl)

        # calculate the transmission map
        _, _, _, _, tm_chop = self.run_socket(s_name='transmission',
                                              method='transmission_map',
                                              map_selection=['tm_chop'],
                                              direct_mode=True,
                                              d_alpha=d_alpha,
                                              d_beta=d_beta,
                                              hfov=self.field_of_search,
                                              image_size=self.image_size
                                              )

        # if applicable, plot the resulting transmission map
        if (self.plot == True):
            self.plot_transmissionmap()

        # normalize the transmission map to instrument performance (required for compatability with functions written
        #   for the old lifesim version
        self.tm_chop = tm_chop * self.single_data_row['int_time'] / self.n_steps * self.data.inst['telescope_area'] \
                       * self.data.inst['eff_tot'] * self.wl_bin_widths[:, np.newaxis, np.newaxis] * 10 ** 6

        return


    def plot_transmissionmap(self):
        '''
        This function plots the differential transmission maps generated at a suitable wavelength, in cartesian and
        polar coordinates, as well as the modulated planet and total signal for one full rotation of the instrument.
        The calculations are similar as in get_transm_map()
        Note that if the transmission map does not fill out the entire plot, this means that the part where nothing
        is shown is outside of the field of view for that wavelength
        '''

        # define the plot resolution
        resolution = 1000

        # show the plot at 12.3 micron (arbitrary, change as desired)
        wl_bin = int(23 / 31 * self.L)
        wl_value = np.round(self.wl_bins[wl_bin] * 10 ** 6, 1)

        # define the cartesian and polar coordinates
        x_lin = np.linspace(-1, 1, resolution, endpoint=False)
        x_lin += 1 * 0.5 / resolution
        x_mat = np.tile(x_lin, (resolution, 1))
        y_mat = np.tile(x_lin, (resolution, 1)).T
        x_mat_wl = np.zeros((self.L, resolution, resolution))
        y_mat_wl = np.zeros((self.L, resolution, resolution))

        phi_lin_plot = np.linspace(0, 2 * np.pi, self.n_steps, endpoint=False)
        phi_mat_plot = np.tile(phi_lin_plot, (resolution, 1))
        theta_lin_plot = np.linspace(0, 1, resolution, endpoint=False)
        theta_lin_plot += 1 * 0.5 / resolution
        theta_mat_plot = np.tile(theta_lin_plot, (self.n_steps, 1)).T
        phi_mat_plot_wl = np.zeros((self.L, resolution, self.n_steps))
        theta_mat_plot_wl = np.zeros((self.L, resolution, self.n_steps))

        # define the required matrices and input parameters
        for i in range(self.L):
            x_mat_wl[i, :, :] = x_mat * self.hfov[i]
            y_mat_wl[i, :, :] = y_mat * self.hfov[i]
            phi_mat_plot_wl[i, :, :] = phi_mat_plot
            theta_mat_plot_wl[i, :, :] = theta_mat_plot * self.hfov[i]

        d_alpha_plot = theta_mat_plot_wl * np.cos(phi_mat_plot_wl)
        d_beta_plot = theta_mat_plot_wl * np.sin(phi_mat_plot_wl)

        # get the transmission maps
        _, _, _, _, tm_chop_plot_lin = self.run_socket(s_name='transmission',
                                                       method='transmission_map',
                                                       map_selection=['tm_chop'],
                                                       direct_mode=True,
                                                       d_alpha=x_mat_wl,
                                                       d_beta=y_mat_wl,
                                                       hfov=self.field_of_search,
                                                       image_size=self.image_size
                                                       )

        _, _, _, _, tm_chop_plot_pol = self.run_socket(s_name='transmission',
                                                       method='transmission_map',
                                                       map_selection=['tm_chop'],
                                                       direct_mode=True,
                                                       d_alpha=d_alpha_plot,
                                                       d_beta=d_beta_plot,
                                                       hfov=self.field_of_search,
                                                       image_size=self.image_size
                                                       )

        # plot the linear map
        X_lin = np.linspace(-self.hfov[wl_bin] * 3600000 * 180 / np.pi, self.hfov[wl_bin] * 3600000 * 180 / np.pi,
                            resolution)
        Y_lin = X_lin
        levels = np.linspace(-1, 1, 100)
        ticks = np.array([-1, -0.5, 0, 0.5, 1])
        contour = plt.contourf(X_lin, Y_lin, tm_chop_plot_lin[wl_bin, :, :], cmap='cividis', levels=levels)
        cbar = plt.colorbar(contour, ticks=ticks)
        cbar.ax.set_ylabel('Transmission')
        circle = plt.Circle((0, 0), radius=self.single_data_row['angsep'] * 1000, fill=False, color='black',
                            linestyle='--', label='planet')
        plt.scatter(0, 0, marker='*', color='white', label='star')
        plt.gca().add_patch(circle)
        plt.title('Linear differential transmission map ($\lambda$=' + str(wl_value) + '$\mu m$)')
        plt.legend(loc='best')
        plt.xlim((-self.single_data_row['angsep'] * 1000 * 1.5, self.single_data_row['angsep'] * 1000 * 1.5))
        plt.ylim((-self.single_data_row['angsep'] * 1000 * 1.5, self.single_data_row['angsep'] * 1000 * 1.5))
        plt.xlabel('$\Delta$RA [mas]')
        plt.ylabel('$\Delta$Dec [mas]')
        plt.show()

        # plot the cartesian map
        X_pol = np.linspace(0, self.n_steps, self.n_steps)
        Y_pol = np.linspace(0, self.hfov[wl_bin] * 3600000 * 180 / np.pi, resolution)
        levels = np.linspace(-1, 1, 100)
        col_ticks = np.array([-1, -0.5, 0, 0.5, 1])
        x_ticks = np.array([0, 90, 180, 270, 360])
        contour = plt.contourf(X_pol, Y_pol, tm_chop_plot_pol[wl_bin, :, :], cmap='cividis', levels=levels)
        cbar = plt.colorbar(contour, ticks=col_ticks)
        cbar.ax.set_ylabel('Transmission')
        plt.title('Polar differential transmission map ($\lambda$=' + str(wl_value) + '$\mu m$)')
        plt.hlines(self.single_data_row['angsep'] * 1000, xmin=0, xmax=360, linestyles='--', color='black',
                   label='planet')
        plt.legend(loc='best')
        plt.ylim((0, self.single_data_row['angsep'] * 1000 * 1.5))
        plt.xlabel('Rotation angle [$^{\circ}$]')
        plt.ylabel('Angular separation $\Theta$ [mas]')
        plt.xticks(x_ticks)
        plt.show()

        # plot the modulated signal
        plt.plot(X_pol, self.signals[wl_bin, :], label='total noisy signal', color='dodgerblue')
        plt.plot(X_pol, self.ideal_signals[wl_bin, :], label='pure planet signal', color='black')
        plt.title('Signal Modulation ($\lambda$=' + str(wl_value) + '$\mu m$)')
        plt.xlabel('Rotation angle [$^{\circ}$]')
        plt.ylabel('Amplitude [photons]')
        plt.xticks(x_ticks)
        plt.grid()
        plt.xlim((0, self.n_steps))
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_powerlimits((0, 2))
        plt.gca().yaxis.set_major_formatter(formatter)
        plt.legend(loc='best')
        plt.show()

        return


    def get_F_estimate(self):
        '''
        This function calculates the estimated flux received by the planet based on the matrices B and C.
        See LIFE II Appendix B for an explanation of the calculations

        :return F: np.ndarray of shape (L,radial_ang_px); total estimated flux (planet plus noise) received for each
                    wl bin at every pixel in the image in [photons]
        :return F_pos: np.ndarray of shape (L,radial_ang_px); positive part of total estimated flux received for each
                    wl bin at every pixel in the image in [photons]
        '''

        # get dimensions of the matrix
        (n_l, n_r, n_p) = self.C.shape

        # if there is no regularization, the calculation can be simplified
        if (self.mu == 0):

            # calculation of the total estimated signal (planet signal plus noise)
            F = self.C / self.B
            F_pos = np.where(F >= 0, F, 0)


        # if regularization is included, F is calculated according to B.6, additionally requiring calculation of D
        else:
            try:
                # calculate D^2
                Dsqr_mat = get_Dsqr_mat(n_l)

                B_diags = np.array([np.diag(self.B[:, i, 0]) for i in range(n_r)])  # B is phi-independent

                # calculate the inverse of (B+mu*D^2)
                S = B_diags + self.mu * Dsqr_mat
                Sinv = np.linalg.inv(S)

                # calculate F and F_pos
                F = np.einsum("rlm,mrp ->lrp", Sinv, self.C)

                F_pos = np.empty_like(F)

                for r in range(n_r):
                    for p in range(n_p):
                        F_pos[:, r, p] = sp.optimize.nnls(S[r], self.C[:, r, p], maxiter=200)[0]

            # for some values of mu, (B+mu*D^2) can be singular. In this case, calculate F without regularization and
            #   print a corresponding warning
            except:
                print('Warning: singular matrix obtained for this value of mu; F is calculated for mu=0')
                F = self.C / self.B
                F_pos = np.where(F >= 0, F, 0)

        return F, F_pos


    def get_F_estimate_whitened(self):
        '''
        This function calculates the estimated flux received by the planet after the signal has been whitened. The
        calculations are analogous to get_F_estimate just using the "whitened" versions of the matrices B&C

        :return F_white: np.ndarray of shape (L,radial_ang_px); estimated planet flux (no noise) received for each
                    wl bin at every pixel in the image in [photons]
        :return F_white_pos: np.ndarray of shape (L,radial_ang_px); positive part of estimated planet flux received
                    for each wl bin at every pixel in the image in [photons]
        :return F_noise: np.ndarray of shape (L,radial_ang_px); estimated noise flux (no signal) received for each
                    wl bin at every pixel in the image in [photons]
        :return F_noise_pos: np.ndarray of shape (L,radial_ang_px); positive part of estimated noise flux received
                    for each wl bin at every pixel in the image in [photons]
        '''

        # get dimensions of the matrix
        (n_l, n_r, n_p) = self.C.shape

        # if no regularization is included, the calculation is simplified
        if (self.mu == 0):
            F = self.C / self.B

            # calculation of the purely noisy signal and the estimated pure planet signal (no noise)
            F_noise = self.C_noise_only / self.B_noise
            F_noise_pos = np.where(F_noise >= 0, F_noise, 0)
            F_white = F - F_noise
            F_white_pos = np.where(F_white >= 0, F_white, 0)


        # if regularization is included, the fluxes are calculated according to B.6, additionally requiring calculation
        #   of the D matrix (all calculations for pure noise and pure planet signal)
        else:
            try:
                # calculate D^2
                Dsqr_mat = get_Dsqr_mat(n_l)

                B_diags = np.array([np.diag(self.B[:, i, 0]) for i in range(n_r)])  # B is phi-independent
                B_diags_noise = np.array([np.diag(self.B_noise[:, i, 0]) for i in range(n_r)])

                # calculate the inverse of (B+mu*D^2)
                S = B_diags + self.mu * Dsqr_mat
                Sinv = np.linalg.inv(S)
                S_noise = B_diags_noise + self.mu * Dsqr_mat
                Sinv_noise = np.linalg.inv(S_noise)

                # calculate F and F_pos for the noisy and pure signal
                F = np.einsum("rlm,mrp ->lrp", Sinv, self.C)
                F_noise = np.einsum("rlm,mrp ->lrp", Sinv_noise, self.C_noise_only)

                F_noise_pos = np.empty_like(F_noise)

                for r in range(n_r):
                    for p in range(n_p):
                        F_noise_pos[:, r, p] = sp.optimize.nnls(S_noise[r], self.C_noise_only[:, r, p], maxiter=200)[0]

                F_white = F - F_noise
                F_white_pos = np.where(F_white >= 0, F_white, 0)

            # for some values of mu, (B+mu*D^2) can be singular. In this case, calculate F without regularization and
            #   print a corresponding warning
            except:
                print('Warning: singular matrix obtained for this value of mu; F is calculated for mu=0')

                F = self.C / self.B

                F_noise = self.C_noise_only / self.B_noise
                F_white = F - F_noise
                F_white_pos = np.where(F_white >= 0, F_white, 0)

        return F_white, F_white_pos, F_noise, F_noise_pos


    def detect_dips(self, spectra, sigmas):
        '''
        This function calculates the t-scores of potential dip detections according to t-statistics for each of the
        molecules considered. Additionally, the curve fit is performed excluding the band range where a dip has
        presumabely been detected

        :param spectra: np.ndarray of size L; contains the extracted planet flux in each wl-bin in [photons]
        :param sigmas: np.ndarray of size L; contains the sigmas of the extracted planet flux in each wl-bin [photons]

        :return t_scores: np.ndarray of size n_molecules; contains the t-score values for each of the
                        considered molecules []
        :return popt_extr: optimal values for (T,R) such that the sum of the squared residuals to the
                        blackbody curve is minimized [K,R_earth]
        :return pcov: np.ndarray of shape (2,2); covariance matrix of the obtained popt, a measure of uncertainty.
                        Dimensions [K^2] and [R_earth^2] on the first and second diagonal element
        :return perr: np.ndarray of size 2; standard deviation of the errors on T and R [K,R_earth]
        '''

        # initialize matrices
        t_scores = np.empty((self.n_molecules))
        dip_extracted = np.zeros((self.n_molecules))
        popts = np.empty((self.n_molecules, 2))
        pcovs = np.empty((self.n_molecules, 2, 2))
        perrs = np.empty((self.n_molecules, 2))

        # loop through all molecules considered
        for i in range(t_scores.size):

            # calculate the bin numbers which would be affected by the dip
            central_dip_bin = np.abs(self.wl_bins - self.atmospheric_scenario.dip_centers[i]).argmin()
            lower_dip_bin = max(0, np.abs(self.wl_bins - (
                        self.atmospheric_scenario.dip_centers[i] - 2 * self.atmospheric_scenario.dip_avg_widths[
                    i])).argmin())
            upper_dip_bin = min(self.L, central_dip_bin + (central_dip_bin - lower_dip_bin))
            dipped_bins = np.arange(lower_dip_bin, upper_dip_bin + 1)

            if (dipped_bins.size <= 1):
                dipped_bins = np.array([central_dip_bin-1,central_dip_bin,central_dip_bin+1])

            truncated_wl_bins_and_distS_array = np.delete(self.wl_bins_and_distS_array, dipped_bins, axis=1)
            truncated_spectra = np.delete(spectra, dipped_bins, axis=0)
            truncated_sigmas = np.delete(sigmas, dipped_bins, axis=0)

            # perform the fit of T and R and calculate the uncertainties
            popts[i], pcovs[i] = sp.optimize.curve_fit(BB_for_fit, xdata=truncated_wl_bins_and_distS_array,
                                                       ydata=truncated_spectra,
                                                       sigma=truncated_sigmas, p0=(300, 1.), absolute_sigma=True,
                                                       maxfev=10000)


            # get the blackbody curve for the estimated parameters
            Fp_fit = BB_for_fit(self.wl_bins_and_distS_array, popts[i][0], popts[i][1])
            residuals = spectra - Fp_fit

            truncated_residuals = np.delete(residuals, dipped_bins, axis=0)
            dip_residuals = residuals[dipped_bins]

            mean_dip_residuals = np.mean(dip_residuals)
            mean_truncated_residuals_data = np.mean(truncated_residuals)
            mean_truncated_residuals_true = 0

            std_truncated_residuals = np.std(truncated_residuals)

            dof = dipped_bins.size - 1

            t_statistic = (mean_dip_residuals - mean_truncated_residuals_true) / (
                    std_truncated_residuals / np.sqrt(dipped_bins.size))

            t_scores[i], _ = sp.integrate.quad(lambda x: t_distribution(x, dof), -np.infty, t_statistic)

            # plot and decide which calculated T to use based on this alpha
            alpha = 0.01

            if (t_scores[i] < alpha):
                dip_extracted[i] = True

            if (self.plot == True):
                critical_value = sp.stats.t.ppf(alpha, dof)

                n_bins = 30
                max_resid = np.max(np.abs(residuals))
                residual_array = np.linspace(-max_resid, max_resid, n_bins)

                counts, bin_edges = np.histogram(truncated_residuals, bins=n_bins)
                area_to_normalize = np.sum(counts * (bin_edges[1] - bin_edges[0]))

                weights_truncated = np.ones_like(truncated_residuals) / area_to_normalize
                weights_dip = np.ones_like(dip_residuals) / area_to_normalize

                x = np.linspace(-max_resid * 1.2, max_resid * 1.2, 1000)

                norm_hist, _, _ = plt.hist(truncated_residuals, label='on-curve residuals', bins=residual_array,
                                           weights=weights_truncated, color='green', alpha=0.8)
                plt.hist(dip_residuals, label='dip residuals', bins=residual_array, weights=weights_dip, color='maroon',
                         alpha=0.8)
                plt.plot(x, sp.stats.norm.pdf(x, loc=0, scale=std_truncated_residuals), label='normal distribution',
                         linewidth=2.5, color='darkblue')
                plt.axvline(x=mean_truncated_residuals_data, color='black', linestyle='--', label='mean res. on-curve')
                plt.axvline(x=mean_dip_residuals, color='red', linestyle='--', label='mean res. dip')
                plt.xlim(-1.2 * max_resid, 1.2 * max_resid)
                plt.ylim(0, 1.2 * max(norm_hist))
                plt.title('Histogram of residual planet flux')
                plt.xlabel(r"residual planet flux $F_\lambda$ [ph $\mathrm{s}^{-1}$m$^{-2}\mu \mathrm{m}^{-1}$]")
                plt.ylabel('normalized counts []')
                plt.legend(loc='best')
                plt.grid()
                plt.show()

                y = np.linspace(-3 * critical_value, 3 * critical_value, 1000)
                t_distr = t_distribution(y, dof)


                plt.plot(y, t_distr, label='t-distribution', linewidth=2, color='navy')
                plt.title('t-distribution')
                plt.axvline(x=critical_value, color='black', linestyle='--', label='critical value')
                plt.text(0.015, 0.95, r'$\alpha$=' + str(alpha), transform=plt.gca().transAxes,
                         bbox=dict(facecolor='white', edgecolor='black', boxstyle='square', pad=0.5))
                plt.legend(loc='best')
                plt.grid()
                plt.show()

        if (np.any(dip_extracted) == True):
            min_t_score = np.argmin(t_scores)

            popt_extr = popts[min_t_score]
            pcov_extr = pcovs[min_t_score]
            perr_extr = perrs[min_t_score]

        else:
            # add the case for no dips
            popt_extr, pcov_extr = sp.optimize.curve_fit(BB_for_fit, xdata=self.wl_bins_and_distS_array,
                                                         ydata=spectra,
                                                         sigma=sigmas, p0=(300, 1.), absolute_sigma=True,
                                                         maxfev=10000)

            perr_extr = np.sqrt(np.diag(pcov_extr))

        return t_scores, popt_extr, pcov_extr, perr_extr


    def get_likelihood(self, depth, width, molecule, Tp, Rp, spectra, sigmas):
        '''
        This function calculates the likelihood function according to Bayesian model statistics, see the written thesis

        :param depth: float; depth of the wavelength dip induced by each molecule as a fraction of the blackbody curve
                                of the planet []
        :param width: float; depth of the wavelength dip induced by each molecule as a fraction of the blackbody curve
                                of the planet []
        :param molecule: str; name of the molecule considered
        :param Tp: float; planet temperature [K]
        :param Rp: float; planet radius [R_earth]
        :param spectra: np.ndarray of size L; contains the extracted planet flux in each wl-bin in [photons]
        :param sigmas: np.ndarray of size L; contains the sigmas of the extracted planet flux in each wl-bin [photons]

        :return likelihood: float; likelihood function as calculated with the input parameters []
        '''

        expected_spectra = self.atmospheric_scenario.expected_spectrum(molecule,
                                                                       self.wl_bins, self.wl_bins_and_distS_array,
                                                                       Tp, Rp, width, depth)

        likelihood = np.prod(
            1 / (np.sqrt(2 * np.pi) * sigmas) * np.exp(-1 / 2 * (spectra - expected_spectra) ** 2 / sigmas ** 2))

        return likelihood


    def get_bayes_factors(self, Tp, Rp, spectra, sigmas):
        '''
        This function calculates the Bayes factor K for each of the molecules considered. See the written thesis for
        the theoretical derivations of the calculations

        :param Tp: float; planet temperature [K]
        :param Rp: float; planet radius [R_earth]
        :param spectra: np.ndarray of size L; contains the extracted planet flux in each wl-bin in [photons]
        :param sigmas: np.ndarray of size L; contains the sigmas of the extracted planet flux in each wl-bin [photons]

        :return bayes_factors: np.ndarray of size self.n_molecules; contains the Bayes factors K []
        '''

        # calculate the evidence for the ground model (no atmosphere)
        expected_spectra_ground = BB_for_fit(self.wl_bins_and_distS_array, Tp, Rp)

        evidence_ground = np.prod(
            1 / (np.sqrt(2 * np.pi) * sigmas) * np.exp(-1 / 2 * (spectra - expected_spectra_ground) ** 2 / sigmas ** 2))


        # calculate the evidences for each of the molecules
        bayes_factors = np.empty((self.n_molecules))

        for i in range(bayes_factors.size):
            evidence, _ = sp.integrate.dblquad(self.get_likelihood, self.atmospheric_scenario.dip_widths[i][0],
                                               self.atmospheric_scenario.dip_widths[i][1],
                                               self.atmospheric_scenario.dip_depths[i][0],
                                               self.atmospheric_scenario.dip_depths[i][1],
                                               args=(
                                               self.atmospheric_scenario.dip_molecules[i], Tp, Rp, spectra, sigmas))
            norm_factor_depth = 1 / (
                        self.atmospheric_scenario.dip_depths[i][1] - self.atmospheric_scenario.dip_depths[i][0])
            norm_factor_width = 1 / (
                        self.atmospheric_scenario.dip_widths[i][1] - self.atmospheric_scenario.dip_widths[i][0])

            normed_evidence = evidence * norm_factor_depth * norm_factor_width

            # for extreme values, set to +/-100 to avoid runtime issues
            with warnings.catch_warnings():
                warnings.filterwarnings("error")
                try:
                    bayes_factors[i] = normed_evidence / evidence_ground
                except RuntimeWarning:
                    bayes_factors[i] = 100
                except RuntimeError:
                    bayes_factors[i] = 100

        return bayes_factors


    def get_T_R_estimate(self, spectra, sigmas):
        '''
        This function obtains the planet temperature and radius by fitting the input spectra and sigmas to blackbody
        curves using scipy.optimize functions

        :param spectra: np.ndarray of size L; contains the extracted planet flux in each wl-bin in [photons]
        :param sigmas: np.ndarray of size L; contains the sigmas of the extracted planet flux in each wl-bin [photons]

        :return popt: np.ndarray of size 2; optimal values for (T,R) such that the sum of the squared residuals to the
                        blackbody curve is minimized. In [T,R_earth]
        :return pcov: np.ndarray of shape (2,2); covariance matrix of the obtained popt, a measure of uncertainty.
                        Dimensions [T^2] and [R_earth^2] on the first and second diagonal element
        :return perr: np.ndarray of size 2; standard deviation of the errors on T and R. In [T,R_earth]
        :return t_score: np.ndarray of size n_molecules; contains the t-score values for each of the
                        considered molecules []
        :return bayes_factors: np.ndarray of size self.n_molecules; contains the Bayes factors K []
        '''

        # combine the input parameters dist_s and wl_bins into one matrix (required to properly call
        #   sp.optimize.curve_fit())
        dist_s_array = np.full(self.L, self.single_data_row['distance_s'])
        self.wl_bins_and_distS_array = np.vstack((self.wl_bins, dist_s_array))

        if (self.include_dips == True):
            t_score, popt, pcov, perr = self.detect_dips(spectra, sigmas)
            bayes_factors = self.get_bayes_factors(popt[0], popt[1], spectra, sigmas)

        else:
            # perform the fit of T and R and calculate the uncertainties
            popt, pcov = sp.optimize.curve_fit(BB_for_fit, xdata=self.wl_bins_and_distS_array, ydata=spectra,
                                               sigma=sigmas, p0=(300, 1.), absolute_sigma=False, maxfev=10000)

            perr = np.sqrt(np.diag(pcov))

            t_score = np.array([1.0])
            bayes_factors = np.array([0])

        # plot the different quantities as defined
        if (self.plot == True):
            # get the blackbody curve for the estimated parameters
            Fp_fit = BB_for_fit(self.wl_bins_and_distS_array, popt[0], popt[1])

            # in the old version, the snr per wavelength according to photon statistics was saved; this is no longer the
            #   case in the new lifesim version, therefor this option is no longer available. Can be re-implemented if
            #   the data is saved again
            snr_photon_stat = None
            Fp_BB = Fp_fit

            # get the ideal blackbody flux of the planet
            Fp = black_body(mode='planet',
                            bins=self.wl_bins,
                            width=self.wl_bin_widths,
                            temp=self.single_data_row['temp_p'],
                            radius=self.single_data_row['radius_p'],
                            distance=self.single_data_row['distance_s']) / self.wl_bin_widths * 10 ** -6


            # plot the measured fluxes as well as the best fit and true blackbody curve
            plot_planet_SED_and_SNR(self.wl_bins, Fp, spectra, sigmas, self.min_wl,
                                    self.max_wl, Fp_BB=Fp_BB, snr_photon_stat=snr_photon_stat, filename=None)

        return popt, pcov, perr, t_score, bayes_factors


    def cost_func_MAP(self):
        '''
        This function calculates the various extracted quantities based on the maximum likelihood method as described
        in LIFE II and the Thesis. In a first step, the auxiliary matrices B&C are calculated, from which the estimated
        flux is calculated. Based in this the cost function is derived, from which the planet position is inferred.
        Also, the planet temperature and radius are extracted by fitting the extracted fluxes to a blackbody curve. To
        calculate the signal-to-noise ratios, the signal is first whitened (if the threshold self.whitening_limit
        is obtained), i.e. the process is repeated but taking the variance only of the noise as opposed to the whole
        signal. Using the whitened signal, the SNR is calculated according to three different methods as described in
        the Thesis

        :return Jmax: float; maximum cost function value (of all pixels) []
        :return r: float; radial coordinate of the maximum cost function value in [radial pixel]
        :return p: int; azimuthal coordinate of the maximum cost function value in [azimuthal pixel]
        :return F_est: np.ndarray of size L; contains the estimated planet flux received at every
                    wavelength at the theta_max pixel coordinates in [photons]
        :return F_est_pos: np.ndarray of size L; contains the positive part of the estimated planet flux received
                    at every wavelength at the theta_max pixel coordinates in [photons]
        :return sigma_est: np.ndarray of size L; contains the sigmas for the extracted spectra per wavelength bin
                    in [photons]
        :return SNR_est: float; total SNR value calculated according to the naive method (see Thesis) over all
                                    wave lengths []
        :return FPR_sigma: float; extracted false positive detection rate (in terms of standard deviations) over the
                                    entire wavelength range as calculated using the true planet position (see Thesis) []
        :return FPR_max_sigma: float; extracted false positive detection rate (in terms of standard deviations) over the
                                    entire wavelength range as calculated using the maximum J'' position (see Thesis) []
        :return T_est: float; extracted planet temperature in [K]
        :return T_sigma: float; uncertainty in the extracted planet temperature [K]
        :return R_est: float; extracted planet radius in [R_earth]
        :return R_sigma: float; uncertainty in the extracted planet radius [R_earth]
        :return t_score: np.ndarray of size n_molecules; contains the t-score values for each of the
                        considered molecules []
        :return bayes_factors: np.ndarray of size self.n_molecules; contains the Bayes factors K []
        '''

        # calculate the matrices B and C for the signal as recorded (unwhitened)
        self.get_B_C()

        # calculate the estimated flux based on the total received signal
        F, F_pos = self.get_F_estimate()

        # calculate the cost function and as well as the position of the maximum value of the cost function
        self.J = (F_pos * self.C).sum(axis=0)
        theta_max = np.unravel_index(np.argmax(self.J, axis=None), self.J.shape)
        (r, p) = theta_max

        # determine the true values of r and phi in order to calculate the value of J at this pixel
        true_r = int(self.single_data_row['angsep'] / 2 / self.field_of_search * self.image_size / 180 / 3600 * np.pi)
        true_phi = 0

        # calculate the estimated flux at the position of Jmax (total and positive part)
        F_est = F[:, r, p]
        F_est_pos = F_pos[:, r, p]

        # calculate sigma at the position of Jmax
        sigma_est = self.B[:, r, p] ** (-1 / 2)

        # calculate the best fit temperature and radius along with the corresponding uncertainties
        try:
            popt, pcov, perr, t_score, bayes_factors = self.get_T_R_estimate(F_est, sigma_est)

        # in very rare instances, the curve fitting does not work; in this case print a warning and continue with the
        #   true values and small uncertainty to avoid the program crashing
        except RuntimeError:
            print('T and R not found')
            popt = np.array([self.single_data_row['temp_p'], self.single_data_row['radius_p']])
            perr = np.array([0.1, 0.1])
            t_score = np.array([1.0])
            bayes_factors = np.array([0.0])

        T_est = popt[0]
        R_est = popt[1]
        T_sigma = perr[0]
        R_sigma = perr[1]

        # determine the SNR before whitening
        SNR_est_before_whitening = np.sum((F_est_pos / sigma_est) ** 2) ** (1 / 2)

        # determine whether the whitening limit is obtained
        if (SNR_est_before_whitening >= self.whitening_limit):
            whitening = True
        else:
            whitening = False

        if (self.plot == True):
            print('Signal whitening: ', whitening)

        # if applicable, the signal is whitened in the following
        if (whitening == True):

            # calculate the blackbody flux of the planet according to the estimated parameters
            est_flux = black_body(mode='planet',
                                  bins=self.wl_bins,
                                  width=self.wl_bin_widths,
                                  temp=T_est,
                                  radius=R_est,
                                  distance=self.single_data_row['distance_s']
                                  )

            # calculate the estimated flux retrieved from the planet with the assumed parameters
            _, self.est_ideal_signal = self.run_socket(s_name='instrument',
                                                       method='get_signal',
                                                       temp_s=self.single_data_row['temp_s'],
                                                       radius_s=self.single_data_row['radius_s'],
                                                       distance_s=self.single_data_row['distance_s'],
                                                       lat_s=self.single_data_row['lat'],
                                                       z=self.single_data_row['z'],
                                                       angsep=2 * r * self.field_of_search / self.image_size
                                                              * 180 * 3600 / np.pi,
                                                       flux_planet_spectrum=[self.wl_bins * u.meter,
                                                                             est_flux / self.wl_bin_widths * u.photon /
                                                                             u.second / (u.meter ** 3)],
                                                       integration_time=self.single_data_row['int_time'],
                                                       phi_n=self.n_steps)

            # calculate the noise by subtracting the (estimated) planet signal from the total signal
            self.noise = self.signals - self.est_ideal_signal

            # recalculate the B and C matrices with the whitened signal
            self.get_B_C_whitened()

            # recalculate the signals using the whitened noise
            F_white, F_white_pos, F_noise, F_noise_pos = self.get_F_estimate_whitened()


        else:
            # this noise is based off ground truth data; however, the quantities derived from it are only used for the
            #   plot that shows the histogram of a theoretical purely noisy signal and do not enter the extracted
            #   quantities in any way
            self.noise = self.signals - self.ideal_signals
            self.get_B_C_whitened()
            _, _, _, F_noise_pos = self.get_F_estimate_whitened()

            # these quantities no longer include ground truth data and are thus set to the un-whitened values
            F_white_pos = F_pos

            self.B_noise = self.B
            self.C_noise_var = self.C

        # take the signal at the derived position
        F_white_pos_max = F_white_pos[:, r, p]

        # recalculate the sigma for the whitened noise
        sigma_est = self.B_noise[:, r, p] ** (-1 / 2)

        # calculate the total SNR over all wavelength bins (method 1)
        SNR_est = np.sum((F_white_pos_max / sigma_est) ** 2) ** (1 / 2)

        # calculate the cost function using the variance only of the noise
        self.J_FPR = (F_pos * self.C_noise_var).sum(axis=0)
        Jmax = np.max(self.J_FPR)

        # if the true position of the planet is outside of the FoS, set J_true_pos to zero
        if (true_r <= self.radial_ang_px-1):
            J_true_pos = self.J_FPR[true_r, true_phi]
        else:
            J_true_pos = 0

        J_max_pos = self.J_FPR[r, p]

        # set the precision to 100 and the false positive rate to one to begin the loop
        precision = 100
        FPR = 1.0
        FPR_max = 1.0

        # calculate the false positive rate at the true planet position. If the calculation fails due to lack of
        #   precision, double the precision and try again
        while ((FPR == 1.0 or FPR_max == 1.0) and precision <= self.precision_limit):
            FPR = cdf_J_precision(self.L, J_true_pos, precision)
            FPR_max = cdf_Jmax_precision(self.L, J_max_pos, precision, self.single_data_row['angsep'])

            precision *= 2

        # if the FPR could not be calculated in the step above, set the FPR_sigmas to 10000 for later reference
        if (precision > self.precision_limit):
            if (self.plot == True):
                print('warning: FPR_sigma could not be calculated')
            FPR_sigma = 10000
            FPR_max_sigma = 10000

        # if SNR_ps is smaller than 20, calculate FPR_sigma directly using the conversion formula (note: this is using
        #   ground truth data, but its only to speed up calculation; in practise, one can just always use one method
        #   at longer runtime cost)
        elif (self.single_data_row['snr_current'] <= 20):
            FPR_sigma = mp.mp.sqrt(2) * mp.erfinv(2 * FPR - 1)
            FPR_sigma = float(FPR_sigma)

            FPR_max_sigma = mp.mp.sqrt(2) * mp.erfinv(2 * FPR_max - 1)
            FPR_max_sigma = float(FPR_max_sigma)

        # calculate FPR_sigma using a lookup table (speed-up for high values)
        else:
            FPR_sigma = alt_sigma_calc(FPR, self.filepath)
            FPR_max_sigma = alt_sigma_calc(FPR_max, self.filepath)

        if (self.plot == True):
            # plot the cost function heatmap
            j_map = pol_to_cart_map(self.J_FPR, self.image_size)
            plot_multi_map(j_map, "Cost Value", self.field_of_search * 3600000 * 180 / np.pi,
                           "inferno", filename_post=None)


            # calculate the J value one would get by a purely noisy signal
            self.J_noise = (F_noise_pos * self.C_noise_only).sum(axis=0)

            # plot a histogram of the values of J (whitened version)
            flat_J = self.J_FPR.flatten()
            n_bins = 100
            j_array = np.linspace(0, 2 * np.max(flat_J), n_bins)

            # plot a histogram of the values of J considering only the noisy part of the signal
            flat_J_noise = self.J_noise.flatten()
            j_array_noise = np.linspace(0, 2 * np.max(flat_J_noise), n_bins)

            counts_to_norm, bin_edges = np.histogram(flat_J_noise, bins=j_array_noise)
            area_to_normalize = np.sum(counts_to_norm * (bin_edges[1] - bin_edges[0]))

            weights_noise = np.ones_like(flat_J_noise) / area_to_normalize
            weights_signal = np.ones_like(flat_J) / area_to_normalize

            # calculate the detections threshold values
            eta_5 = get_detection_threshold(self.L, 5)
            eta_max_5 = get_detection_threshold_max(self.L, 5, self.single_data_row['angsep'])

            # calculate the theoretical pdf and cdf of the J'' function
            pdf_J2prime = pdf_J(self.L, j_array_noise)
            cdf_J2prime = np.empty_like(j_array_noise)
            cdf_Jmax = np.empty_like(j_array_noise)

            for i in range(j_array_noise.size):
                cdf_J2prime[i] = cdf_J_precision(self.L, j_array_noise[i], 100)
                cdf_Jmax[i] = cdf_Jmax_precision(self.L, j_array_noise[i], 100, self.single_data_row['angsep'])

            # plot with the noisy values, p(J'') and the detection threshold
            counts_noise, bins_noise, _ = plt.hist(flat_J_noise, j_array_noise, weights=weights_noise,
                                                   label='measured J\u2032\u2032 values noise',
                                                   color='darkblue', rwidth=0.75)
            plt.plot(j_array_noise, pdf_J2prime, label='p(J\u2032\u2032) theoretical', color='maroon')
            plt.title('J\u2032\u2032 pure noisy signal')
            plt.axvline(x=eta_5, color='red', linestyle='--', label='det. thres.')
            plt.legend(loc='best')
            plt.xlabel('Cost function J\u2032\u2032')
            plt.ylabel('Normalized probability density')
            plt.grid()
            # Uncomment the following line to save plot
            # plt.savefig(self.filepath+'J_plot1.pdf')
            plt.show()

            # plot with additionally including the cdfs and the detection threshold for Jmax
            counts_noise, bins_noise, _ = plt.hist(flat_J_noise, j_array_noise, weights=weights_noise,
                                                   label='measured J\u2032\u2032 values noise',
                                                   color='darkblue', rwidth=0.75)
            plt.plot(j_array_noise, pdf_J2prime, label='p(J\u2032\u2032) theoretical', color='maroon')
            plt.plot(j_array_noise, cdf_J2prime, label=r'''$\Phi(J'')$''', color='red')
            plt.plot(j_array_noise, cdf_Jmax, label=r'''$\Phi_{max}(J'')$''', color='black')
            plt.title('J\u2032\u2032 pure noisy signal')
            plt.axvline(x=eta_5, color='red', linestyle='--', label='det. thres.')
            plt.axvline(x=eta_max_5, color='black', linestyle='--', label='det. thres. max')
            plt.legend(loc='best')
            plt.xlabel('Cost function J\u2032\u2032')
            plt.ylabel('Normalized probability density')
            plt.grid()
            # Uncomment the following line to save plot
            # plt.savefig(self.filepath+'J_plot2.pdf')
            plt.show()

            # plot additionally showing the histogram of the J'' values of the whitened signal
            counts, bins, _ = plt.hist(flat_J, j_array, weights=weights_signal,
                                       label='measured J\u2032\u2032 values signal', color='gold', rwidth=0.75)
            counts_noise, bins_noise, _ = plt.hist(flat_J_noise, j_array_noise, weights=weights_noise,
                                                   label='measured J\u2032\u2032 values noise', color='darkblue',
                                                   rwidth=0.75)
            plt.plot(j_array_noise, pdf_J2prime, label='p(J\u2032\u2032) theoretical', color='maroon')
            plt.title('J\u2032\u2032 signal including planet')
            plt.axvline(x=eta_5, color='red', linestyle='--', label='det. thres.')
            # plt.axvline(x=eta_max_5, color='black', linestyle='--', label='det. thres. max')
            plt.legend(loc='best')
            plt.xlabel('Cost function J\u2032\u2032')
            plt.ylabel('Normalized probability density')
            plt.grid()
            plt.xlim((0, 1.5 * np.max(flat_J)))
            plt.ylim((0, 20 * 1 / (self.image_size / 2) ** 2))
            # Uncomment the following line to save plot
            # plt.savefig(self.filepath+'J_plot3.pdf')
            plt.show()

        return Jmax, r, p, F_est, F_est_pos, sigma_est, SNR_est, FPR_sigma, FPR_max_sigma, T_est, T_sigma, \
            R_est, R_sigma, t_score, bayes_factors


    def single_spectrum_extraction(self, n_run=1):
        '''
        This is the core function of the ML_Extraction class that executes the signal extraction for a single planet.
        In each run, the extracted spectrum, snr, sigma, planet positions and parameters are calculated based
        on the parameters in the bus catalog and random noise generation. A single run goes through the following
        process:

        In a first step, the signal (including noise) is generated by calling the instrument.get_signals function.
        Next, get_transm_map() is called to generate the transmission maps and the matrices required to calculate the
        estimated planet flux and the cost function J.
        The third step consists of calling cost_func_MAP to get the estimated planet flux as well as the cost functions,
        the positions corresponding to the maxima of the cost function and the extracted planet temperature and radius.
        The signal is then whitened and the signal-to-noise ratio is calculated with three different methods (as
        described in the Thesis). All the quantities of the run are finally saved in a list

        :param n_run: integer; number of runs to perform []

        :return extracted spectra: np.ndarray of shape (n_run,L); contains the extracted spectra per wavelength bin for
                        each of the runs in [photons]
        :return extracted snrs: np.ndarray of size n_run; contains the extracted snrs for each of the runs over the
                        entire wavelength range as calculated using the naive approach (see Thesis) []
        :return extracted sigmas: np.ndarray of shape (n_run,L); contains the sigmas for the extracted spectra
                        per wavelength bin for each of the runs in [photons]
        :return extracted Jmaxs: np.ndarray of size n_run containing the maximum values of J over all wavelengths
                        for every run []
        :return rss: np.ndarray of size n_run; contains the extracted angular separations for each of the runs
                        in [arcsec]
        :return: phis: np.ndarray of size n_run; contains the extracted azimuthal positions on the image for each of
                        the runs in [degrees]
        :return Ts: np.ndarray of size n_run; extracted planet temperatures of every run in [K]
        :return Ts_sigma: np.ndarray of size n_run; uncertainties of the extracted planet temperatures of every run [K]
        :return Rs: np.ndarray of size n_run; extracted planet radii of every run in [R_earth]
        :return Rs_sigma: np.ndarray of size n_run; uncertainties of the extracted planet radii of every run [R_earth]
        :return FPRs: np.ndarray of size n_run; contains the extracted false positive detection rate for the extracted
                        planet for every run in [number of sigmas] when taking the true position
        :return FPR_maxs: np.ndarray of size n_run; contains the extracted false positive detection rate for the
                        extracted planet for every run in [number of sigmas] when taking the maximum position
        :return induced_dipss: np.ndarray of shape (n_runs,n_molecules); denotes for each molecule whether a
                        corresponding dip was induced in each run []
        :return t_scores: np.ndarray of shape (n_runs,n_molecules); contains the t-score values for each of the
                        considered molecules in each run []
        :return SNR_ps_news: np.ndarray of size n_runs; new SNR_ps values for the planet after a potential dip
                        was induced for each run []
        :return bayes_factors: np.ndarray of shape (n_runs,n_molecules); contains the Bayes factors K for each
                        molecule in each run []
        '''

        extracted_spectra = []
        extracted_snrs = []
        extracted_sigmas = []
        extracted_Jmaxs = []
        rss = []
        phiss = []
        Ts = []
        Ts_sigma = []
        Rs = []
        Rs_sigma = []
        FPRs = []
        FPR_maxs = []
        induced_dipss = []
        t_scores = []
        SNR_ps_news = []
        bayes_factors = []

        for n in range(n_run):
            if (self.plot == True):
                print('run:', n)

            if (self.include_dips == True):
                induced_dips, SNR_ps_new = self.create_dips()

            else:
                induced_dips = np.array([0])
                SNR_ps_new = 0

            # get the signals (with and without noise) for the specified constellation
            self.signals, self.ideal_signals = self.run_socket(s_name='instrument',
                                                               method='get_signal',
                                                               temp_s=self.single_data_row['temp_s'],
                                                               radius_s=self.single_data_row['radius_s'],
                                                               distance_s=self.single_data_row['distance_s'],
                                                               lat_s=self.single_data_row['lat'],
                                                               z=self.single_data_row['z'],
                                                               angsep=self.single_data_row['angsep'],
                                                               angle_p=self.single_data_row['angle_p'],
                                                               flux_planet_spectrum=[self.wl_bins * u.meter,
                                                                                     self.single_data_row[
                                                                                         'planet_flux_use'][
                                                                                         0] / 3600 /
                                                                                     self.data.inst[
                                                                                         'telescope_area'] /
                                                                                     self.data.inst[
                                                                                         'eff_tot'] / self.wl_bin_widths
                                                                                     * u.photon / u.second /
                                                                                     (u.meter ** 3)],
                                                               integration_time=self.single_data_row['int_time'],
                                                               phi_n=self.n_steps)

            # create the transmission maps
            self.get_transm_map()

            # extract all the quantities
            Jmax, r, p, Fp_est, Fp_est_pos, sigma_est, SNR_est, FPR_extr, FPR_max_extr, T_est, T_sigma, R_est, \
                R_sigma, t_score, bayes_factor = self.cost_func_MAP()

            # add the quantities from this run to the final list
            extracted_spectra.append(Fp_est)
            extracted_snrs.append(SNR_est)
            extracted_sigmas.append(sigma_est)
            extracted_Jmaxs.append(Jmax)
            rss.append(2 * r * self.field_of_search / self.image_size * 180 * 3600 / np.pi)
            phiss.append(p * self.n_steps / 360)
            Ts.append(T_est)
            Rs.append(R_est)
            Ts_sigma.append(T_sigma)
            Rs_sigma.append(R_sigma)
            FPRs.append(FPR_extr)
            FPR_maxs.append(FPR_max_extr)
            induced_dipss.append(induced_dips)
            t_scores.append(t_score)
            SNR_ps_news.append(SNR_ps_new)
            bayes_factors.append(bayes_factor)

        # convert the final lists to arrays
        extracted_spectra = np.array(extracted_spectra)
        extracted_snrs = np.array(extracted_snrs)
        extracted_sigmas = np.array(extracted_sigmas)
        extracted_Jmaxs = np.array(extracted_Jmaxs)
        rss = np.array(rss)
        phiss = np.array(phiss)
        Ts = np.array(Ts)
        Ts_sigma = np.array(Ts_sigma)
        Rs = np.array(Rs)
        Rs_sigma = np.array(Rs_sigma)
        FPRs = np.array(FPRs)
        FPR_maxs = np.array(FPR_maxs)
        induced_dipss = np.array(induced_dipss)
        t_scores = np.array(t_scores)
        SNR_ps_news = np.array(SNR_ps_news)
        bayes_factors = np.array(bayes_factors)

        return extracted_spectra, extracted_snrs, extracted_sigmas, extracted_Jmaxs, rss, phiss, Ts, Ts_sigma, Rs, \
            Rs_sigma, FPRs, FPR_maxs, induced_dipss, t_scores, SNR_ps_news, bayes_factors


    def main_parameter_extraction(self, n_run=1, mu=0, whitening_limit=0, n_processes=1, precision_limit=1600,
                                  max_FoS=False, plot=False, single_planet_mode=False, planet_number=0,
                                  include_dips=False, atmospheric_scenario=None, save_mode=True,
                                  filepath=None):
        '''
        The main_parameter_extraction function is the function that should get called by other files. In defines all the
        required parameters and then runs single_spectrum_extraction for either one specified planet
        (single_planet_mode=True) or for all of the planets in the catalog. In single_planet_mode=True, the extracted
        quantities are returned. If =False, then a number of processes equal to n_processes are created to run the
        extraction of the planets in the catalog of the bus in parallel. The output is directly adjusted and saved
        (provided save_mode=True) in path/changeme.csv

        :param n_run: int; number of runs to perform for each planet []
        :param mu: float; regularization parameter for the calculation of the cost function J as described in
                    section 3.1 of LIFE II []
        :param whitening_limit: float; threshold above which a signal with SNR=whitening_limit should be whitened []
        :param n_processes: int; number of processes to run in parallel for multi-planet extraction (ignored if
                    single_planet_mode=True) []
        :param precision_limit: int; determines the maximum precision the calculation of the false positive rate
                    will use before deeming that it is high enough and set to 10000 []
        :param max_FoS: boolean; if True, the extraction is performed using the maximum FoS, i.e. 10x the outer
                    habitable zone angular separation. Otherwise, the outer habitable zone is used
        :param plot: boolean; determines whether to show plots throughout the runs (only advised for small n_run)
        :param single_planet_mode: boolean; if True, signal extraction is performed for only one planet in the catalog
                    as defined by the following parameter
        :param planet_number: integer; index of the planet in the catalog of the bus of which the signal is to be
                    extracted. Only applicable if single_planet_mode=True
        :param include_dips: boolean; if True, atmospheric dips are included, otherwise perfect blackbodies are assumed
        :param atmospheric scenario: class object; determines which atmospheric scenario initialized in
                    Extraction_auxiliary.py is used. Only applicable if include_dips = True
        :param save_mode: boolean; if True, the catalog with added columns for the extracted parameters will be saved
                    in path/changeme.csv. Only applicable if single_planet_mode=False
        :param filepath: str; path to where the the file will be saved (only applicable if save_mode=True and
                    single_planet_mode=False)

        :return spectra, snrs, sigmas, Jmaxs, rss, phiss, Ts, Ts_sigma, Rs, Rs_sigma, FPRs, FPR_maxs, induced_dips,
                    t_scores, SNR_ps_news, bayes_factors; Only applicable if single_planet_mode=True, else nothing
                    is returned. See single_spectrum_extraction for description

        Attributes
        -----------
        Most of the attributes are simply abbreviations of the respective objects inherited from the bus. These are:
        - 'n_planets': int; Number of planets in the catalog []
        - 'L': int; Number of wavelength bins as given by the wavelength range and the resolution parameter R []
        - 'min_wl': float; minimum wavelength captured by the instrument [m]
        - 'max_wl': float; maximum wavelength captured by the instrument [m]
        - 'wl_bins': np.ndarray of size L; consists of all the wavelength bins (middle wavelength of each bin) in [m]
        - 'wl_bin_widths': np.ndarray of size L; consists of all the widths of the wavelength bins in [m]
        - 'image_size': int; precision that the image can be resolved to in one dimension in [number of pixels]
        - 'radial_ang_px': int; precision that the image can be resolved to in the radial coordinate in
                                [number of pixels]
        - 'hfov': np.ndarray of size L; consists of the half field of views of the instrument in each wavelength bin
                                in [radians]
        - 'mu': float; regularization parameter for the calculation of the cost function J as described above []
        - 'whitening_limit': float; threshold above which a signal with SNR=whitening_limit should be whitened []
        - 'n_run': int; number of runs to perform for each planet []
        - 'precision_limit': int; maximum precision the calculation of the false positive rate []
        - 'max_FoS': boolean; determines whether to use the maximum FoS or not
        - 'filepath': str; path to where the files should be saved
        - 'plot': boolean; determines whether to generate plots during the extraction
        - 'include_dips': boolean; determines whether to include dips in the blackbody spectra caused by molecules
        - 'atmospheric scenario': class object; determines the atmospheric scenario to be used if required
        - 'n_molecules': int; number of molecular dips considered

        Two attributes concerning the angular dimensioning of the image:
        - 'planet_azimuth': Angular position of the planet on the image. Contrary to the radial position, which is
                            given by the angular separation, this attribute is not known a priori from the simulation
                            and thus set to 0 [degrees] for simplicity. Changing this or making it adaptable would
                            require some modification of the transmission functions to comply with the rest of the code
        - 'n_steps':        Resolution of the angular coordinate. Given as an integer, e.g. 360 corresponds to the
                            angular resolution being split into 360 parts, i.e. one part has a size of 1 degree []

        The following attributes are planet-specific and must be redefined for each new planet:
        - 'single_data_row': pd.series; catalog data from the bus corresponding to the planet described by planet_number
        - 'field_of_search': float; Half field of view used for the calculation of the cost map in [radian]. This value
                            is set based upon the outer habitable zone angular separation and the value of max_FoS
        '''

        self.n_planets = len(self.data.catalog.index)
        self.L = self.data.inst['wl_bins'].size
        self.min_wl = self.data.inst['wl_bin_edges'][0]
        self.max_wl = self.data.inst['wl_bin_edges'][-1]
        self.wl_bins = self.data.inst['wl_bins']
        self.wl_bin_widths = self.data.inst['wl_bin_widths']
        self.image_size = self.data.options.other['image_size']
        self.radial_ang_px = int(self.image_size / 2)
        self.hfov = self.data.inst['hfov']
        self.mu = mu
        self.whitening_limit = whitening_limit
        self.n_run = n_run
        self.precision_limit = precision_limit
        self.max_FoS = max_FoS
        self.filepath = filepath
        self.plot = plot
        self.include_dips = include_dips
        self.atmospheric_scenario = atmospheric_scenario
        self.n_molecules = len(self.atmospheric_scenario.dip_molecules)

        self.planet_azimuth = 0
        self.n_steps = 360

        # if in single_planet_mode, define the planet-specific parameters, run single_spectrum_extraction once and
        #   return the extracted parameters
        if (single_planet_mode == True):
            self.single_data_row = self.data.catalog.iloc[planet_number]
            if (self.max_FoS == True):
                #self.field_of_search = 0.3 * np.pi / 3600 / 180
                self.field_of_search = 10 * self.single_data_row['hz_out'] / self.single_data_row['distance_s'] * \
                                       np.pi / 3600 / 180
            else:
                self.field_of_search = self.single_data_row['hz_out'] / self.single_data_row['distance_s'] * np.pi / \
                                       3600 / 180


            spectra, snrs, sigmas, Jmaxs, rss, phiss, Ts, Ts_sigma, Rs, Rs_sigma, FPRs, FPR_maxs, induced_dips, \
                t_scores, SNR_ps_news, bayes_factors \
                = self.single_spectrum_extraction(n_run=n_run)

            return spectra, snrs, sigmas, Jmaxs, rss, phiss, Ts, Ts_sigma, Rs, Rs_sigma, FPRs, FPR_maxs, \
                induced_dips, t_scores, SNR_ps_news, bayes_factors


        # if single_planet_mode=False, proceed here
        else:

            # create lists to store extracted data
            extracted_spectra_tot = []
            extracted_snrs_tot = []
            extracted_sigmas_tot = []
            extracted_Jmaxs_tot = []
            extracted_rss_tot = []
            extracted_phiss_tot = []
            extracted_Ts_tot = []
            extracted_Ts_sigma_tot = []
            extracted_Rs_tot = []
            extracted_Rs_sigma_tot = []
            extracted_FPRs_tot = []
            extracted_FPR_maxs_tot = []
            extracted_induced_dips_tot = []
            extracted_t_scores_tot = []
            extracted_SNR_ps_news_tot = []
            extracted_bayes_factors_tot = []

            # divide the planets into equal ranges for parallel processing
            n_processes = n_processes
            planet_indices = []

            for i in range(n_processes):
                lower_index = int(np.floor(self.n_planets / n_processes * i))
                if (i == n_processes):
                    upper_index = int(self.n_planets)
                else:
                    upper_index = int(np.floor(self.n_planets / n_processes * (i + 1)))
                index_range = np.arange(lower_index, upper_index, 1)
                planet_indices.append(index_range)

            # define processes, queues and events:
            #   the processes are the objects which will run the extraction in parallel
            #   the queues are where the processes store their output
            #   the process has an event attributed to it; it's function is to be 'set' when the process is completed so
            #       that the main code can continue as soon as all events are set, i.e. all processes are finished
            processes = []
            events = []
            res_queue = multiprocessing.Queue()
            num_queue = multiprocessing.Queue()

            for i in range(n_processes):
                e = multiprocessing.Event()
                p = multiprocessing.Process(target=self.execute_multiprocessing, args=[planet_indices[i],
                                                                                       res_queue, num_queue, e, i])
                p.start()
                events.append(e)
                processes.append(p)

            # wait for all processes to finish
            for event in events:
                event.wait()

            # get the results from the queues and store them in lists. The numbers list is used to keep track of what
            #   order the queues finished in
            results = []
            numbers = []

            for i in range(n_processes):
                result = res_queue.get()
                number = num_queue.get()
                results.append(result)
                numbers.append(number)

            # add the results to the main list in the correct order
            for i in range(n_processes):
                place_in_queue = int(numbers.index(i))

                extracted_spectra_tot.append(results[place_in_queue][0])
                extracted_snrs_tot.append(results[place_in_queue][1])
                extracted_sigmas_tot.append(results[place_in_queue][2])
                extracted_Jmaxs_tot.append(results[place_in_queue][3])
                extracted_rss_tot.append(results[place_in_queue][4])
                extracted_phiss_tot.append(results[place_in_queue][5])
                extracted_Ts_tot.append(results[place_in_queue][6])
                extracted_Ts_sigma_tot.append(results[place_in_queue][7])
                extracted_Rs_tot.append(results[place_in_queue][8])
                extracted_Rs_sigma_tot.append(results[place_in_queue][9])
                extracted_FPRs_tot.append(results[place_in_queue][10])
                extracted_FPR_maxs_tot.append(results[place_in_queue][11])
                extracted_induced_dips_tot.append(results[place_in_queue][12])
                extracted_t_scores_tot.append(results[place_in_queue][13])
                extracted_SNR_ps_news_tot.append(results[place_in_queue][14])
                extracted_bayes_factors_tot.append(results[place_in_queue][15])

            # add the data to the bus catalog
            self.data.catalog['extracted_spectra'] = sum(extracted_spectra_tot, [])
            self.data.catalog['extracted_snrs'] = sum(extracted_snrs_tot, [])
            self.data.catalog['extracted_sigmas'] = sum(extracted_sigmas_tot, [])
            self.data.catalog['extracted_Jmaxs'] = sum(extracted_Jmaxs_tot, [])
            self.data.catalog['extracted_rss'] = sum(extracted_rss_tot, [])
            self.data.catalog['extracted_phiss'] = sum(extracted_phiss_tot, [])
            self.data.catalog['extracted_Ts'] = sum(extracted_Ts_tot, [])
            self.data.catalog['extracted_Ts_sigma'] = sum(extracted_Ts_sigma_tot, [])
            self.data.catalog['extracted_Rs'] = sum(extracted_Rs_tot, [])
            self.data.catalog['extracted_Rs_sigma'] = sum(extracted_Rs_sigma_tot, [])
            self.data.catalog['extracted_FPRs'] = sum(extracted_FPRs_tot, [])
            self.data.catalog['extracted_FPR_maxs'] = sum(extracted_FPR_maxs_tot, [])
            self.data.catalog['extracted_induced_dips'] = sum(extracted_induced_dips_tot, [])
            self.data.catalog['extracted_t_scores'] = sum(extracted_t_scores_tot, [])
            self.data.catalog['extracted_SNR_ps_new'] = sum(extracted_SNR_ps_news_tot, [])
            self.data.catalog['extracted_bayes_factors'] = sum(extracted_bayes_factors_tot, [])

            # save the catalog
            if (save_mode == True):
                self.data.catalog.to_csv(self.filepath + 'changeme.csv')

            print('main_parameter_extraction completed')

            return


    def execute_multiprocessing(self, process_range, res, num, event, n_process):
        '''
        This function is called by each of the parallel running processes and executes the signal extraction for the
        planets in its range

        :param process_range: np.array of size self.n_planets; contains the indices of the planets in the catalog that
                                    the process should extract
        :param res: multiprocessing.Queue; this is where the results are stored
        :param num: multiprocessing.Queue; keeps track of the order in which the processes fill up the results queue
        :param event: multiprocessing.Event; this object is 'set' as soon as the process is complete
        :param n_process: int; indicates the process number [0, n_processes] (label-like)
        '''

        # create lists to store extracted data
        extracted_spectra = []
        extracted_snrs = []
        extracted_sigmas = []
        extracted_Jmaxs = []
        extracted_rss = []
        extracted_phiss = []
        extracted_Ts = []
        extracted_Ts_sigma = []
        extracted_Rs = []
        extracted_Rs_sigma = []
        extracted_FPRs = []
        extracted_FPR_maxs = []
        extracted_induced_dips = []
        extracted_t_scores = []
        extracted_SNR_ps_news = []
        extracted_bayes_factors = []

        print('Process #', n_process, ' started')

        # loop through all of the planets in the process range
        for j in tqdm(process_range):
            self.single_data_row = self.data.catalog.iloc[j]

            if (self.max_FoS == True):
                # self.field_of_search = 0.3 * np.pi / 3600 / 180
                self.field_of_search = 10 * self.single_data_row['hz_out'] / self.single_data_row[
                    'distance_s'] * np.pi / 3600 / 180
            else:
                self.field_of_search = self.single_data_row['hz_out'] / self.single_data_row['distance_s'] * np.pi / 3600 / 180

            # call the extraction function for a single planet
            spectra, snrs, sigmas, Jmaxs, rss, phiss, Ts, Ts_sigma, Rs, Rs_sigma, FPRs, FPR_maxs, induced_dips, \
                t_scores, SNR_ps_news, bayes_factors \
                = self.single_spectrum_extraction(n_run=self.n_run)

            # store the data in the lists
            extracted_spectra.append(spectra.tolist())
            extracted_snrs.append(snrs.tolist())
            extracted_sigmas.append(sigmas.tolist())
            extracted_Jmaxs.append(Jmaxs.tolist())
            extracted_rss.append(rss.tolist())
            extracted_phiss.append(phiss.tolist())
            extracted_Ts.append(Ts.tolist())
            extracted_Ts_sigma.append(Ts_sigma.tolist())
            extracted_Rs.append(Rs.tolist())
            extracted_Rs_sigma.append(Rs_sigma.tolist())
            extracted_FPRs.append(FPRs.tolist())
            extracted_FPR_maxs.append(FPR_maxs.tolist())
            extracted_induced_dips.append(induced_dips.tolist())
            extracted_t_scores.append(t_scores.tolist())
            extracted_SNR_ps_news.append(SNR_ps_news.tolist())
            extracted_bayes_factors.append(bayes_factors.tolist())

        # add the process number and the results to the queue and set the event
        num.put(n_process)
        res.put([extracted_spectra, extracted_snrs, extracted_sigmas, extracted_Jmaxs, extracted_rss,
                 extracted_phiss, extracted_Ts, extracted_Ts_sigma, extracted_Rs, extracted_Rs_sigma, extracted_FPRs,
                 extracted_FPR_maxs, extracted_induced_dips, extracted_t_scores,
                 extracted_SNR_ps_news, extracted_bayes_factors])
        event.set()
        print('Process #', n_process, ' finished')

        return