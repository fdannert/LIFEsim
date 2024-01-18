import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sp
from scipy.special import factorial as fact
import matplotlib.ticker as mtick
from mpl_toolkits.axes_grid1 import AxesGrid
from lifesim.util.constants import c, h, k, radius_earth, m_per_pc
import mpmath as mp


### Code Author: Philipp Binkert

class Atmospheric_scenarios:
    '''
    This class defines different atmospheric scenarios. At the moment, the molecules included are CO_2, O_3 and H_20
    at different levels and probabilities of existance. Their strength and existance is not defined by concentrations
    but rather by how deep a Gaussian dip they induce into the planet atmosphere
    '''

    def __init__(self, name, likelihood_CO2, likelihood_O3, likelihood_H2O, spread):
        '''
        This function initializes an instance of the class

        :param name: str; name of the instance
        :param likelihood_CO2: float; assumed fraction of the planets in the universe that contain CO_2
        :param likelihood_O3: float; assumed fraction of the planets in the universe that contain O_3
        :param likelihood_H2O: float; assumed fraction of the planets in the universe that contain H_2O
        :param spread: float; determines the fraction with which the values for dip_widths and dip_depths can deviate
                        from the mean

        Attributes:
        ----------------
        - 'name': str; name of the atmospheric scenario
        - 'dip_molecules': list; all molecules that could be included in the scenario (as strings)
        - 'dip_centers': list; center of the wavelength dip induced by each molecule in [micron]
        - 'dip_avg_widths': list; average widths of the wavelength dip induced by each molecule in [microns]
        - 'dip_widths': list; width of the wavelength dip induced by each molecule in [micron]. Given as a np.array with
                                minimum and maximum wavelength over which the width is uniformly selected
        - 'dip_avg_depths': list; average depth of the wavelength dip induced by each molecule as a fraction []
        - 'dip_depths': list; depth of the wavelength dip induced by each molecule as a fraction of the blackbody curve
                                of the planet []. Given as a np.array with minimum and maximum fraction over which
                                the depth is uniformly selected
        - 'dip_likelihoods': list; assumed fractions of the planets containing each molecule []
        '''

        self.name = name
        self.dip_molecules = [r'$CO_2$', r'$O_3$', r'$H_2O$']
        self.dip_centers = [15*10**-6, 9.7*10**-6, 6.5*10**-6]
        self.dip_avg_widths = [1.0 * 10 ** -6, 0.2 * 10 ** -6, 0.7 * 10 ** -6]
        self.dip_widths = [np.array([(1.0-spread*1.0)*10**-6,(1.0+spread*1.0)*10**-6]),
                           np.array([(0.2-spread*0.2)*10**-6,(0.2+spread*0.2)*10**-6]),
                           np.array([(0.7-spread*0.7)*10**-6,(0.7+spread*0.7)*10**-6])]
        self.dip_avg_depths = [0.45, 0.45, 0.45]
        self.dip_depths = [np.array([0.45-spread*0.45,0.45+spread*0.45]),
                           np.array([0.45-spread*0.45,0.45+spread*0.45]),
                           np.array([0.45-spread*0.45,0.45+spread*0.45])]
        self.dip_likelihoods = [likelihood_CO2, likelihood_O3, likelihood_H2O]

        pass


    def get_width(self, molecule):
        '''
        Function that selects the dip-width of a molecule for a planet based on a random uniform selection between the
        minimum and maximum assumed widths

        :param molecule: str; molecule in question of which to determine the width

        :return width: float; selected width the wavelength dip induced by the molecule for the planet in [micron]
        '''

        index = self.dip_molecules.index(molecule)
        width = np.random.uniform(self.dip_widths[index][0],self.dip_widths[index][1])

        return width


    def get_depth(self, molecule):
        '''
        Function that selects the dip-depth of a molecule for a planet based on a random uniform selection between the
        minimum and maximum assumed depths

        :param molecule: str; molecule in question of which to determine the depth

        :return depth: float; selected depth of the wavelength dip induced by the molecule for the planet as a fraction
                            of the blackbody curve of the planet []
        '''

        index = self.dip_molecules.index(molecule)
        depth = np.random.uniform(self.dip_depths[index][0], self.dip_depths[index][1])

        return depth


    def expected_spectrum(self, molecule, wl_bins, wl_bins_and_distS_array, Tp, Rp, dip_width, dip_depth):
        '''
        This function calculates the blackbody spectrum of a planet with Gaussian dips arising from a molecule

        :param molecule: str; molecule in question that induced the dip in the spectrum
        :param wl_bins: np.ndarray of size L; consists of all the wavelength bins (middle wavelength of each bin) in [m]
        :param wl_bins_and_distS_array: np.ndarray of size (L,2); matrix containing the wavelength bins in [m] and the
                        distance of Earth to the host star in [parsec]
        :param Tp: float; planet temperature in [K]
        :param Rp: float; planet radius in [R_earth]
        :param dip_width: float; standard deviation of the Gaussian dip induced by the molecule in [m]
        :param dip_depth: float; depth defined as the fraction the planet blackbody curve is reduced by at the center
                        of the induced Gaussian dip []

        :return expected spectrum: np.ndarray of size L; photon flux in each of the L wavelength bin in [photons]
        '''

        # planet blackbody spectrum before taking the dip into account
        expected_blackbody = BB_for_fit(wl_bins_and_distS_array, Tp, Rp)

        index = self.dip_molecules.index(molecule)

        # calculate the Gaussian dip induced by the molecule
        gaussian_dip = sp.stats.norm.pdf(wl_bins, loc=self.dip_centers[index], scale=dip_width)
        gaussian_dip_normalized = gaussian_dip / np.max(gaussian_dip) * expected_blackbody[np.abs(wl_bins -
                                                                                    self.dip_centers[index]).argmin()]
        # calculate the planet blackbody spectrum minus the Gaussian dip
        expected_spectrum = expected_blackbody - dip_depth * gaussian_dip_normalized

        return expected_spectrum


# define the atmospheric scenarios
only_CO2 = Atmospheric_scenarios('only CO2', 1.0, 0.0, 0.0, 0.5)
only_O3 = Atmospheric_scenarios('only O3', 0.0, 1.0, 0.0, 0.5)
only_H2O = Atmospheric_scenarios('only H2O', 0.0, 0.0, 1.0, 0.5)
mix_40 = Atmospheric_scenarios('mix 40', 0.4, 0.4, 0.4, 0.5)

Earth_like = Atmospheric_scenarios('Earth-like', 0.8, 0.5, 0.5, 0.001)
Earth_like_loose = Atmospheric_scenarios('Earth-like loose', 0.8, 0.5, 0.5, 0.25)
Earth_like_veryloose = Atmospheric_scenarios('Earth-like very loose', 0.8, 0.5, 0.5, 0.5)

Jupiter_like = Atmospheric_scenarios('Jupiter-like', 0.6, 0.3, 0.6, 0.5)



def get_ratio_safe(list):
    '''
    This function takes lists with "0"s and "1"s as input and calculates the fraction of "1"s

    :param list: list; contains elements equal to "0" or "1"

    :return: float; fraction of "1"s in the list
    '''

    if (len(list) == 0):
        return 0
    else:
        return sum(list)/len(list)


def cdf_J(L, J):
    '''
    This function calculates the cumulative density function for the cost function J'' as described in LIFE II
    equation (31). It uses floats and is therefore only useable for values not too close to 1

    :param L: int; number of wavelength bins []
    :param J: float; value of the cost function []

    :return cdf: float; Cumulative probability density of sum of L bins of the chi-squared distribution []
    '''

    fact = sp.special.factorial

    # calculate the cdf as in equation (31) LIFE II
    cdf = 1 / 2 ** L
    for l in range(0, L):
        cdf += fact(L) / (2 ** L * fact(l) * fact(L - l)) * sp.special.gammainc((L - l) / 2, J / 2)

    return cdf


def cdf_Jmax(L, J, angsep):
    '''
    This function calculates the cumulative distribution function for the cost function J'' when factoring in that the
    maximum pixel is selected as described in Thesis. It calls upon the cdf_J and raises it to the power of the number
    of total pixels in the image. It uses floats and is therefore only useable for values not too close to 1

    :param L: int; number of wavelength bins []
    :param J: float; value of the cost function []
    :param angsep: float; angular separation of the field of search [arcsec]

    :return cdf_Jmax: float; Cumulative probability density of sum of L bins of the chi-squared distribution when
                                selecting the maximum position cluster []
    '''

    # get the number of clusters from the empirical fit
    n_clusters = get_clusters_read(angsep)

    # calculate the cdf suing this number of clusters
    cdf_Jmax = cdf_J(L, J)**n_clusters

    return cdf_Jmax


def pdf_J(L, J):
    '''
    This function calculates the probability density function for the cost function J'' as described in LIFE II
    equation (30). Note that the dirac delta function is not implemented

    :param L: int; number of wavelength bins []
    :param J: float; value of the cost function []

    :return pdf: float; Normalized probability density of sum of L bins of the chi-squared distribution []
    '''

    fact = sp.special.factorial

    # calculate the pdf as in equation (30) LIFE II
    pdf = 0
    for l in range(0, L):
        pdf += fact(L) / (2 ** L * fact(l) * fact(L - l)) * sp.stats.chi2.pdf(J, L - l)

    return pdf


def cdf_J_precision(L, J, precision):
    '''
    This function calculates the cumulative density function for the cost function J'' as described in LIFE II
    equation (31). It uses the mpmath package to achieve better performance for values close to 1

    :param L: int; number of wavelength bins []
    :param J: mpmath_object; value of the cost function []
    :param precision: int; degree of precision used by mpmath []

    :return cdf: mpmath_object; Cumulative probability density of sum of L bins of the chi-squared distribution []
    '''

    # set the precision
    mp.mp.dps = precision

    # calculate the cdf as in equation (31) LIFE II
    cdf = mp.power(mp.mp.mpf('0.5'), L)
    for l in range(0, L):
        cdf += mp.fac(L) / (mp.mp.power(mp.mp.mpf('2'), L) * mp.fac(l) * mp.fac(L - l)) *\
                    mp.gammainc((L - l) / 2, 0, J / 2, regularized=True)

    return cdf


def cdf_Jmax_precision(L, J, precision, angsep):
    '''
    This function calculates the cumulative distribution function for the cost function J'' when factoring in that the
    maximum pixel is selected as described in Thesis. It calls upon the cdf_J_precision and raises it to the power of
    the number of total pixels in the image. It uses the mpmath package to achieve better performance for values close
    to 1

    :param L: int; number of wavelength bins []
    :param J: mpmath_object; value of the cost function []
    :param precision: int; degree of precision used by mpmath []
    :param angsep: float; angular separation of the field of search [arcsec]

    :return cdf_Jmax: mpmath_object; Cumulative probability density of sum of L bins of the chi-squared distribution
                                        when selecting the maximum pixel []
    '''

    # get the number of clusters from the empirical fit
    n_clusters = get_clusters_read(angsep)

    # calculate the cdf suing this number of clusters
    cdf_Jmax = cdf_J_precision(L, J, precision)**n_clusters

    return cdf_Jmax


def alt_sigma_calc(FPR, filepath):
    '''
    This function is used to calculate the number of standard deviations a given input value of a cdf function
    corresponds to. It does this by calling upon a lookup table, which is much faster than performing the direct
    conversion for large sigmas.

    :param FPR: mpmath_object; value of the cdf function that is to be converted []
    :param filepath: str; path to the lookup table

    :return: float; number of sigmas that are inside the confidence interval of the given cdf output value
    '''

    # load the lookup table; the two approaches are for the server/local data structure
    try:
        lookup_table = pd.read_csv('/Users/fdannert/Documents/projects/LIFEsim/working/development/signal_extraction/laura/prob2sigma_conversiontable_minimal.csv',dtype={'prob': str})
    except FileNotFoundError:
        lookup_table = pd.read_csv(filepath+'Auxiliary/'+'prob2sigma_conversiontable_minimal.csv', dtype={'prob': str})

    # create a series with all the probability values as mpmath objects
    mp_series = lookup_table['prob'].map(lambda x: mp.mpf(str(x)))

    # find the probability value closest to the input value
    abs_diff = np.abs(mp_series - FPR)
    min_index = np.argmin(abs_diff)

    # take the sigma corresponding to the closest value
    sigma = lookup_table['sigma'][min_index]

    return sigma


def get_detection_threshold(L, sigma):
    '''
    This function calculates the threshold above which one can be certain to 'sigma' sigmas that a detection is not
    a false positive. See LIFE II section 3.2 for a detailed description

    :param L: int; Number of wavelength bins as given by the wavelength range and the resolution parameter R []
    :param sigma: float; # of sigmas outside of which a false positive must lie []

    :return eta_threshold_sigma: float; threshold is terms of the cost function J []
    '''

    # create the input linspace
    eta = np.linspace(0, 300, int(10 ** 4))

    # calculate the cdf values for each element in the eta-linspace
    cdf = cdf_J(L, eta)

    # find the threshold value eta
    eta_ind_sigma = np.searchsorted(cdf, sp.stats.norm.cdf(sigma))
    eta_threshold_sigma = eta[eta_ind_sigma]

    return eta_threshold_sigma


def get_detection_threshold_max(L, sigma, angsep):
    '''
    This function calculates the threshold above which one can be certain to 'sigma' sigmas that a detection is not
    a false positive. See LIFE II section 3.2 for a detailed description

    :param L: int; Number of wavelength bins as given by the wavelength range and the resolution parameter R []
    :param sigma: float; # of sigmas outside of which a false positive must lie []
    :param angsep: float; angular separation of the field of search [arcsec]

    :return eta_threshold_sigma: float; threshold is terms of the cost function J []
    '''

    # create the input linspace
    eta = np.linspace(0, 300, int(10 ** 3))

    # calculate the cdf values for each element in the eta-linspace
    cdf = np.empty_like(eta)
    for i in range(eta.size):
        cdf[i] = cdf_Jmax(L, eta[i], angsep)

    # find the threshold value eta
    eta_ind_sigma = np.searchsorted(cdf, sp.stats.norm.cdf(sigma))
    eta_threshold_sigma = eta[eta_ind_sigma]

    return eta_threshold_sigma


def BB_for_fit(wl_and_distS, Tp, Rp):
    '''
    This function calculates the flux received at Earth from an object with radius Rp radiating at temperature T,
    at distance dist_s and at wavelengths wl

    :param wl_and_distS: np.ndarray of shape (L,2); wavelength bins and distance to host star (constant) in [m,pc]
    :param Tp: float; planet temperature in [K]
    :param Rp: float; planet radius in [R_earth]

    :return fgamma: np.ndarray of size L; contains the total blackbody fluxes in each of the L wavelength bins
                    in [photons]
    '''

    # unpack the input array (required in this format for the scipy.optimizer function)
    wl = wl_and_distS[0]
    dist_s = wl_and_distS[1]

    fact1 = 2 * c / (wl ** 4)
    fact2 = (h * c) / (k * Tp * wl)

    # calculate the standard planck function
    fgamma = np.array(fact1 / (np.exp(fact2) - 1.0)) * 10 ** -6 * np.pi * (
            (Rp * radius_earth) / (dist_s * m_per_pc)) ** 2

    return fgamma


def t_distribution(x,dof):
    '''
    Calculates the Student's t-distribution for a given input

    :param x: np.dnarray of variable size; input quantity (here usually J'') []
    :param dof: int; degrees of freedom of th distribution []

    :return t_distr: np.dnarray of variable size; t-distribution values for at each position x
    '''
    t_distr = sp.stats.t.pdf(x=x,df=dof)

    return t_distr


def get_Dsqr_mat(L):
    '''
    This function calcualtes the D^2 matrix as described in LIFE II Appendix B used for the calculation of the
    estimated planet flux. See this paper for the explanations of the exact calculations

    :param L : int; number of wavelegth bins []

    :return Dsqr_mat: matrix of dimensions (L,radial_ang_px,n_steps); required for calculation of the estimated
            planet flux []
    '''

    dif_diag = -2 * np.diag(np.ones(L))
    dif_pre = 1 * np.diag(np.ones(L - 1), 1)
    dif_post = 1 * np.diag(np.ones(L - 1), -1)

    D_mat = dif_diag + dif_pre + dif_post
    D_mat[0, :2] = np.array([-1, 1])
    D_mat[-1, -2:] = np.array([1, -1])

    Dsqr_mat = np.matmul(D_mat, D_mat)

    return Dsqr_mat


def cartesian2polar_for_map(outcoords, inputshape):
    '''
    Transforms cartesian coordinates into polar coordinates; auxiliary function for pol_to_cart_map (adapted from
    the old LIFEsim version)

    :param outcoords: format of the output coordinates
    :param inputshape: shape of the input coordinates

    :return (r,phi_index): tuple; indices r and phi of the output polar map
    '''

    y, x = outcoords
    x = x - (inputshape[0] - 0.5)
    y = y - (inputshape[0] - 0.5)

    r = np.sqrt(x ** 2 + y ** 2)
    phi = np.arctan2(-y, -x)
    phi_index = (phi + np.pi) * inputshape[1] / (2 * np.pi)

    return (r, phi_index)


def pol_to_cart_map(image, image_size):
    '''
    Produces a cartesian map from a given input image (adapted from the old LIFEsim version)

    :param image: input image
    :param image_size: int; size of the input image

    :return cartesian map: transformed image in cartesian coordinates
    '''

    # create new column at end (360 deg) with same values as first column (0 deg) to get complete image
    image_new = np.empty((image.shape[0], image.shape[1] + 1))
    image_new[:, :-1] = image
    image_new[:, -1] = image[:, 0]

    cartesian_map = sp.ndimage.geometric_transform(image_new, cartesian2polar_for_map, order=1,
                                                   output_shape=(image_size, image_size), mode="constant",
                                                   extra_keywords={'inputshape': (image.shape)})

    return cartesian_map


def plot_planet_SED_and_SNR(wl_bins, Fp, Fp_est, sigma, wl_min, wl_max, Fp_BB=None, snr_photon_stat=None,
                            filename=None):
    '''
    Plots the true and extracted blackbody curves (from the old lifesim version: modules/plotting/plotter.py)

    :param wl_bins: np.ndarray of size L; consists of all the wavelength bins (middle wavelength of each bin) in [m]
    :param Fp: np.ndarray of size L; true blackbody curve of the planet [photons]
    :param Fp_est: np.ndarray of size L; extracted fluxes for each wl bin [photons]
    :param sigma: np.ndarray of size L; uncertainties of the extracted fluxes for each wl bin [photons]
    :param wl_min: float; minimum wavelength captured by the instrument [m]
    :param wl_max: float; maximum wavelength captured by the instrument [m]
    :param Fp_BB: np.ndarray of size L; fitted blackbody curve of the planet [photons]
    :param snr_photon_stat: np.ndarray of size L; snr per wavelength bin as calculated using photon statistics
                                (no longer calculated in the new lifesim version, always set to 'None')
    :param filename: str; name of the file if the plot should be saved. If 'None', no plot is saved
    '''

    fig, ax = plt.subplots()

    ax.plot(wl_bins * 1e6, Fp, color="mediumblue", linestyle="-", label="True spectrum")


    if snr_photon_stat is not None:
        ax.fill_between(wl_bins * 1e6, Fp * (1 - 1 / snr_photon_stat), Fp * (1 + 1 / snr_photon_stat),
                        color="mediumblue", alpha=0.1, label=r"Photon noise $\sigma$")

    if Fp_BB is not None:
        ax.plot(wl_bins * 1e6, Fp_BB, color="green", linestyle="-", label="Fit BB Spectrum")

    ax.plot(wl_bins * 1e6, Fp_est, color="red", marker=".", linestyle="")

    with np.errstate(divide='ignore'):
        ax.errorbar(wl_bins * 1e6, Fp_est, yerr=sigma, capsize=3, marker=".",
                    color="red", ecolor="red", linestyle="none",
                    elinewidth=1, alpha=0.4,
                    label="Estimated spectrum")

    ax.set_xlabel(r"$\lambda$ [$\mu$m]", fontsize=12)
    ax.set_ylabel(r"Planet flux $F_\lambda$ [ph $\mathrm{s}^{-1}$m$^{-2}\mu \mathrm{m}^{-1}$]", fontsize=12)
    ax.set_xlim(wl_min * 10 ** 6, wl_max * 10 ** 6)
    ax.set_ylim(0, 1.6 * np.max(Fp))
    ax.grid()
    ax.set_title('Blackbody curve fit')
    ax.legend(fontsize=10)

    if filename is not None:
        plt.savefig("C:\\Users\\Philipp Binkert\\OneDrive\\ETH\\Master_Thesis\\06_plots\\" + filename + ".pdf",
                        bbox_inches='tight')

    plt.show()

    return


def plot_multi_map(maps, map_type, hfov_mas, colormap="inferno", vmin=None, vmax=None,
                   show_detections=False, filename_post=None):
    '''
    Plots the cost function J'' for each of the pixels in the image (from the old
    lifesim version: modules/plotting/plotanalysis.py)

    :param maps: np.ndarray; input map data
    :param map_type: str or list of str; labels for the maps
    :param hfov_mas: float; half field of view used to create the map
    :param colormap: str; style of the colormap
    :param vmin: float; minimum value of the colorbar to be shown
    :param vmax: float; maximum value of the colorbar to be shown
    :param show_detections: boolean; if True, mark the position of the detected planet in the plot
    :param filename_post: str; name of the file if the plot should be saved. If 'None', no plot is saved
    '''

    if len(np.shape(maps)) < 3:
        maps = [maps]

    n = len(maps)

    if n > 1:
        cbar_location = "bottom"
    else:
        cbar_location = "right"

    fig = plt.figure(figsize=(n * 4.8 + 1, 4.8), dpi=300)
    grid = AxesGrid(fig, 111,
                    nrows_ncols=(1, n),
                    axes_pad=0.2,
                    cbar_mode='each',
                    cbar_location=cbar_location,
                    cbar_pad=0.05,
                    cbar_size="5%"
                    )

    sf_mas = hfov_mas

    i = 0
    for ax in grid:
        im = ax.imshow(maps[i], cmap=colormap, origin="lower",
                       extent=[sf_mas, -sf_mas, -sf_mas, sf_mas], vmin=vmin, vmax=vmax)

        ax.set_title('Heatmap of the Cost function J\u2032\u2032')

        ax.set_xticks([-np.round(sf_mas / 2,0), 0, np.round(sf_mas / 2,0)])
        ax.set_yticks([-np.round(sf_mas / 2, 0), 0, np.round(sf_mas / 2, 0)])

        ax.tick_params(pad=1)
        plt.setp(ax.get_yticklabels(), rotation=90, va='center')

        ax.set_xlabel(r"$\Delta$RA [mas]")
        ax.set_ylabel(r"$\Delta$Dec [mas]")

        if show_detections:
            theta_max = np.unravel_index(np.argmax(maps[i], axis=None), maps[i].shape)
            t = tuple(ti / len(maps[i]) for ti in theta_max)

            ax.annotate(str(i + 1), xy=(t[1], t[0] - 0.02), xytext=(t[1], t[0] - 0.2), va="center", ha='center',
                        color="white",
                        arrowprops=dict(color='white', width=1, headwidth=10, shrink=0.1), xycoords=ax.transAxes,
                        fontsize=18)

        sfmt = mtick.ScalarFormatter(useMathText=True)
        sfmt.set_powerlimits((-3, 3))
        sfmt.set_scientific(False)


        cbar = grid.cbar_axes[i].colorbar(im, format=sfmt)

        if n > 1:
            if type(map_type) == list:
                label = map_type[i]
            else:
                label = map_type

            ax.xaxis.set_label_position('top')
            ax.xaxis.tick_top()
            cbar.ax.set_xlabel(label, fontsize=12)
        else:
            cbar.ax.set_ylabel(map_type, fontsize=12)

        if maps[i].max() == 1:
            cbar.set_ticks([maps[i].min(), 0, maps[i].max()])

        i += 1

    if filename_post is not None:
        plt.savefig(filename_post + '.pdf', bbox_inches='tight')

    plt.show()

    return


def get_rates(retrieved_values, true_values, threshold, bound_below):
    '''
    This function calculates the confusion matrix for a given set of retrieved and actual values

    :param retrieved_values: np.ndarray of size corresponding to the number of planets; confidence levels in the measure
                                of the threshold of the retrieved values (t-score, log of bayesian factor) []
    :param true_values: np.ndarray of size corresponding to the number of planets; confidence levels in the measure
                                of the threshold of the true values (t-score, log of bayesian factor) []
    :param threshold: float; threshold below/above which a retreived value is considered a detection []
    :param bound_below: boolean; if True, retrieved values below the threshold are considered detections
                                (alpha, t-score), otherwise above (jeffrey, log of bayes factor)

    :return TP: float; True Positive, number of instances correctly predicted as positive by the model used []
    :return FP: float; False Positive, number of instances incorrectly predicted as positive by the model used []
    :return FN: float; False Negative, number of instances incorrectly predicted as negative by the model used []
    :return TN: float; True Negative, number of instances correctly predicted as negative by the model used []
    :return TPR: float; True Positive Rate, proportion of actual positive instances that are correctly identified as
                            positive by the model used []
    :return FPR: float; False Positive Rate, proportion of actual negative instances that are incorrectly classified as
                            positive by the model used []
    '''

    TP = np.zeros((true_values.size))
    FP = np.zeros((true_values.size))
    FN = np.zeros((true_values.size))
    TN = np.zeros((true_values.size))

    # determine which category each planet falls into
    for j in range(true_values.size):
        if (bound_below == True):
            if (retrieved_values[j] <= threshold and true_values[j] == 1):
                TP[j] = 1
            elif (retrieved_values[j] <= threshold and true_values[j] == 0):
                FP[j] = 1
            elif (retrieved_values[j] > threshold and true_values[j] == 1):
                FN[j] = 1
            elif (retrieved_values[j] > threshold and true_values[j] == 0):
                TN[j] = 1

        else:
            if (retrieved_values[j] >= threshold and true_values[j] == 1):
                TP[j] = 1
            elif (retrieved_values[j] >= threshold and true_values[j] == 0):
                FP[j] = 1
            elif (retrieved_values[j] < threshold and true_values[j] == 1):
                FN[j] = 1
            elif (retrieved_values[j] < threshold and true_values[j] == 0):
                TN[j] = 1

    # calculate the TPR and the FPR from the existing arrays
    if (sum(TP) == 0 and sum(FN) == 0):
        TPR = 0
    else:
        TPR = sum(TP)/(sum(TP)+sum(FN))

    if(sum(FP) == 0 and sum(TN) == 0):
        FPR = 0
    else:
        FPR = sum(FP)/(sum(FP)+sum(TN))

    return TP, FP, FN, TN, TPR, FPR


def get_clusters(j_map):
    '''
    This function calculates the number of clusters obtained defined as the number of minima and maxima. No longer
    actually used in the final version, as the parameter n was fitted directly

    :param j_map: np.ndarray of shape (radial_ang_pixels,n_steps); cost function values in each element []

    :return clusters: int; number of clusters calculated []
    '''

    local_maxima = 0
    local_minima = 0

    # loop through each pixel and determine if it is a local maximum or minimum by comparing to its four neighbors
    for l in range(1, j_map.shape[0] - 1):
        for k in range(1, j_map.shape[1] - 1):
            neighborhood = j_map[l - 1:l + 2, k - 1:k + 2]

            if (j_map[l, k] == np.max(neighborhood) and j_map[l, k] > 0 and np.all(neighborhood != 0)):
                local_maxima += 1
            if (j_map[l, k] == np.min(neighborhood) and j_map[l, k] > 0 and np.all(neighborhood != 0)):
                local_minima += 1

    clusters = local_maxima + local_minima

    return clusters


def get_clusters_read(angsep):
    '''
    This function returns the value n according to the empirical fit performed in evaluate_Jmax_runs.py. See the
    written thesis for a detailed explanation

    :param angsep: float; angular separation of the field of search [arcsec]

    :return clusters: float; number of clusters (abstractly) n as determined by the fit
    '''

    # fitted parameters
    fit_parameters_below = [1054093, -7422, 38]
    fit_parameters_above = [118019, -1735]

    threshold = 0.05

    # determine the value of n based on the angular separation
    if (angsep < threshold):
        clusters = fit_parameters_below[0] * angsep**2 + fit_parameters_below[1] * angsep + fit_parameters_below[2]
    else:
        clusters = fit_parameters_above[0] * angsep + fit_parameters_above[1]

    return clusters