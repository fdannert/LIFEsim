import numpy as np
from tqdm import tqdm
from spectres import spectres

from lifesim.core.modules import InstrumentModule
from lifesim.util.habitable import single_habitable_zone
from lifesim.util.radiation import black_body


class Instrument(InstrumentModule):
    """
    The Instrument class represents the central module for simulating the LIFE array. It connects
    to 'plugin' modules which calculate signal and noise terms and distributes tasks and data
    between them. The instrument class features two socket types:
        a)  For calculation of the instrument transmission map a single socket of f_type
            'transmission'
        b)  For simulation of the photon noise sources a number (set in the options class) of
            sockets of f_type 'photon_noise'

    Attributes
    ----------
    options : lifesim.Options
        The options class containing all setting for the array and computations
    bl : float
        Length of the shorter, nulling baseline in [m]
    telescope_area : float
        Area of all array apertures combined in [m^2]
    eff_tot : float
        Total efficiency of the telescope as ratio of generated counts over incoming photons
        (dimensionless)
    wl_bins : np.ndarray
        Central values of the spectral bins in the wavelength regime in [m]
    wl_bin_widths : np.ndarray
        Widths of the spectral wavelength bins in [m]
    wl_bin_edges : np.ndarray
        Edges of the spectral wavelength bins in [m]. For N bins, this array will contain N+1 edges
    hfov : np.ndarray
        Contains the half field of view of the observatory in [rad] for each of the spectral bins
    hfov_mas : np.ndarray
        Contains the half field of view of the observatory in [milliarcseconds] for each of the
        spectral bins
    rad_pix : np.ndarray
        Contains the size of each pixel projected to the sky in [rad]
    mas_pix : np.ndarray
        Contains the size of each pixel projected to the sky in [milliarcseconds]
    apertures : np.ndarray
        Positions of the collector spacecraft relative to the beam combiner in [m]
    x_map : np.ndarray
        A map used for speeding up calculations. Contains the indices of the pixels on the detector
        in x-direction (dimensionless)
    y_map : np.ndarray
        A map used for speeding up calculations. Contains the indices of the pixels on the detector
        in x-direction (dimensionless)
    r_square_map : np.ndarray
        A map used for speeding up calculations. Contains the square of the distance of a pixel
        from the center of the detector in [pix]
    r_map : np.ndarray
        A map used for speeding up calculations. Contains the distance of a pixel
        from the center of the detector in [pix]
    """

    def __init__(self,
                 name: str):
        """
        Parameters
        ----------
        name : str
            Name of the instrument module
        options : lifesim.Options
            The options class containing all setting for the array and computations
        """

        super().__init__(name=name)

    def apply_options(self):
        """
        Applies the options given to the instrument module and recalculates all necessary values

        Parameters
        ----------
        options : lifesim.Options
            The options class containing all setting for the array and computations
        """

        # Get array parameters from options for faster calculation
        self.data.inst['bl'] = self.data.options.array['baseline']

        self.data.inst['telescope_area'] = 4. * np.pi \
                                           * (self.data.options.array['diameter'] / 2.) ** 2
        self.data.inst['eff_tot'] = self.data.options.array['quantum_eff'] \
                                    * self.data.options.array['throughput']

        self.data.inst['wl_bins'], \
        self.data.inst['wl_bin_widths'], \
        self.data.inst['wl_bin_edges'] = self.get_wl_bins_const_spec_res()

        # fov = wl / D -> hfov=wl/(2*D)
        # TODO remove the double usage of mas and rad, stick to only one
        self.data.inst['hfov'] = self.data.inst['wl_bins'] \
                                 / (2. * self.data.options.array['diameter'])

        self.data.inst['hfov_mas'] = self.data.inst['hfov'] * (3600000. * 180.) / np.pi
        self.data.inst['rad_pix'] = (2 * self.data.inst['hfov']) \
                                    / self.data.options.other['image_size']  # Radians per pixel
        self.data.inst['mas_pix'] = (2 * self.data.inst['hfov_mas']) \
                                    / self.data.options.other['image_size']  # mas per pixel

        # apertures defines the telescope positions (and *relative* radius)
        self.data.inst['apertures'] = np.array([
            [-self.data.inst['bl'] / 2,
             -self.data.options.array['ratio'] * self.data.inst['bl'] / 2., 1.],
            [self.data.inst['bl'] / 2,
             -self.data.options.array['ratio'] * self.data.inst['bl'] / 2., 1.],
            [self.data.inst['bl'] / 2,
             self.data.options.array['ratio'] * self.data.inst['bl'] / 2., 1.],
            [-self.data.inst['bl'] / 2,
             self.data.options.array['ratio'] * self.data.inst['bl'] / 2., 1.]
        ])

        # coordinate maps for faster calculations
        x_map = np.tile(np.array(range(0, self.data.options.other['image_size'])),
                        (self.data.options.other['image_size'], 1))
        y_map = x_map.T
        r_square_map = ((x_map - (self.data.options.other['image_size'] - 1) / 2) ** 2
                        + (y_map - (self.data.options.other['image_size'] - 1) / 2) ** 2)
        self.data.inst['radius_map'] = np.sqrt(r_square_map)

    def get_wl_bins_const_spec_res(self):
        """
        Create the wavelength bins for the given spectral resolution and wavelength limits
        """
        wl_edge = self.data.options.array['wl_min']
        wl_bins = []
        wl_bin_widths = []
        wl_bin_edges = [wl_edge]

        while wl_edge < self.data.options.array['wl_max']:

            # set the wavelength bin width according to the spectral resolution
            wl_bin_width = wl_edge / self.data.options.array['spec_res'] / \
                           (1 - 1 / self.data.options.array['spec_res'] / 2)

            # make the last bin shorter when it hits the wavelength limit
            if wl_edge + wl_bin_width > self.data.options.array['wl_max']:
                wl_bin_width = self.data.options.array['wl_max'] - wl_edge

            # calculate the center and edges of the bins
            wl_center = wl_edge + wl_bin_width / 2
            wl_edge += wl_bin_width

            wl_bins.append(wl_center)
            wl_bin_widths.append(wl_bin_width)
            wl_bin_edges.append(wl_edge)

        # convert everything to [m]
        wl_bins = np.array(wl_bins) * 1e-6  # in m
        wl_bin_widths = np.array(wl_bin_widths) * 1e-6  # in m
        wl_bin_edges = np.array(wl_bin_edges) * 1e-6  # in m

        return wl_bins, wl_bin_widths, wl_bin_edges

    # TODO does not take the inclination into account!
    def adjust_bl_to_hz(self,
                        hz_center: float,
                        distance_s: float):
        """
        Adjusts the baseline of the array to be optimal for observations in the habitable zone of
        the target star for the selected optimal wavelength

        Parameters
        ----------
        hz_center : float
            Separation of the center of the habitable zone in [AU]
        distance_s : float
            Distance between the observed star and the LIFE array in [pc]
        """

        # convert the habitable zone to radians
        hz_center_rad = hz_center / distance_s / (3600 * 180) * np.pi  # in rad

        # put first transmission peak of optimal wl on center of HZ
        self.data.inst['bl'] = 0.589645 / hz_center_rad * self.data.options.other['wl_optimal']*10**(-6)

        # make sure that the baseline does not exeed the set baseline limits
        self.data.inst['bl'] = np.maximum(self.data.inst['bl'], self.data.options.array['bl_min'])
        self.data.inst['bl'] = np.minimum(self.data.inst['bl'], self.data.options.array['bl_max'])

        # update the position of the apertures
        self.data.inst['apertures'] = np.array([
            [-self.data.inst['bl'] / 2,
             -self.data.options.array['ratio'] * self.data.inst['bl'] / 2., 1.],
            [self.data.inst['bl'] / 2,
             -self.data.options.array['ratio'] * self.data.inst['bl'] / 2., 1.],
            [self.data.inst['bl'] / 2,
             self.data.options.array['ratio'] * self.data.inst['bl'] / 2., 1.],
            [-self.data.inst['bl'] / 2,
             self.data.options.array['ratio'] * self.data.inst['bl'] / 2., 1.]
        ])

    # TODO Re-add functionality for calculating the SNR without certain noise term
    def get_snr(self,
                safe_mode:bool = False):
        """
        Calculates the signal-to-noise ration for all planets within the catalog if the are
        observed by the LIFE array for the given observing time

        Parameters
        ----------
        c : lifesim.Catalog
            The catalog for which the SNRs will be calculated and added to
        time : float
            The observation time per planet in [s]
        """

        self.apply_options()
        integration_time = 60 * 60

        self.data.catalog['snr_1h'] = np.zeros_like(self.data.catalog.nstar, dtype=float)
        if safe_mode:
            self.data.catalog['noise_astro'] = None
            self.data.catalog['planet_flux_use'] = None

        # create mask returning only unique stars
        _, temp = np.unique(self.data.catalog.nstar, return_index=True)
        star_mask = np.zeros_like(self.data.catalog.nstar, dtype=bool)
        star_mask[temp] = True

        # iterate over all stars
        for i, n in enumerate(tqdm(np.where(star_mask)[0])):
            # if i == 10:
            #     break
            nstar = self.data.catalog.nstar.iloc[n]

            # adjust baseline of array and give new baseline to transmission generator plugin
            self.adjust_bl_to_hz(hz_center=float(self.data.catalog.hz_center.iloc[n]),
                                 distance_s=float(self.data.catalog.distance_s.iloc[n]))

            # get transmission map
            _, _, self.data.inst['t_map'], _, _ = self.run_socket(s_name='transmission',
                                                                  method='transmission_map',
                                                                  map_selection='tm3')

            # calculate the noise from the background sources
            noise_bg_list = self.run_socket(s_name='photon_noise',
                                            method='noise',
                                            index=n)

            if type(noise_bg_list) == list:
                noise_bg = np.zeros_like(noise_bg_list[0])
                for _, noise in enumerate(noise_bg_list):
                    noise_bg += noise
            else:
                noise_bg = noise_bg_list

            noise_bg = noise_bg * integration_time * self.data.inst['eff_tot'] * 2

            # go through all planets for the chosen star
            for _, n_p in enumerate(np.argwhere(
                    self.data.catalog.nstar.to_numpy() == nstar)[:, 0]):

                # calculate the photon flux originating from the planet
                flux_planet_thermal = black_body(mode='planet',
                                                 bins=self.data.inst['wl_bins'],
                                                 width=self.data.inst['wl_bin_widths'],
                                                 temp=self.data.catalog['temp_p'].iloc[n_p],
                                                 radius=self.data.catalog['radius_p'].iloc[n_p],
                                                 distance=self.data.catalog['distance_s'].iloc[n_p]
                                                 )

                transm_eff, transm_noise = self.run_socket(s_name='transmission',
                                                           method='transmission_efficiency',
                                                           index=n_p)

                # calculate the signal and photon noise flux received from the planet
                flux_planet = flux_planet_thermal \
                              * transm_eff \
                              * integration_time \
                              * self.data.inst['eff_tot'] \
                              * self.data.inst['telescope_area']
                noise_planet = flux_planet_thermal \
                               * transm_noise \
                               * integration_time \
                               * self.data.inst['eff_tot'] \
                               * self.data.inst['telescope_area'] \
                               * 2

                # Add up the noise and caluclate the SNR
                noise = noise_bg + noise_planet
                self.data.catalog.snr_1h.iat[n_p] = np.sqrt((flux_planet ** 2 / noise).sum())

                if safe_mode:
                    self.data.catalog.noise_astro.iat[n_p] = [noise_bg]
                    self.data.catalog.planet_flux_use.iat[n_p] = [flux_planet_thermal \
                                                                  * integration_time \
                                                                  * self.data.inst['eff_tot'] \
                                                                  * self.data.inst['telescope_area']]

    def get_snr_t(self,
                safe_mode:bool = False):
        """
        Calculates the signal-to-noise ration for all planets within the catalog if the are
        observed by the LIFE array for the given observing time

        Parameters
        ----------
        c : lifesim.Catalog
            The catalog for which the SNRs will be calculated and added to
        time : float
            The observation time per planet in [s]
        """

        self.apply_options()
        
        self.data.inst["rotation_period"] = 60 * 60 * 1
        self.data.inst["rotations"] = 1
        self.data.inst["integration_time"] = self.data.inst["rotations"] * self.data.inst["rotation_period"]
        integration_time = self.data.inst["integration_time"]
        # timestep option; goal is to calculate angular seperation and the transmission coefficient there
        # then sum noise and singal over all steps
        # angular seperation calculation starts from true anomaly which is random


        self.data.catalog['snr_1h'] = np.zeros_like(self.data.catalog.nstar, dtype=float)
        if safe_mode:
            self.data.catalog['noise_astro'] = None
            self.data.catalog['planet_flux_use'] = None

        # create mask returning only unique stars
        _, temp = np.unique(self.data.catalog.nstar, return_index=True)
        star_mask = np.zeros_like(self.data.catalog.nstar, dtype=bool)
        star_mask[temp] = True

        # iterate over all stars
        for i, n in enumerate(tqdm(np.where(star_mask)[0])):
            # if i == 10:
            #     break
            nstar = self.data.catalog.nstar.iloc[n]

            # adjust baseline of array and give new baseline to transmission generator plugin
            self.adjust_bl_to_hz(hz_center=float(self.data.catalog.hz_center.iloc[n]),
                                 distance_s=float(self.data.catalog.distance_s.iloc[n]))

            # get transmission map
            _, _, self.data.inst['t_map'], _, _ = self.run_socket(s_name='transmission',
                                                                  method='transmission_map',
                                                                  map_selection='tm3')

            # calculate the noise from the background sources
            noise_bg_list = self.run_socket(s_name='photon_noise',
                                            method='noise',
                                            index=n)

            if type(noise_bg_list) == list:
                noise_bg = np.zeros_like(noise_bg_list[0])
                for _, noise in enumerate(noise_bg_list):
                    noise_bg += noise
            else:
                noise_bg = noise_bg_list

            noise_bg = noise_bg * integration_time * self.data.inst['eff_tot'] * 2

            # go through all planets for the chosen star
            for _, n_p in enumerate(np.argwhere(
                    self.data.catalog.nstar.to_numpy() == nstar)[:, 0]):

                # calculate the photon flux originating from the planet
                flux_planet_thermal = black_body(mode='planet',
                                                 bins=self.data.inst['wl_bins'],
                                                 width=self.data.inst['wl_bin_widths'],
                                                 temp=self.data.catalog['temp_p'].iloc[n_p],
                                                 radius=self.data.catalog['radius_p'].iloc[n_p],
                                                 distance=self.data.catalog['distance_s'].iloc[n_p]
                                                 )

                transm_eff, transm_noise = self.run_socket(s_name='transmission',
                                                           method='transmission_efficiency_t',
                                                           index=n_p)

                # calculate the signal and photon noise flux received from the planet
                flux_planet = flux_planet_thermal \
                              * transm_eff \
                              * integration_time \
                              * self.data.inst['eff_tot'] \
                              * self.data.inst['telescope_area']
                noise_planet = flux_planet_thermal \
                               * transm_noise \
                               * integration_time \
                               * self.data.inst['eff_tot'] \
                               * self.data.inst['telescope_area'] \
                               * 2

                # Add up the noise and caluclate the SNR
                noise = noise_bg + noise_planet
                self.data.catalog.snr_1h.iat[n_p] = np.sqrt((flux_planet ** 2 / noise).sum())

                if safe_mode:
                    self.data.catalog.noise_astro.iat[n_p] = [noise_bg]
                    self.data.catalog.planet_flux_use.iat[n_p] = [flux_planet_thermal \
                                                                  * integration_time \
                                                                  * self.data.inst['eff_tot'] \
                                                                  * self.data.inst['telescope_area']]

    # TODO: fix units in documentation
    def get_spectrum(self,
                     temp_s: float,  # in K
                     radius_s: float,  # in R_sun
                     distance_s: float,  # in pc
                     lat_s: float,  # in radians
                     z: float,  # in zodis
                     angsep: float,  # in arcsec
                     flux_planet_spectrum: list,  # in ph m-3 s-1 over m
                     integration_time: float):
        """
        Simulates the spectrum of an exoplanet as seen by LIFE array, i.e. with the respective
        noise terms

        Parameters
        ----------
        pathtofile : str
            Path to the .txt file holding the input spectrum. The input spectrum need to be in the
            units of a photon flux density per micron, i.e. [s^-1 m^-2 microns^-1]
        temp_s : float
            Temperature of the observed star in [K]
        radius_s : float
            Radius of the observed star in [sun radii]
        distance_s : float
            Distance between the observed star and the LIFE array in [pc]
        lat_s : float
            Ecliptic latitude of the observed star in [rad]
        z : float
            Zodi level in the observed system in [zodis]
        angsep : float
            Angular separation between the observed star and the observed exoplanet in [arcsec]
        radius_p : float
            Radius of the observed exoplanet in [earth radii]
        radius_spec : float
            Radius of the observed exoplanet that was assumed in the creation of the input spectrum
            in [earth radii]
        distance_spec : float
            Distance between the observed star and the observer that was assumed in the creation of
            the input spectrum in [pc]
        integration_time : float
            Time that the LIFE array spends for observing the observed planet in [s]

        Returns
        -------
        Tuple[wl_bins, snr_spec]
            Some text
        flux_planet
            Some text
        noise
            Some text
        """

        self.apply_options()

        self.data.single['temp_s'] = temp_s
        self.data.single['radius_s'] = radius_s
        self.data.single['distance_s'] = distance_s
        self.data.single['lat'] = lat_s
        self.data.single['z'] = z
        self.data.single['angsep'] = angsep

        # calculate the habitable zone of the specified star
        s_in, s_out, l_sun, \
            hz_in, hz_out, \
            hz_center = single_habitable_zone(model=self.data.options.models['habitable'],
                                              temp_s=temp_s,
                                              radius_s=radius_s)

        self.data.single['l_sun'] = l_sun

        # adjust the baseline of the array to the habitable zone
        self.adjust_bl_to_hz(hz_center=hz_center,
                             distance_s=distance_s)

        _, _, self.data.inst['t_map'], _, _ = self.run_socket(s_name='transmission',
                                                              method='transmission_map',
                                                              map_selection='tm3')

        # update and run the photon noise plugins
        flux_planet_spectrum = spectres(new_wavs=self.data.inst['wl_bin_edges'],
                                        spec_wavs=flux_planet_spectrum[0].value,
                                        spec_fluxes=flux_planet_spectrum[1].value,
                                        edge_mode=True)

        transm_eff, transm_noise = self.run_socket(s_name='transmission',
                                                   method='transmission_efficiency',
                                                   index=None)

        # calculate the signal and photon noise flux received from the planet
        # TODO: to be consistent with get_snr, make it such that bin_width is multiplied elsewhere
        flux_planet = flux_planet_spectrum \
                      * transm_eff \
                      * integration_time \
                      * self.data.inst['eff_tot'] \
                      * self.data.inst['telescope_area'] \
                      * self.data.inst['wl_bin_widths']
        noise_planet = flux_planet_spectrum \
                       * transm_noise \
                       * integration_time \
                       * self.data.inst['eff_tot'] \
                       * self.data.inst['telescope_area'] \
                       * self.data.inst['wl_bin_widths']

        # calculate the noise from the background sources
        noise_bg_list = self.run_socket(s_name='photon_noise',
                                        method='noise',
                                        index=None)

        if type(noise_bg_list) == list:
            noise_bg = np.zeros_like(noise_bg_list[0])
            for _, noise in enumerate(noise_bg_list):
                noise_bg += noise
        else:
            noise_bg = noise_bg_list

        noise_bg = noise_bg * integration_time * self.data.inst['eff_tot']

        # Add up the noise and caluclate the SNR
        noise = (noise_bg + noise_planet) * 2
        snr_spec = np.sqrt((flux_planet ** 2 / noise))

        return ([self.data.inst['wl_bins'], snr_spec],
                flux_planet,
                noise)



