import numpy as np
from tqdm import tqdm
from spectres import spectres

from lifesim.core.modules import InstrumentModule
from lifesim.util.habitable import single_habitable_zone
from lifesim.util.radiation import black_body


class Instrument(InstrumentModule):
    """
    The Instrument class represents the central module for simulating the LIFE array. It connects
    to other modules which calculate signal and noise terms and distributes tasks and data
    between them. The instrument class features two socket types:
        a)  For calculation of the instrument transmission map a single socket of f_type
            'transmission'.
        b)  For simulation of the photon noise sources a number (set in the options class) of
            sockets of f_type 'photon_noise'.

    Notes
    -----
    Note, that all attributes are saved in the data class.

    Attributes
    ----------
    data.options : lifesim.Options
        The options class containing all setting for the array and computations.
    data.inst['bl'] : float
        Length of the shorter, nulling baseline in [m].
    data.inst['telescope_area'] : float
        Area of all array apertures combined in [m^2].
    data.inst['eff_tot'] : float
        Total efficiency of the telescope as ratio of generated counts over incoming photons
        (dimensionless).
    data.inst['wl_bins'] : np.ndarray
        Central values of the spectral bins in the wavelength regime in [m].
    data.inst['wl_bin_widths'] : np.ndarray
        Widths of the spectral wavelength bins in [m].
    data.inst['wl_bin_edges'] : np.ndarray
        Edges of the spectral wavelength bins in [m]. For N bins, this array will contain N+1
        edges.
    data.inst['hfov'] : np.ndarray
        Contains the half field of view of the observatory in [rad] for each of the spectral bins.
    data.inst['hfov_mas'] : np.ndarray
        Contains the half field of view of the observatory in [milliarcseconds] for each of the
        spectral bins.
    data.inst['rad_pix'] : np.ndarray
        Contains the size of each pixel projected to the sky in [rad].
    data.inst['mas_pix'] : np.ndarray
        Contains the size of each pixel projected to the sky in [milliarcseconds].
    data.inst['apertures'] : np.ndarray
        Positions of the collector spacecraft relative to the beam combiner in [m].
    data.inst['radius_map'] : np.ndarray
        A map used for speeding up calculations. Contains the distance of a pixel
        from the center of the detector in [pix].
    """

    def __init__(self,
                 name: str):
        """
        Parameters
        ----------
        name : str
            Name of the instrument module.
        """

        super().__init__(name=name)

    def apply_options(self):
        """
        Applies the options given to the instrument module and recalculates all necessary values.
        """

        # Get array parameters from options for faster calculation
        self.data.inst['bl'] = self.data.options.array['baseline']

        self.data.inst['telescope_area'] = 4. * np.pi \
                                           * (self.data.options.array['diameter'] / 2.) ** 2
        self.data.inst['eff_tot'] = self.data.options.array['quantum_eff'] \
                                    * self.data.options.array['throughput']

        # Calculate the spectral channels with a constant spectral resolution
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
        Create the wavelength bins for the given spectral resolution and wavelength limits.
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
        the target star for the selected optimal wavelength.

        Parameters
        ----------
        hz_center : float
            Separation of the center of the habitable zone in [AU].
        distance_s : float
            Distance between the observed star and the LIFE array in [pc].
        """

        # convert the habitable zone to radians
        hz_center_rad = hz_center / distance_s / (3600 * 180) * np.pi  # in rad

        # put first transmission peak of optimal wl on center of HZ
        # for the origin of the value 0.5.. see Ottiger+2021
        self.data.inst['bl'] = (0.589645 / hz_center_rad
                                * self.data.options.other['wl_optimal'] * 10 ** (-6))

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

    def get_snr(self,
                safe_mode: bool = False):
        """
        Calculates the one-hour signal-to-noise ration for all planets in the catalog.

        Parameters
        ----------
        safe_mode : bool
            If safe mode is enables, the individual photon counts of the planet and noise sources
            are written to the catalog.
        """

        # options are applied before the simulation run
        self.apply_options()

        # currently, the choice of integration time here is arbitrary. Since the background limited
        # case is assumed, the SNR scales with sqrt(integration time) and through this, the SNR
        # for any integration time can be calculated by knowing the SNR of a specific integration
        # time
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

            # TODO: Reinstate the method in which the noise list is keyed by the name of the
            #  producing noise module
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

                # calculate the transmission efficiency of the planets separation
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
                    self.data.catalog.planet_flux_use.iat[n_p] = (
                        [flux_planet_thermal
                         * integration_time
                         * self.data.inst['eff_tot']
                         * self.data.inst['telescope_area']])

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
        Calculate the signal-to-noise ratio per spectral bin of a given spectrum of a single
        planet.

        Parameters
        ----------
        temp_s : float
            Temperature of the observed star in [K].
        radius_s : float
            Radius of the observed star in [sun radii].
        distance_s : float
            Distance between the observed star and the LIFE array in [pc].
        lat_s : float
            Ecliptic latitude of the observed star in [rad].
        z : float
            Zodi level in the observed system in [zodis], i.e. the dust surface density of the
            observed system is z-times as high as in the solar system.
        angsep : float
            Angular separation between the observed star and the observed exoplanet in [arcsec].
        flux_planet_spectrum : list
            Spectrum of the planet. In the first element of the list `flux_planet_spectrum[0]`, the
            wavelength bins of the spectrum must be given in [m]. In the second element
            `flux_planet_spectrum[1]`, the photon count within the spectral bin must be given in
            [photons m-3 s-1].
        integration_time : float
            Time that the LIFE array spends for integrating on the observed planet in [s].

        Returns
        -------
        Tuple[wl_bins, snr_spec]
            Returns the wavelength bins in [m] in the first element and the SNR per wavelength bin
            in the second element.
        flux_planet
            Returns the flux of the planet as used in the simulation in [photons]
        noise
            Returns the noise contribution in [photons]
        """

        # options are applied before the simulation run
        self.apply_options()

        # write the given parameters to the single planet data in the bus. If the connected modules
        # are given an empty index to specify the star, they will use the data saved in this single
        # planet location
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

        # calculate the transmission map
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



