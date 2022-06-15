from warnings import warn

import numpy as np
from tqdm import tqdm
from spectres import spectres
from PyQt5.QtGui import QGuiApplication

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
        baseline = (0.589645 / hz_center_rad
                                * self.data.options.other['wl_optimal'] * 10 ** (-6))

        self.apply_baseline(baseline=baseline)

    def apply_baseline(self,
                       baseline: float,
                       print_warning: bool = False):
        """
        Adjusts the nulling baseline of the array to the specified value.

        Parameters
        ----------
        baseline : float
            Length of the nulling baseline in [m].
        print_warning : bool
            If set to true, function will print a warning if the specified baseline lies outside
            the allow baseline range.
        """
        # make sure that the baseline does not exeed the set baseline limits
        self.data.inst['bl'] = np.maximum(baseline,
                                          self.data.options.array['bl_min'])
        self.data.inst['bl'] = np.minimum(baseline,
                                          self.data.options.array['bl_max'])
        if (self.data.inst['bl'] != baseline) and print_warning:
            warn('Specified baseline exceeded baseline limits. Baseline fixed to '
                 'respective limit')

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
                save_mode: bool = False):
        # options are applied before the simulation run
        self.apply_options()

        # currently, the choice of integration time here is arbitrary. Since the background limited
        # case is assumed, the SNR scales with sqrt(integration time) and through this, the SNR
        # for any integration time can be calculated by knowing the SNR of a specific integration
        # time
        integration_time = 60 * 60

        self.data.catalog['snr_1h'] = np.zeros_like(self.data.catalog.nstar, dtype=float)
        self.data.catalog['baseline'] = np.zeros_like(self.data.catalog.nstar, dtype=float)
        if save_mode:
            self.data.catalog['noise_astro'] = None
            self.data.catalog['planet_flux_use'] = None
            self.data.catalog['photon_rate_planet'] = None
            self.data.catalog['photon_rate_noise'] = None

        # create mask returning only unique stars
        _, temp = np.unique(self.data.catalog.nstar, return_index=True)
        star_mask = np.zeros_like(self.data.catalog.nstar, dtype=bool)
        star_mask[temp] = True

        # iterate over all stars to calculate noise specific to stars
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

            # calculate the noise from the background sources specific to star
            noise_bg_list_star = self.run_socket(s_name='photon_noise_star',
                                                 method='noise',
                                                 index=n)

            if type(noise_bg_list_star) == list:
                noise_bg_star = np.zeros_like(noise_bg_list_star[0])
                for _, noise in enumerate(noise_bg_list_star):
                    noise_bg_star += noise
            else:
                noise_bg_star = noise_bg_list_star

            # calculate the noise from the background sources specific to universe
            noise_bg_list_universe = self.run_socket(s_name='photon_noise_universe',
                                                     method='noise',
                                                     index=n)

            if type(noise_bg_list_universe) == list:
                noise_bg_universe = np.zeros_like(noise_bg_list_universe[0])
                for _, noise in enumerate(noise_bg_list_universe):
                    noise_bg_universe += noise
            else:
                noise_bg_universe = noise_bg_list_star

            # iterate throgh all universes
            universes = np.unique(self.data.catalog.nuniverse[self.data.catalog.nstar == nstar])
            for nuniverse in universes:
                n_u = np.where(np.logical_and(self.data.catalog.nstar == nstar,
                                              self.data.catalog.nuniverse == nuniverse))[0][0]

                noise_bg_universe_temp = noise_bg_universe * self.data.catalog.z.iloc[n_u] / self.data.catalog.z.iloc[n]

                noise_bg = (noise_bg_star + noise_bg_universe_temp) * integration_time * self.data.inst['eff_tot'] * 2

                # go through all planets for the chosen star
                for _, n_p in enumerate(np.argwhere(
                        np.logical_and(self.data.catalog.nstar.to_numpy() == nstar,
                                       self.data.catalog.nuniverse.to_numpy() == nuniverse))[:, 0]):

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
                    flux_planet = (flux_planet_thermal
                                   * transm_eff
                                   * integration_time
                                   * self.data.inst['eff_tot']
                                   * self.data.inst['telescope_area'])
                    noise_planet = (flux_planet_thermal
                                    * transm_noise
                                    * integration_time
                                    * self.data.inst['eff_tot']
                                    * self.data.inst['telescope_area']
                                    * 2)

                    # Add up the noise and caluclate the SNR
                    noise = noise_bg + noise_planet
                    self.data.catalog.snr_1h.iat[n_p] = np.sqrt((flux_planet ** 2 / noise).sum())

                    # save baseline
                    self.data.catalog['baseline'].iat[n_p] = self.data.inst['bl']

                    if save_mode:
                        self.data.catalog.noise_astro.iat[n_p] = [noise_bg]
                        self.data.catalog.planet_flux_use.iat[n_p] = (
                            [flux_planet_thermal
                             * integration_time
                             * self.data.inst['eff_tot']
                             * self.data.inst['telescope_area']])
                        self.data.catalog['photon_rate_planet'].iat[n_p] = (flux_planet
                                                                            / integration_time
                                                                            / self.data.inst['eff_tot']).sum()
                        self.data.catalog['photon_rate_noise'].iat[n_p] = (noise
                                                                           / integration_time
                                                                           / self.data.inst['eff_tot']).sum()

    # TODO: fix units in documentation
    def get_spectrum(self,
                     temp_s: float,  # in K
                     radius_s: float,  # in R_sun
                     distance_s: float,  # in pc
                     lat_s: float,  # in radians
                     z: float,  # in zodis
                     angsep: float,  # in arcsec
                     flux_planet_spectrum: list,  # in ph m-3 s-1 over m
                     integration_time: float,  # in s
                     pbar = None,
                     baseline_to_planet: bool = False,
                     baseline: float = None,
                     safe_mode: bool = False):
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
        pbar
            Takes a PyQt5 QProgressBar to display the progress of the baseline optimization.
        baseline_to_planet : bool
            If set to True, the baseline will be optimized to the position of the planet. If set
            to False, the baseline will be optimized to the center of the habitable zone of the
            host star.
        baseline : float
            Specifies a custom baseline. To have an effect, baseline_to_planet must be set to
            False.

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

        # use spectres to rescale the spectrum onto the correct wl bins
        flux_planet_spectrum_input = flux_planet_spectrum
        flux_planet_spectrum = spectres(new_wavs=self.data.inst['wl_bin_edges'],
                                        spec_wavs=flux_planet_spectrum[0].value,
                                        spec_fluxes=flux_planet_spectrum[1].value,
                                        edge_mode=True)

        # adjust the baseline
        if baseline_to_planet:
            # adjust baseline to planet
            bl = np.linspace(self.data.options.array['bl_min'],
                             self.data.options.array['bl_max'],
                             20)
            snr_analog = np.zeros_like(bl)
            for i in range(len(bl)):
                spec_snr, _, _ = self.get_spectrum(temp_s=temp_s,
                                                   radius_s=radius_s,
                                                   distance_s=distance_s,
                                                   lat_s=lat_s,
                                                   z=z,
                                                   angsep=angsep,
                                                   flux_planet_spectrum=flux_planet_spectrum_input,
                                                   integration_time=integration_time,
                                                   baseline=bl[i])
                snr_analog[i] = np.sqrt((spec_snr[1]**2).sum())
                if pbar is not None:
                    pbar.setValue(30+i/20*30)
                    QGuiApplication.processEvents()
            max_int = np.argmax(snr_analog)

            bl = np.linspace(bl[np.amax((max_int-1, 0))],
                             bl[np.amin((max_int+1, len(bl)-1))],
                             20)
            snr_analog = np.zeros_like(bl)
            for i in range(len(bl)):
                spec_snr, _, _ = self.get_spectrum(temp_s=temp_s,
                                                   radius_s=radius_s,
                                                   distance_s=distance_s,
                                                   lat_s=lat_s,
                                                   z=z,
                                                   angsep=angsep,
                                                   flux_planet_spectrum=flux_planet_spectrum_input,
                                                   integration_time=integration_time,
                                                   baseline=bl[i])
                snr_analog[i] = np.sqrt((spec_snr[1] ** 2).sum())
                if pbar is not None:
                    pbar.setValue(60+i/20*30)
                    QGuiApplication.processEvents()
            self.apply_baseline(baseline=bl[np.argmax(snr_analog)])

        else:
            if baseline is not None:
                # set baseline manually
                self.apply_baseline(baseline=baseline,
                                    print_warning=True)
            else:
                # adjust baseline to HZ
                self.adjust_bl_to_hz(hz_center=hz_center,
                                     distance_s=distance_s)

        # calculate the transmission map
        _, _, self.data.inst['t_map'], _, _ = self.run_socket(s_name='transmission',
                                                              method='transmission_map',
                                                              map_selection='tm3')

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
        # noise_bg_list = self.run_socket(s_name='photon_noise',
        #                                 method='noise',
        #                                 index=None)
        #
        # if type(noise_bg_list) == list:
        #     noise_bg = np.zeros_like(noise_bg_list[0])
        #     for _, noise in enumerate(noise_bg_list):
        #         noise_bg += noise
        # else:
        #     noise_bg = noise_bg_list
        #
        # noise_bg = noise_bg * integration_time * self.data.inst['eff_tot']

        # calculate the noise from the background sources specific to star
        noise_bg_list_star = self.run_socket(s_name='photon_noise_star',
                                             method='noise',
                                             index=None)

        if type(noise_bg_list_star) == list:
            if not noise_bg_list_star:
                noise_bg_star = np.zeros_like(self.data.inst['wl_bins'])
            else:
                noise_bg_star = np.zeros_like(noise_bg_list_star[0])
                for _, noise in enumerate(noise_bg_list_star):
                    noise_bg_star += noise
        else:
            noise_bg_star = noise_bg_list_star

        # calculate the noise from the background sources specific to universe
        noise_bg_list_universe = self.run_socket(s_name='photon_noise_universe',
                                                 method='noise',
                                                 index=None)

        if type(noise_bg_list_universe) == list:
            if not noise_bg_list_universe:
                noise_bg_universe = np.zeros_like(self.data.inst['wl_bins'])
            else:
                noise_bg_universe = np.zeros_like(noise_bg_list_universe[0])
                for _, noise in enumerate(noise_bg_list_universe):
                    noise_bg_universe += noise
        else:
            noise_bg_universe = noise_bg_list_star

        noise_bg = (noise_bg_star + noise_bg_universe) * integration_time * self.data.inst['eff_tot'] * 2

        # Add up the noise and caluclate the SNR
        noise = (noise_bg + noise_planet) * 2
        snr_spec = np.sqrt((flux_planet ** 2 / noise))

        if not safe_mode:
            return ([self.data.inst['wl_bins'], snr_spec],
                    flux_planet,
                    noise)
        else:
            return ([self.data.inst['wl_bins'], snr_spec],
                    flux_planet,
                    [noise, noise_bg_list_star, noise_bg_list_universe])

def multiprocessing_runner(input_dict: dict):
    # create mask returning only unique stars
    universes = np.unique(
        input_dict['catalog'].nuniverse[input_dict['catalog'].nstar == input_dict['nstar']], return_index=False)

    # get transmission map
    _, _, self.data.inst['t_map'], _, _ = self.run_socket(s_name='transmission',
                                                          method='transmission_map',
                                                          map_selection='tm3')

    for nuniverse in universes:
        inst.z = input_dict['catalog'][np.logical_and(input_dict['catalog'].nstar == input_dict['nstar'],
                                                  input_dict['catalog'].nuniverse == nuniverse)].z.iloc[0]

        # redo calculation for exozodi
        inst.create_exozodi()
        inst.sensitivity_coefficients(exozodi_only=True)
        inst.fundamental_noise(exozodi_only=True)

        # go through all planets for the chosen star
        for _, n_p in enumerate(np.argwhere(
                np.logical_and(input_dict['catalog'].nstar.to_numpy() == input_dict['nstar'],
                               input_dict['catalog'].nuniverse.to_numpy() == nuniverse))[:, 0]):
            inst.temp_planet = input_dict['catalog']['temp_p'].iloc[n_p]
            inst.radius_planet = input_dict['catalog']['radius_p'].iloc[n_p]
            inst.separation_planet = (input_dict['catalog']['angsep'].iloc[n_p]
                                      * input_dict['catalog']['distance_s'].iloc[n_p])

            # ----- must be repeated for every planet -----
            inst.create_planet(force=True)
            inst.planet_signal()

            if (inst.chopping == 'nchop'):
                inst.sn_nchop()
            else:
                inst.sn_chop()

            # save baseline
            input_dict['catalog']['baseline'].iat[n_p] = input_dict['baseline']

            # save snr results
            if (inst.chopping == 'nchop'):
                input_dict['catalog'].t_rot.iat[n_p] = input_dict['integration_time']
                input_dict['catalog'].signal.iat[n_p] = inst.photon_rates.loc['signal', 'nchop'].sum()
                input_dict['catalog'].photon_noise.iat[n_p] = (
                    np.sqrt((inst.photon_rates.loc['pn', 'nchop'] ** 2).sum()))
                input_dict['catalog'].systematic_noise.iat[n_p] = (
                    np.sqrt((inst.photon_rates.loc['sn', 'nchop'] ** 2).sum()))
            else:
                input_dict['catalog'].t_rot.iat[n_p] = input_dict['integration_time']
                input_dict['catalog'].signal.iat[n_p] = inst.photon_rates.loc['signal', 'chop'].sum()
                input_dict['catalog'].photon_noise.iat[n_p] = (
                    np.sqrt((inst.photon_rates.loc['pn', 'chop'] ** 2).sum()))
                input_dict['catalog'].systematic_noise.iat[n_p] = (
                    np.sqrt((inst.photon_rates.loc['sn', 'chop'] ** 2).sum()))

            if input_dict['safe_mode']:
                if (inst.chopping == 'nchop'):
                    input_dict['noise_catalog'].loc[input_dict['catalog']['id'].iloc[n_p]] = (
                        inst.photon_rates.nchop)
                else:
                    input_dict['noise_catalog'].loc[input_dict['catalog']['id'].iloc[n_p]] = (
                        inst.photon_rates.chop)

    return_dict = {'catalog': input_dict['catalog']}
    if input_dict['safe_mode']:
        return_dict['noise_catalog'] = input_dict['noise_catalog']

    return return_dict



