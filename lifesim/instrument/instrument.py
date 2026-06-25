from warnings import warn
from copy import deepcopy
from typing import List, Optional, Sequence, Any

import numpy as np
import pandas as pd
from tqdm import tqdm
from spectres import spectres
from PyQt5.QtGui import QGuiApplication
from joblib import Parallel, delayed, parallel_config
from joblib_progress import joblib_progress

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

        self.data.inst['telescope_area'] = self.data.options.array['num_apertures'] * np.pi \
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

        # set size of the integrated image, adjusted to a threshold measured on the FoV taper
        if self.data.options.models['fov_taper'] == 'gaussian':
            self.data.inst['image_angle'] = self.data.inst['hfov'] * 4 / np.pi * np.sqrt(
                -np.log(self.data.options.other['fov_threshold'])
            )
        elif self.data.options.models['fov_taper'] == 'none':
            self.data.inst['image_angle'] = self.data.inst['hfov']
        else:
            raise ValueError('Nonexistent fov taper model')

        self.data.inst['image_angle_mas'] = self.data.inst['image_angle'] * (3600000. * 180.) / np.pi

        self.data.inst['rad_pix'] = (2 * self.data.inst['image_angle']) \
                                    / self.data.options.other['image_size']  # Radians per pixel
        self.data.inst['mas_pix'] = (2 * self.data.inst['image_angle_mas']) \
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
        # for the origin of the value 0.5.. see Dannert+2022
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
        if self.data.options.array['fixed_baseline']:
            self.data.inst['bl'] = self.data.options.array['baseline']

        else:
        # make sure that the baseline does not exeed the set baseline limits
            self.data.inst['bl'] = np.maximum(baseline,
                                              self.data.options.array['bl_min'])
            self.data.inst['bl'] = np.minimum(self.data.inst['bl'],
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

    def _unpack_socket(self,
                       socket_return):
        """
        Checks if the object returned from the socket is a list. If it is, add all values in the list.
        """

        if type(socket_return) == list:
            if not socket_return:
                output = np.zeros_like(self.data.inst['wl_bins'])
            else:
                output = np.zeros_like(socket_return[0])
                for _, noise in enumerate(socket_return):
                    output += noise
        else:
            output = socket_return

        return output

    def get_snr(self,
                save_mode: bool = False):
        """
        Calculates the one-hour signal-to-noise ration for all planets in the catalog. Switches between single and
        multi processing depending on the n_cpu option.

        Parameters
        ----------
        safe_mode : bool
            If save mode is enables, the individual photon counts of the planet and noise sources
            are written to the catalog.
        """

        # options are applied before the simulation run
        self.apply_options()

        self.data.catalog['snr_1h'] = np.zeros_like(self.data.catalog.nstar, dtype=float)
        self.data.catalog['baseline'] = np.zeros_like(self.data.catalog.nstar, dtype=float)
        if save_mode:
            self.data.catalog['noise_astro'] = None
            self.data.catalog['planet_flux_use'] = None
            self.data.catalog['photon_rate_planet'] = None
            self.data.catalog['photon_rate_noise'] = None

        if self.data.options.other['n_cpu'] == 1:
            self.get_snr_single_processing(save_mode=save_mode)
        elif self.data.options.other['n_cpu'] > 1:
            self.get_snr_multi_processing(save_mode=save_mode)
        else:
            raise ValueError('n_cpu option must be >= 1')

    def get_snr_single_processing(self,
                                  save_mode: bool = False,
                                  verbose : bool = True):
        """
        Calculates the one-hour signal-to-noise ration for all planets in the catalog.
        Parameters
        ----------
        safe_mode : bool
            If save mode is enables, the individual photon counts of the planet and noise sources
            are written to the catalog.
        """

        # currently, the choice of integration time here is arbitrary. Since the background limited
        # case is assumed, the SNR scales with sqrt(integration time) and through this, the SNR
        # for any integration time can be calculated by knowing the SNR of a specific integration
        # time
        integration_time = 60 * 60

        # calculate instrument noise once, since it is the same for all planets
        noise_list_thermal = self.run_socket(s_name='photon_noise_instrument',
                                             method='noise',
                                             index=None)

        noise_thermal = self._unpack_socket(noise_list_thermal)

        noise_inst = (
                np.sum(noise_thermal, axis=0)
                * integration_time
                * self.data.options.array['quantum_eff']
                * self.data.options.array['num_outputs']
        )

        # if type(noise_list_thermal) == list:
        #     if not noise_list_thermal:
        #         noise_thermal = np.zeros_like(self.data.inst['wl_bins'])
        #     else:
        #         noise_thermal = np.zeros_like(noise_list_thermal[0])
        #         for _, noise in enumerate(noise_list_thermal):
        #             noise_thermal += noise
        # else:
        #     noise_thermal = noise_list_thermal
        
        # noise_inst = (noise_thermal[0] * integration_time * self.data.inst['eff_tot'] * self.data.options.array['num_outputs']) \
        #                  + (noise_thermal[1] * integration_time * self.data.options.array['quantum_eff'] * self.data.options.array['num_outputs'])

        # calculate the dark current noise from the detector once, since it is the same for all planets
        noise_dc_list = self.run_socket(s_name='electron_noise_detector',
                                        method='noise',
                                        index=None)
        
        if type(noise_dc_list) == list:
            if not noise_dc_list:
                noise_dc_d = np.zeros_like(self.data.inst['wl_bins'])
            else:
                noise_dc_d = np.zeros_like(noise_dc_list[0])
                for _, noise in enumerate(noise_dc_list):
                    noise_dc_d += noise
        else:
            noise_dc_d = noise_dc_list
        
        noise_dc = noise_dc_d * integration_time

        # create mask returning only unique stars
        _, temp = np.unique(self.data.catalog.nstar, return_index=True)
        star_mask = np.zeros_like(self.data.catalog.nstar, dtype=bool)
        star_mask[temp] = True

        # iterate over all stars to calculate noise specific to stars
        for i, n in enumerate(tqdm(np.where(star_mask)[0], disable=not verbose)):
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

                noise_bg_universe_temp = (noise_bg_universe * self.data.catalog.z.iloc[n_u]
                                      / self.data.catalog.z.iloc[n])

                noise_bg = ((noise_bg_star + noise_bg_universe_temp)
                        * integration_time * self.data.inst['eff_tot'] * self.data.options.array['num_outputs'])

                # go through all planets for the chosen star
                for _, n_p in enumerate(np.argwhere(
                        np.logical_and(self.data.catalog.nstar.to_numpy() == nstar,
                                       self.data.catalog.nuniverse.to_numpy() == nuniverse))[:, 0]):

                    # calculate the photon flux originating from the planet
                    flux_planet_thermal = black_body(
                        mode='planet',
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
                                    * self.data.options.array['num_outputs'])

                    # Add up the noise and calculate the SNR
                    noise = noise_bg + noise_planet + noise_inst + noise_dc

                    # use index label to avoid chained assignment / view-copy problems
                    idx_label = self.data.catalog.index[n_p]
                    self.data.catalog.loc[idx_label, 'snr_1h'] = np.sqrt((flux_planet ** 2 / noise).sum())

                    if self.data.options.optimization['iwa_cut'] is not None:
                        curve_chop, _ = self.run_socket(s_name='transmission',
                                                        method='transmission_curve',
                                                        angsep=self.data.catalog.angsep.iloc[n_p])
                        if np.min(np.max(curve_chop[:, 0, :], axis=1)) < self.data.options.optimization['iwa_cut']:
                            self.data.catalog.loc[idx_label, 'snr_1h'] = 0.

                    # save baseline
                    self.data.catalog.loc[idx_label, 'baseline'] = self.data.inst['bl']

                    if save_mode:
                        self.data.catalog.loc[idx_label, 'noise_astro'] = [noise_bg]
                        self.data.catalog.loc[idx_label, 'planet_flux_use'] = (
                            [flux_planet_thermal
                             * integration_time
                             * self.data.inst['eff_tot']
                             * self.data.inst['telescope_area']])
                        self.data.catalog.loc[idx_label, 'photon_rate_planet'] = (
                                flux_planet
                                / integration_time
                                / self.data.inst['eff_tot']
                        ).sum()
                        self.data.catalog.loc[idx_label, 'photon_rate_noise'] = (
                                noise
                                / integration_time
                                / self.data.inst['eff_tot']
                        ).sum()

        # ...existing code...

    def get_snr_multi_processing(self,
                                    save_mode: bool = False):

        # divide the catalog into roughly equal chunks for each cpu
        n_star, occ_star = np.unique(self.data.catalog.nstar, return_counts=True)
        star_groups = balanced_partition_greedy(occ=occ_star, items=n_star, n_groups=self.data.options.other['n_cpu']*10)

        sub_catalogs = [self.data.catalog[np.isin(self.data.catalog.nstar, sg)] for sg in star_groups]
        group_sizes = [len(sc) for sc in sub_catalogs]
        per_dev = (np.max(group_sizes) - np.min(group_sizes)) / np.mean(group_sizes)
        print(f'Maximum deviation in group sizes: {per_dev*100:.1f}%')

        reference_bus = deepcopy(self)
        del reference_bus.data.catalog

        with parallel_config(
                backend="loky", inner_max_num_threads=1
        ), joblib_progress(
            description="Running SNR calculation in parallel ...",
            total=len(star_groups),
        ):
            results = Parallel(n_jobs=self.data.options.other['n_cpu'])(
                delayed(mp_runner)(
                    bus=reference_bus,
                    catalog=sc,
                    save_mode=save_mode
                )
                for sc in sub_catalogs)

        # combine results back into main catalog
        self.data.catalog = pd.concat(results)
        self.data.catalog.sort_values('id', inplace=True)


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
        flux_planet = (flux_planet_spectrum
                       * transm_eff
                       * integration_time
                       * self.data.inst['eff_tot']
                       * self.data.inst['telescope_area']
                       * self.data.inst['wl_bin_widths'])
        noise_planet = (flux_planet_spectrum
                        * transm_noise
                        * integration_time
                        * self.data.inst['eff_tot']
                        * self.data.inst['telescope_area']
                        * self.data.inst['wl_bin_widths']
                        * self.data.options.array['num_outputs'])

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
            noise_bg_universe = noise_bg_list_universe

        noise_bg = ((noise_bg_star + noise_bg_universe)
                    * integration_time * self.data.inst['eff_tot'] * self.data.options.array['num_outputs'])
        
        # calculate the thermal noise from the instrument
        noise_list_thermal = self.run_socket(s_name='photon_noise_instrument',
                                             method='noise',
                                             index=None)

        noise_thermal = self._unpack_socket(noise_list_thermal)

        noise_inst = (
                np.sum(noise_thermal, axis=0)
                * integration_time
                * self.data.options.array['quantum_eff']
                * self.data.options.array['num_outputs']
        )

        # if type(noise_list_thermal) == list:
        #     if not noise_list_thermal:
        #         noise_thermal = np.zeros_like(self.data.inst['wl_bins'])
        #     else:
        #         noise_thermal = np.zeros_like(noise_list_thermal[0])
        #         for _, noise in enumerate(noise_list_thermal):
        #             noise_thermal += noise
        # else:
        #     noise_thermal = noise_list_thermal
        #
        # # output is two arrays (due to mirror and detector leakage) so combine like this
        # noise_inst = (noise_thermal[0] * integration_time * self.data.inst['eff_tot'] * self.data.options.array['num_outputs']) \
        #                  + (noise_thermal[1] * integration_time * self.data.options.array['quantum_eff'] * self.data.options.array['num_outputs'])

        # calculate the dark current noise from the detector
        noise_dc_list = self.run_socket(s_name='electron_noise_detector',
                                                 method='noise',
                                                 index=None)
        if type(noise_dc_list) == list:
            if not noise_dc_list:
                noise_dc_d = np.zeros_like(self.data.inst['wl_bins'])
            else:
                noise_dc_d = np.zeros_like(noise_dc_list[0])
                for _, noise in enumerate(noise_dc_list):
                    noise_dc_d += noise
        else:
            noise_dc_d = noise_dc_list
        
        noise_dc = noise_dc_d * integration_time

        # Add up the noise and calculate the SNR
        noise = (noise_bg + noise_planet + noise_inst + noise_dc)
        snr_spec = np.sqrt((flux_planet ** 2 / noise))

        if not safe_mode:
            return ([self.data.inst['wl_bins'], snr_spec],
                    flux_planet,
                    noise)
        else:
            return ([self.data.inst['wl_bins'], snr_spec],
                    flux_planet,
                    [noise, noise_bg_list_star, noise_bg_list_universe, noise_thermal, noise_dc_list, noise_inst])


    def get_signal(self,
                   temp_s: float,  # in K
                   radius_s: float,  # in R_sun
                   distance_s: float,  # in pc
                   lat_s: float,  # in radians
                   z: float,  # in zodis
                   angsep: float,  # in arcsec
                   flux_planet_spectrum: list,  # in ph m-3 s-1 over m
                   integration_time: float,  # in s
                   phi_n: int = 360):
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

        # TODO: remove by 2024
        warn('The get_spectrum function was implemented with a major bug between versions 0.2.16 '
             'and 0.2.24 in which the noise level was twice as large as the correct value. If '
             'you created results with the versions in question, please validate them with the '
             'latest version of LIFEsim.')

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

        # adjust baseline to habitable zone
        self.adjust_bl_to_hz(hz_center=hz_center,
                             distance_s=distance_s)

        # use spectres to rescale the spectrum onto the correct wl bins
        flux_planet_spectrum_input = flux_planet_spectrum
        flux_planet_spectrum = spectres(new_wavs=self.data.inst['wl_bin_edges'],
                                        spec_wavs=flux_planet_spectrum[0].value,
                                        spec_fluxes=flux_planet_spectrum[1].value,
                                        edge_mode=True)

        # calculate the transmission map
        _, _, self.data.inst['t_map'], _, _ = self.run_socket(s_name='transmission',
                                                              method='transmission_map',
                                                              map_selection='tm3')

        curve_chop, curve_tm4 = self.run_socket(s_name='transmission',
                                                method='transmission_curve',
                                                angsep=angsep,
                                                phi_n=phi_n)

        # calculate the signal and photon noise flux received from the planet per time bin
        flux_planet = (flux_planet_spectrum[:, np.newaxis]
                      * np.squeeze(curve_chop, axis=1)
                      * integration_time
                      / phi_n
                      * self.data.inst['eff_tot']
                      * self.data.inst['telescope_area']
                      * self.data.inst['wl_bin_widths'][:, np.newaxis])
        noise_planet = (flux_planet_spectrum[:, np.newaxis]
                       * np.squeeze(curve_tm4, axis=1)
                       / phi_n
                       * integration_time
                       * self.data.inst['eff_tot']
                       * self.data.inst['telescope_area']
                       * self.data.inst['wl_bin_widths'][:, np.newaxis])


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

        noise_bg = ((noise_bg_star + noise_bg_universe)
                    * integration_time / phi_n * self.data.inst['eff_tot'] * self.data.options.array['num_outputs'])

        noise = (noise_bg[:, np.newaxis] + noise_planet)

        # draw noise
        noise_drawn = np.random.poisson(lam=noise) - np.random.poisson(lam=noise)

        # add up to noisy signal
        signal = flux_planet + noise_drawn

        return signal, flux_planet

def mp_runner(bus,
              catalog,
              save_mode: bool = False):
    bus = deepcopy(bus)
    bus.data.catalog = catalog
    bus.get_snr_single_processing(save_mode=save_mode,
                                  verbose=False)
    return bus.data.catalog

def balanced_partition_greedy(occ: Sequence[int],
                                    n_groups: int,
                                    items: Optional[Sequence[Any]] = None) -> List[List[Any]]:
    """
    Greedy multi-way partition using NumPy for faster min selection.
    Assigns the largest items first to the current smallest-sum group.

    Parameters
    ----------
    occ : Sequence[int]
        Sequence containing the weight (e.g. number of sub-objects) for each item.
    n_groups : int
        Number of groups to partition into. Must be >= 1.
    items : Optional[Sequence[Any]]
        Sequence of items corresponding to `occ`. If ``None``, the function will use the
        indices ``range(len(occ))`` as the items.

    Returns
    -------
    List[List[Any]]
        List of length ``n_groups`` where each element is a list of the assigned items
        (or indices if ``items`` was ``None``). Group totals are balanced using a greedy
        heuristic; some groups may remain empty if ``n_groups`` > ``len(occ)``.

    Raises
    ------
    ValueError
        If ``n_groups`` < 1 or if ``items`` is provided but its length does not match ``occ``.
    """
    if n_groups <= 0:
        raise ValueError("n_groups must be >= 1")
    n = len(occ)
    if items is None:
        items = list(range(n))
    if len(items) != n:
        raise ValueError("items and occ must have same length")

    occ_arr = np.asarray(occ, dtype=np.int64)
    # sort indices by descending weight (largest first)
    idxs = np.argsort(-occ_arr)

    groups: List[List[Any]] = [[] for _ in range(n_groups)]
    sums = np.zeros(n_groups, dtype=np.int64)

    # assign each item to the group with the smallest current sum
    for i in idxs:
        g = int(np.argmin(sums))  # fast C-level operation
        groups[g].append(items[int(i)])
        sums[g] += int(occ_arr[int(i)])

    return groups
