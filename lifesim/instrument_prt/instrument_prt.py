import numpy as np
from tqdm import tqdm
from spectres import spectres

from lifesim.core.modules import InstrumentModule
import inlifesim as ils
from lifesim.instrument.instrument import Instrument
from lifesim.util.habitable import single_habitable_zone


class InstrumentPrt(InstrumentModule):
    """
        XXX

        Notes
        -----
        Note, that all attributes are saved in the data class.

        Attributes
        ----------
        data.options : lifesim.Options
            The options class containing all setting for the array and computations.
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
        self.add_socket(s_name='instrument',
                        s_type=InstrumentModule,
                        s_number=1)

    def apply_options(self,
                      hz_center,
                      distance_s):
        self.run_socket(s_name='instrument',
                        method='apply_options')
        self.run_socket(s_name='instrument',
                        method='adjust_bl_to_hz',
                        hz_center=hz_center,
                        distance_s=distance_s)

    def get_snr(self,
                safe_mode=True):
        # currently, the choice of integration time here is arbitrary. Since the background limited
        # case is assumed, the SNR scales with sqrt(integration time) and through this, the SNR
        # for any integration time can be calculated by knowing the SNR of a specific integration
        # time
        integration_time = 60 * 60

        self.data.catalog['snr_1h'] = np.zeros_like(self.data.catalog.nstar, dtype=float)
        self.data.catalog['snr_1h_fundamental'] = np.zeros_like(self.data.catalog.nstar, dtype=float)
        self.data.catalog['baseline'] = np.zeros_like(self.data.catalog.nstar, dtype=float)
        if safe_mode:
            self.data.noise_catalog_from_catalog()

        # create mask returning only unique stars
        _, temp = np.unique(self.data.catalog.nstar, return_index=True)
        star_mask = np.zeros_like(self.data.catalog.nstar, dtype=bool)
        star_mask[temp] = True

        # iterate over all stars
        for i, n in enumerate(tqdm(np.where(star_mask)[0])):
            nstar = self.data.catalog.nstar.iloc[n]

            # adjust baseline of array and give new baseline to transmission generator plugin
            self.apply_options(hz_center=float(self.data.catalog.hz_center.iloc[n]),
                               distance_s=float(self.data.catalog.distance_s.iloc[n]))

            col_pos = np.array(((-self.data.inst['bl'] / 2, -self.data.inst['bl']*self.data.options.array['ratio'] / 2),
                                (-self.data.inst['bl'] / 2, self.data.inst['bl']*self.data.options.array['ratio'] / 2),
                                (self.data.inst['bl'] / 2, -self.data.inst['bl']*self.data.options.array['ratio'] / 2),
                                (self.data.inst['bl'] / 2, self.data.inst['bl']*self.data.options.array['ratio'] / 2)))

            # TODO: Correct treatment of quantum efficiency
            inst = ils.Instrument(
                # ----- static parameters -----
                wl_bins=self.data.inst['wl_bins'],  # wavelength bins center position in m
                wl_bin_widths=self.data.inst['wl_bin_widths'],  # wavelength bin widhts in m
                integration_time=integration_time,
                image_size=self.data.options.other['image_size'],  # size of image used to simulate exozodi in pix
                diameter_ap=self.data.options.array['diameter'],  # diameter of the primary mirrors in m
                flux_division=self.data.options.array['flux_division'],
                # division of the flux between the primary mirrors, e.g. in baseline case [0.25, 0.25, 0.25, 0.25]
                throughput=self.data.options.array['throughput']*self.data.options.array['quantum_eff'],
                # fraction of light that is sustained through the optical train
                phase_response=self.data.options.array['phase_response'],  # phase response of each collector arm in rad
                phase_response_chop=self.data.options.array['phase_response_chop'],
                # phase response of each collector arm in the chopped state in rad
                t_rot=self.data.options.array['t_rot'],  # rotation period of the array in seconds
                chopping=self.data.options.array['chopping'],
                # run calculation with or without chopping, 'chop', 'nchop', 'both'
                pix_per_wl=self.data.options.array['pix_per_wl'],  # pixels on detector used per wavelength channel
                n_sampling_rot=self.data.options.array['n_sampling_rot'],
                # number of sampling points per array rotation
                detector_dark_current='manual',
                # detector type, 'MIRI' or 'manual'. Specify dark_current_pix in 'manual'
                dark_current_pix=0.,  # detector dark current in electrons s-1 px-1
                detector_thermal='MIRI',  # detector type, 'MIRI'
                det_temp=0.,  # temperature of the detector environment in K
                magnification=15.73,  # telescope magnification
                f_number=20.21,  # telescope f-number, i.e. ratio of focal length to aperture size
                secondary_primary_ratio=0.114,  # ratio of secondary to primary mirror sizes
                primary_emissivity=0.,  # emissivity epsilon of the primary mirror
                primary_temp=0.,  # temperature of the primary mirror in K
                pink_noise_co=10000,  # cutoff frequency for the pink noise spectra
                n_cpu=1,  # number of cores used in the simulation
                rms_mode='lay',  # mode for rms values, 'lay', 'static', 'wavelength'
                agnostic_mode=True,  # derive instrumental photon noise from agnostic mode
                eps_cold=0.,  # scaling constant for cold agnostic photon noise spectrum
                eps_hot=0.,  # scaling constant for hot agnostic photon noise spectrum
                eps_white=0.,  # scaling constant white agnostic photon noise spectrum
                agnostic_spacecraft_temp=0.,  # cold-side spacecraft temperature in the
                # agnostic case
                n_sampling_max=10000,  # largest fourier mode used in noise sampling
                d_a_rms=None,  # relative amplitude error rms
                d_phi_rms=None,  # phase error rms
                d_pol_rms=None,  # polarization error rms
                d_x_rms=None,  # collector position rms, x-direction
                d_y_rms=None,  # collector position rms, y-direction
                wl_resolution=200,  # number of wavelength bins simulated for the thermal background
                flux_planet=None,  # substitute flux input in ph m-2 s-1
                simultaneous_chopping=True,
                # ----- parameters change with star -----
                dist_star=self.data.catalog.distance_s.iloc[n],  # distance to the target system in pc
                radius_star=self.data.catalog.radius_s.iloc[n],  # radius of the star in stellar radii
                temp_star=self.data.catalog.temp_s.iloc[n],  # temperature of the host star in Kelvin
                lat_star=self.data.catalog.lat.iloc[n],  # ecliptic latitude of the target star
                l_sun=self.data.catalog.l_sun.iloc[n],  # stellar luminosity in solar luminosities
                z=self.data.catalog.z.iloc[n],
                # zodi level: the exozodi dust is z-times denser than the localzodi dust
                col_pos=col_pos,  # collector position in m
                # ----- parameters change with planet -----
                temp_planet=0.,  # planet temperature in Kelvin
                radius_planet=0.,  # planet radius in earth radii
                separation_planet=0.,  # separation of target planet from host star in AU
            )

            # ----- same for every star -----
            inst.instrumental_parameters()
            inst.create_star()
            inst.create_localzodi()
            inst.create_exozodi()
            inst.sensitivity_coefficients()
            inst.fundamental_noise()

            if inst.agnostic_mode:
                inst.pn_agnostic()
            else:
                inst.pn_dark_current()
                inst.pn_thermal_background_detector()
                inst.pn_thermal_primary_mirror()

            # create mask returning only unique stars
            universes = np.unique(self.data.catalog.nuniverse[self.data.catalog.nstar == nstar], return_index=False)

            for nuniverse in universes:
                inst.z = self.data.catalog[np.logical_and(self.data.catalog.nstar == nstar,
                                                          self.data.catalog.nuniverse == nuniverse)].z.iloc[0]

                # redo calculation for exozodi
                inst.create_exozodi()
                inst.sensitivity_coefficients(exozodi_only=True)
                inst.fundamental_noise(exozodi_only=True)

                # go through all planets for the chosen star
                for _, n_p in enumerate(np.argwhere(
                        np.logical_and(self.data.catalog.nstar.to_numpy() == nstar,
                                       self.data.catalog.nuniverse.to_numpy() == nuniverse))[:, 0]):
                    inst.temp_planet = self.data.catalog['temp_p'].iloc[n_p]
                    inst.radius_planet = self.data.catalog['radius_p'].iloc[n_p]
                    inst.separation_planet = (self.data.catalog['angsep'].iloc[n_p]
                                              * self.data.catalog['distance_s'].iloc[n_p])

                    # ----- must be repeated for every planet -----
                    inst.create_planet(force=True)
                    inst.planet_signal()

                    if (inst.chopping == 'nchop'):
                        inst.sn_nchop()
                    else:
                        inst.sn_chop()

                    # save baseline
                    self.data.catalog['baseline'].iat[n_p] = self.data.inst['bl']

                    # save snr results
                    if (inst.chopping == 'nchop'):
                        self.data.catalog.snr_1h.iat[n_p] = (np.sqrt(
                            (inst.photon_rates.loc['snr', 'nchop'] ** 2).sum()))
                        self.data.catalog.snr_1h_fundamental.iat[n_p] = (
                                np.sqrt(((inst.photon_rates.loc['signal', 'nchop']
                                          / inst.photon_rates.loc['fundamental', 'nchop']) ** 2).sum()))
                    else:
                        self.data.catalog.snr_1h.iat[n_p] = (np.sqrt((inst.photon_rates.loc['snr', 'chop'] ** 2).sum()))
                        self.data.catalog.snr_1h_fundamental.iat[n_p] = (
                                np.sqrt(((inst.photon_rates.loc['signal', 'chop']
                                          / inst.photon_rates.loc['fundamental', 'chop']) ** 2).sum()))

                    if safe_mode:
                        if (inst.chopping == 'nchop'):
                            self.data.noise_catalog.loc[self.data.catalog['id'].iloc[n_p]] = inst.photon_rates.nchop
                        else:
                            self.data.noise_catalog.loc[self.data.catalog['id'].iloc[n_p]] = inst.photon_rates.chop

    def get_spectrum(self,
                     temp_s: float,  # in K
                     radius_s: float,  # in R_sun
                     distance_s: float,  # in pc
                     lat_s: float,  # in radians
                     z: float,  # in zodis
                     angsep: float,  # in arcsec
                     flux_planet_spectrum: list,  # in ph m-3 s-1 over m
                     integration_time: float,  # in s
                     pbar: bool = None,
                     baseline_to_planet: bool = False,
                     baseline: float = None,
                     safe_mode: bool = False):

        # TODO: Implement baseline_to_planet option

        # calculate the habitable zone of the specified star
        s_in, s_out, l_sun, \
        hz_in, hz_out, \
        hz_center = single_habitable_zone(model=self.data.options.models['habitable'],
                                          temp_s=temp_s,
                                          radius_s=radius_s)

        flux_planet_spectrum_input = flux_planet_spectrum
        flux_planet_spectrum = spectres(new_wavs=self.data.inst['wl_bin_edges'],
                                        spec_wavs=flux_planet_spectrum[0].value,
                                        spec_fluxes=flux_planet_spectrum[1].value,
                                        edge_mode=True)

        flux_planet_spectrum *= self.data.inst['wl_bin_widths']

        self.run_socket(s_name='instrument',
                        method='apply_options')

        if baseline is not None:
            # set baseline manually
            self.run_socket(s_name='instrument',
                            method='apply_baseline',
                            baseline=baseline,
                            print_warning=True)
        else:
            # adjust baseline to HZ
            self.run_socket(s_name='instrument',
                            method='adjust_bl_to_hz',
                            hz_center=hz_center,
                            distance_s=distance_s)

        col_pos = np.array(((-self.data.inst['bl'] / 2, -self.data.inst['bl'] * self.data.options.array['ratio'] / 2),
                            (-self.data.inst['bl'] / 2, self.data.inst['bl'] * self.data.options.array['ratio'] / 2),
                            (self.data.inst['bl'] / 2, -self.data.inst['bl'] * self.data.options.array['ratio'] / 2),
                            (self.data.inst['bl'] / 2, self.data.inst['bl'] * self.data.options.array['ratio'] / 2)))

        inst = ils.Instrument(
            # ----- static parameters -----
            wl_bins=self.data.inst['wl_bins'],  # wavelength bins center position in m
            wl_bin_widths=self.data.inst['wl_bin_widths'],  # wavelength bin widhts in m
            integration_time=integration_time,
            image_size=self.data.options.other['image_size'],  # size of image used to simulate exozodi in pix
            diameter_ap=self.data.options.array['diameter'],  # diameter of the primary mirrors in m
            flux_division=self.data.options.array['flux_division'],
            # division of the flux between the primary mirrors, e.g. in baseline case [0.25, 0.25, 0.25, 0.25]
            throughput=self.data.options.array['throughput']*self.data.options.array['quantum_eff'],
            # fraction of light that is sustained through the optical train
            phase_response=self.data.options.array['phase_response'],  # phase response of each collector arm in rad
            phase_response_chop=self.data.options.array['phase_response_chop'],
            # phase response of each collector arm in the chopped state in rad
            t_rot=self.data.options.array['t_rot'],  # rotation period of the array in seconds
            chopping=self.data.options.array['chopping'],
            # run calculation with or without chopping, 'chop', 'nchop', 'both'
            pix_per_wl=self.data.options.array['pix_per_wl'],  # pixels on detector used per wavelength channel
            n_sampling_rot=self.data.options.array['n_sampling_rot'],
            # number of sampling points per array rotation
            detector_dark_current='manual',
            # detector type, 'MIRI' or 'manual'. Specify dark_current_pix in 'manual'
            dark_current_pix=1e-4,  # detector dark current in electrons s-1 px-1
            detector_thermal='MIRI',  # detector type, 'MIRI'
            det_temp=11.,  # temperature of the detector environment in K
            magnification=15.73,  # telescope magnification
            f_number=20.21,  # telescope f-number, i.e. ratio of focal length to aperture size
            secondary_primary_ratio=0.114,  # ratio of secondary to primary mirror sizes
            primary_emissivity=0.,  # emissivity epsilon of the primary mirror
            primary_temp=0.,  # temperature of the primary mirror in K
            pink_noise_co=10000,  # cutoff frequency for the pink noise spectra
            n_cpu=1,  # number of cores used in the simulation
            rms_mode='lay',  # mode for rms values, 'lay', 'static', 'wavelength'
            agnostic_mode=False,  # derive instrumental photon noise from agnostic mode
            eps_cold=0.,  # scaling constant for cold agnostic photon noise spectrum
            eps_hot=0.,  # scaling constant for hot agnostic photon noise spectrum
            eps_white=0.,  # scaling constant white agnostic photon noise spectrum
            agnostic_spacecraft_temp=0.,  # cold-side spacecraft temperature in the
            # agnostic case
            n_sampling_max=self.data.options.other['n_sampling_max'],  # largest fourier mode used in noise sampling
            d_a_rms=None,  # relative amplitude error rms
            d_phi_rms=None,  # phase error rms
            d_pol_rms=None,  # polarization error rms
            d_x_rms=None,  # collector position rms, x-direction
            d_y_rms=None,  # collector position rms, y-direction
            wl_resolution=200,  # number of wavelength bins simulated for the thermal background
            flux_planet=flux_planet_spectrum,  # substitute flux input in ph m-2 s-1
            simultaneous_chopping=True,
            # ----- parameters change with star -----
            dist_star=distance_s,  # distance to the target system in pc
            radius_star=radius_s,  # radius of the star in stellar radii
            temp_star=temp_s,  # temperature of the host star in Kelvin
            lat_star=lat_s,  # ecliptic latitude of the target star
            l_sun=l_sun,  # stellar luminosity in solar luminosities
            z=z,
            # zodi level: the exozodi dust is z-times denser than the localzodi dust
            col_pos=col_pos,  # collector position in m
            # ----- parameters change with planet -----
            temp_planet=0.,  # planet temperature in Kelvin
            radius_planet=0.,  # planet radius in earth radii
            separation_planet=angsep * distance_s,  # separation of target planet from host star in AU
        )

        inst.run()

        return inst.photon_rates

