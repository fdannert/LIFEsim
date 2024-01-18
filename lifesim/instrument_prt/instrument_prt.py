import multiprocessing as mp
from copy import deepcopy
from typing import Union

import numpy as np
from tqdm import tqdm
from spectres import spectres
import pandas as pd
import xarray as xr

from lifesim.core.modules import InstrumentModule
from inlifesim.instrument import Instrument
#from lifesim.instrument.instrument import Instrument
from lifesim.util.habitable import single_habitable_zone
from lifesim.instrument.instrument import adjust_sampling


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
        self.inst_prt = None

    def apply_options(self,
                      hz_center: float = 0.,
                      distance_s: float = 0.,
                      run_baseline: bool = True):
        self.run_socket(s_name='instrument',
                        method='apply_options')
        if run_baseline:
            self.run_socket(s_name='instrument',
                            method='adjust_bl_to_hz',
                            hz_center=hz_center,
                            distance_s=distance_s)

    def get_snr(self,
                safe_mode:bool = True,
                lookup_table:str = 'none'):
        '''
        Calculate the SNR for all stars in the catalog.
        Parameters
        ----------
        safe_mode : bool
        lookup_table : str, None
            If set to 'none', no lookup table is used. If set to 'input:path' the lookup
            table is read from the specified path. If set to 'output:path' the lookup table is saved
            to the specified path.

        Returns
        -------

        '''

        if safe_mode and (lookup_table.split(':')[0] == 'output'):
            raise ValueError('Save mode cannot be used when creating a lookup table.')

        # currently, the choice of integration time here is arbitrary. Since the background limited
        # case is assumed, the SNR scales with sqrt(integration time) and through this, the SNR
        # for any integration time can be calculated by knowing the SNR of a specific integration
        # time
        integration_time = self.data.options.array['t_rot']

        self.data.catalog['t_rot'] = np.zeros_like(self.data.catalog.nstar, dtype=float)
        self.data.catalog['signal'] = np.zeros_like(self.data.catalog.nstar, dtype=float)
        self.data.catalog['photon_noise'] = np.zeros_like(self.data.catalog.nstar, dtype=float)
        self.data.catalog['systematic_noise'] = np.zeros_like(self.data.catalog.nstar, dtype=float)
        self.data.catalog['baseline'] = np.zeros_like(self.data.catalog.nstar, dtype=float)
        self.data.catalog['n_sampling_rot'] = np.zeros_like(self.data.catalog.nstar, dtype=int)
        self.data.catalog['image_size'] = np.zeros_like(self.data.catalog.nstar, dtype=int)

        if safe_mode:
            self.data.noise_catalog = {}

        # load lookup table
        if lookup_table.split(':')[0] == 'input':
            lookup_table_in = pd.read_hdf(lookup_table.split(':')[1]).to_dict()

        # create mask returning only unique stars
        _, temp = np.unique(self.data.catalog.nstar, return_index=True)
        star_mask = np.zeros_like(self.data.catalog.nstar, dtype=bool)
        star_mask[temp] = True

        # create list of multiprocessing input dictionaries
        input_dict_list = []

        # iterate over all stars
        print('\nPreparing sample...')
        for i, n in enumerate(tqdm(np.where(star_mask)[0])):
            nstar = self.data.catalog.nstar.iloc[n]

            # adjust baseline of array and give new baseline to transmission generator plugin
            self.apply_options(hz_center=float(self.data.catalog.hz_center.iloc[n]),
                               distance_s=float(self.data.catalog.distance_s.iloc[n]))

            col_pos = np.array(((-self.data.inst['bl'] / 2,
                                 -self.data.inst['bl']*self.data.options.array['ratio'] / 2),
                                (-self.data.inst['bl'] / 2,
                                 self.data.inst['bl']*self.data.options.array['ratio'] / 2),
                                (self.data.inst['bl'] / 2,
                                 -self.data.inst['bl']*self.data.options.array['ratio'] / 2),
                                (self.data.inst['bl'] / 2,
                                 self.data.inst['bl']*self.data.options.array['ratio'] / 2)))

            self.run_socket(s_name='instrument',
                            method='adjust_image_size')

            if lookup_table.split(':')[0] == 'input':
                lt_in = lookup_table_in[nstar]
            else:
                lt_in = None

            # create single input dictionary

            input_dict = {'catalog': self.data.catalog[self.data.catalog.nstar == nstar],
                          'wl_bins': self.data.inst['wl_bins'],
                          'wl_bin_widths': self.data.inst['wl_bin_widths'],
                          'wl_min': self.data.options.array['wl_min'],
                          'integration_time': integration_time,
                          'image_size': self.data.inst['image_size'],
                          'diameter_ap': self.data.options.array['diameter'],
                          'flux_division': self.data.options.array['flux_division'],
                          'throughput': self.data.options.array['throughput']
                                        * self.data.options.array['quantum_eff'],
                          'phase_response': self.data.options.array['phase_response'],
                          'phase_response_chop': self.data.options.array['phase_response_chop'],
                          't_rot': self.data.options.array['t_rot'],
                          'chopping': self.data.options.array['chopping'],
                          'pix_per_wl': self.data.options.array['pix_per_wl'],
                          'col_pos': col_pos,
                          'bl': self.data.inst['bl'],
                          'ratio': self.data.options.array['ratio'],
                          'n_sampling_multiplier': self.data.options.array['n_sampling_multiplier'],
                          'nstar': nstar,
                          'baseline': self.data.inst['bl'],
                          'safe_mode': safe_mode,
                          'd_a_rms': self.data.options.array['d_a_rms'],
                          'd_phi_rms': self.data.options.array['d_phi_rms'],
                          'd_x_rms': self.data.options.array['d_x_rms'],
                          'd_y_rms': self.data.options.array['d_y_rms'],
                          'd_pol_rms': self.data.options.array['d_pol_rms'],
                          'agn_phot_hot': self.data.options.array['agn_phot_hot'],
                          'agn_phot_cold': self.data.options.array['agn_phot_cold'],
                          'agn_phot_white': self.data.options.array['agn_phot_white'],
                          'agn_spacecraft_temp': self.data.options.array['agn_spacecraft_temp'],
                          'rms_mode': self.data.options.array['rms_mode'],
                          'lookup_table': lookup_table.split(':')[0],
                          'lt_in': lt_in,
                          }

            # if safe_mode:
            #     input_dict['noise_catalog'] = self.data.noise_catalog.loc[
            #         self.data.catalog.id[self.data.catalog.nstar == nstar]]

            input_dict_list.append(input_dict)

        self.data.catalog = None

        # if safe_mode:
        #     store = pd.HDFStore(self.data.options.other['output_path']
        #                         + self.data.options.other['output_filename'] + '.hdf5')

        output_dict_list = []

        if self.data.options.other['n_cpu'] == 1:
            print('\nRunning in single processing...')
            for input_dict in tqdm(input_dict_list):
                output_dict_list.append(multiprocessing_runner(input_dict=input_dict))

        else:
            print('\nRunning in multiprocessing...')
            pool = mp.Pool(self.data.options.other['n_cpu'])
            output_dict_list = []
            for result in tqdm(pool.imap_unordered(multiprocessing_runner, input_dict_list),
                               total=len(input_dict_list)):
                output_dict_list.append(result)

        self.data.catalog = pd.concat([output_dict['catalog'] for output_dict in output_dict_list])
        # if safe_mode:
        #     self.data.noise_catalog = pd.concat([output_dict['noise_catalog'] for output_dict in
        #     output_dict_list])

        # if lookup table is in output mode, collect all lookup data and save it to a file
        if lookup_table.split(':')[0] == 'output':
            lookup_table_out = pd.DataFrame.from_dict(
                {output_dict['lookup_table']['nstar']:output_dict['lookup_table']
                                for output_dict in output_dict_list}
            )
            lookup_table_out.to_hdf(lookup_table.split(':')[1],
                                    key='lookup_table', mode='a')

        if safe_mode:
            for output_dict in output_dict_list:
                self.data.noise_catalog.update(output_dict['noise_catalog'])

            self.data.noise_catalog = xr.Dataset(self.data.noise_catalog).to_array()
            self.data.noise_catalog = self.data.noise_catalog.rename({'dim_0': 'wl_bins',
                                                                      'dim_1': 'params',
                                                                      'variable': 'ids'})
            self.data.noise_catalog = self.data.noise_catalog.assign_coords(
                wl_bins=self.data.inst['wl_bins'],
                params=self.data.noise_catalog.coords['params'].values.astype(str),
                ids=self.data.noise_catalog.coords['ids'].values.astype(int)
            )

            self.data.noise_catalog = self.data.noise_catalog.astype(float)
            # self.data.pivot_noise_catalog(to_wavelength=True)

        # if safe_mode:
        #     store.close()


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
                     safe_mode: bool = False,
                     run: bool = True):

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

        col_pos = np.array(((-self.data.inst['bl'] / 2,
                             -self.data.inst['bl'] * self.data.options.array['ratio'] / 2),
                            (-self.data.inst['bl'] / 2,
                             self.data.inst['bl'] * self.data.options.array['ratio'] / 2),
                            (self.data.inst['bl'] / 2,
                             -self.data.inst['bl'] * self.data.options.array['ratio'] / 2),
                            (self.data.inst['bl'] / 2,
                             self.data.inst['bl'] * self.data.options.array['ratio'] / 2)))

        self.run_socket(s_name='instrument',
                        method='adjust_sampling_rate',
                        angsep=angsep)

        print('Sampling rate: ', self.data.inst['n_sampling_rot'])

        self.inst_prt = Instrument(
            # ----- static parameters -----
            wl_bins=self.data.inst['wl_bins'],  # wavelength bins center position in m
            wl_bin_widths=self.data.inst['wl_bin_widths'],  # wavelength bin widhts in m
            integration_time=integration_time,
            image_size=self.data.inst['image_size'],
            # size of image used to simulate exozodi in pix
            diameter_ap=self.data.options.array['diameter'],
            # diameter of the primary mirrors in m
            flux_division=self.data.options.array['flux_division'],
            # division of the flux between the primary mirrors, e.g. in baseline case
            # [0.25, 0.25, 0.25, 0.25]
            throughput=self.data.options.array['throughput']*self.data.options.array['quantum_eff'],
            # fraction of light that is sustained through the optical train
            phase_response=self.data.options.array['phase_response'],
            # phase response of each collector arm in rad
            phase_response_chop=self.data.options.array['phase_response_chop'],
            # phase response of each collector arm in the chopped state in rad
            t_rot=self.data.options.array['t_rot'],  # rotation period of the array in seconds
            chopping=self.data.options.array['chopping'],
            # run calculation with or without chopping, 'chop', 'nchop', 'both'
            pix_per_wl=self.data.options.array['pix_per_wl'],
            # pixels on detector used per wavelength channel
            n_sampling_rot=self.data.inst['n_sampling_rot'],
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
            agnostic_mode=True,  # derive instrumental photon noise from agnostic mode
            eps_cold=0.,  # scaling constant for cold agnostic photon noise spectrum
            eps_hot=0.,  # scaling constant for hot agnostic photon noise spectrum
            eps_white=0.,  # scaling constant white agnostic photon noise spectrum
            agnostic_spacecraft_temp=0.,  # cold-side spacecraft temperature in the
            # agnostic case
            n_sampling_max=self.data.options.other['n_sampling_max'],
            # largest fourier mode used in noise sampling
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
            separation_planet=angsep * distance_s,
            # separation of target planet from host star in AU
        )

        if run:
            self.inst_prt.run()

            if self.data.options.array['chopping'] == 'nchop':
                return self.inst_prt.photon_rates_nchop
            else:
                return self.inst_prt.photon_rates_chop

def multiprocessing_runner(input_dict: dict):
    # TODO: Correct treatment of quantum efficiency
    inst = ils.Instrument(
        # ----- static parameters -----
        wl_bins=input_dict['wl_bins'],  # wavelength bins center position in m
        wl_bin_widths=input_dict['wl_bin_widths'],  # wavelength bin widhts in m
        integration_time=input_dict['integration_time'],
        image_size=input_dict['image_size'],  # size of image used to simulate exozodi in pix
        diameter_ap=input_dict['diameter_ap'],  # diameter of the primary mirrors in m
        flux_division=input_dict['flux_division'],
        # division of the flux between the primary mirrors, e.g. in baseline case
        # [0.25, 0.25, 0.25, 0.25]
        throughput=input_dict['throughput'],
        # fraction of light that is sustained through the optical train
        phase_response=input_dict['phase_response'],  # phase response of each collector arm in rad
        phase_response_chop=input_dict['phase_response_chop'],
        # phase response of each collector arm in the chopped state in rad
        t_rot=input_dict['t_rot'],  # rotation period of the array in seconds
        chopping=input_dict['chopping'],
        # run calculation with or without chopping, 'chop', 'nchop', 'both'
        pix_per_wl=input_dict['pix_per_wl'],  # pixels on detector used per wavelength channel
        n_sampling_rot=360,
        # number of sampling points per array rotation
        detector_dark_current='manual',
        # detector type, 'MIRI' or 'manual'. Specify dark_current_pix in 'manual'
        dark_current_pix=0.,  # detector dark current in electrons s-1 px-1
        detector_thermal='MIRI',  # detector type, 'MIRI'
        det_temp=0.,  # temperature of the detector environment in K
        magnification=15.73,  # tele# scope magnification
        f_number=20.21,  # telescope f-number, i.e. ratio of focal length to aperture size
        secondary_primary_ratio=0.114,  # ratio of secondary to primary mirror sizes
        primary_emissivity=0.,  # emissivity epsilon of the primary mirror
        primary_temp=0.,  # temperature of the primary mirror in K
        pink_noise_co=10000,  # cutoff frequency for the pink noise spectra
        n_cpu=1,  # number of cores used in the simulation
        rms_mode=input_dict['rms_mode'],  # mode for rms values, 'lay', 'static', 'wavelength'
        agnostic_mode=True,  # derive instrumental photon noise from agnostic mode
        eps_cold=input_dict['agn_phot_cold'],
        # scaling constant for cold agnostic photon noise spectrum
        eps_hot=input_dict['agn_phot_hot'],
        # scaling constant for hot agnostic photon noise spectrum
        eps_white=input_dict['agn_phot_white'],
        # scaling constant white agnostic photon noise spectrum
        agnostic_spacecraft_temp=input_dict['agn_spacecraft_temp'],
        # cold-side spacecraft temperature in the agnostic case
        n_sampling_max=10000,  # largest fourier mode used in noise sampling
        d_a_rms=input_dict['d_a_rms'],  # relative amplitude error rms
        d_phi_rms=input_dict['d_phi_rms'],  # phase error rms
        d_pol_rms=input_dict['d_pol_rms'],  # polarization error rms
        d_x_rms=input_dict['d_x_rms'],  # collector position rms, x-direction
        d_y_rms=input_dict['d_y_rms'],  # collector position rms, y-direction
        wl_resolution=200,  # number of wavelength bins simulated for the thermal background
        flux_planet=None,  # substitute flux input in ph m-2 s-1
        simultaneous_chopping=True,
        # ----- parameters change with star -----
        dist_star=input_dict['catalog'].distance_s.iloc[0],  # distance to the target system in pc
        radius_star=input_dict['catalog'].radius_s.iloc[0],  # radius of the star in stellar radii
        temp_star=input_dict['catalog'].temp_s.iloc[0],  # temperature of the host star in Kelvin
        lat_star=input_dict['catalog'].lat.iloc[0],  # ecliptic latitude of the target star
        l_sun=input_dict['catalog'].l_sun.iloc[0],  # stellar luminosity in solar luminosities
        z=input_dict['catalog'].z.iloc[0],
        # zodi level: the exozodi dust is z-times denser than the localzodi dust
        col_pos=input_dict['col_pos'],  # collector position in m
        # ----- parameters change with planet -----
        temp_planet=0.,  # planet temperature in Kelvin
        radius_planet=0.,  # planet radius in earth radii
        separation_planet=0.,  # separation of target planet from host star in AU
    )
    return_dict = {'noise_catalog': {}}

    if input_dict['lookup_table'] != 'input':
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

    if input_dict['lookup_table'] == 'output':
        return_dict['lookup_table'] = {'nstar': input_dict['nstar'],
                                       'A': inst.A,
                                       'wl_bins': inst.wl_bins,
                                       'num_a': inst.num_a,
                                       'rms_mode': inst.rms_mode,
                                       'n_sampling_max': inst.n_sampling_max,
                                       't_rot': inst.t_rot,
                                       't_int': inst.t_int,
                                       'flux_star': inst.flux_star,
                                       'universe': {}}
    elif input_dict['lookup_table'] == 'input':
        inst.flux_star = input_dict['lt_in']['flux_star']
        inst.instrumental_parameters()

    # create mask returning only unique stars
    universes = np.unique(
        input_dict['catalog'].nuniverse[input_dict['catalog'].nstar == input_dict['nstar']],
        return_index=False
    )

    for nuniverse in universes:
        if input_dict['lookup_table'] != 'input':
            inst.z = input_dict['catalog'][np.logical_and(
                input_dict['catalog'].nstar == input_dict['nstar'],
                input_dict['catalog'].nuniverse == nuniverse
            )].z.iloc[0]

            # redo calculation for exozodi
            inst.create_exozodi()
            inst.sensitivity_coefficients(exozodi_only=True)
            inst.fundamental_noise(exozodi_only=True)

        if input_dict['lookup_table'] == 'output':
            return_dict['lookup_table']['universe'][nuniverse] = {'c_a': inst.c_a,
                                                                  'c_phi': inst.c_phi,
                                                                  'c_x': inst.c_x,
                                                                  'c_y': inst.c_y,
                                                                  'c_aa': inst.c_aa,
                                                                  'c_phiphi': inst.c_phiphi,
                                                                  'c_aphi': inst.c_aphi,
                                                                  'c_thetatheta':
                                                                      inst.c_thetatheta}

            copy_params = ['pn_sgl', 'pn_ez', 'pn_lz', 'pn_dc', 'pn_tbd', 'pn_tbpm', 'pn_ag_ht',
                           'pn_ag_cld', 'pn_ag_wht']
            for param in copy_params:
                return_dict['lookup_table']['universe'][nuniverse][param] = (
                    inst.photon_rates_nchop[param]
                )

            return_dict['lookup_table']['universe'][nuniverse]['planet'] = {}

        # load parameters from lookup table
        elif input_dict['lookup_table'] == 'input':
            inst.c_a = input_dict['lt_in']['universe'][nuniverse]['c_a']
            inst.c_phi = input_dict['lt_in']['universe'][nuniverse]['c_phi']
            inst.c_x = input_dict['lt_in']['universe'][nuniverse]['c_x']
            inst.c_y = input_dict['lt_in']['universe'][nuniverse]['c_y']
            inst.c_aa = input_dict['lt_in']['universe'][nuniverse]['c_aa']
            inst.c_phiphi = input_dict['lt_in']['universe'][nuniverse]['c_phiphi']
            inst.c_aphi = input_dict['lt_in']['universe'][nuniverse]['c_aphi']
            inst.c_thetatheta = input_dict['lt_in']['universe'][nuniverse]['c_thetatheta']

            copy_params = ['pn_sgl', 'pn_ez', 'pn_lz', 'pn_dc', 'pn_tbd', 'pn_tbpm', 'pn_ag_ht',
                           'pn_ag_cld', 'pn_ag_wht']
            for param in copy_params:
                inst.photon_rates_nchop[param] = input_dict['lt_in']['universe'][nuniverse][
                    param
                ]

        # go through all planets for the chosen star
        for _, n_p in enumerate(np.argwhere(
                np.logical_and(input_dict['catalog'].nstar.to_numpy() == input_dict['nstar'],
                               input_dict['catalog'].nuniverse.to_numpy() == nuniverse))[:, 0]):

            # ----- must be repeated for every planet -----

            if input_dict['lookup_table'] != 'input':
                # adjust the temporal sampling rate to the baseline and planet separation
                inst.n_sampling_rot = adjust_sampling(
                    angsep=input_dict['catalog']['angsep'].iloc[n_p],
                    baseline=input_dict['bl'],
                    baseline_ratio=input_dict['ratio'],
                    n_sampling_multiplier=input_dict['n_sampling_multiplier'],
                    wl_min=input_dict['wl_min']
                )

                inst.temp_planet = input_dict['catalog']['temp_p'].iloc[n_p]
                inst.radius_planet = input_dict['catalog']['radius_p'].iloc[n_p]
                inst.separation_planet = (input_dict['catalog']['angsep'].iloc[n_p]
                                          * input_dict['catalog']['distance_s'].iloc[n_p])

                # create the planet signal and template function
                inst.create_planet(force=True)
                inst.planet_signal()

            # create lookup table for planets if requested
            if input_dict['lookup_table'] == 'output':
                if inst.chopping == 'chop':
                    return_dict['lookup_table']['universe'][nuniverse]['planet'][
                        input_dict['catalog']['id'].iloc[n_p]
                    ] = {'planet_template_chop': inst.planet_template_chop}
                else:
                    return_dict['lookup_table']['universe'][nuniverse]['planet'][
                        input_dict['catalog']['id'].iloc[n_p]
                    ] = {'planet_template_nchop': inst.planet_template_nchop}

            else:
                # load planet template from lookup table if requested
                if input_dict['lookup_table'] == 'input':
                    if inst.chopping == 'nchop':
                        inst.planet_template_nchop = input_dict['lt_in']['universe'][nuniverse][
                            'planet'
                        ][input_dict['catalog']['id'].iloc[n_p]]['planet_template_nchop']

                    else:
                        inst.planet_template_chop = input_dict['lt_in']['universe'][nuniverse][
                            'planet'
                        ][input_dict['catalog']['id'].iloc[n_p]]['planet_template_chop']

                if (inst.chopping == 'nchop'):
                    inst.sn_nchop()
                else:
                    inst.sn_chop()

                # save baseline
                input_dict['catalog']['baseline'].iat[n_p] = deepcopy(input_dict['baseline'])

                # save sampling rates
                input_dict['catalog']['n_sampling_rot'].iat[n_p] = deepcopy(inst.n_sampling_rot)
                input_dict['catalog']['image_size'].iat[n_p] = deepcopy(inst.image_size)

                # save snr results
                if (inst.chopping == 'nchop'):
                    input_dict['catalog'].t_rot.iat[n_p] = deepcopy(input_dict['integration_time'])
                    input_dict['catalog'].signal.iat[n_p] = inst.photon_rates_nchop['signal'].sum()
                    input_dict['catalog'].photon_noise.iat[n_p] = (
                        np.sqrt((inst.photon_rates_nchop['pn'] ** 2).sum()))
                    input_dict['catalog'].systematic_noise.iat[n_p] = (
                        np.sqrt((inst.photon_rates_nchop['sn'] ** 2).sum()))
                else:
                    input_dict['catalog'].t_rot.iat[n_p] = deepcopy(input_dict['integration_time'])
                    input_dict['catalog'].signal.iat[n_p] = inst.photon_rates_chop['signal'].sum()
                    input_dict['catalog'].photon_noise.iat[n_p] = (
                        np.sqrt((inst.photon_rates_chop['pn'] ** 2).sum()))
                    input_dict['catalog'].systematic_noise.iat[n_p] = (
                        np.sqrt((inst.photon_rates_chop['sn'] ** 2).sum()))

                if input_dict['safe_mode']:
                    if (inst.chopping == 'nchop'):
                        return_dict['noise_catalog'][str(input_dict['catalog'].id.iat[n_p])] = (
                            deepcopy(inst.photon_rates_nchop)
                        )
                    else:
                        return_dict['noise_catalog'][str(input_dict['catalog'].id.iat[n_p])] = (
                            deepcopy(inst.photon_rates_chop)
                        )

    return_dict['catalog'] = input_dict['catalog']
    return_dict['nstar'] = input_dict['nstar']

    return return_dict
