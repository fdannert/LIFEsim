import multiprocessing as mp
from copy import deepcopy
from typing import Union

import numpy as np
from tqdm import tqdm
from spectres import spectres
import pandas as pd
import xarray as xr
import h5py
from joblib import Parallel, delayed, parallel_config
from joblib_progress import joblib_progress

from lifesim.core.modules import InstrumentModule
from inlifesim.observatory import Instrument
#from lifesim.instrument.instrument import Instrument
from lifesim.util.habitable import single_habitable_zone
from lifesim.instrument.instrument import adjust_sampling
from lifesim.core.data import save_to_hdf5, load_from_hdf5, invert_coefficients


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
        self.data.catalog['t_exp'] = np.zeros_like(self.data.catalog.nstar, dtype=float)
        self.data.catalog['signal'] = np.zeros_like(self.data.catalog.nstar, dtype=float)
        self.data.catalog['photon_noise'] = np.zeros_like(self.data.catalog.nstar, dtype=float)
        self.data.catalog['systematic_noise'] = np.zeros_like(self.data.catalog.nstar, dtype=float)
        self.data.catalog['baseline'] = np.zeros_like(self.data.catalog.nstar, dtype=float)
        self.data.catalog['n_sampling_rot'] = np.zeros_like(self.data.catalog.nstar, dtype=int)
        self.data.catalog['image_size'] = np.zeros_like(self.data.catalog.nstar, dtype=int)
        self.data.catalog['fundamental_noise'] = np.zeros_like(self.data.catalog.nstar, dtype=float)
        self.data.catalog['fundamental_snr_1h'] = np.zeros_like(self.data.catalog.nstar, dtype=float)
        self.data.catalog['snr_1h'] = np.zeros_like(self.data.catalog.nstar, dtype=float)
        self.data.catalog['pn_ez'] = np.zeros_like(self.data.catalog.nstar, dtype=float)
        self.data.catalog['pn_lz'] = np.zeros_like(self.data.catalog.nstar, dtype=float)
        self.data.catalog['pn_sgl'] = np.zeros_like(self.data.catalog.nstar, dtype=float)

        if safe_mode:
            self.data.noise_catalog = {}

        # load lookup table
        # if lookup_table.split(':')[0] == 'input':
        #     lookup_table_in = pd.read_hdf(lookup_table.split(':')[1]).to_dict()
        # Load lookup table
        if lookup_table.split(':')[0] == 'input':
            # lookup_table_in = []
            #
            # if lookup_table.split(':')[0] == 'input':
            #     with h5py.File(lookup_table.split(':')[1], "r") as h5file:
            #         for group_name in h5file.keys():
            #             group = h5file[group_name]
            #             lookup_table_in.append({"lookup_table": load_from_hdf5(group)})
            with h5py.File(
                    lookup_table.split(':')[1],
                    'r') as h5file:
                lookup_table_in = load_from_hdf5(h5file)

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
                          'zodi_reference': self.data.options.other['zodi_reference'],
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
                          't_exp': self.data.options.array['t_exp'],
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
                          'd_a_co': self.data.options.array['d_a_co'],
                          'd_phi_co': self.data.options.array['d_phi_co'],
                          'd_x_co': self.data.options.array['d_x_co'],
                          'd_y_co': self.data.options.array['d_y_co'],
                          'd_pol_co': self.data.options.array['d_pol_co'],
                          'd_a_period_bin': self.data.options.array['d_a_period_bin'],
                          'd_phi_period_bin': self.data.options.array['d_phi_period_bin'],
                          'd_x_period_bin': self.data.options.array['d_x_period_bin'],
                          'd_y_period_bin': self.data.options.array['d_y_period_bin'],
                          'd_pol_period_bin': self.data.options.array['d_pol_period_bin'],
                          'agn_phot_hot': self.data.options.array['agn_phot_hot'],
                          'agn_phot_cold': self.data.options.array['agn_phot_cold'],
                          'agn_phot_white': self.data.options.array['agn_phot_white'],
                          'agn_spacecraft_temp': self.data.options.array['agn_spacecraft_temp'],
                          'rms_mode': self.data.options.array['rms_mode'],
                          'hyperrot_noise': self.data.options.array['hyperrot_noise'],
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
            with parallel_config(
                    backend="loky", inner_max_num_threads=1
            ), joblib_progress(
                description="Running stars in parallel ...",
                total=int(len(input_dict_list)),
            ):
                output_dict_list = Parallel(n_jobs=self.data.options.other['n_cpu'],
                                            verbose=10)(
                    delayed(safe_function)(
                        input_dict
                    )
                    for input_dict in input_dict_list
                )
            # pool = mp.Pool(self.data.options.other['n_cpu'])
            # output_dict_list = []
            # for result in tqdm(pool.map(multiprocessing_runner, input_dict_list),
            #                    total=len(input_dict_list)):
            #     output_dict_list.append(result)

        # with parallel_config(
        #         backend="loky", inner_max_num_threads=1
        # ), joblib_progress(
        #     description="Calculating time series ...",
        #     total=int(self.n_draws / self.n_draws_per_run),
        # ):
        #     results = Parallel(n_jobs=self.n_cpu)(
        #         delayed(draw_sample)(
        #             params=params,
        #             return_variables=self.time_samples_return_values,
        #         )
        #         for _ in range(int(self.n_draws / self.n_draws_per_run))
        #     )

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

            # with h5py.File(lookup_table.split(':')[1], "w") as h5file:
            #     for i, output_dict in enumerate(output_dict_list):
            #         lt_out = output_dict['lookup_table']
            #         nstar = lt_out['nstar']
            #
            #
            #         group = h5file.create_group(f"entry_{i}")  # Create a group for each entry
            #         save_to_hdf5(group, output_dict["lookup_table"])  # Save the "lookup_table" contents

            # with h5py.File(
            #         lookup_table.split(':')[1],
            #         'w') as h5file:
            #     save_to_hdf5(output_dict_list, h5file)

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
                     exposure_time: float, # in s
                     n_rot: int,
                     hyperrot_noise: str,
                     pbar: bool = None,
                     baseline_to_planet: bool = False,
                     baseline: float = None,
                     safe_mode: bool = False,
                     run: bool = True,
                     single_bw: bool = False,
                     draw_samples: bool = False,
                     n_draws: int = int(1e2),
                     n_draws_per_run: int = int(1e1),
                     get_single_bracewell: bool = False,
                     wl_bin: Union[np.ndarray, type(None)] = None,
                     wl_bin_width: Union[np.ndarray, type(None)] = None,
                     verbose: bool = True,
                     instrumental_source: Union[str, type(None)] = 'None'
                     ):

        # TODO: Implement baseline_to_planet option

        if wl_bin is not None:
            self.data.inst['wl_bins'] = np.array([wl_bin])
            self.data.inst['wl_bin_widths'] = np.array([wl_bin_width])
            self.data.inst['wl_bin_edges'] = np.array((
                wl_bin - wl_bin_width / 2,
                wl_bin + wl_bin_width / 2
            ))

        # calculate the habitable zone of the specified star
        s_in, s_out, l_sun, \
        hz_in, hz_out, \
            hz_center = single_habitable_zone(
            model=self.data.options.models['habitable'],
            temp_s=temp_s,
            radius_s=radius_s
        )

        flux_planet_spectrum = spectres(
            new_wavs=self.data.inst['wl_bin_edges'],
            spec_wavs=flux_planet_spectrum[0].value,
            spec_fluxes=flux_planet_spectrum[1].value,
            edge_mode=True
        )

        flux_planet_spectrum *= self.data.inst['wl_bin_widths']

        self.run_socket(s_name='instrument',
                        method='apply_options')

        if wl_bin is not None:
            self.data.inst['wl_bins'] = np.array([wl_bin])
            self.data.inst['wl_bin_widths'] = np.array([wl_bin_width])

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
        if not single_bw:
            col_pos = np.array((
                (-self.data.inst['bl'] / 2,
                 -self.data.inst['bl'] * self.data.options.array['ratio'] / 2),
                (-self.data.inst['bl'] / 2,
                 self.data.inst['bl'] * self.data.options.array['ratio'] / 2),
                (self.data.inst['bl'] / 2,
                 -self.data.inst['bl'] * self.data.options.array['ratio'] / 2),
                (self.data.inst['bl'] / 2,
                 self.data.inst['bl'] * self.data.options.array['ratio'] / 2)
            ))
        else:
            col_pos = np.array((
                (-self.data.inst['bl'] / 2,
                 -self.data.inst['bl'] * self.data.options.array['ratio'] / 2),
                (self.data.inst['bl'] / 2,
                 -self.data.inst['bl'] * self.data.options.array['ratio'] / 2)
            ))

        self.run_socket(s_name='instrument',
                        method='adjust_sampling_rate',
                        angsep=angsep)

        if single_bw and (self.data.inst['wl_bins'].shape[0] > 1):
            null = []
            for i in range(self.data.inst['wl_bins'].shape[0]):
                wl_bins = np.array([self.data.inst['wl_bins'][i]])
                wl_bin_widths = np.array([self.data.inst['wl_bin_widths'][i]])
                fp_spec = np.array([flux_planet_spectrum[i]])

                self.inst_prt = Instrument(
                    # ----- static parameters -----
                    wl_bins=wl_bins,
                    # wavelength bins center position in m
                    wl_bin_widths=wl_bin_widths,
                    # wavelength bin widhts in m
                    t_total=integration_time,
                    # total integration time in s
                    t_exp=exposure_time,
                    # time of a single exposure in s
                    n_rot=n_rot,
                    # number of array rotations
                    get_single_bracewell=get_single_bracewell,
                    hyperrot_noise=hyperrot_noise,
                    image_size=self.data.inst['image_size'],
                    # size of image used to simulate exozodi in pix
                    diameter_ap=self.data.options.array['diameter'],
                    # diameter of the primary mirrors in m
                    flux_division=self.data.options.array['flux_division'],
                    # division of the flux between the primary mirrors, e.g. in
                    # baseline case [0.25, 0.25, 0.25, 0.25]
                    throughput=self.data.options.array['throughput']
                               * self.data.options.array['quantum_eff'],
                    # fraction of light that is sustained through the optical train
                    phase_response=self.data.options.array['phase_response'],
                    # phase response of each collector arm in rad
                    phase_response_chop=self.data.options.array['phase_response_chop'],
                    # phase response of each collector arm in the chopped state in rad
                    d_a_co=self.data.options.array['d_a_co'],
                    d_phi_co=self.data.options.array['d_phi_co'],
                    d_pol_co=self.data.options.array['d_pol_co'],
                    d_x_co=self.data.options.array['d_x_co'],
                    d_y_co=self.data.options.array['d_y_co'],
                    n_cpu=1,  # number of cores used in the simulation
                    rms_mode=self.data.options.array['rms_mode'],
                    # mode for rms values, 'lay', 'static', 'wavelength'
                    n_sampling_max=self.data.options.other['n_sampling_max'],
                    # largest fourier mode used in noise sampling
                    d_a_rms=self.data.options.array['d_a_rms'],
                    # relative amplitude error rms
                    d_phi_rms=self.data.options.array['d_phi_rms'],  # phase error rms
                    d_pol_rms=self.data.options.array['d_pol_rms'],
                    # polarization error rms
                    d_x_rms=self.data.options.array['d_x_rms'],
                    # collector position rms, x-direction
                    d_y_rms=self.data.options.array['d_y_rms'],
                    # collector position rms, y-direction
                    simultaneous_chopping=True,
                    draw_samples=draw_samples,
                    n_draws=n_draws,
                    n_draws_per_run=n_draws_per_run,
                    verbose=verbose,
                    # ----- parameters change with star -----
                    dist_star=distance_s,  # distance to the target system in pc
                    radius_star=radius_s,  # radius of the star in stellar radii
                    temp_star=temp_s,  # temperature of the host star in Kelvin
                    lat_star=lat_s,  # ecliptic latitude of the target star
                    l_sun=l_sun,  # stellar luminosity in solar luminosities
                    z=z,
                    # zodi level: the exozodi dust is z-times denser than the
                    # localzodi dust
                    col_pos=col_pos,  # collector position in m
                    # ----- parameters change with planet -----
                    temp_planet=0.,  # planet temperature in Kelvin
                    radius_planet=0.,  # planet radius in earth radii
                    separation_planet=angsep * distance_s,
                    # separation of target planet from host star in AU
                    flux_planet=fp_spec,
                    # substitute flux input in ph m-2 s-1
                    instrumental_source=instrumental_source,
                )
                self.inst_prt.run()
                print(instrumental_source)

                null.append({'pn_timeseries': self.inst_prt.time_samples['pn_timeseries'],
                             'sys_timeseries': self.inst_prt.time_samples['sys_timeseries'],
                             'noise_timeseries_singlebw': self.inst_prt.time_samples['noise_timeseries_singlebw'],
                             'timeseries_singlebw': self.inst_prt.time_samples['timeseries_singlebw'],
                             'planet_signal_nchop': self.inst_prt.planet_signal_nchop,
                             'star_signal': self.inst_prt.flux_star * self.inst_prt.A**2 * self.inst_prt.t_exp,})

            return null

        else:
            self.inst_prt = Instrument(
                # ----- static parameters -----
                wl_bins=self.data.inst['wl_bins'],
                # wavelength bins center position in m
                wl_bin_widths=self.data.inst['wl_bin_widths'],
                # wavelength bin widhts in m
                t_total=integration_time,
                # total integration time in s
                t_exp=exposure_time,
                # time of a single exposure in s
                n_rot=n_rot,
                # number of array rotations
                image_size=self.data.inst['image_size'],
                # size of image used to simulate exozodi in pix
                diameter_ap=self.data.options.array['diameter'],
                # diameter of the primary mirrors in m
                flux_division=self.data.options.array['flux_division'],
                # division of the flux between the primary mirrors, e.g. in
                # baseline case [0.25, 0.25, 0.25, 0.25]
                hyperrot_noise=hyperrot_noise,
                throughput=self.data.options.array['throughput']
                           *self.data.options.array['quantum_eff'],
                # fraction of light that is sustained through the optical train
                phase_response=self.data.options.array['phase_response'],
                # phase response of each collector arm in rad
                phase_response_chop=self.data.options.array['phase_response_chop'],
                # phase response of each collector arm in the chopped state in rad
                d_a_co=self.data.options.array['d_a_co'],
                d_phi_co=self.data.options.array['d_phi_co'],
                d_pol_co=self.data.options.array['d_pol_co'],
                d_x_co=self.data.options.array['d_x_co'],
                d_y_co=self.data.options.array['d_y_co'],
                n_cpu=1,  # number of cores used in the simulation
                rms_mode=self.data.options.array['rms_mode'],
                # mode for rms values, 'lay', 'static', 'wavelength'
                n_sampling_max=self.data.options.other['n_sampling_max'],
                # largest fourier mode used in noise sampling
                d_a_rms=self.data.options.array['d_a_rms'],
                # relative amplitude error rms
                d_phi_rms=self.data.options.array['d_phi_rms'],  # phase error rms
                d_pol_rms=self.data.options.array['d_pol_rms'],
                # polarization error rms
                d_x_rms=self.data.options.array['d_x_rms'],
                # collector position rms, x-direction
                d_y_rms=self.data.options.array['d_y_rms'],
                # collector position rms, y-direction
                simultaneous_chopping=True,
                draw_samples=draw_samples,
                n_draws=n_draws,
                n_draws_per_run=n_draws_per_run,
                verbose=verbose,
                # ----- parameters change with star -----
                dist_star=distance_s,  # distance to the target system in pc
                radius_star=radius_s,  # radius of the star in stellar radii
                temp_star=temp_s,  # temperature of the host star in Kelvin
                lat_star=lat_s,  # ecliptic latitude of the target star
                l_sun=l_sun,  # stellar luminosity in solar luminosities
                z=z,
                # zodi level: the exozodi dust is z-times denser than the
                # localzodi dust
                col_pos=col_pos,  # collector position in m
                # ----- parameters change with planet -----
                temp_planet=0.,  # planet temperature in Kelvin
                radius_planet=0.,  # planet radius in earth radii
                separation_planet=angsep * distance_s,
                # separation of target planet from host star in AU
                flux_planet=flux_planet_spectrum,
                # substitute flux input in ph m-2 s-1
            )

            if run:
                self.inst_prt.run()
                return self.inst_prt.photon_rates_chop

def multiprocessing_runner(input_dict: dict):
    # TODO: Correct treatment of quantum efficiency
    # inst = ils.Instrument(
    #     # ----- static parameters -----
    #     wl_bins=input_dict['wl_bins'],  # wavelength bins center position in m
    #     wl_bin_widths=input_dict['wl_bin_widths'],  # wavelength bin widhts in m
    #     integration_time=input_dict['integration_time'],
    #     image_size=input_dict['image_size'],  # size of image used to simulate exozodi in pix
    #     diameter_ap=input_dict['diameter_ap'],  # diameter of the primary mirrors in m
    #     flux_division=input_dict['flux_division'],
    #     # division of the flux between the primary mirrors, e.g. in baseline case
    #     # [0.25, 0.25, 0.25, 0.25]
    #     throughput=input_dict['throughput'],
    #     # fraction of light that is sustained through the optical train
    #     phase_response=input_dict['phase_response'],  # phase response of each collector arm in rad
    #     phase_response_chop=input_dict['phase_response_chop'],
    #     # phase response of each collector arm in the chopped state in rad
    #     t_rot=input_dict['t_rot'],  # rotation period of the array in seconds
    #     chopping=input_dict['chopping'],
    #     # run calculation with or without chopping, 'chop', 'nchop', 'both'
    #     pix_per_wl=input_dict['pix_per_wl'],  # pixels on detector used per wavelength channel
    #     n_sampling_rot=360,
    #     # number of sampling points per array rotation
    #     detector_dark_current='manual',
    #     # detector type, 'MIRI' or 'manual'. Specify dark_current_pix in 'manual'
    #     dark_current_pix=0.,  # detector dark current in electrons s-1 px-1
    #     detector_thermal='MIRI',  # detector type, 'MIRI'
    #     det_temp=0.,  # temperature of the detector environment in K
    #     magnification=15.73,  # tele# scope magnification
    #     f_number=20.21,  # telescope f-number, i.e. ratio of focal length to aperture size
    #     secondary_primary_ratio=0.114,  # ratio of secondary to primary mirror sizes
    #     primary_emissivity=0.,  # emissivity epsilon of the primary mirror
    #     primary_temp=0.,  # temperature of the primary mirror in K
    #     pink_noise_co=10000,  # cutoff frequency for the pink noise spectra
    #     n_cpu=1,  # number of cores used in the simulation
    #     rms_mode=input_dict['rms_mode'],  # mode for rms values, 'lay', 'static', 'wavelength'
    #     agnostic_mode=True,  # derive instrumental photon noise from agnostic mode
    #     eps_cold=input_dict['agn_phot_cold'],
    #     # scaling constant for cold agnostic photon noise spectrum
    #     eps_hot=input_dict['agn_phot_hot'],
    #     # scaling constant for hot agnostic photon noise spectrum
    #     eps_white=input_dict['agn_phot_white'],
    #     # scaling constant white agnostic photon noise spectrum
    #     agnostic_spacecraft_temp=input_dict['agn_spacecraft_temp'],
    #     # cold-side spacecraft temperature in the agnostic case
    #     n_sampling_max=10000,  # largest fourier mode used in noise sampling
    #     d_a_rms=input_dict['d_a_rms'],  # relative amplitude error rms
    #     d_phi_rms=input_dict['d_phi_rms'],  # phase error rms
    #     d_pol_rms=input_dict['d_pol_rms'],  # polarization error rms
    #     d_x_rms=input_dict['d_x_rms'],  # collector position rms, x-direction
    #     d_y_rms=input_dict['d_y_rms'],  # collector position rms, y-direction
    #     wl_resolution=200,  # number of wavelength bins simulated for the thermal background
    #     flux_planet=None,  # substitute flux input in ph m-2 s-1
    #     simultaneous_chopping=True,
    #     # ----- parameters change with star -----
    #     dist_star=input_dict['catalog'].distance_s.iloc[0],  # distance to the target system in pc
    #     radius_star=input_dict['catalog'].radius_s.iloc[0],  # radius of the star in stellar radii
    #     temp_star=input_dict['catalog'].temp_s.iloc[0],  # temperature of the host star in Kelvin
    #     lat_star=input_dict['catalog'].lat.iloc[0],  # ecliptic latitude of the target star
    #     l_sun=input_dict['catalog'].l_sun.iloc[0],  # stellar luminosity in solar luminosities
    #     z=input_dict['catalog'].z.iloc[0],
    #     # zodi level: the exozodi dust is z-times denser than the localzodi dust
    #     col_pos=input_dict['col_pos'],  # collector position in m
    #     # ----- parameters change with planet -----
    #     temp_planet=0.,  # planet temperature in Kelvin
    #     radius_planet=0.,  # planet radius in earth radii
    #     separation_planet=0.,  # separation of target planet from host star in AU
    # )
    inst = Instrument(
        wl_bins=input_dict['wl_bins'],  # wavelength bins center position in m
        wl_bin_widths=input_dict['wl_bin_widths'],  # wavelength bin widths in m
        image_size=input_dict['image_size'],  # size of the image used to simulate exozodi in pix
        diameter_ap=input_dict['diameter_ap'],  # diameter of the primary mirrors in m
        flux_division=input_dict['flux_division'],  # division of the flux between primary mirrors
        throughput=input_dict['throughput'],  # fraction of light sustained through the optical train
        dist_star=input_dict['catalog'].distance_s.iloc[0],  # distance to the target system in pc
        radius_star=input_dict['catalog'].radius_s.iloc[0],  # radius of the star in stellar radii
        temp_star=input_dict['catalog'].temp_s.iloc[0],  # temperature of the host star in Kelvin
        lat_star=input_dict['catalog'].lat.iloc[0],  # ecliptic latitude of the target star
        l_sun=input_dict['catalog'].l_sun.iloc[0],  # stellar luminosity in solar luminosities
        z=input_dict['zodi_reference'],  # zodi level
        temp_planet=0.,  # planet temperature in Kelvin
        radius_planet=0.,  # planet radius in Earth radii
        separation_planet=0.,  # separation of target planet from host star in AU
        col_pos=input_dict['col_pos'],  # collector position in m
        phase_response=input_dict['phase_response'],  # phase response of each collector arm in rad
        phase_response_chop=input_dict['phase_response_chop'],  # phase response in the chopped state in rad
        n_rot=1,  # NEW: total number of rotations over the observation time
        t_total=input_dict['t_rot'],  # NEW: total observation time in seconds
        t_exp=input_dict['t_exp'],  # NEW: exposure time per sampling in seconds
        n_cpu=1,  # number of cores used in the simulation
        rms_mode=input_dict['rms_mode'],  # mode for RMS values: 'lay', 'static', or 'wavelength'
        hyperrot_noise=input_dict['hyperrot_noise'],  # NEW: hyperrotation noise source, e.g., "pink" or None
        d_a_rms=input_dict['d_a_rms'],  # relative amplitude error RMS
        d_phi_rms=input_dict['d_phi_rms'],  # phase error RMS
        d_pol_rms=input_dict['d_pol_rms'],  # polarization error RMS
        d_x_rms=input_dict['d_x_rms'],  # collector position RMS, x-direction
        d_y_rms=input_dict['d_y_rms'],  # collector position RMS, y-direction
        d_a_co=input_dict['d_a_co'],  # NEW: amplitude error cutoff frequency
        d_phi_co=input_dict['d_phi_co'],  # NEW: phase error cutoff frequency
        d_pol_co=input_dict['d_pol_co'],  # NEW: polarization error cutoff frequency
        d_x_co=input_dict['d_x_co'],  # NEW: position error cutoff frequency, x-direction
        d_y_co=input_dict['d_y_co'],  # NEW: position error cutoff frequency, y-direction
        d_a_period_bin=input_dict['d_a_period_bin'],  # NEW: amplitude periodic error in binning mode
        d_phi_period_bin=input_dict['d_phi_period_bin'],  # NEW: phase periodic error in binning mode
        d_pol_period_bin=input_dict['d_pol_period_bin'],  # NEW: polarization periodic error in binning mode
        d_x_period_bin=input_dict['d_x_period_bin'],  # NEW: position periodic error in binning mode, x-direction
        d_y_period_bin=input_dict['d_y_period_bin'],  # NEW: position periodic error in binning mode, y-direction
        simultaneous_chopping=True)
    return_dict = {'noise_catalog': {}}

    if input_dict['lookup_table'] != 'input':
        # ----- same for every star -----
        inst.run(run_method=['star'])

        b_ez = deepcopy(inst.b_ez)
        # inst.instrumental_parameters()
        # inst.create_star()
        # inst.create_localzodi()
        # inst.create_exozodi()
        # inst.sensitivity_coefficients()
        # inst.fundamental_noise()
        #
        # if inst.agnostic_mode:
        #     inst.pn_agnostic()
        # else:
        #     inst.pn_dark_current()
        #     inst.pn_thermal_background_detector()
        #     inst.pn_thermal_primary_mirror()

    if input_dict['lookup_table'] == 'output':
        return_dict['lookup_table'] = {'nstar': int(input_dict['nstar']),
                                       'A': inst.A,
                                       'wl_bins': inst.wl_bins,
                                       'num_a': inst.num_a,
                                       'rms_mode': inst.rms_mode,
                                       'n_sampling_total': inst.n_sampling_total,
                                       't_total': inst.t_total,
                                       'n_rot': inst.n_rot,
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

            inst.b_ez = b_ez * inst.z / input_dict['zodi_reference']

            # redo calculation for exozodi
            inst.run(run_method=['exozodi'])
            # inst.create_exozodi()
            # inst.sensitivity_coefficients(exozodi_only=True)
            # inst.fundamental_noise(exozodi_only=True)

        if input_dict['lookup_table'] == 'output':
            return_dict['lookup_table']['universe'][nuniverse] = {
                'grad_n_coeff': invert_coefficients(inst.grad_n_coeff),
                'hess_n_coeff': invert_coefficients(inst.hess_n_coeff),
                'grad_n_coeff_chop': invert_coefficients(inst.grad_n_coeff_chop),
                'hess_n_coeff_chop': invert_coefficients(inst.hess_n_coeff_chop)
            }

            copy_params = ['pn_sgl', 'pn_ez', 'pn_lz']
            for param in copy_params:
                return_dict['lookup_table']['universe'][nuniverse][param] = (
                    inst.photon_rates_nchop[param].to_numpy()
                )

            return_dict['lookup_table']['universe'][nuniverse]['planet'] = {}

        # load parameters from lookup table
        elif input_dict['lookup_table'] == 'input':
            inst.grad_n_coeff = invert_coefficients(input_dict['lt_in']['universe'][nuniverse]['grad_n_coeff'])
            inst.hess_n_coeff = invert_coefficients(input_dict['lt_in']['universe'][nuniverse]['hess_n_coeff'])
            inst.grad_n_coeff_chop = invert_coefficients(input_dict['lt_in']['universe'][nuniverse]['grad_n_coeff_chop'])
            inst.hess_n_coeff_chop = invert_coefficients(input_dict['lt_in']['universe'][nuniverse]['hess_n_coeff_chop'])

            copy_params = ['pn_sgl', 'pn_ez', 'pn_lz']
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
                n_sampling_rot = adjust_sampling(
                    angsep=input_dict['catalog']['angsep'].iloc[n_p],
                    baseline=input_dict['bl'],
                    baseline_ratio=input_dict['ratio'],
                    n_sampling_multiplier=input_dict['n_sampling_multiplier'],
                    wl_min=input_dict['wl_min']
                )

                if n_sampling_rot % 2 == 0:
                    n_sampling_rot += 1

                inst.t_exp = inst.t_rot / n_sampling_rot
                inst.n_sampling_total = int(np.round(inst.t_total / inst.t_exp))
                inst.n_sampling_rot = int(np.round(inst.t_rot / inst.t_exp))


                inst.temp_planet = input_dict['catalog']['temp_p'].iloc[n_p]
                inst.radius_planet = input_dict['catalog']['radius_p'].iloc[n_p]
                inst.separation_planet = (input_dict['catalog']['angsep'].iloc[n_p]
                                          * input_dict['catalog']['distance_s'].iloc[n_p])

                # create the planet signal and template function

                inst.run(run_method=['planet'])
                # inst.create_planet(force=True)
                # inst.planet_signal()

            # create lookup table for planets if requested
            if input_dict['lookup_table'] == 'output':
                return_dict['lookup_table']['universe'][nuniverse]['planet'][
                    input_dict['catalog']['id'].iloc[n_p]
                ] = {'planet_template_chop': inst.planet_template_chop,
                     't_exp': inst.t_exp,
                     'n_sampling_total': inst.n_sampling_total,
                     'n_sampling_rot': inst.n_sampling_rot,
                     'signal_nchop': inst.photon_rates_nchop['signal'].to_numpy(),
                     'signal_chop': inst.photon_rates_chop['signal'].to_numpy(),}

            else:
                # load planet template from lookup table if requested
                if input_dict['lookup_table'] == 'input':
                    inst.planet_template_chop = input_dict['lt_in']['universe'][nuniverse][
                        'planet'
                    ][input_dict['catalog']['id'].iloc[n_p]]['planet_template_chop']

                    inst.t_exp = input_dict['lt_in']['universe'][nuniverse][
                        'planet'
                    ][input_dict['catalog']['id'].iloc[n_p]]['t_exp']

                    inst.n_sampling_total = input_dict['lt_in']['universe'][nuniverse][
                        'planet'
                    ][input_dict['catalog']['id'].iloc[n_p]]['n_sampling_total']

                    inst.n_sampling_rot = input_dict['lt_in']['universe'][nuniverse][
                        'planet'
                    ][input_dict['catalog']['id'].iloc[n_p]]['n_sampling_rot']

                    inst.photon_rates_nchop['signal'] = input_dict['lt_in']['universe'][nuniverse][
                        'planet'
                    ][input_dict['catalog']['id'].iloc[n_p]]['signal_nchop']

                    inst.photon_rates_chop['signal'] = input_dict['lt_in']['universe'][nuniverse][
                        'planet'
                    ][input_dict['catalog']['id'].iloc[n_p]]['signal_chop']


                # if (inst.chopping == 'nchop'):
                #     inst.sn_nchop()
                # else:
                #     inst.sn_chop()

                inst.run(run_method=['systematic'])

                # save baseline
                input_dict['catalog']['baseline'].iat[n_p] = deepcopy(input_dict['baseline'])

                # save sampling rates
                input_dict['catalog']['n_sampling_rot'].iat[n_p] = deepcopy(inst.n_sampling_rot)
                input_dict['catalog']['image_size'].iat[n_p] = deepcopy(inst.image_size)


                input_dict['catalog'].t_rot.iat[n_p] = deepcopy(input_dict['integration_time'])
                input_dict['catalog'].t_exp.iat[n_p] = deepcopy(inst.t_exp)
                input_dict['catalog'].signal.iat[n_p] = inst.photon_rates_chop['signal'].sum()
                input_dict['catalog'].photon_noise.iat[n_p] = (
                    np.sqrt((inst.photon_rates_chop['pn'] ** 2).sum()))
                input_dict['catalog'].systematic_noise.iat[n_p] = (
                    np.sqrt((inst.photon_rates_chop['sn'] ** 2).sum()))

                input_dict['catalog'].fundamental_snr_1h.iat[n_p] = np.sqrt(
                    np.sum(
                        (inst.photon_rates_chop['signal'] / inst.photon_rates_chop['fundamental'])**2
                    )
                ) * np.sqrt(60 * 60 / input_dict['t_rot'])

                input_dict['catalog'].pn_ez.iat[n_p] = np.sqrt(np.sum(inst.photon_rates_chop['pn_ez'] ** 2))
                input_dict['catalog'].pn_lz.iat[n_p] = np.sqrt(np.sum(inst.photon_rates_chop['pn_lz'] ** 2))
                input_dict['catalog'].pn_sgl.iat[n_p] = np.sqrt(np.sum(inst.photon_rates_chop['pn_sgl'] ** 2))

                input_dict['catalog'].snr_1h.at[n_p] = np.sqrt(
                    np.sum(
                        (inst.photon_rates_chop['signal'] / inst.photon_rates_chop['noise'])**2
                    )
                ) * np.sqrt(60 * 60 / input_dict['t_rot'])



                if input_dict['safe_mode']:
                    # if (inst.chopping == 'nchop'):
                    #     return_dict['noise_catalog'][str(input_dict['catalog'].id.iat[n_p])] = (
                    #         deepcopy(inst.photon_rates_nchop)
                    #     )
                    # else:
                    return_dict['noise_catalog'][str(input_dict['catalog'].id.iat[n_p])] = (
                        deepcopy(inst.photon_rates_chop)
                    )

    return_dict['catalog'] = input_dict['catalog']
    return_dict['nstar'] = input_dict['nstar']

    return return_dict

def safe_function(arg):
    try:
        result = multiprocessing_runner(arg)
        return result
    except Exception as e:
        print(f"Worker failed with input {arg['nstar']} and error: {e}")
        return None

