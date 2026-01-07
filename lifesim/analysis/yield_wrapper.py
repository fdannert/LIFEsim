import time
import shutil
import os
import contextlib
from copy import deepcopy
import json
from typing import Union
import logging

import numpy as np
from joblib import Parallel, delayed, parallel_config
from joblib_progress import joblib_progress
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt

import lifesim

class ScienceYield:
    def __init__(self,
                 config_path,
                 catalog_path,
                 output_path,
                 n_cpu: int = 1,
                 cat_from_ppop: bool = True):
        self.config_path = config_path
        self.catalog_path = catalog_path
        self.output_path = output_path
        self.n_cpu = n_cpu
        self.cat_from_ppop = cat_from_ppop

        # Create logger that writes to `yield_wrapper_log.txt` in append mode.
        os.makedirs(self.output_path, exist_ok=True)
        log_file_path = os.path.join(self.output_path, 'yield_wrapper.log')

        self.logger = logging.getLogger('lifesim.yield_wrapper')
        self.logger.setLevel(logging.INFO)

        # Avoid duplicating handlers for the same file when multiple instances are created
        file_path_abs = os.path.abspath(log_file_path)
        existing_handler = None
        for h in list(self.logger.handlers):
            if isinstance(h, logging.FileHandler) and getattr(h, 'baseFilename', None) == file_path_abs:
                existing_handler = h
                break

        if existing_handler is None:
            fh = logging.FileHandler(log_file_path, mode='a')
            fh.setLevel(logging.INFO)
            fh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
            self.logger.addHandler(fh)
            self.logger.propagate = False

        # Log instance creation and parameters
        self.logger.info('%s', '=' * 80)
        self.logger.info('ScienceYield instance created at: %s', time.ctime())
        self.logger.info('Parameters:')
        self.logger.info('  config_path: %s', self.config_path)
        self.logger.info('  catalog_path: %s', self.catalog_path)
        self.logger.info('  output_path: %s', self.output_path)
        self.logger.info('  n_cpu: %s', self.n_cpu)

    def _compute_snrs(self,
                     output_path,
                     output_filename,
                     run_maxsep,
                     diameter=None):

        print('START OF RUN: ', time.ctime())
        print('RUN NAME: ', output_filename)
        print('----------------------------')

        t = time.time()
        # create bus
        bus = lifesim.Bus()

        # setting the options
        bus.build_from_config(filename=self.config_path)
        bus.data.options.set_manual(n_cpu=self.n_cpu) # speed up calculation

        if diameter is not None:
            bus.data.options.set_manual(diameter=diameter)

        bus.data.options.set_manual(
            output_path=output_path)
        bus.data.options.set_manual(output_filename=output_filename)


        # ---------- Loading the Catalog ----------
        if self.cat_from_ppop:
            bus.data.catalog_from_ppop(input_path=self.catalog_path)
        else:
            bus.data.import_catalog(input_path=self.catalog_path)

        # ---------- Creating the Instrument ----------

        # create modules and add to bus
        instrument = lifesim.Instrument(name='inst')
        bus.add_module(instrument)

        transm = lifesim.TransmissionMap(name='transm')
        bus.add_module(transm)

        exo = lifesim.PhotonNoiseExozodi(name='exo')
        bus.add_module(exo)
        local = lifesim.PhotonNoiseLocalzodi(name='local')
        bus.add_module(local)
        star = lifesim.PhotonNoiseStar(name='star')
        bus.add_module(star)

        # connect all modules
        bus.connect(('inst', 'transm'))
        bus.connect(('inst', 'exo'))
        bus.connect(('inst', 'local'))
        bus.connect(('inst', 'star'))

        bus.connect(('star', 'transm'))

        if run_maxsep:
            # only run maxsep SNR for planets that are part of an experiment
            bus.data.catalog['is_interesting'] = False
            for exp in bus.data.options.optimization['experiments'].keys():
                mask_exp = ((bus.data.catalog.radius_p
                             >= bus.data.options.optimization['experiments'][exp]['radius_p_min'])
                            & (bus.data.catalog.radius_p
                               <= bus.data.options.optimization['experiments'][exp]['radius_p_max'])
                            & (bus.data.catalog.temp_s
                               >= bus.data.options.optimization['experiments'][exp]['temp_s_min'])
                            & (bus.data.catalog.temp_s
                               <= bus.data.options.optimization['experiments'][exp]['temp_s_max']))

                if bus.data.options.optimization['experiments'][exp]['in_HZ']:
                    mask_exp = (mask_exp
                                & (bus.data.catalog['habitable']))

                bus.data.catalog['exp_' + exp] = mask_exp

                bus.data.catalog['is_interesting'] = np.logical_or(mask_exp, bus.data.catalog['is_interesting'])
            
            bus.data.catalog = bus.data.catalog[bus.data.catalog.is_interesting]
            bus.data.catalog['angsep'] = bus.data.catalog['maxangsep']

        instrument.get_snr()
        bus.save()

        del bus

        print('Generation took ', (time.time() - t) / 60, ' minutes to complete.')

    def run_aperture_sweep_snr(self,
                               mirror_diameters,
                               run_name):

        final_output_path = os.path.join(self.output_path, run_name)
        if not os.path.exists(final_output_path):
            os.makedirs(final_output_path)
        else:
            raise ValueError('Directory already exists: ' + final_output_path)

        for ndim, diameter in enumerate(mirror_diameters):
            print('')
            print('')
            print('STARTING RUN FOR MIRROR DIAMETER: ', diameter)
            print('AT TIME: ', time.ctime())

            print('Preparing directories... ', end='')

            output_directory = os.path.join(final_output_path, 'diam_' + str(np.round(diameter, 2)).replace('.', '_'))
            if not os.path.exists(output_directory):
                os.makedirs(output_directory)

            print('[Done]')

            print('Commencing base run... ')
            self._compute_snrs(output_path=f'{output_directory}/',
                              output_filename='sweep_diam_' + str(np.round(diameter, 2)).replace('.', '_'),
                              run_maxsep=False,
                              diameter=float(diameter))
            print('[Done]')

            print('Commencing maxsep run... ')
            self._compute_snrs(output_path=f'{output_directory}/',
                              output_filename='sweep_diam_' + str(np.round(diameter, 2)).replace('.', '_') + '_maxsep',
                              run_maxsep=True,
                              diameter=float(diameter))
            print('[Done]')

        # Append single multi-line log summary for this method
        t_end = time.time()
        try:
            diam_dirs = sorted([d for d in os.listdir(final_output_path)
                                if os.path.isdir(os.path.join(final_output_path, d)) and d.startswith('diam_')])
            diameters_created = [d.replace('diam_', '').replace('_', '.') for d in diam_dirs]
        except Exception:
            diam_dirs = []
            diameters_created = []
        msg = f"""run_aperture_sweep_snr summary:
                  run_name: {run_name}
                  requested_mirror_diameters: {mirror_diameters}
                  mirror_directories_found: {diam_dirs}
                  mirror_diameters_reported: {diameters_created}
                  final_output_path: {final_output_path}
                  start_time: {time.ctime(t_end - (t_end - t_end))}  # placeholder, exact start not stored in this scope
                  end_time: {time.ctime(t_end)}
                  elapsed_seconds: {round(0.0, 2)}  # elapsed not measured here to avoid changing existing code flow
                """
        self.logger.info(msg)

    def run_optimizer_sweep(self,
                            run_name,
                            source_name,
                            characterization: bool = False,
                            opt_limit_factor: Union[None, float] = None):
        source_path = os.path.join(self.output_path, source_name)
        # get the names of all subdirectories in output_path (which contain subdirectories for different mirror diameters)
        subdirs = [d for d in os.listdir(source_path) if os.path.isdir(os.path.join(source_path, d))]

        # check a directory with the name run_name already exists in output_path, otherwise create it
        final_output_path = os.path.join(self.output_path, run_name)
        if not os.path.exists(final_output_path):
            os.makedirs(final_output_path)
        else:
            raise ValueError('Directory already exists: ' + final_output_path)

        run_configs = []

        for subdir in subdirs:
            # make a list of all files ending in .hdf5 in subdir, then keep only the part of the filename before '_catalog.hdf5' and only if the sting does not contain 'maxsep'
            catalog_files = [f.split('_catalog.hdf5')[0] for f in os.listdir(os.path.join(source_path, subdir))
                                if f.endswith('_catalog.hdf5') and 'maxsep' not in f]

            # create a subdir of the same name in final_output_path, no existence check required
            output_directory = os.path.join(final_output_path, subdir)
            os.makedirs(output_directory)

            for catalog_file in catalog_files:
                if self.n_cpu == 1:
                    compute_yields_mp(
                        output_filename=catalog_file,
                        output_path=output_directory + '/',
                        catalog_path=os.path.join(source_path, subdir, catalog_file + '_catalog.hdf5'),
                        config_path=self.config_path,
                        characterization=characterization,
                        opt_limit_factor=opt_limit_factor
                    )
                else:
                    run_configs.append({'output_filename': catalog_file,
                                        'output_path': output_directory + '/',
                                        'catalog_path': os.path.join(source_path, subdir, catalog_file + '_catalog.hdf5'),
                                        })

        if self.n_cpu > 1:
            with parallel_config(
                    backend="loky", inner_max_num_threads=1
            ), joblib_progress(
                description="Running SNR calculation in parallel ...",
                total=len(run_configs),
            ):
                Parallel(n_jobs=self.n_cpu)(
                    delayed(compute_yields_mp)(
                        output_filename=rc['output_filename'],
                        output_path=rc['output_path'],
                        catalog_path=rc['catalog_path'],
                        config_path=self.config_path,
                        characterization=characterization,
                        opt_limit_factor=opt_limit_factor
                    )
                    for rc in run_configs)

        # Append single multi-line log summary for this method
        t_end = time.time()
        # derive counts from filesystem to avoid modifying existing logic
        catalog_counts = 0
        catalog_list = {}
        for d in subdirs:
            files = [f for f in os.listdir(os.path.join(source_path, d)) if f.endswith('_catalog.hdf5') and 'maxsep' not in f]
            catalog_counts += len(files)
            catalog_list[d] = files
        msg = f"""run_optimizer_sweep summary:
                  run_name: {run_name}
                  source_name: {source_name}
                  source_path: {source_path}
                  n_cpu: {self.n_cpu}
                  subdirs_found: {subdirs}
                  total_subdirs_count: {len(subdirs)}
                  total_catalog_files_count: {catalog_counts}
                  catalog_files_by_subdir: {catalog_list}
                  run_configs_queued_parallel: {len(run_configs)}
                  final_output_path: {final_output_path}
                  characterization: {characterization}
                  opt_limit_factor: {opt_limit_factor}
                  end_time: {time.ctime(t_end)}
                """
        self.logger.info(msg)

    def combine_catalog_maxsep(self,
                               source_name,):
        source_path = os.path.join(self.output_path, source_name)

        print('START OF combine_catalog_maxsep: ', time.ctime())
        print('SOURCE NAME: ', source_name)
        print('SOURCE PATH: ', source_path)
        t_all = time.time()

        subdirs = [d for d in os.listdir(source_path) if os.path.isdir(os.path.join(source_path, d))]

        for subdir in subdirs:
            print('--------------------------------------------------')
            print('Processing subdir: ', subdir)
            t_sub = time.time()

            # make a list of all files ending in .hdf5 in subdir, then keep only the part of the filename before '_catalog.hdf5' and only if the sting does not contain 'maxsep'
            base_catalog_files = [f.split('_catalog.hdf5')[0] for f in os.listdir(os.path.join(source_path, subdir))
                                if f.endswith('_catalog.hdf5') and 'maxsep' not in f]

            if len(base_catalog_files) == 0:
                print('  No base catalog files found in', os.path.join(source_path, subdir))
                print('  Skipping subdir.')
                print('Subdir processing took', round(time.time() - t_sub, 2), 's')
                continue

            for catalog_file in base_catalog_files:
                print('  ----------------------------------------------')
                print('  Processing catalog:', catalog_file)
                t_cat = time.time()

                base_path = os.path.join(source_path, subdir, catalog_file + '_catalog.hdf5')
                maxsep_path = os.path.join(source_path, subdir, catalog_file + '_maxsep_catalog.hdf5')

                print('    Reading base catalog from:', base_path)
                base_catalog = pd.read_hdf(base_path)
                print('    Base catalog entries:', len(base_catalog))

                print('    Reading maxsep catalog from:', maxsep_path)
                maxsep_catalog = pd.read_hdf(maxsep_path)
                print('    Maxsep catalog entries:', len(maxsep_catalog))

                # map and fill missing values
                base_catalog['maxsep_snr_1h'] = base_catalog['id'].map(maxsep_catalog.set_index('id')['snr_1h'])
                missing_before = int(base_catalog['maxsep_snr_1h'].isna().sum())
                base_catalog['maxsep_snr_1h'].fillna(0, inplace=True)
                print(f'    Mapped maxsep snr_1h, filled {missing_before} missing values with 0')

                print('    Saving combined catalog to:', base_path)
                base_catalog.to_hdf(base_path, key='catalog', mode='w')

                print('  Done processing', catalog_file, '- took', round(time.time() - t_cat, 2), 's')

            print('Finished subdir:', subdir, '- took', round(time.time() - t_sub, 2), 's')

        print('ALL combine_catalog_maxsep finished. Total time:', round(time.time() - t_all, 2), 's')

        # Append single multi-line log summary for this method
        t_end = time.time()
        # derive file counts from filesystem (no change to processing logic)
        subdirs_list = [d for d in os.listdir(source_path) if os.path.isdir(os.path.join(source_path, d))]
        total_base_catalogs = sum(len([f for f in os.listdir(os.path.join(source_path, d))
                                       if f.endswith('_catalog.hdf5') and 'maxsep' not in f]) for d in subdirs_list)
        total_maxsep_catalogs = sum(len([f for f in os.listdir(os.path.join(source_path, d))
                                         if f.endswith('_maxsep_catalog.hdf5')]) for d in subdirs_list)
        msg = f"""combine_catalog_maxsep summary:
                  source_name: {source_name}
                  source_path: {source_path}
                  subdirs_found: {subdirs_list}
                  subdirs_total: {len(subdirs_list)}
                  total_base_catalog_files: {total_base_catalogs}
                  total_maxsep_catalog_files: {total_maxsep_catalogs}
                  start_time: {time.ctime(t_all)}
                  end_time: {time.ctime(t_end)}
                  elapsed_seconds: {round(t_end - t_all, 2)}
                """
        self.logger.info(msg)

    def get_mission_time(self,
                         catalog_path,
                         config_path):
        bus = lifesim.Bus()

        # loading the config and catalog, make sure that the catalog is already combined with maxsep SNRs
        bus.build_from_config(filename=config_path)
        bus.data.import_catalog(input_path=catalog_path)

        # collect the experiments
        exps = [col[4:] for col in bus.data.catalog.columns if col.startswith('exp_')]
        
        # recalculate the interesting flag
        bus.data.catalog['is_interesting'] = False
        for exp in exps:
            mask_exp = ((bus.data.catalog.radius_p
                         >= bus.data.options.optimization['experiments'][exp]['radius_p_min'])
                        & (bus.data.catalog.radius_p
                           <= bus.data.options.optimization['experiments'][exp]['radius_p_max'])
                        & (bus.data.catalog.temp_s
                           >= bus.data.options.optimization['experiments'][exp]['temp_s_min'])
                        & (bus.data.catalog.temp_s
                           <= bus.data.options.optimization['experiments'][exp]['temp_s_max']))

            if bus.data.options.optimization['experiments'][exp]['in_HZ']:
                mask_exp = (mask_exp
                            & (bus.data.catalog['habitable']))

            bus.data.catalog['exp_' + exp] = mask_exp

            bus.data.catalog['is_interesting'] = np.logical_or(mask_exp, bus.data.catalog['is_interesting'])

        cat_det = bus.data.catalog[np.logical_and(
            bus.data.catalog.is_interesting,
            bus.data.catalog.detected
        )].sort_values('t_detected')

        cat_det['follow_up'] = False
        cat_det['t_orbit'] = 0.
        cat_det['t_char'] = 0.

        bus.data.catalog['follow_up'] = False

        # set up the time sheet that records the mission time per universe and per experiment
        columns = ['detection', 'orbit', 'characterization', 'total']
        for exp in exps:
            columns.append('n_' + exp)
        time_sheet = pd.DataFrame(index=np.unique(cat_det.nuniverse),
                                           columns=columns)

        # 1. get total time for detection campaign from every universe
        t_det = bus.data.catalog.t_detected.max()
        time_sheet['detection'] = t_det

        # 2. identify the follow_up targets for each universe
        cat_det.sort_values('maxsep_snr_1h', ascending=False, inplace=True)
        for nu in np.unique(cat_det.nuniverse):
            for exp in exps:
                mask_det = np.logical_and.reduce((cat_det.nuniverse == nu,
                                                  cat_det['exp_' + exp],
                                                  cat_det.t_detected <= t_det))
                # only set the top N targets to follow up where N is the sample size for the experiment
                mask_det_indices = cat_det[mask_det].index[:bus.data.options.optimization['experiments'][exp]['sample_size']]
                cat_det.loc[mask_det_indices, 'follow_up'] = True
                time_sheet.loc[nu, 'n_' + exp] = len(mask_det_indices)

        # 3. calculate follow-up time for orbit and characterization
        cat_det.loc[cat_det.follow_up, 't_orbit'] = (
                (((bus.data.options.optimization['snr_target'] / cat_det[cat_det.follow_up].maxsep_snr_1h) ** 2)
                 * 60 * 60
                 + bus.data.options.array['t_slew'])
                * (bus.data.options.optimization['n_orbits'] - 1)
        )

        cat_det.loc[cat_det.follow_up, 't_char'] = (
                (((bus.data.options.optimization['snr_char'] / cat_det[cat_det.follow_up].maxsep_snr_1h) ** 2)
                 * 60 * 60
                 + bus.data.options.array['t_slew'])
        )

        for nu in np.unique(cat_det.nuniverse):
            mask_followup = np.logical_and.reduce(
                (cat_det.nuniverse == nu, cat_det.follow_up,
                 ))
            time_sheet.loc[nu, 'orbit'] = cat_det.loc[mask_followup, 't_orbit'].sum()
            time_sheet.loc[nu, 'characterization'] = cat_det.loc[mask_followup, 't_char'].sum()

        # 4. get total time
        time_sheet['total'] = time_sheet['detection'] + time_sheet['orbit'] + time_sheet['characterization']

        # 5. copy to original catalog
        cols = ['follow_up', 't_orbit', 't_char']
        mapping_df = cat_det.set_index('id')[cols]

        for col, fill_value, out_type in [
            ('follow_up', False, bool),
            ('t_orbit', 0.0, float),
            ('t_char', 0.0, float),
        ]:
            # map and infer object dtypes first
            s = bus.data.catalog['id'].map(mapping_df[col]).infer_objects(copy=False)
            # replace missing values without using .fillna
            s_filled = s.where(s.notna(), other=fill_value)
            # assign with desired type
            bus.data.catalog[col] = s_filled.astype(out_type)

        bus.save()

        return time_sheet

    def sweep_mission_time(self,
                           source_name):
        source_path = os.path.join(self.output_path, source_name)
        subdirs = [d for d in os.listdir(source_path) if os.path.isdir(os.path.join(source_path, d))]

        for subdir in subdirs:
            print('--------------------------------------------------')
            print('Processing subdir: ', subdir)
            t_sub = time.time()

            # make a list of all files ending in .hdf5 in subdir, then keep only the part of the filename before '_catalog.hdf5' and only if the sting does not contain 'maxsep'
            run_names = [f.split('_catalog.hdf5')[0] for f in os.listdir(os.path.join(source_path, subdir))
                                  if f.endswith('_catalog.hdf5')]

            for run_name in run_names:
                print('  ----------------------------------------------')
                print('  Processing run: ', run_name)
                t_run = time.time()

                catalog_path = os.path.join(source_path, subdir, run_name + '_catalog.hdf5')
                config_path = os.path.join(source_path, subdir, run_name + '.yaml')

                time_sheet = self.get_mission_time(catalog_path=catalog_path,
                                                  config_path=config_path)

                # save time_sheet to csv files, one per experiment
                output_csv_path = os.path.join(source_path, subdir, run_name + '_mission_time.csv')
                time_sheet.to_csv(output_csv_path)
                print('    Saved mission time to ', output_csv_path)

                print('  Done processing run:', run_name, '- took', round(time.time() - t_run, 2), 's')

        # Append single multi-line log summary for this method
        t_end = time.time()
        subdirs_list = [d for d in os.listdir(source_path) if os.path.isdir(os.path.join(source_path, d))]
        total_runs = sum(len([f for f in os.listdir(os.path.join(source_path, d)) if f.endswith('_catalog.hdf5')]) for d in subdirs_list)
        mission_csvs = sum(len([f for f in os.listdir(os.path.join(source_path, d)) if f.endswith('_mission_time.csv')]) for d in subdirs_list)
        msg = f"""sweep_mission_time summary:
                  source_name: {source_name}
                  source_path: {source_path}
                  subdirs_found: {subdirs_list}
                  subdirs_total: {len(subdirs_list)}
                  total_runs_catalogs_found: {total_runs}
                  mission_time_csvs_found: {mission_csvs}
                  end_time: {time.ctime(t_end)}
                """
        self.logger.info(msg)

    def process_mission_time(self,
                             source_name):
        # Start logging for this processing run (style consistent with other functions)
        print('START OF process_mission_time: ', time.ctime())
        source_path = os.path.join(self.output_path, source_name)
        print('SOURCE NAME: ', source_name)
        print('SOURCE PATH: ', source_path)
        t_all = time.time()

        subdirs = [d for d in os.listdir(source_path) if os.path.isdir(os.path.join(source_path, d))]
        diams_float = [float('.'.join(d.split('_')[1:])) for d in subdirs]
        diams = ['_'.join(d.split('_')[1:]) for d in subdirs]

        # -- 1. CREATE TYPETABLE --
        timetable = {}

        for subdir, d in zip(subdirs, diams):
            csv_file = [f for f in os.listdir(os.path.join(source_path, subdir))
                        if f.endswith('.csv')]
            if len(csv_file) != 1:
                raise ValueError('More than one csv file found in subdir ' + subdir)
            timetable[d] = pd.read_csv(os.path.join(source_path, subdir, csv_file[0]))

        exps = [c.split('n_')[1] for c in timetable[diams[0]].columns if c.startswith('n_')]

        typetable = pd.DataFrame(columns=['detection', 'orbit', 'char', 'total'], index=diams_float)

        for d, df in zip(diams, diams_float):
            typetable.loc[df, 'detection'] = np.max(timetable[d]['detection']) * 1.25
            typetable.loc[df, 'orbit'] = np.mean(timetable[d]['orbit']) * 1.25
            typetable.loc[df, 'char'] = np.mean(timetable[d]['characterization']) * 1.25
            typetable.loc[df, 'total'] = typetable.loc[df, 'detection'] + typetable.loc[df, 'orbit'] + typetable.loc[
                df, 'char']

        typetable.sort_index(inplace=True, ascending=True)

        # -- 2. PLOT TYPETABLE --
        # fill under stepped line
        columns = ['detection', 'orbit', 'char']
        columns_label = ['Detection', 'Orbit', 'Characterization']
        x = np.asarray(typetable.index, dtype=float)
        y0 = np.zeros_like(x, dtype=float)

        colors = ['#2066a8', '#3594cc', '#8cc5e3']

        fig, ax = plt.subplots()

        for i in range(len(columns)):
            y = typetable[columns[i]].to_numpy(dtype=float) / (365.25 * 24 * 60 * 60)
            # ax.step(x, y+y0, label=columns[i], where='mid', color=colors[i])
            ax.fill_between(x, y0, y + y0, step='mid', alpha=1, color=colors[i], label=columns_label[i], edgecolor=None)

            y0 += typetable[columns[i]].to_numpy(dtype=float) / (365.25 * 24 * 60 * 60)

        ax.axhline(y=5, color='tab:blue', linestyle='--')
        ax.set_xlabel('Mirror Diameter (m)')
        ax.set_ylabel('Time (years)')
        ax.legend()
        type_fig_path = os.path.join(self.output_path, source_name, source_name + '_by_type_time.pdf')
        fig.savefig(type_fig_path)
        plt.close()

        # Save typetable CSV and report as a single cohesive step
        typetable_csv_path = os.path.join(self.output_path, source_name, source_name + '_typetable.csv')

        # -- 3. SAVE TYPETABLE --
        typetable /= (365.25 * 24 * 60 * 60)
        typetable.to_csv(typetable_csv_path)
        print('--------------------------------------------------')
        print('Step 1/2: Typetable generated and saved.')
        print('  Figure: ', type_fig_path)
        print('  CSV:    ', typetable_csv_path)

        # -- 4. CREATE EXPTABLE --
        exptable = {}

        for exp in exps:
            exptable[exp] = pd.DataFrame(
                columns=['detection', 'orbit', 'orbit-1s', 'orbit+1s', 'char', 'char-1s', 'char+1s', 'total',
                         'total-1s', 'total+1s'], index=diams_float)
            exptable[exp].sort_index(inplace=True)
        for subdir in subdirs:
            df = float('.'.join(subdir.split('_')[1:]))
            catalog_file = [f for f in os.listdir(os.path.join(source_path, subdir))
                        if f.endswith('.hdf5') and 'maxsep' not in f]
            if len(catalog_file) != 1:
                raise ValueError('More than one catalog file found in subdir ' + subdir)
            config_file = [f for f in os.listdir(os.path.join(source_path, subdir))
                        if f.endswith('.yaml')]

            bus = lifesim.Bus()

            # loading the config and catalog, make sure that the catalog is already combined with maxsep SNRs
            bus.build_from_config(
                filename=os.path.join(source_path, subdir, config_file[0]))
            bus.data.import_catalog(
                input_path=os.path.join(source_path, subdir, catalog_file[0]))
            cat = bus.data.catalog

            for exp in exps:
                mask = np.logical_and.reduce((cat.detected, cat['exp_' + exp]))
                exptable[exp].loc[df, 'detection'] = (np.sum(np.unique(cat[mask].int_time)) + len(
                    np.unique(cat[mask].int_time)) * bus.data.options.array['t_slew']) * 1.25

                mask = np.logical_and.reduce((cat.detected, cat['exp_' + exp], cat.follow_up))
                temp_t_orbit = []
                temp_t_char = []
                temp_t_total = []
                for nu in np.unique(cat[mask].nuniverse):
                    temp_t_orbit.append(np.sum(cat[np.logical_and(mask, cat.nuniverse == nu)].t_orbit))
                    temp_t_char.append(np.sum(cat[np.logical_and(mask, cat.nuniverse == nu)].t_char))
                    temp_t_total.append(np.sum(cat[np.logical_and(mask, cat.nuniverse == nu)].t_orbit) + np.sum(
                        cat[np.logical_and(mask, cat.nuniverse == nu)].t_char) + exptable[exp].loc[
                                            df, 'detection'] * 0.8)
                exptable[exp].loc[df, 'orbit'] = np.mean(temp_t_orbit) * 1.25
                exptable[exp].loc[df, 'orbit+1s'] = np.quantile(temp_t_orbit, 0.841) * 1.25
                exptable[exp].loc[df, 'orbit-1s'] = np.quantile(temp_t_orbit, 0.159) * 1.25
                exptable[exp].loc[df, 'char'] = np.mean(temp_t_char) * 1.25
                exptable[exp].loc[df, 'char+1s'] = np.quantile(temp_t_char, 0.841) * 1.25
                exptable[exp].loc[df, 'char-1s'] = np.quantile(temp_t_char, 0.159) * 1.25
                exptable[exp].loc[df, 'total'] = np.mean(temp_t_total) * 1.25
                exptable[exp].loc[df, 'total+1s'] = np.quantile(temp_t_total, 0.841) * 1.25
                exptable[exp].loc[df, 'total-1s'] = np.quantile(temp_t_total, 0.159) * 1.25

            del bus
            del cat

        # -- 5. PLOT EXPTABLE --
        colors = [['#2066a8', '#3594cc', '#8cc5e3'],
                  ['#a00000', '#c46666', '#d8a6a6'],
                  ['#1f6f6f', '#54a1a1', '#9fc8c8'], ]

        columns = ['detection', 'orbit', 'char']
        x = np.asarray(exptable[exps[0]].index, dtype=float)
        y0 = np.zeros_like(x, dtype=float)

        fig, ax = plt.subplots()

        for j, exp in enumerate(exps):
            for i in range(len(columns)):
                y = exptable[exp][columns[i]].to_numpy(dtype=float) / (365.25 * 24 * 60 * 60)
                # ax.step(x, y+y0, label=columns[i], where='mid', color=colors[i])
                if i == 0:
                    ax.fill_between(x, y0, y + y0, step='mid', alpha=1, color=colors[j][i], edgecolor=None, label=exp)
                else:
                    ax.fill_between(x, y0, y + y0, step='mid', alpha=1, color=colors[j][i], edgecolor=None)
                y0 += y

        ax.set_xlabel('Mirror Diameter (m)')
        ax.set_ylabel('Time (years)')
        ax.axhline(y=5, color='tab:blue', linestyle='--')
        ax.legend()

        exptable_fig_path = os.path.join(self.output_path, source_name, source_name + '_by_exp_time.pdf')
        fig.savefig(exptable_fig_path)
        plt.close()

        # -- 6. SAVE EXPTABLE --
        for exp in exps:
            exptable[exp] /= (365.25 * 24 * 60 * 60)
            exptable[exp].to_csv(os.path.join(self.output_path, source_name, source_name + '_' + exp + '_exptable.csv'))

        print('--------------------------------------------------')
        print('Step 2/2: Exptable plotted and CSVs saved.')
        print('  Figure: ', exptable_fig_path)
        print('  CSVs:   ', os.path.join(self.output_path, source_name))

        # Final log with elapsed time
        print('ALL process_mission_time finished. Total time:', round(time.time() - t_all, 2), 's')


    def run_covergence_test(self,
                            run_name,
                            catalog_path,
                            min_universes,
                            num_steps,
                            plot=False):

        # check a directory with the name run_name already exists in output_path, otherwise create it
        output_path = os.path.join(self.output_path, run_name)
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        else:
            raise ValueError('Directory already exists: ' + output_path)

        # determine number universes
        catalog = pd.read_hdf(catalog_path, key='catalog')
        total_universes = np.unique(catalog.nuniverse).shape[0]
        del catalog

        if min_universes > total_universes-1:
            raise ValueError('min_universes must be less than the total number of universes in the catalog.')

        if num_steps is None:
            num_steps = total_universes - min_universes + 1
        elif num_steps > total_universes - min_universes + 1:
            raise ValueError('num_steps is too large for the given min_universes and total universes in the catalog.')

        uni_sel = np.linspace(min_universes, total_universes, num_steps, dtype=int)

        if self.n_cpu > 1:
            with parallel_config(
                    backend="loky", inner_max_num_threads=1
            ), joblib_progress(
                description="Running SNR calculation in parallel ...",
                total=len(uni_sel),
            ):
                results = Parallel(n_jobs=self.n_cpu)(
                    delayed(compute_yields_mp)(
                        output_filename=run_name + '_nuni_' + str(us),
                        output_path=output_path,
                        catalog_path=catalog_path,
                        config_path=self.config_path,
                        uni_sel=us,
                        return_yields=True
                    )
                    for us in uni_sel)
        else:
            results = []
            for us in tqdm(uni_sel):
                res = compute_yields_mp(
                    output_filename=run_name + '_nuni_' + str(us),
                    output_path=output_path,
                    catalog_path=catalog_path,
                    config_path=self.config_path,
                    uni_sel=us,
                    return_yields=True
                )

                results.append(res)

        experiments = list(results[0][1].keys())
        yield_convergence = {exp: {'n_universes': [],
                                   'mean_yield': [],
                                   'p16_yield': [],
                                   'p84_yield': []}
                             for exp in experiments}
        for res in results:
            n_uni, yields = res
            for exp in experiments:
                yield_convergence[exp]['n_universes'].append(n_uni)
                yield_convergence[exp]['mean_yield'].append(yields[exp][0])
                yield_convergence[exp]['p16_yield'].append(yields[exp][1])
                yield_convergence[exp]['p84_yield'].append(yields[exp][2])

        # save yield_convergence to a json file
        with open(os.path.join(output_path, run_name + '_yield_convergence.json'), 'w') as f:
            json.dump(yield_convergence, f, indent=4)

        if plot:
            fig, ax = plt.subplots(figsize=(8, 6))
            for exp in experiments:
                ax.plot(yield_convergence[exp]['n_universes'],
                        yield_convergence[exp]['mean_yield'],
                        label=exp)
                ax.fill_between(yield_convergence[exp]['n_universes'],
                                yield_convergence[exp]['p16_yield'],
                                yield_convergence[exp]['p84_yield'],
                                alpha=0.3)
            ax.set_xlabel('Number of Universes')
            ax.set_ylabel('Yield')
            ax.set_title('Yield Convergence Test')
            ax.legend()
            plt.savefig(os.path.join(output_path, run_name + '_yield_convergence.pdf'))
            plt.close()


def get_yields(bus,
               return_yields=False):

    yields = {}

    for exp in bus.data.options.optimization['experiments'].keys():
        mask = np.logical_and.reduce((bus.data.catalog.temp_s
                                      >= bus.data.options.optimization['experiments'][exp]['temp_s_min'],
                                      bus.data.catalog.temp_s 
                                      <= bus.data.options.optimization['experiments'][exp]['temp_s_max'],
                                      bus.data.catalog.radius_p
                                      >= bus.data.options.optimization['experiments'][exp]['radius_p_min'],
                                      bus.data.catalog.radius_p
                                      <= bus.data.options.optimization['experiments'][exp]['radius_p_max'],
                                      bus.data.catalog.detected
                                     ))
        if bus.data.options.optimization['experiments'][exp]['in_HZ']:
            mask = np.logical_and(mask,
                                  bus.data.catalog.habitable)

        _, tyield = np.unique(bus.data.catalog[mask].nuniverse, return_counts=True)

        if tyield.size != 0:
            yields[exp] = [float(np.mean(tyield)), float(np.percentile(tyield, 15.9)), float(np.percentile(tyield, 84.1))]
        else:
            yields[exp] = [0.0, 0.0, 0.0]

    return yields

def compute_yields_mp(output_filename,
                      output_path,
                      catalog_path,
                      config_path,
                      uni_sel=None,
                      return_yields=False,
                      characterization: bool = False,
                      opt_limit_factor: Union[None, float] = None):

    t = time.time()
    # create bus
    bus = lifesim.Bus()

    # setting the options
    bus.build_from_config(filename=config_path)
    bus.data.options.set_manual(n_cpu=1)  # speed up calculation

    if opt_limit_factor is not None:
        bus.data.options.optimization['opt_limit_factor'] = opt_limit_factor

    bus.data.options.set_manual(output_path=output_path)
    bus.data.options.set_manual(output_filename=output_filename)
    bus.data.options.optimization['characterization'] = characterization
    # ---------- Loading the Catalog ----------
    bus.data.import_catalog(input_path=catalog_path)

    if uni_sel is not None:
        selected_universes = np.random.choice(np.unique(bus.data.catalog.nuniverse), replace=False, size=uni_sel)
        bus.data.catalog = bus.data.catalog[np.isin(bus.data.catalog.nuniverse, selected_universes)]

    # ---------- Creating the Instrument ----------

    # create modules and add to bus
    instrument = lifesim.Instrument(name='inst')
    bus.add_module(instrument)

    transm = lifesim.TransmissionMap(name='transm')
    bus.add_module(transm)

    # connect all modules
    bus.connect(('inst', 'transm'))

    # optimizing the result
    opt = lifesim.Optimizer(name='opt')
    bus.add_module(opt)
    ahgs = lifesim.AhgsModule(name='ahgs')
    bus.add_module(ahgs)

    bus.connect(('transm', 'opt'))
    bus.connect(('inst', 'opt'))
    bus.connect(('opt', 'ahgs'))

    # with contextlib.redirect_stdout(None):
    opt.ahgs()
    bus.save()

    if return_yields:
        yields = get_yields(bus=bus,
                            return_yields=return_yields)
        return int(uni_sel), yields
