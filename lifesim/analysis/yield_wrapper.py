import time
import shutil
import os
import contextlib
from copy import deepcopy
import json

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

    def compute_snrs(self,
                     output_path,
                     output_filename,
                     run_maxsep):

        print('START OF RUN: ', time.ctime())
        print('RUN NAME: ', output_filename)
        print('----------------------------')

        t = time.time()
        # create bus
        bus = lifesim.Bus()

        # setting the options
        bus.build_from_config(filename=self.config_path)
        bus.data.options.set_manual(n_cpu=self.n_cpu) # speed up calculation

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
            bus.data.catalog = bus.data.catalog[bus.data.catalog.habitable]
            bus.data.catalog['angsep'] = bus.data.catalog['maxangsep']

        instrument.get_snr()
        bus.save()

        del bus

        print('Generation took ', (time.time() - t) / 60, ' minutes to complete.')

    def run_aperture_sweep_snr(self,
                               mirror_diameters):

        for ndim, diameter in enumerate(mirror_diameters):
            print('')
            print('')
            print('STARTING RUN FOR MIRROR DIAMETER: ', diameter)
            print('AT TIME: ', time.ctime())

            print('Preparing directories... ', end='')

            output_directory = os.path.join(self.output_path, 'diam_' + str(np.round(diameter, 2)).replace('.', '_'))
            if not os.path.exists(output_directory):
                os.makedirs(output_directory)

            print('[Done]')

            print('Commencing base run... ')
            self.compute_snrs(output_path=f'{output_directory}/',
                              output_filename='sweep_diam_' + str(np.round(diameter, 2)).replace('.', '_'),
                              run_maxsep=False)
            print('[Done]')

            print('Commencing maxsep run... ')
            self.compute_snrs(output_path=f'{output_directory}/',
                              output_filename='sweep_diam_maxsep_' + str(np.round(diameter, 2)).replace('.', '_'),
                              run_maxsep=True)
            print('[Done]')

    def run_optimizer_sweep(self,
                            run_name):
        # get the names of all subdirectories in output_path
        subdirs = [d for d in os.listdir(self.output_path) if os.path.isdir(os.path.join(self.output_path, d))]

        # check a directory with the name run_name already exists in output_path, otherwise create it
        final_output_path = os.path.join(self.output_path, run_name)
        if not os.path.exists(final_output_path):
            os.makedirs(final_output_path)
        else:
            raise ValueError('Directory already exists: ' + final_output_path)

        run_configs = []

        for subdir in subdirs:
            # make a list of all files ending in .hdf5 in subdir, then keep only the part of the filename before '_catalog.hdf5' and only if the sting does not contain 'maxsep'
            catalog_files = [f.split('_catalog.hdf5')[0] for f in os.listdir(os.path.join(self.output_path, subdir))
                                if f.endswith('_catalog.hdf5') and 'maxsep' not in f]

            # create a subdir of the same name in final_output_path, no existence check required
            output_directory = os.path.join(final_output_path, subdir)
            os.makedirs(output_directory)

            for catalog_file in catalog_files:
                if self.n_cpu == 1:
                    compute_yields_mp(
                        output_filename=catalog_file,
                        output_path=output_directory + '/',
                        catalog_path=os.path.join(self.output_path, subdir, catalog_file + '_catalog.hdf5'),
                        config_path=self.config_path
                    )
                else:
                    run_configs.append({'output_filename': catalog_file,
                                        'output_path': output_directory + '/',
                                        'catalog_path': os.path.join(self.output_path, subdir, catalog_file + '_catalog.hdf5')
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
                        config_path=self.config_path
                    )
                    for rc in run_configs)

    def run_covergence_test(self,
                            run_name,
                            catalog_path,
                            output_path,
                            min_universes,
                            num_steps,
                            plot=False):
        if not os.path.exists(output_path):
            os.makedirs(output_path)

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
                      return_yields=False):

    t = time.time()
    # create bus
    bus = lifesim.Bus()

    # setting the options
    bus.build_from_config(filename=config_path)
    bus.data.options.set_manual(n_cpu=1)  # speed up calculation

    bus.data.options.set_manual(
        output_path=output_path)
    bus.data.options.set_manual(output_filename=output_filename)

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

    with contextlib.redirect_stdout(None):
        opt.ahgs()

    bus.save()

    if return_yields:
        yields = get_yields(bus=bus,
                            return_yields=return_yields)
        return int(uni_sel), yields
