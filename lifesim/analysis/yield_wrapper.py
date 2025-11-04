import time
import shutil
import os

import numpy as np

import lifesim

class ScienceYield:
    def __init__(self,
                 config_path,
                 catalog_path,
                 output_path,
                 n_cpu: int = 1):
        self.config_path = config_path
        self.catalog_path = catalog_path
        self.output_path = output_path
        self.n_cpu = n_cpu

    def compute_yield(self,
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

    def run_aperture_sweep(self,
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
            self.compute_yield(output_path=f'{output_directory}/',
                      output_filename='sweep_diam_' + str(np.round(diameter, 2)).replace('.', '_'),
                      run_maxsep=False)
            print('[Done]')

            print('Commencing maxsep run... ')
            self.compute_yield(output_path=f'{output_directory}/',
                               output_filename='sweep_diam_maxsep_' + str(np.round(diameter, 2)).replace('.', '_'),
                               run_maxsep=True)
            print('[Done]')
