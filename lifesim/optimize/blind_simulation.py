import numpy as np
import pandas as pd
from copy import deepcopy

from lifesim.core.modules import RealisticModule


class BlindSimulationModule(RealisticModule):
    def __init__(self,
                 name: str):
        super().__init__(name=name)

    def run_simulation(self):
        self.data.optm['t_current'] = 0.
        self.data.optm['array_ang_current'] = 0.


        # calculate cartesian vectors for all stars
        for i in range(self.data.stars.shape[0]):
            self.data.stars.vector.iat[i] = [np.sin(np.pi/2 - self.data.stars.lat.iloc[i])
                                                * np.cos(self.data.stars.lon.iloc[i]),
                                                np.sin(np.pi / 2 - self.data.stars.lat.iloc[i])
                                                * np.sin(self.data.stars.lon.iloc[i]),
                                                np.cos(np.pi / 2 - self.data.stars.lat.iloc[i])]

        # add time to detection for each planet
        self.data.catalog['t_detection'] = (60 * 60
                                            * self.data.options.optimization['snr_target'] ** 2
                                            / self.data.catalog.snr_1h ** 2)

        # do horrible frankenstein switching of catalogs
        self.data.statistic = deepcopy(self.data.catalog)
        self.data.catalog = self.data.universe
        self.data.catalog['snr_current'] = 0.
        self.data.catalog['detected'] = False
        self.data.catalog['t_to_dect'] = np.inf

        # precalculate stuff
        # self.run_socket(s_name='instrument',
        #                 method='get_snr',
        #                 safe_mode=True)
        # self.run_socket(s_name='interest',
        #                 method='do_precalc')
        # self.data.export_catalog(output_path='temp1')
        # self.data.stars.to_hdf(path_or_buf='temp2', key='stars', mode='w')
        self.run_socket(s_name='instrument',
                        method='apply_options')
        self.data.stars = pd.read_hdf(path_or_buf='temp2',
                                      key='stars')
        self.data.catalog = None
        self.data.import_catalog(input_path='temp')

        obs_last = (None, None)

        while (self.data.optm['t_current']
               < (self.data.options.optimization['t_search']
                  * self.data.options.array['t_efficiency'])):
            print(self.data.optm['t_current'] / (60 * 60 * 24 * 365))
            self.run_socket(s_name='interest',
                            method='in_for')
            for i, nstar in enumerate(self.data.stars.nstar):
                self.run_socket(s_name='interest',
                                method='get_star',
                                nstar=nstar)

            obs_nstar = self.data.stars.nstar.iloc[self.data.stars.interest_current.argmax()]

            print(obs_nstar)
            t = self.run_socket(s_name='observation',
                                method='observe',
                                nstar=obs_nstar)
            self.data.optm['t_current'] += t
            self.data.optm['array_ang_current'] += 2 * np.pi / (60 * 60 * 24 * 365) * t
            if self.data.options.optimization['multi_visit']:
                mask = self.data.stars.nstar == obs_nstar
                if self.data.stars.loc[mask, 'frustration'].iloc[0] != 0:
                    self.data.stars.loc[mask, 'i_number'].iat[0] = (
                        self.data.stars.loc[mask, 'frustration'])
                    self.run_socket(s_name='interest',
                                    method='precalc_single',
                                    nstar=obs_nstar)

            if obs_nstar == obs_last[1]:
                break
            obs_last = (obs_nstar, obs_last[0])

