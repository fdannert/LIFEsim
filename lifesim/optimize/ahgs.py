from copy import deepcopy

import numpy as np

from lifesim.core.modules import SlopeModule

class AhgsModule(SlopeModule):
    def __init__(self,
                 name: str):
        super().__init__(name=name)

        self.tot_time = 0

    def obs_array_star(self, nstar):
        mask = np.logical_and.reduce((self.data.catalog.nstar == nstar,
                                np.invert(self.data.catalog.detected),
                                      self.data.catalog.is_interesting))

        if not np.any(mask):
            return np.array((np.inf, np.inf))
        else:
            obs = (60 * 60 *
                   (self.data.options.optimization['snr_target'] ** 2
                    - self.data.catalog['snr_current'].loc[mask] ** 2)
                   / self.data.catalog.snr_1h.loc[mask] ** 2)
            obs -= self.data.catalog.t_slew.loc[mask]
            obs = np.sort(obs) / np.arange(1, np.count_nonzero(mask) + 1, 1)
            return obs

    def observe_star(self,
                     nstar,
                     int_time,
                     delete=False):
        mask = self.data.catalog.nstar == nstar

        if not delete:
            self.data.optm['tot_time'] += int_time

            slew_time = self.data.catalog.loc[mask, 't_slew'].iloc[0]
            if not (slew_time == 0):
                if (slew_time + int_time) < 0:
                    self.data.catalog.loc[mask, 't_slew'] += int_time
                    int_actual = 0
                else:
                    self.data.catalog.loc[mask, 't_slew'] = 0
                    self.data.catalog.loc[mask, 'int_time'] += (slew_time + int_time)
                    int_actual = slew_time + int_time
            else:
                self.data.catalog.loc[mask, 'int_time'] += int_time
                int_actual = int_time

            self.data.catalog.loc[mask, 'snr_current'] = np.sqrt(
                self.data.catalog.loc[mask, 'snr_current'] ** 2
                + (self.data.catalog.loc[mask, 'snr_1h']
                   * np.sqrt(int_actual
                             / (60 * 60)))**2)

            for _, i in enumerate(np.where(mask)[0]):
                if (not self.data.catalog.detected.iloc[i]) and \
                        (self.data.catalog.snr_current.iloc[i]
                         >= self.data.options.optimization['snr_target']):
                    self.data.catalog.detected.iat[i] = True
                    self.data.catalog.t_detected.iat[i] = deepcopy(self.tot_time)
                    exp_cols = [col for col in self.data.catalog.columns if col.startswith('exp_')]
                    true_experiments = [col[4:] for col in exp_cols if self.data.catalog.at[i, col]]

                    for exp in true_experiments:
                        self.data.optm['exp_detected'][exp] += 1
        else:
            raise ValueError('Delete mode not implemented for AHGS optimizer.')

    def distribute_time(self):
        stars, n = np.unique(ar=self.data.catalog.nstar,
                             return_counts=True)
        obs = np.zeros((stars.shape[0], np.max(n))) + np.inf

        # fill the observation time array
        for i, nstar in enumerate(stars):
            temp = self.obs_array_star(nstar=nstar)
            obs[i, :temp.shape[0]] = temp

        obs_time = (self.data.options.optimization['t_search']
                    * self.data.options.array['t_efficiency'])

        self.tot_time = 0

        print('Number of planets detected for each experiment:')

        while self.tot_time < obs_time:
            # find the best global slope and observe star
            no_star, ind_t = np.unravel_index(np.argmin(obs), obs.shape)
            if (self.tot_time + obs[no_star, ind_t] * (ind_t + 1) + 0.01) > obs_time:
                rem_time = obs_time - self.tot_time
                self.observe_star(nstar=stars[no_star],
                                  int_time=rem_time)
                self.tot_time += rem_time
            else:
                self.observe_star(nstar=stars[no_star],
                                  int_time=obs[no_star, ind_t] * (ind_t + 1) + 0.01)
                temp = self.obs_array_star(nstar=stars[no_star])
                self.tot_time += obs[no_star, ind_t] * (ind_t + 1) + 0.01
                obs[no_star, :] = np.inf
                obs[no_star, :temp.shape[0]] = temp

            out_string = ''
            for key in self.data.optm['exp_detected'].keys():
                out_string += (key + ': '
                               + str(self.data.optm['exp_detected'][key] / self.data.optm['num_universe'])
                               + '  ')
            out_string += ('-  (' + str(np.round(self.tot_time / 60 / 60 / 24 / 365.25, decimals=1)) + ' / '
                                          + str(np.round(obs_time / 60 / 60 / 24 / 365.25, decimals=1))
                           + ') yrs observed')
            print('\r' + out_string, end='')

            if any(
                    self.data.optm['exp_detected'][exp] / self.data.optm['num_universe']
                    > self.data.options.optimization['experiments'][exp]['sample_size']
                    and not self.data.optm['hit_limit'][exp]
                    for exp in self.data.optm['exp_detected']
            ):
                over_limit_experiments = [
                    exp for exp in self.data.optm['exp_detected']
                    if
                    self.data.optm['exp_detected'][exp] > self.data.options.optimization['experiments'][exp][
                        'sample_size']
                    and not self.data.optm['hit_limit'][exp]
                ]
                for exp in over_limit_experiments:
                    self.data.optm['hit_limit'][exp] = True
                    self.data.catalog['is_interesting'] = False
                    for exp_interesting in [exp for exp, hit in self.data.optm['hit_limit'].items() if not hit]:
                        self.data.catalog['is_interesting'] = np.logical_or(self.data.catalog['is_interesting'],
                                                                           self.data.catalog['exp_' + exp_interesting])

                if self.data.catalog['is_interesting'].sum() == 0:
                    print('\nAll experiments have been completed, spending remaining mission time on all HZ planets.')
                    self.data.catalog['is_interesting'] = self.data.catalog['habitable']

                else:
                    print('\nCompleted experiments: ' + ', '.join(over_limit_experiments)
                          + ', RECOUNTING -------------------')

                obs = np.zeros((stars.shape[0], np.max(n))) + np.inf

                # fill the observation time array
                for i, nstar in enumerate(stars):
                    temp = self.obs_array_star(nstar=nstar)
                    obs[i, :temp.shape[0]] = temp
