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

            if self.data.options.optimization['characterization']:
                return np.array(((np.inf, np.inf), (np.inf, np.inf)))
            else:
                return np.array((np.inf, np.inf))
        else:
            obs = (60 * 60 *
                   (self.data.options.optimization['snr_target'] ** 2
                    - self.data.catalog['snr_current'].loc[mask] ** 2)
                   / self.data.catalog.snr_1h.loc[mask] ** 2)

            # included slew time
            # needs to be divided by the number of universes for proper optimization, since the characterization happens
            # 'per universe'
            if self.data.options.optimization['characterization']:
                obs_char = (60 * 60 *
                            self.data.options.optimization['snr_char'] ** 2
                            / self.data.catalog.maxsep_snr_1h.loc[mask] ** 2
                            + self.data.options.array['t_slew']) / self.data.optm['num_universe']

                # characterization time is needed for optimization, detection time is needed to understand how much
                # mission time was used for one observation
                met = np.stack((obs, obs_char), axis=0)

                met = met[:, np.argsort(met[0])]
                met[1] = np.cumsum(met[1])
                met[0] -= self.data.catalog.t_slew.loc[mask]
                met[1] += met[0]
                met = met / np.arange(1, np.count_nonzero(mask) + 1, 1)[np.newaxis, :]

                return met
            else:
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

            # Changed: use label-based .loc to avoid chained-assignment warnings
            for _, i in enumerate(np.where(mask)[0]):
                idx = self.data.catalog.index[i]
                if (not self.data.catalog.loc[idx, 'detected']) and \
                        (self.data.catalog.loc[idx, 'snr_current']
                         >= self.data.options.optimization['snr_target']):
                    self.data.catalog.loc[idx, 'detected'] = True
                    self.data.catalog.loc[idx, 't_detected'] = deepcopy(self.tot_time + int_time)
                    exp_cols = [col for col in self.data.catalog.columns if col.startswith('exp_')]
                    true_experiments = [col[4:] for col in exp_cols if self.data.catalog.loc[idx, col]]

                    for exp in true_experiments[1:]:
                        self.data.optm['exp_detected'][exp] += 1
                        self.data.optm['exp_detected_uni'][exp][
                            1,
                            self.data.optm['exp_detected_uni'][exp][0, :] == self.data.catalog.loc[idx, 'nuniverse']
                        ] += 1
        else:
            raise ValueError('Delete mode not implemented for AHGS optimizer.')

    def distribute_time(self):
        stars, n = np.unique(ar=self.data.catalog.nstar,
                             return_counts=True)
        if self.data.options.optimization['characterization']:
            obs = np.zeros((stars.shape[0], 2, np.max(n))) + np.inf
            # fill the observation time array
            for i, nstar in enumerate(stars):
                temp = self.obs_array_star(nstar=nstar)
                obs[i, :, :temp.shape[1]] = temp
        else:
            obs = np.zeros((stars.shape[0], np.max(n))) + np.inf
            # fill the observation time array
            for i, nstar in enumerate(stars):
                temp = self.obs_array_star(nstar=nstar)
                obs[i, :temp.shape[0]] = temp

        obs_time = (self.data.options.optimization['t_search']
                    * self.data.options.array['t_efficiency'])

        self.tot_time = 0

        print('Number of planets detected for each experiment:')

        run_bool = True

        while run_bool:
            if self.data.options.optimization['characterization']:
                # find the best global slope and observe star
                no_star, ind_t = np.unravel_index(np.argmin(obs[:, 1, :]), obs[:, 1, :].shape)
                if not np.isfinite(obs[no_star, 1, ind_t]):
                    print('Not sufficient targets remaining to continue characterization optimization.')
                    break
                if (((self.tot_time + obs[no_star, 0, ind_t] * (ind_t + 1) + 0.01) > obs_time)
                        and (self.data.options.optimization['opt_limit'] == 'time')):
                    rem_time = obs_time - self.tot_time
                    self.observe_star(nstar=stars[no_star],
                                      int_time=rem_time)
                    self.tot_time += rem_time
                else:
                    if ((obs[no_star, 1, ind_t] - obs[no_star, 0, ind_t]) * (ind_t + 1) + 0.01) < 0:
                        raise ValueError('Negative time difference encountered.')
                    self.observe_star(nstar=stars[no_star],
                                      int_time=obs[no_star, 0, ind_t] * (ind_t + 1) + 0.01)
                    temp = self.obs_array_star(nstar=stars[no_star])
                    self.tot_time += obs[no_star, 0, ind_t] * (ind_t + 1) + 0.01
                    obs[no_star, :, :] = np.inf
                    obs[no_star, :, :temp.shape[1]] = temp
            else:
                # find the best global slope and observe star
                no_star, ind_t = np.unravel_index(np.argmin(obs), obs.shape)
                if (((self.tot_time + obs[no_star, ind_t] * (ind_t + 1) + 0.01) > obs_time)
                        and (self.data.options.optimization['opt_limit'] == 'time')):
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

            if self.data.options.optimization['opt_limit'] == 'experiments':
                out_string += ('-  ' + str(np.round(self.tot_time / 60 / 60 / 24 / 365.25, decimals=1))
                               + ' yrs observed')
            else:
                out_string += ('-  (' + str(np.round(self.tot_time / 60 / 60 / 24 / 365.25, decimals=1)) + ' / '
                                              + str(np.round(obs_time / 60 / 60 / 24 / 365.25, decimals=1))
                               + ') yrs observed')
            print('\r' + out_string, end='')

            # if any(
            #         self.data.optm['exp_detected'][exp] / self.data.optm['num_universe']
            #         > self.data.options.optimization['experiments'][exp]['sample_size']
            #         and not self.data.optm['hit_limit'][exp]
            #         for exp in self.data.optm['exp_detected']
            # ):
            if any([
                ((self.data.optm['exp_detected_uni'][exp][1, :]
                 > self.data.options.optimization['experiments'][exp]['sample_size']).sum() >
                (self.data.options.optimization['opt_limit_factor']
                        * self.data.optm['num_universe'])) and not self.data.optm['hit_limit'][exp]
                for exp in self.data.optm['exp_detected_uni']]):
                # over_limit_experiments = [
                #     exp for exp in self.data.optm['exp_detected']
                #     if
                #     self.data.optm['exp_detected'][exp] > self.data.options.optimization['experiments'][exp][
                #         'sample_size']
                #     and not self.data.optm['hit_limit'][exp]
                # ]
                over_limit_experiments = [
                    exp for exp in self.data.optm['exp_detected_uni']
                    if
                    (self.data.optm['exp_detected_uni'][exp][1, :]
                     > self.data.options.optimization['experiments'][exp]['sample_size']).sum() >
                    (self.data.options.optimization['opt_limit_factor']
                            * self.data.optm['num_universe'])
                    and not self.data.optm['hit_limit'][exp]
                ]

                for exp in over_limit_experiments:
                    self.data.optm['hit_limit'][exp] = True
                    self.data.catalog['is_interesting'] = False
                    for exp_interesting in [exp for exp, hit in self.data.optm['hit_limit'].items() if not hit]:
                        self.data.catalog['is_interesting'] = np.logical_or(self.data.catalog['is_interesting'],
                                                                           self.data.catalog['exp_' + exp_interesting])

                if self.data.catalog['is_interesting'].sum() == 0:
                    if self.data.options.optimization['opt_limit'] == 'time':
                        print('\nAll experiments have been completed, spending remaining mission time on all HZ planets.')
                    self.data.catalog['is_interesting'] = self.data.catalog['habitable']

                else:
                    print('\nCompleted experiments: ' + ', '.join(over_limit_experiments)
                          + ', RECOUNTING -------------------')

                # fill the observation time array
                if self.data.options.optimization['characterization']:
                    obs = np.zeros((stars.shape[0], 2, np.max(n))) + np.inf
                    # fill the observation time array
                    for i, nstar in enumerate(stars):
                        temp = self.obs_array_star(nstar=nstar)
                        obs[i, :, :temp.shape[1]] = temp
                else:
                    obs = np.zeros((stars.shape[0], np.max(n))) + np.inf
                    # fill the observation time array
                    for i, nstar in enumerate(stars):
                        temp = self.obs_array_star(nstar=nstar)
                        obs[i, :temp.shape[0]] = temp

            if self.data.options.optimization['opt_limit'] == 'time':
                run_bool = self.tot_time < obs_time
            elif self.data.options.optimization['opt_limit'] == 'experiments':
                run_bool = not all(self.data.optm['hit_limit'].values())
            else:
                raise ValueError('Optimization limit not recognized.')