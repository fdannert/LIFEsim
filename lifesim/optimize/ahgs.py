import numpy as np

from lifesim.core.modules import SlopeModule

class AhgsModule(SlopeModule):
    def __init__(self,
                 name: str):
        super().__init__(name=name)

    def obs_array_star(self, nstar):
        mask = self.data.catalog.nstar == nstar

        # return infinity if the detection limit is reached for this stype
        if not bool(np.isin(element=self.data.catalog.stype.loc[mask].iloc[0],
                            test_elements=np.array(list(self.data.options.optimization['limit'].keys()))[np.invert(
                                self.data.optm['hit_limit'])])):
            return np.array((np.inf, np.inf))
        else:
            if self.data.options.optimization['habitable']:
                mask = np.logical_and.reduce((mask,
                                              self.data.catalog.habitable,
                                              np.invert(self.data.catalog.detected)))
            else:
                mask = np.logical_and(mask, np.invert(self.data.catalog.detected))

            if self.data.options.optimization['instrumental_opt']:
                obs = (self.data.catalog['photon_noise'].loc[mask] ** 2 * self.data.catalog['t_rot'].loc[mask]
                       / ((self.data.catalog['signal'].loc[mask]
                           / self.data.options.optimization['snr_target']) ** 2
                          - self.data.catalog['systematic_noise'].loc[mask] ** 2))
                obs[obs <= 0] = np.inf
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

            self.data.catalog.loc[mask, 'int_time_actual'] += int_actual

            if self.data.options.optimization['instrumental_opt']:
                self.data.catalog.loc[mask, 'snr_current'] = (
                        self.data.catalog.loc[mask, 'signal']
                        / np.sqrt(self.data.catalog.loc[mask, 'photon_noise'] ** 2
                                  * self.data.catalog.loc[mask, 't_rot']
                                  / self.data.catalog.loc[mask, 'int_time_actual']
                                  + self.data.catalog.loc[mask, 'systematic_noise'] ** 2))
            else:
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
                    if self.data.catalog.habitable.iloc[i]:
                        self.data.optm['sum_detected'][
                            np.where(np.array(list(self.data.options.optimization['limit'].keys()))
                                     == self.data.catalog.stype.iloc[i])] += 1
        else:
            pass
            # self.planets.slew_time[mask_star] = -self.t_slew
            # self.planets.int_time[mask_star] = 0
            # self.planets.SNR_prediction_sum[mask_star] = 0
            # for _, i in enumerate(np.where(mask_star)[0]):
            #     if self.planets.detected[i] and self.planets.habitable[i]:
            #         self.sum_detected[np.where(self.limits[0][:] == self.planets.Stype[i])] \
            #             -= 1
            # self.planets.detected[mask_star] = 0

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

        tot_time = 0

        print('Number of planets detected by stellar type:')

        while tot_time < obs_time:
            # find the best global slope and observe star
            no_star, ind_t = np.unravel_index(np.argmin(obs), obs.shape)
            if (tot_time + obs[no_star, ind_t] * (ind_t + 1) + 0.01) > obs_time:
                rem_time = obs_time - tot_time
                self.observe_star(nstar=stars[no_star],
                                  int_time=rem_time)
                tot_time += rem_time
                # if False:
                #     tot_time += rem_time
                #     return rem_time
                # else:
                #     self.observe_star(nstar=stars[no_star],
                #                       int_time=rem_time)
                #     tot_time += rem_time
            else:
                self.observe_star(nstar=stars[no_star],
                                  int_time=obs[no_star, ind_t] * (ind_t + 1) + 0.01)
                temp = self.obs_array_star(nstar=stars[no_star])
                tot_time += obs[no_star, ind_t] * (ind_t + 1) + 0.01
                obs[no_star, :] = np.inf
                obs[no_star, :temp.shape[0]] = temp

            out_string = ''
            for key in self.data.options.optimization['limit'].keys():
                out_string += (key + ': '
                               + str((self.data.optm['sum_detected'] / self.data.optm['num_universe'])[
                                         np.where(np.array(list(self.data.options.optimization['limit'].keys()))
                                                  == key)][0])
                               + '  ')
            out_string += ('-  (' + str(np.round(tot_time/60/60/24/365.25, 1)) + ' / '
                           + str(np.round(obs_time/60/60/24/365.25, 1)) + ') yrs observed')
            print('\r' + out_string, end='')

            if np.any(
                    np.logical_and(
                        (self.data.optm['sum_detected'] / self.data.optm['num_universe'])
                        > np.array(list(self.data.options.optimization['limit'].values())),
                        np.invert(self.data.optm['hit_limit']))):
                print('\n')
                print('HIT LIMIT, RECOUNTING -------------------')
                self.data.optm['hit_limit'] = ((self.data.optm['sum_detected']
                                                / (self.data.optm['num_universe']))
                                               >= np.array(
                            list(self.data.options.optimization['limit'].values())
                        ))
                obs = np.zeros((stars.shape[0], np.max(n))) + np.inf

                # fill the observation time array
                for i, nstar in enumerate(stars):
                    temp = self.obs_array_star(nstar=nstar)
                    obs[i, :temp.shape[0]] = temp
