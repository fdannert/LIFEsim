import numpy as np

from lifesim.core.modules import SlopeModule

class AhgsCharModule(SlopeModule):
    def __init__(self,
                 name: str):
        super().__init__(name=name)
        self.t_series = []
        self.yield_series = {'A': [],
                             'F': [],
                             'G': [],
                             'K': [],
                             'M': []}
        self.id_series = []

    def obs_array_star(self, nstar):
        mask = self.data.catalog.nstar == nstar

        # return infinity if the detection limit is reached for this stype
        if not bool(np.isin(element=self.data.catalog.stype.loc[mask].iloc[0],
                            test_elements=np.array(list(self.data.options.optimization['limit'].keys()))[np.invert(
                                self.data.optm['hit_limit'])])):
            return np.array((np.inf, np.inf))
        else:
            mask = np.logical_and.reduce((mask,
                                          self.data.catalog.habitable,
                                          np.invert(self.data.catalog.detected)))

            obs = (60 * 60 *
                   (self.data.options.optimization['snr_target'] ** 2
                    - self.data.catalog['snr_current'].loc[mask] ** 2)
                   / self.data.catalog.snr_1h.loc[mask] ** 2)

            obs_char = (60 * 60 *
                   self.data.options.optimization['snr_char'] ** 2
                   / self.data.catalog.maxsep_snr_1h.loc[mask] ** 2)

            met = np.stack((obs, obs_char), axis=0)
            met = met[:, np.argsort(met[0])]
            met[1] = np.cumsum(met[1])
            met[0] -= self.data.catalog.t_slew.loc[mask]
            met[1] += met[0]
            met = met / np.arange(1, np.count_nonzero(mask) + 1, 1)[np.newaxis, :]

            return met

    def observe_star(self,
                     nstar,
                     int_time):
        mask = self.data.catalog.nstar == nstar

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

        self.data.catalog.loc[mask, 'snr_current'] = np.sqrt(
            self.data.catalog.loc[mask, 'snr_current'] ** 2
            + (self.data.catalog.loc[mask, 'snr_1h']
               * np.sqrt(int_actual
                         / (60 * 60)))**2)

        # self.data.catalog.loc[mask, 'snr_char_current'] = np.sqrt(
        #     self.data.catalog.loc[mask, 'snr_char_current'] ** 2
        #     + (self.data.catalog.loc[mask, 'maxsep_snr_1h']
        #        * np.sqrt(int_char_time
        #                  / (60 * 60))) ** 2)

        det_ids = []

        for _, i in enumerate(np.where(mask)[0]):
            if (not self.data.catalog.detected.iloc[i]) and \
                    (self.data.catalog.snr_current.iloc[i]
                     >= self.data.options.optimization['snr_target']):
                self.data.catalog.detected.iat[i] = True
                det_ids.append(self.data.catalog.id.iat[i])
                if self.data.catalog.habitable.iloc[i]:
                    self.data.optm['sum_detected'][
                        np.where(np.array(list(self.data.options.optimization['limit'].keys()))
                                 == self.data.catalog.stype.iloc[i])] += 1
        self.id_series.append(det_ids)

    def distribute_time(self):
        if not self.data.options.optimization['habitable']:
            raise ValueError('Characterization optimization only implemented for habitable planets.')
        stars, n = np.unique(ar=self.data.catalog.nstar,
                             return_counts=True)
        obs = np.zeros((stars.shape[0], 2, np.max(n))) + np.inf

        # fill the observation time array
        for i, nstar in enumerate(stars):
            temp = self.obs_array_star(nstar=nstar)
            obs[i, :, :temp.shape[1]] = temp

        obs_time = (self.data.options.optimization['t_search']
                    * self.data.options.array['t_efficiency'])

        tot_time = 0
        if self.data.options.optimization['verbose']:
            print('Number of planets detected by stellar type:')

        while tot_time < obs_time:
            # find the best global slope and observe star
            no_star, ind_t = np.unravel_index(np.argmin(obs[:, 1, :]), obs[:, 1, :].shape)
            if (tot_time + obs[no_star, 0, ind_t] * (ind_t + 1) + 0.01) > obs_time:
                rem_time = obs_time - tot_time
                self.observe_star(nstar=stars[no_star],
                                  int_time=rem_time)
                tot_time += rem_time

            else:
                if ((obs[no_star, 1, ind_t] - obs[no_star, 0, ind_t]) * (ind_t + 1) + 0.01) < 0:
                    raise ValueError('Negative time difference encountered.')
                self.observe_star(nstar=stars[no_star],
                                  int_time=obs[no_star, 0, ind_t] * (ind_t + 1) + 0.01)
                temp = self.obs_array_star(nstar=stars[no_star])
                tot_time += obs[no_star, 0, ind_t] * (ind_t + 1) + 0.01
                obs[no_star, :, :] = np.inf
                obs[no_star, :, :temp.shape[1]] = temp

            out_string = ''
            for key in self.data.options.optimization['limit'].keys():
                temp_yield = (self.data.optm['sum_detected'] / self.data.optm['num_universe'])[
                    np.where(np.array(list(self.data.options.optimization['limit'].keys()))
                             == key)][0]
                out_string += (key + ': '
                               + str(temp_yield)
                               + '  ')
                self.yield_series[key].append(temp_yield)

            out_string += ('-  (' + str(np.round(tot_time/60/60/24/365.25, 1)) + ' / '
                           + str(np.round(obs_time/60/60/24/365.25, 1)) + ') yrs observed')

            self.t_series.append(tot_time)

            if self.data.options.optimization['verbose']:
                print('\r' + out_string, end='')

            if np.any(
                    np.logical_and(
                        (self.data.optm['sum_detected'] / self.data.optm['num_universe'])
                        > np.array(list(self.data.options.optimization['limit'].values())),
                        np.invert(self.data.optm['hit_limit']))):
                if self.data.options.optimization['verbose']:
                    print('\n')
                    print('HIT LIMIT, RECOUNTING -------------------')
                self.data.optm['hit_limit'] = ((self.data.optm['sum_detected']
                                                / (self.data.optm['num_universe']))
                                               >= np.array(
                            list(self.data.options.optimization['limit'].values())
                        ))
                obs = np.zeros((stars.shape[0], 2, np.max(n))) + np.inf

                # fill the observation time array
                for i, nstar in enumerate(stars):
                    temp = self.obs_array_star(nstar=nstar)
                    obs[i, :, :temp.shape[1]] = temp
