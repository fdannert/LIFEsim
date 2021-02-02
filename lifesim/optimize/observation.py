import warnings

import numpy as np
import pandas as pd

from lifesim.core.modules import ObservationModule


class BLAdjustObservationModule(ObservationModule):
    def __init__(self,
                 name: str):
        super().__init__(name=name)

    def observe(self,
                nstar: int):
        pass


class SimpleObservationModule(ObservationModule):
    def __init__(self,
                 name: str):
        super().__init__(name=name)

    def observe(self,
                nstar: int):
        # just to be sure, a for check
        mask_s = self.data.stars.nstar == nstar
        mask_c = self.data.catalog.nstar == nstar
        pointing = np.array((np.cos(self.data.optm['array_ang_current']),
                             np.sin(self.data.optm['array_ang_current']),
                             0))
        angle = np.arccos(np.dot(pointing, self.data.stars.loc[mask_s, 'vector'].iloc[0]))
        if not ((angle > self.data.options.array['ang_low'])
                and (angle < self.data.options.array['ang_high'])):
            warnings.warn('Star not in field of regard')
            return 0

        target_n = self.data.stars.loc[mask_s, 'i_number'].iloc[0]
        target_t = self.data.stars.loc[mask_s, 'time_p'].iloc[0]

        if mask_c.sum() == 0:
            self.data.stars.loc[mask_s, 'i_efficiency'] = 0.
            return target_t

        self.data.stars.loc[mask_s, 'lon'].iat[0] = (self.data.stars.loc[mask_s, 'lon'].iloc[0]
                                                     - target_t / (60 * 60 * 24 * 365) * 2 * np.pi)

        for _, n_p in enumerate(np.argwhere(
                self.data.catalog.nstar.to_numpy() == nstar)[:, 0]):
            theta_p = (self.data.catalog.theta_p.iloc[n_p]
                       + 2 * np.pi / self.data.catalog.p_orb.iloc[n_p]) % (2 * np.pi)
            self.data.catalog.angsep.iat[n_p] = (self.data.catalog.semimajor_p.iloc[n_p]
                                                 / self.data.catalog.distance_s.iloc[n_p]
                                                 * np.sqrt(
                        np.cos(self.data.catalog.small_omega_p.iloc[n_p] + theta_p) ** 2
                        + np.cos(self.data.catalog.inc_p.iloc[n_p]) ** 2
                        * np.sin(self.data.catalog.small_omega_p.iloc[n_p] + theta_p) ** 2))

        self.run_socket(s_name='instrument',
                        method='get_snr_single_star',
                        index=np.where(mask_c.to_numpy())[0][0],
                        integration_time=60)

        temp = pd.DataFrame(columns=['n_p', 'time'])
        for _, n_p in enumerate(np.argwhere(
                self.data.catalog.nstar.to_numpy() == nstar)[:, 0]):
            self.data.catalog.t_to_dect.iat[n_p] = (
                    60 * 60
                    * (self.data.options.optimization['snr_target'] ** 2
                       - self.data.catalog.snr_current.iloc[n_p] ** 2)
                    / self.data.catalog.snr_1h.iloc[n_p] ** 2)
            temp = temp.append({'n_p': n_p,
                                'time': self.data.catalog.t_to_dect.iloc[n_p]},
                               ignore_index=True)

        temp = temp.sort_values(by='time')
        if not (len(temp) < target_n) and (temp.time.iloc[target_n - 1] < target_t):
            t_integration = temp.time.iloc[target_n - 1]
        else:
            t_integration = target_t

        for _, n_p in enumerate(np.argwhere(
                self.data.catalog.nstar.to_numpy() == nstar)[:, 0]):
            self.data.catalog.snr_current.iat[n_p] = np.sqrt(
                self.data.catalog.snr_current.iat[n_p] ** 2
                + t_integration/(60 * 60) * self.data.catalog.snr_1h.iloc[n_p] ** 2)
            if ((self.data.catalog.snr_current.iloc[n_p]
                >= self.data.options.optimization['snr_target'])
                    and not self.data.catalog.detected.iloc[n_p]):
                self.data.catalog.detected.iat[n_p] = True

        self.data.stars.loc[mask_s, 'i_efficiency'] = 0.
        pass
                # n = self.data.stars.loc[mask_s, 'number_dect'].iloc[0] + 1
                # self.data.stars.loc[self.data.stars.nstar == nstar, 'number_dect'].iat[0] = n
                # self.data.discover.append({'name': ,
                #                           'discovery_epoch':, 't_detection', 'snr_current'})

        return t_integration

            # noise_bg = 0.
            # for name in self.data.catalog.noise_diff.iloc[n_p].keys():
            #     if name != "<class 'lifesim.instrument.pn_localzodi.PhotonNoiseLocalzodi'>":
            #         noise_bg += self.data.catalog.noise_diff.iloc[n_p][name]
            #
            # self.adjust_bl_to_hz(hz_center=float(self.data.catalog.hz_center.iloc[n_p]),
            #                      distance_s=float(self.data.catalog.distance_s.iloc[n_p]))
            #
            # _, _, self.data.inst['t_map'], _, _ = self.run_socket(s_name='transmission',
            #                                                       method='transmission_map',
            #                                                       map_selection='tm3')
            #
            # noise_lz = self.run_socket(s_name='localzodi',
            #                            method='noise',
            #                            index=n_p)




