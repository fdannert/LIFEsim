import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from tqdm import tqdm

from lifesim.core.modules import InterestModule
import lifesim.util.constants as c

matplotlib.use('Qt5Agg')


class MultiInterestModule(InterestModule):

    def __init__(self,
                 name: str):
        super().__init__(name=name)

    def get_star(self,
                 nstar: int):
        mask = self.data.stars.nstar == nstar
        if self.data.stars.loc[mask, 'in_for'].iloc[0]:
            self.get_efficiency(nstar=nstar)
            self.run_socket(s_name='int_time',
                            method='times_for_single',
                            nstar=nstar,
                            order=self.data.stars.loc[mask, 'i_number'].iloc[0])
            self.get_start_time(nstar=nstar)
            if self.data.stars.loc[mask, 'time_p_0'].iloc[0] is not None:
                self.get_time_interest(nstar=nstar)
                self.get_localzodi(nstar=nstar)
                self.data.stars.loc[mask, 'interest_current'] = (
                        self.data.stars.loc[mask, 'i_efficiency'].iloc[0]
                        * self.data.stars.loc[mask, 'i_time'].iloc[0]
                        * self.data.stars.loc[mask, 'i_localzodi'].iloc[0])
                self.data.stars.loc[mask, 'interest'].append(
                    self.data.stars.loc[mask, 'interest_current'].iloc[0])
            else:
                self.data.stars.loc[mask, 'interest_current'] = 0.
                self.data.stars.loc[mask, 'interest'].append(0.)
        else:
            self.data.stars.loc[mask, 'interest_current'] = 0.
            self.data.stars.loc[mask, 'interest'].append(0.)

    def in_for(self):
        self.data.stars.in_for = False
        pointing = np.array((np.cos(self.data.optm['array_ang_current']),
                             np.sin(self.data.optm['array_ang_current']),
                             0))
        for i in range(self.data.stars.shape[0]):
            angle = np.arccos(np.dot(pointing, self.data.stars.vector.iloc[i]))
            self.data.stars.in_for.iat[i] = ((angle > self.data.options.array['ang_low'])
                                             and (angle < self.data.options.array['ang_high']))

    def get_start_time(self,
                       nstar: int):
        mask = self.data.stars.nstar == nstar
        ang_tp = self.data.stars.loc[mask, 'time_p'].iloc[0] * 2 * np.pi / (60*60*24*365)

        pointing = np.array((np.cos(self.data.optm['array_ang_current']),
                             np.sin(self.data.optm['array_ang_current']),
                             0))

        # beta_large check
        if self.between(low=self.data.stars.loc[mask, 'lon'].iloc[0] - ang_tp,
                        high=self.data.stars.loc[mask, 'lon'].iloc[0],
                        mid=self.data.optm['array_ang_current']
                            - np.sqrt(self.data.options.array['ang_high'] ** 2
                                      - self.data.stars.loc[mask, 'lat'].iloc[0] ** 2)):
            # CASE .3
            print('.3')
            tp_0 = None
        else:
            s_lon = np.array((np.cos(self.data.stars.loc[mask, 'lon'].iloc[0]),
                              np.sin(self.data.stars.loc[mask, 'lon'].iloc[0]),
                              0))

            # lon 'rot' check
            if np.cross(s_lon, pointing)[-1] > 0:
                # CASE D.1
                print('D.1')
                tp_0 = self.data.optm['t_current']

            else:
                # lat check
                if ((abs(self.data.stars.loc[mask, 'lat'])
                     > self.data.options.array['ang_low']).iloc[0]):
                    # CASE F
                    if self.between(low=self.data.stars.loc[mask, 'lon'].iloc[0] - ang_tp / 2,
                                    high=self.data.stars.loc[mask, 'lon'].iloc[0],
                                    mid=self.data.optm['array_ang_current']):
                        # CASE F.1
                        print('F.1')
                        tp_0 = self.data.optm['t_current']

                    else:
                        # CASE F.2
                        print('F.2')
                        tp_0 = (self.data.optm['t_current']
                                + np.arccos(np.dot(pointing, s_lon))
                                * 60 * 60 * 24 * 365 / (2 * np.pi)
                                - self.data.stars.loc[mask, 'time_p'].iloc[0] / 2)

                else:
                    # CASE E
                    # beta small check
                    if self.between(low=self.data.stars.loc[mask, 'lon'].iloc[0] - ang_tp,
                                    high=self.data.stars.loc[mask, 'lon'].iloc[0],
                                    mid=self.data.optm['array_ang_current']
                                        + np.sqrt(self.data.options.array['ang_low'] ** 2
                                                  - self.data.stars.loc[mask, 'lat'].iloc[0] ** 2)
                                    ):
                        # CASE E.3
                        print('E.3')
                        tp_0 = None

                    else:
                        # CASE E.1
                        print('E.1')
                        d_ang = (np.arccos(np.dot(pointing, s_lon))
                                 - ang_tp
                                 - np.sqrt(self.data.options.array['ang_low'] ** 2
                                           - self.data.stars.loc[mask, 'lat'].iloc[0] ** 2))
                        tp_0 = (self.data.optm['t_current']
                                + d_ang * 60 * 60 * 24 * 365 / (2 * np.pi))
        self.data.stars.loc[mask, 'time_p_0'] = tp_0
        if tp_0 is not None:
            self.data.stars.loc[mask, 'time_p_mid'] = tp_0 + self.data.stars.loc[mask, 'time_p']
        # print(tp_0)
        # self.make_plot(nstar=nstar, ang=ang_tp, t_0=tp_0, t_p=self.data.stars.loc[mask, 'time_p'].iloc[0])
        # pass

    def get_efficiency(self,
                       nstar: int):
        mask = self.data.statistic.nstar == nstar
        self.data.stars.loc[mask, 'i_number'] = int(np.round((
                self.data.statistic.loc[mask, 'detected'].sum()
                / (self.data.statistic.nuniverse.max() + 1)
        )))
        if self.data.stars.loc[mask, 'i_number'].iloc[0] == 0:
            self.data.stars.loc[mask, 'i_number'] = 1
        self.data.stars.loc[mask, 'i_efficiency'] = (
                self.data.statistic.loc[mask, 'detected'].sum()
                / self.data.statistic.loc[mask, 'int_time'].sum()
                * 60 * 60 * 24
        )

    def get_time_interest(self,
                          nstar: int):
        mask = self.data.stars.nstar == nstar
        self.data.stars.loc[mask, 'i_time'] = (
                (1 +
                 np.exp(
                     (-8 / self.data.options.optimization['time_scaler']
                      / self.data.stars.loc[mask, 'time_p'].iloc[0])
                     * (self.data.optm['t_current']
                        + self.data.options.optimization['time_scaler'] / 2
                        * self.data.stars.loc[mask, 'time_p'].iloc[0]
                        - self.data.stars.loc[mask, 'time_p_0'].iloc[0]))) ** (-1)
        )

    def get_localzodi(self,
                      nstar: int):
        mask = self.data.stars.nstar == nstar
        G = ((
                     (np.pi
                      / np.arccos(np.cos(self.data.stars.loc[mask, 'lon'].iloc[0]
                                         - self.data.optm['array_ang_current'])
                                  * np.cos(self.data.stars.loc[mask, 'lat'].iloc[0])))
                     /
                     (np.sin(self.data.stars.loc[mask, 'lat'].iloc[0]) ** 2
                      + 0.36
                      * (self.data.options.optimization['wl_optimal_lz'] / (11 * 1e-6)) ** (-0.8)
                      * np.cos(self.data.stars.loc[mask, 'lat'].iloc[0]))
             ) ** 0.5)

        G_max = ((np.pi / self.data.options.array['ang_high'])
                 / (0.36
                    * (self.data.options.optimization['wl_optimal_lz']
                       / (11 * 1e-6)) ** (-0.8))) ** 0.5

        self.data.stars.loc[mask, 'i_localzodi'] = (
                1
                - self.data.options.optimization['localzodi_scaler']
                * ((G - G_max) / (np.sqrt(2) - G_max))
        )


    def between(self,
                low: float,
                high: float,
                mid: float):
        end = (high - low) % (2 * np.pi)
        mid = (mid - low) % (2 * np.pi)

        return mid < end

    def make_plot(self, nstar, ang, t_0, t_p):
        mask = self.data.stars.nstar == nstar
        plt.subplot(111, projection='aitoff')
        # plt.subplot(111)
        plt.grid(True)
        # plt.scatter(gal.l.wrap_at('180d').radian, gal.b.radian)
        a = np.arange(start=0, stop=2*np.pi, step=0.1)
        lon = self.data.options.array['ang_low'] * np.cos(a) + self.data.optm['array_ang_current']
        lat = self.data.options.array['ang_low'] * np.sin(a)
        plt.plot(lon, lat, c='gray')
        lon = self.data.options.array['ang_high'] * np.cos(a) + self.data.optm['array_ang_current']
        lat = self.data.options.array['ang_high'] * np.sin(a)
        plt.plot(lon, lat, c='gray')
        plt.scatter(self.data.stars.loc[mask, 'lon'].iloc[0], self.data.stars.loc[mask, 'lat'].iloc[0])
        plt.scatter(self.data.stars.loc[mask, 'lon'].iloc[0] - 2 * np.pi,
                    self.data.stars.loc[mask, 'lat'].iloc[0])
        plt.scatter(self.data.stars.loc[mask, 'lon'].iloc[0] - ang,
                    self.data.stars.loc[mask, 'lat'].iloc[0])
        plt.scatter(self.data.stars.loc[mask, 'lon'].iloc[0] - 2 * np.pi - ang,
                    self.data.stars.loc[mask, 'lat'].iloc[0])
        if t_0 is not None:
            # plt.scatter(self.data.stars.loc[mask, 'lon'].iloc[0] - t_0 / (60*60*24*365) * 2 * np.pi,
            #             self.data.stars.loc[mask, 'lat'].iloc[0], marker='*')
            # plt.scatter(self.data.stars.loc[mask, 'lon'].iloc[0] - 2 * np.pi - t_0 / (60*60*24*365) * 2 * np.pi,
            #             self.data.stars.loc[mask, 'lat'].iloc[0], marker='*')
            plt.scatter([self.data.stars.loc[mask, 'lon'].iloc[0] - t_0 / (60*60*24*365) * 2 * np.pi,
                      self.data.stars.loc[mask, 'lon'].iloc[0] - (t_0 + t_p) / (60*60*24*365) * 2 * np.pi],
                     [self.data.stars.loc[mask, 'lat'].iloc[0],
                      self.data.stars.loc[mask, 'lat'].iloc[0]], marker='*')
            plt.scatter(
                [self.data.stars.loc[mask, 'lon'].iloc[0] - t_0 / (60 * 60 * 24 * 365) * 2 * np.pi - 2 * np.pi,
                 self.data.stars.loc[mask, 'lon'].iloc[0] - (t_0 + t_p) / (
                         60 * 60 * 24 * 365) * 2 * np.pi - 2 * np.pi],
                [self.data.stars.loc[mask, 'lat'].iloc[0],
                 self.data.stars.loc[mask, 'lat'].iloc[0]], marker='*')
        plt.show()


class SimpleInterestModule(InterestModule):

    def __init__(self,
                 name: str):
        super().__init__(name=name)

    def do_precalc(self):
        for _, nstar in enumerate(tqdm(self.data.stars.nstar)):
            self.precalc_single(nstar=nstar)

    def precalc_single(self,
                       nstar: int):
        mask = self.data.stars.nstar == nstar
        self.get_efficiency(nstar=nstar)
        self.run_socket(s_name='int_time',
                        method='times_for_single',
                        nstar=nstar,
                        order=self.data.stars.loc[mask, 'i_number'].iloc[0])

    def get_star(self,
                 nstar: int):
        mask = self.data.stars.nstar == nstar
        if self.data.stars.loc[mask, 'in_for'].iloc[0]:
            self.get_start_time(nstar=nstar)
            if self.data.stars.loc[mask, 'time_p_0'].iloc[0] is not None:
                self.get_localzodi(nstar=nstar)
                interest = (
                        self.data.stars.loc[mask, 'i_efficiency'].iloc[0]
                        * self.data.stars.loc[mask, 'i_localzodi'].iloc[0])
                if self.data.options.optimization['multi_visit']:
                    self.get_multi_visit(nstar=nstar)
                    interest *= self.data.stars.loc[mask, 'i_multi_visit'].iloc[0]
            else:
                interest = 0
        else:
            interest = 0
        self.data.stars.loc[mask, 'interest_current'] = interest
        self.data.stars.loc[mask, 'interest'].iat[0].append(interest)

    def in_for(self):
        self.data.stars.in_for = False
        pointing = np.array((np.cos(self.data.optm['array_ang_current']),
                             np.sin(self.data.optm['array_ang_current']),
                             0))
        for i in range(self.data.stars.shape[0]):
            angle = np.arccos(np.dot(pointing, self.data.stars.vector.iloc[i]))
            self.data.stars.in_for.iat[i] = ((angle > self.data.options.array['ang_low'])
                                             and (angle < self.data.options.array['ang_high']))

    def get_start_time(self,
                       nstar: int):
        mask = self.data.stars.nstar == nstar
        ang_tp = self.data.stars.loc[mask, 'time_p'].iloc[0] * 2 * np.pi / (60*60*24*365)

        # beta_large check
        if self.between(low=self.data.stars.loc[mask, 'lon'].iloc[0] - ang_tp,
                        high=self.data.stars.loc[mask, 'lon'].iloc[0],
                        mid=self.data.optm['array_ang_current']
                            - np.sqrt(self.data.options.array['ang_high'] ** 2
                                      - self.data.stars.loc[mask, 'lat'].iloc[0] ** 2)
                        ):
            # CASE .3
            tp_0 = None
        else:
            # lat check
            if ((abs(self.data.stars.loc[mask, 'lat'])
                 > self.data.options.array['ang_low']).iloc[0]):
                # CASE D/F.1
                tp_0 = self.data.optm['t_current']
            else:
                # beta small test
                if self.between(low=self.data.stars.loc[mask, 'lon'].iloc[0] - ang_tp,
                                high=self.data.stars.loc[mask, 'lon'].iloc[0],
                                mid=self.data.optm['array_ang_current']
                                    + np.sqrt(self.data.options.array['ang_low'] ** 2
                                              - self.data.stars.loc[mask, 'lat'].iloc[0] ** 2)
                                ):
                    # CASE E.3
                    tp_0 = None
                else:
                    # CASE E.1
                    tp_0 = self.data.optm['t_current']
        self.data.stars.loc[mask, 'time_p_0'] = tp_0
        if tp_0 is not None:
            self.data.stars.loc[mask, 'time_p_mid'] = tp_0 + self.data.stars.loc[mask, 'time_p']

    def get_efficiency(self,
                       nstar: int):
        mask = self.data.statistic.nstar == nstar
        self.data.stars.loc[mask, 'i_number'] = int(np.round((
                self.data.statistic.loc[mask, 'detected'].sum()
                / (self.data.statistic.nuniverse.max() + 1)
        )))
        if self.data.stars.loc[mask, 'i_number'].iloc[0] == 0:
            self.data.stars.loc[mask, 'i_number'] = 1
        if self.data.statistic.loc[mask, 'detected'].sum() != 0:
            self.data.stars.loc[mask, 'i_efficiency'] = (
                    self.data.statistic.loc[mask, 'detected'].sum()
                    / self.data.statistic.loc[mask, 'int_time'].sum()
                    * 60 * 60 * 24
            )
        else:
            self.data.stars.loc[mask, 'i_efficiency'] = 0

    def get_localzodi(self,
                      nstar: int):
        mask = self.data.stars.nstar == nstar
        ang_tp = self.data.stars.loc[mask, 'time_p'].iloc[0] * 2 * np.pi / (60 * 60 * 24 * 365)
        G = ((
                     (np.pi
                      / np.arccos(np.cos(self.data.stars.loc[mask, 'lon'].iloc[0] - ang_tp / 2
                                         - self.data.optm['array_ang_current'])
                                  * np.cos(self.data.stars.loc[mask, 'lat'].iloc[0])))
                     /
                     (np.sin(self.data.stars.loc[mask, 'lat'].iloc[0]) ** 2
                      + 0.36
                      * (self.data.options.optimization['wl_optimal_lz'] / (11 * 1e-6)) ** (-0.8)
                      * np.cos(self.data.stars.loc[mask, 'lat'].iloc[0]))
             ) ** 0.5)

        G_max = ((np.pi / self.data.options.array['ang_high'])
                 / (0.36
                    * (self.data.options.optimization['wl_optimal_lz']
                       / (11 * 1e-6)) ** (-0.8))) ** 0.5

        self.data.stars.loc[mask, 'i_localzodi'] = (
                1
                - self.data.options.optimization['localzodi_scaler']
                * ((G - G_max) / (np.sqrt(2) - G_max))
        )

    def get_multi_visit(self,
                        nstar: int):
        mask = self.data.stars.nstar == nstar
        if self.data.stars.loc[mask, 'frustration'].iloc[0] is None:
            I = 1
        elif self.data.stars.loc[mask, 'frustration'].iloc[0] != 0:
            p_orb_hz = 2 * np.pi * np.sqrt(
                (self.data.stars.loc[mask, 'hz_center'].iloc[0] * c.m_per_au) ** 3
                / (c.grav_const * self.data.stars.loc[mask, 'mass_s'].iloc[0] * c.mass_sun))

            I = (self.data.stars.loc[mask, 'frustration'].iloc[0]
                 * self.data.options.optimization['multi_scaler']
                 * np.sin(2 * np.pi / p_orb_hz
                          * (self.data.optm['t_current']
                             + self.data.stars.loc[mask, 'time_p'].iloc[0] / 2
                             - self.data.stars.loc[mask, 't_first'].iloc[0])
                          ) ** 2)
        else:
            I = 0

        self.data.stars.loc[mask, 'i_multi_visit'] = I

    def between(self,
                low: float,
                high: float,
                mid: float):
        end = (high - low) % (2 * np.pi)
        mid = (mid - low) % (2 * np.pi)

        return mid < end