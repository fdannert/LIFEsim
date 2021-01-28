import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from lifesim.core.modules import InterestModule

matplotlib.use('Qt5Agg')


class MultiInterestModule(InterestModule):

    def __init__(self,
                 name: str):
        super().__init__(name=name)

    def get_star(self):
        pass

    def in_for(self):
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

        s_tp = np.array((np.sin(np.pi/2 - self.data.stars.loc[mask, 'lat'].iloc[0])
                         * np.cos(self.data.stars.loc[mask, 'lon'].iloc[0] - ang_tp),
                         np.sin(np.pi / 2 - self.data.stars.loc[mask, 'lat'].iloc[0])
                         * np.sin(self.data.stars.loc[mask, 'lon'].iloc[0] - ang_tp),
                         np.cos(np.pi / 2 - self.data.stars.loc[mask, 'lat'].iloc[0])))

        if (abs(self.data.stars.loc[mask, 'lat']) > self.data.options.array['ang_low']).iloc[0]:
            # CASE B
            #
            # decide on optimization if angular distance between star after tp and anti-sun
            # direction is less than the outer FOR angle and the star after tp is 'on the right'
            # of the anti-sun vector
            if ((np.dot(s_tp, pointing) < self.data.options.array['ang_high'])
                    and not self.between(low=self.data.stars.loc[mask, 'lon'].iloc[0] - ang_tp,
                                         high=self.data.stars.loc[mask, 'lon'].iloc[0],
                                         mid=self.data.optm['array_ang_current']
                                             - np.sqrt(self.data.options.array['ang_high'] ** 2
                                                       - self.data.stars.loc[mask, 'lat'].iloc[0] ** 2)
                                         )):
                # CASE B.1 or B.2

                # if self.between(low=self.data.optm['array_ang_current'] - ang_tp / 2,
                #                 high=self.data.optm['array_ang_current'] + ang_tp / 2,
                #                 mid=self.data.stars.loc[mask, 'lon'].iloc[0]):
                if self.between(low=self.data.stars.loc[mask, 'lon'].iloc[0] - ang_tp / 2,
                                high=self.data.stars.loc[mask, 'lon'].iloc[0],
                                mid=self.data.optm['array_ang_current']):
                    # CASE B.2
                    print('B.2')
                    tp_0 = self.data.optm['t_current']
                else:
                    # CASE B.1
                    print('B.1')
                    ang_s_p = np.arccos(
                        np.dot(np.array((np.cos(self.data.optm['array_ang_current']),
                                         np.sin(self.data.optm['array_ang_current']),
                                         0)),
                               np.array((np.cos(self.data.stars.loc[mask, 'lon'].iloc[0]),
                                         np.sin(self.data.stars.loc[mask, 'lon'].iloc[0]),
                                         0))))
                    tp_0 = self.data.optm['t_current'] + ang_s_p * 60 * 60 * 24 * 365 / (2 * np.pi)

            else:
                # CASE B.3
                print('B.3')
                tp_0 = None

        else:
            # case A or case C
            if np.cross(np.array((np.cos(self.data.stars.loc[mask, 'lon'].iloc[0]),
                                  np.sin(self.data.stars.loc[mask, 'lon'].iloc[0]),
                                  0)),
                        np.array((np.cos(self.data.optm['array_ang_current']),
                                  np.sin(self.data.optm['array_ang_current']),
                                  0)))[-1] > 0:
                # CASE C
                #
                # set to t_current if angular distance between star after tp and anti-sun direction
                # is less than the outer FOR angle and the star after tp is 'on the right' of the
                # anti-sun vector
                if ((np.dot(s_tp, pointing) < self.data.options.array['ang_high'])
                        and not self.between(low=self.data.stars.loc[mask, 'lon'].iloc[0] - ang_tp,
                                         high=self.data.stars.loc[mask, 'lon'].iloc[0],
                                         mid=self.data.optm['array_ang_current']
                                             - np.sqrt(self.data.options.array['ang_high'] ** 2
                                                       - self.data.stars.loc[mask, 'lat'].iloc[
                                                           0] ** 2)
                                         )):
                    # CASE C.1
                    print('C.1')
                    tp_0 = self.data.optm['t_current']
                else:
                    # CASE C.3
                    print('C.3')
                    tp_0 = None

            else:
                # CASE A
                #
                # set to t_current if angular distance between star after tp and anti-sun direction
                # is more than the inner FOR angle and the star after tp is 'on the left' of the
                # anti-sun vector
                if ((np.dot(s_tp, pointing) > self.data.options.array['ang_low'])
                        and not self.between(low=self.data.stars.loc[mask, 'lon'].iloc[0] - ang_tp,
                                         high=self.data.stars.loc[mask, 'lon'].iloc[0],
                                         mid=self.data.optm['array_ang_current']
                                             - np.sqrt(self.data.options.array['ang_high'] ** 2
                                                       - self.data.stars.loc[mask, 'lat'].iloc[
                                                           0] ** 2)
                                         )):
                    # CASE A.1
                    print('A.1')
                    kappa = np.arccos(self.data.options.array['ang_low']
                                      / np.sin(np.pi / 2 - self.data.stars.loc[mask, 'lat']))
                    ang_s_p = np.arccos(np.dot(
                        np.array((np.cos(self.data.optm['array_ang_current']),
                                  np.sin(self.data.optm['array_ang_current']),
                                  0)),
                        np.array((np.cos(self.data.stars.loc[mask, 'lon']),
                                  np.sin(self.data.stars.loc[mask, 'lon']),
                                  0))
                    ))
                    ang_s_start = self.data.optm['array_ang_current'] + kappa + ang_tp
                    tp_0 = ((ang_s_p - ang_s_start) * 60 * 60 * 24 * 365 / (2 * np.pi)
                            + self.data.optm['t_current'])
                else:
                    # CASE A.3
                    print('A.3')
                    tp_0 = None
        print(tp_0)
        self.make_plot(nstar=nstar, ang=ang_tp, t_0=tp_0)
        pass

    def between(self,
                low: float,
                high: float,
                mid: float):
        end = (high - low) % (2 * np.pi)
        mid = (mid - low) % (2 * np.pi)

        return mid < end

    def make_plot(self, nstar, ang, t_0):
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
            plt.scatter(self.data.stars.loc[mask, 'lon'].iloc[0] - t_0 / (60*60*24*365) * 2 * np.pi,
                        self.data.stars.loc[mask, 'lat'].iloc[0], marker='*')
            plt.scatter(self.data.stars.loc[mask, 'lon'].iloc[0] - 2 * np.pi - t_0 / (60*60*24*365) * 2 * np.pi,
                        self.data.stars.loc[mask, 'lat'].iloc[0], marker='*')
        plt.show()
