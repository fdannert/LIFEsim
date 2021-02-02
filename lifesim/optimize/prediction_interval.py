import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from lifesim.core.modules import IntegrationTimeModule

matplotlib.use('Qt5Agg')


class PredictionIntervalModule(IntegrationTimeModule):
    def __init__(self,
                 name: str):
        super().__init__(name=name)

    def times_for_all(self):
        for i in range(self.data.stars.shape[0]):
            pass

    def times_for_single(self,
                         nstar: int,
                         order: int):
        times = []
        not_present = 0
        for i in range(self.data.statistic.nuniverse.max()+1):
            temp = self.data.statistic.loc[
                np.logical_and.reduce((self.data.statistic.nstar == nstar, self.data.statistic.nuniverse == i)),
                't_detection'].sort_values()
            if len(temp) < order:
                not_present += 1
            else:
                times.append(temp.iloc[order - 1])
        n_1 = int(len(times) + 1)
        coeff_matrix = np.array(((1, 1), (1, -1)))
        ord_values = np.array((n_1, self.data.options.optimization['stat_size']*n_1))
        res = np.linalg.solve(a=coeff_matrix,
                              b=ord_values)
        res = np.sort(np.round(res, decimals=0).astype(int))
        times = np.sort(np.array(times))
        median_time_p = np.median(times)
        time_p = times[res[1]]

        self.data.stars.loc[self.data.stars.nstar == nstar, 'median_time_p'] = median_time_p
        self.data.stars.loc[self.data.stars.nstar == nstar, 'time_p'] = time_p