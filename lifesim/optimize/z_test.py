import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from lifesim.core.modules import IntegrationTimeModule

matplotlib.use('Qt5Agg')


class ZTestModule(IntegrationTimeModule):
    def __init__(self,
                 name: str):
        super().__init__(name=name)

    def times_for_all(self):
        for i in range(self.data.stars.shape[0]):
            pass


    def times_for_single(self,
                         nstar: int,
                         order: int):
        snrs = []
        not_present = 0
        for i in range(self.data.statistic.nuniverse.max()+1):
            temp = self.data.statistic.loc[
                np.logical_and.reduce((self.data.statistic.nstar == nstar, self.data.statistic.nuniverse == i)),
                't_detection'].sort_values()
            # temp = self.data.statistic.loc[
            #     np.logical_and.reduce((self.data.statistic.nstar == nstar, self.data.statistic.nuniverse == i)),
            #     'snr_1h'].sort_values()
            if len(temp) < order:
                not_present += 1
            else:
                # snrs.append(temp.iloc[-order])
                snrs.append(temp.iloc[order - 1])
        pass