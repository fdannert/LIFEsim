import numpy as np

from lifesim.core.modules import RealisticModule


class BlindSimulationModule(RealisticModule):
    def __init__(self,
                 name: str):
        super().__init__(name=name)

    def run_simulation(self):
        self.data.optm['t_current'] = 0.
        self.data.optm['array_ang_current'] = 0.
        self.data.stars['vector'] = None

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
