import numpy as np

from lifesim.core.modules import AnalysisModule

class SampleAnalysisModule(AnalysisModule):
    def __init__(self,
                 name: str):
        """
        Parameters
        ----------
        name : str
            Name of the instrument module.
        """

        super().__init__(name=name)

    def get_fundamental_snr(self):
        if (self.data.catalog is None) or (self.data.noise_catalog is None):
            raise ValueError('Catalog and noise catalog need to be specified')

        self.data.catalog['snr_1h'] = np.zeros_like(self.data.catalog.nstar, dtype=float)

        for id in self.data.noise_catalog.keys():
            self.data.catalog.loc[self.data.catalog.id == int(id), 'snr_1h'] = np.sqrt(
                np.sum((self.data.noise_catalog[id]['signal']
                        /self.data.noise_catalog[id]['fundamental'])**2) *
                60 * 60 / self.data.options.array['t_rot'])


    def planet_count(self):
        pass