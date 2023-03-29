import numpy as np
from tqdm import tqdm
import pandas as pd

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

        self.data.catalog = pd.merge(left=self.data.catalog,
                                     right=np.sqrt((
                                                           (self.data.noise_catalog.sel(params='signal')
                                                            / self.data.noise_catalog.sel(params='fundamental')) ** 2
                                                   ).sum(axis=1) * 3600 / self.data.options.array['t_rot']
                                                   ).to_pandas().rename('snr_1h_prt'),
                                     left_on='id',
                                     right_on='ids')

    def planet_count(self):
        pass