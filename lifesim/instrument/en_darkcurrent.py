import numpy as np
from typing import Union

from lifesim.core.modules import ElectronNoiseDetectorModule


class ElectronNoiseDarkCurrent(ElectronNoiseDetectorModule):
    """
    This class simulates the dark current noise contribution of the detector to the interferometric
    measurement of LIFE.
    """

    def __init__(self,
                 name: str):
        super().__init__(name=name)
        """
        Parameters
        ----------
        name : str
            Name of the module.
        """

    def noise(self,
              index: Union[int, type(None)]):
        """
        Simulates the amount of electron noise originating from the dark current of the detector
        leaking into the LIFE array measurement.

        Parameters
        ----------
        index: Union[int, type(None)]
            Specifies the planet for which to calculate the noise contribution. If an integer n is
            given, the noise will be calculated for the n-th row in the `data.catalog`. If `None`
            is given, the noise is caluculated for the parameters located in `data.single`.

        Returns
        -------
        dc_leak
            Dark current leakage of the detector in [electron s-1] per wavelength bin.

        Notes
        -----
        All of the following parameters are needed for the calculation of the dark current noise
        contribution and should be specified either in `data.catalog` or `data.single` or 'data.inst' or 'data.options'.

        data.inst['wl_bins'] : np.ndarray
            Central values of the spectral bins in the wavelength regime in [m].
        data.options.array['dc_per_pix'] : float
            Dark current per pixel in [electron s-1 px-1].
        """

        # read data on detector
        dc_per_pix = self.data.options.array['dc_per_pix']
        total_pixels = self.data.options.array['pix_per_wl'] * len(self.data.inst['wl_bins']) # minimum number of detector pixels (nyquist rate)

        # calculate total dark current noise
        dc_leak = np.full(self.data.inst['wl_bins'].shape, dc_per_pix * total_pixels)

        return dc_leak