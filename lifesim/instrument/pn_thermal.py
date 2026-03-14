import numpy as np
from typing import Union
from lifesim.util import constants

from lifesim.core.modules import PhotonNoiseInstrumentModule
from lifesim.util.radiation import black_body


class PhotonNoiseThermal(PhotonNoiseInstrumentModule):
    """
    This class simulates the thermal noise contribution of the mirror to the interferometric
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
        Simulates the amount of photon noise originating from the thermal emission of the mirror and detector
        leaking into the LIFE array measurement.

        Parameters
        ----------
        index: Union[int, type(None)]
            Specifies the planet for which to calculate the noise contribution. If an integer n is
            given, the noise will be calculated for the n-th row in the `data.catalog`. If `None`
            is given, the noise is caluculated for the parameters located in `data.single`.

        Returns
        -------
        tm_leak
            Thermal leakage of the mirror in [photon s-1] per wavelength bin.
        td_leak
            Thermal leakage of the detector in [photon s-1] per wavelength bin.

        Notes
        -----
        All of the following parameters are needed for the calculation of the thermal mirror and detector noise
        contribution and should be specified either in `data.catalog` or `data.single` or 'data.inst'.

        data.inst['wl_bins'] : np.ndarray
            Central values of the spectral bins in the wavelength regime in [m].
        data.inst['wl_widths'] : np.ndarray
            Widths of the spectral wavelength bins in [m].
        data.inst['telescope_area'] : float
            Area of all array apertures combined in [m^2].
    
        """
        # read data on mirror
        mirror_temp = 40
        mirror_area = self.data.inst['telescope_area']
        emissivity = 0.9

        # calculate noise from the mirror
        # calculate the black body radiation emitted by the mirror
        # emissivity per wavelength bin?
        mirror_bb = emissivity * black_body(mode='wavelength',
                                            bins=self.data.inst['wl_bins'],
                                            width=self.data.inst['wl_bin_widths'],
                                            temp=mirror_temp)

        # integrate over area and solid angle
        tm_leak = mirror_bb * mirror_area * np.pi

        # calculate noise from the detector WIP
        detector_temp = 23
        td_leak = np.full(self.data.inst['wl_bins'].shape, detector_temp, dtype=float)


        return tm_leak
