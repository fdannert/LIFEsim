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
        data.options.array['diameter'] : float
            Diameter of the array in [m].
        data.options.array['m_temp'] : float
            Temperature of the mirror in [K].
        data.options.array['m_emissivity'] : float
            Emissivity of the mirror (dimensionless).
        data.options.array['d_temp'] : float
            Temperature of the detector in [K].
        data.options.array['pixel_size'] : float
            Size of the pixels in [m]. (length of one side of the square pixel)
        data.options.array['beam_size'] : float
            Diameter of the beam in [m].
        """

        # read data on mirror
        mirror_emissivity = self.data.options.array['m_emissivity']
        mirror_temp = self.data.options.array['m_temp']
        mirror_area = np.pi *(self.data.options.array['diameter'] / 2.) ** 2
        beam_radius = self.data.options.array['beam_size']/2
        distance = 2.5 * self.data.options.array['diameter']
        solid_angle = (np.pi * beam_radius ** 2) / (distance ** 2)
        
        # calculate noise from the mirror
        mirror_bb = black_body(mode='wavelength',
                                            bins=self.data.inst['wl_bins'],
                                            width=self.data.inst['wl_bin_widths'],
                                            temp=mirror_temp)

        tm_leak = mirror_emissivity * mirror_bb * mirror_area * solid_angle


        # read data on detector
        detector_temp = self.data.options.array['d_temp']
        pixel_area = self.data.options.array['pixel_size'] ** 2
        total_area = pixel_area * 2 * len(self.data.inst['wl_bins']) # minimum number of detector pixels (nyquist rate)

        # calculate noise from the detector
        detector_bb = black_body(mode='wavelength',
                                   bins=self.data.inst['wl_bins'],
                                   width=self.data.inst['wl_bin_widths'],
                                   temp=detector_temp)
        
        td_leak = np.pi * total_area * detector_bb

        return tm_leak, td_leak
