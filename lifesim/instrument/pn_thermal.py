import numpy as np
from typing import Union

from lifesim.core.modules import PhotonNoiseInstrumentModule
from lifesim.util.radiation import black_body


class PhotonNoiseThermal(PhotonNoiseInstrumentModule):
    """
    This class simulates the thermal noise contribution of the mirror and detector to the interferometric
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
        contribution and should be specified either in 'data.inst' or 'data.options'.

        data.inst['hfov'] : np.ndarray
            Contains the half field of view of the observatory in [rad] for each of the spectral bins.
        data.inst['wl_bins'] : np.ndarray
            Central values of the spectral bins in the wavelength regime in [m].
        data.inst['wl_widths'] : np.ndarray
            Widths of the spectral wavelength bins in [m].
        data.options.array['primary_temp'] : float
            Temperature of the mirror in [K].
        data.options.array['primary_emissivity'] : float
            Emissivity of the mirror (dimensionless).
        data.inst['telescope_area'] : float
            Area of all array apertures combined in [m^2].
        data.options.array['num_apertures'] : int
            Number of apertures in the array.
        data.options.array['pixel_size'] : float
            Size of the pixels in [m]. (length of one side of the square pixel)
        data.options.array['pix_per_wl'] : int
            Number of pixels per wavelength bin (Nyquist rate).
        data.options.array['detector_wl_min'] : float
            Minimum wavelength of the detector sensitivity range in [m].
        data.options.array['detector_wl_max'] : float
            Maximum wavelength of the detector sensitivity range in [m].
        data.options.array['d_temp'] : float
            Temperature of the detector environment in [K].
        """

        # solid angle is governed by the fiber pick-up, which for single mode is lambda / D
        solid_angle = np.pi * (self.data.inst['hfov'])**2
        
        # calculate noise from the mirror
        mirror_bb = black_body(mode='wavelength',
                                            bins=self.data.inst['wl_bins'],
                                            width=self.data.inst['wl_bin_widths'],
                                            temp=self.data.options.array['primary_temp'])

        tm_leak = (solid_angle
                   * self.data.options.array['primary_emissivity']
                   * self.data.inst['telescope_area'] / self.data.options.array['num_apertures']
                   * mirror_bb)

        # detector collects thermal noise photons across its whole sensitivity range (at least from the detector
        # housing). Define temporary wl bins. Delta_wl is chosen to be small enough to capture the shape of the black
        # body curve and does not need to be adjusted

        delta_wl = 1e-7
        total_area = self.data.options.array['pixel_size'] ** 2 * self.data.options.array['pix_per_wl'] # minimum number of detector pixels (nyquist rate)
        solid_angle = np.pi  # half sphere, considering angle relevant to the normal of the detector surface for the irradiance and spherical coordinates

        wl_bins = np.arange(self.data.options.array['detector_wl_min'],
                            self.data.options.array['detector_wl_max'],
                            step=delta_wl)
        wl_bin_widths = np.full_like(wl_bins, delta_wl)

        # calculate noise from the detector
        detector_bb = black_body(mode='wavelength',
                                 bins=wl_bins,
                                 width=wl_bin_widths,
                                 temp=self.data.options.array['d_temp']) / wl_bin_widths

        # integral over all wavelengths
        detector_bb_int = np.trapz(y=detector_bb, x=wl_bins)

        td_leak = solid_angle * total_area * detector_bb_int * np.ones_like(self.data.inst['wl_bins'])

        return tm_leak, td_leak
