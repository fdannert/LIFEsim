import numpy as np
from typing import Union

from lifesim.core.modules import PhotonNoiseUniverseModule
from lifesim.util import constants
from lifesim.util.radiation import black_body


class PhotonNoiseExozodi(PhotonNoiseUniverseModule):
    """
    This class simulates the noise contribution of an exozodi disk to the interferometric
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
        Simulates the amount of photon noise originating from the exozodi of the observed system
        leaking into the LIFE array measurement.

        Parameters
        ----------
        index: Union[int, type(None)]
            Specifies the planet for which to calculate the noise contribution. If an integer n is
            given, the noise will be calculated for the n-th row in the `data.catalog`. If `None`
            is given, the noise is caluculated for the parameters located in `data.single`.

        Returns
        -------
        ez_leak
            Exozodi leakage in [photon s-1] per wavelength bin.

        Notes
        -----
        All of the following parameters are needed for the calculation of the exozodi noise
        contribution and should be specified either in `data.catalog` or in `data.single`.

        l_sun : float
            Luminosity of the observed star in [solar luminosities].
        distance_s : float
            Distance between the observed star and the LIFE array in [pc].
        z : float
            Zodi level in the observed system in [zodis].
        mas_pix : np.ndarray
            Contains the size of each pixel projected to the sky in [milliarcseconds].
        rad_pix : np.ndarray
            Contains the size of each pixel projected to the sky in [radians].
        data.inst['radius_map'] : np.ndarray
            Contains the distance of a pixel from the center of the detector in [pix].
        data.options.other['image_size']
            Number of pixels on one axis of a square detector (dimensionless). I.e. for a 512x512
            detector this value is 512.
        wl_bins : np.ndarray
            Central values of the spectral bins in the wavelength regime in [m].
        wl_widths : np.ndarray
            Widths of the spectral wavelength bins in [m].
        data.inst['telescope_area'] : float
            Area of all array apertures combined in [m^2].
        data.inst['t_map'] : np.ndarray
            Transmission map of the TM3 mode of the array created by the
            lifesim.TransmissionMap module.
        """

        # read from catalog or single data depending on index specification
        if index is None:
            l_sun = self.data.single['l_sun']
            distance_s = self.data.single['distance_s']
            z = self.data.single['z']
        else:
            l_sun = self.data.catalog.l_sun.iloc[index]
            distance_s = self.data.catalog.distance_s.iloc[index]
            z = 1

        # calculate the parameters required by Kennedy2015
        alpha = 0.34
        r_in = 0.034422617777777775 * np.sqrt(l_sun)
        r_0 = np.sqrt(l_sun)
        sigma_zero = 7.11889e-8  # Sigma_{m,0} from Kennedy+2015 (doi:10.1088/0067-0049/216/2/23)

        # reshape the mas per pixel array for calculation (to (n, 1, 1))
        mas_pix = np.array([self.data.inst['mas_pix']])
        if mas_pix.shape[-1] > 1:
            mas_pix = np.reshape(mas_pix, (mas_pix.shape[-1], 1, 1))
        rad_pix = np.array([self.data.inst['rad_pix']])
        if rad_pix.shape[-1] > 1:
            rad_pix = np.reshape(rad_pix, (rad_pix.shape[-1], 1, 1))

        au_pix = mas_pix / 1e3 * distance_s

        # the radius as measured from the central star for every pixel in [AU]
        r_au = self.data.inst['radius_map'] * au_pix

        # identify all pixels where the radius is larges than the inner radius by Kennedy+2015
        r_cond = ((r_au >= r_in)
                  & (r_au <= self.data.options.other['image_size'] / 2 * au_pix))

        # calculate the temperature at all pixel positions according to Kennedy2015 Eq. 2
        temp_map = np.where(r_cond,
                            278.3 * (l_sun ** 0.25) / np.sqrt(r_au), 0)

        # calculate the Sigma (Eq. 3) in Kennedy2015 and set everything inside the inner radius to 0
        sigma = np.where(r_cond,
                         sigma_zero * z *
                         (r_au / r_0) ** (-alpha), 0)

        wl_bins = np.array([self.data.inst['wl_bins']])
        if wl_bins.shape[-1] > 1:
            wl_bins = np.reshape(wl_bins, (wl_bins.shape[-1], 1, 1))

        wl_bin_widths = np.array([self.data.inst['wl_bin_widths']])
        if wl_bin_widths.shape[-1] > 1:
            wl_bin_widths = np.reshape(wl_bin_widths, (wl_bin_widths.shape[-1], 1, 1))

        # get the black body radiation emitted by the interexoplanetary dust
        f_nu_disk = black_body(bins=wl_bins,
                               width=wl_bin_widths,
                               temp=temp_map,
                               mode='wavelength') \
                    * sigma * rad_pix ** 2 * self.data.inst['telescope_area']

        ap = np.where(self.data.inst['radius_map']
                      <= self.data.options.other['image_size'] / 2, 1, 0)
        # add the transmission map
        ez_leak = (f_nu_disk * self.data.inst['t_map'] * ap).sum(axis=(-2, -1))

        return ez_leak
