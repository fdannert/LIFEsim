from typing import Union

import numpy as np

from lifesim.core.modules import PhotonNoiseModule
from lifesim.util.radiation import black_body


class PhotonNoiseLocalzodi(PhotonNoiseModule):
    def __init__(self,
                 name: str):
        super().__init__(name=name)

    def noise(self,
              index: Union[int, type(None)]):
        """
        Simulates the amount of photon noise originating from the exozodi of the observed system
        leaking into the LIFE array measurement

        Parameters
        ----------
        image_size : int
            Number of pixels on one axis of a square detector (dimensionless). I.e. for a 512x512
            detector this value is 512
        telescope_area : float
            Area of all array apertures combined in [m^2]
        radius_map : np.ndarray
            Contains the distance of a pixel from the center of the detector in [pix]
        wl_bins : np.ndarray
            Central values of the spectral bins in the wavelength regime in [m]
        wl_bin_widths : np.ndarray
            Widths of the spectral wavelength bins in [m]
        t_map : np.ndarray
            Transmission map of the TM3 mode of the array created by the
            lifesim.TransmissionMap module
        lz_model : str
            Specifies which localzodi model will be used
        lat_s : str
            Ecliptic latitude of the observed star in [rad]
        hfov : np.ndarray
            Contains the half field of view of the observatory in [rad] for each of the spectral bins

        Returns
        -------
        lz_leak
            Localzodi leakage in [s^-1] per wavelength bin

        Raises
        ______

        ValueError
            If the specified localzodi model does not exits
        """
        # TODO Implement longitude dependence of localzodi
        # TODO Find model after which this is calculated and reference
        # TODO complete comments

        if index is None:
            lat_s = self.data.single['lat']
        else:
            lat_s = self.data.catalog.lat.iloc[index]

        # check if the model exists
        if not ((self.data.options.models['localzodi'] == 'glasse')
                or (self.data.options.models['localzodi'] == 'darwinsim')):
            raise ValueError('Specified model does not exist')

        long = 3 / 4 * np.pi
        lat = lat_s

        ap = np.where(self.data.inst['radius_map']
                      <= self.data.options.other['image_size'] / 2, 1, 0)

        if self.data.options.models['localzodi'] == 'glasse':
            temp = 270
            epsilon = 4.30e-8
            lz_flux_sr = epsilon * black_body(mode='wavelength',
                                              bins=self.data.inst['wl_bins'],
                                              width=self.data.inst['wl_bin_widths'],
                                              temp=temp)

        elif self.data.options.models['localzodi'] == 'darwinsim':
            radius_sun_au = 0.00465047  # in AU
            tau = 4e-8
            temp_eff = 265
            temp_sun = 5777
            a = 0.22

            b_tot = black_body(mode='wavelength',
                               bins=self.data.inst['wl_bins'],
                               width=self.data.inst['wl_bin_widths'],
                               temp=temp_eff) + a \
                    * black_body(mode='wavelength',
                                 bins=self.data.inst['wl_bins'],
                                 width=self.data.inst['wl_bin_widths'],
                                 temp=temp_sun) \
                    * (radius_sun_au / 1.5) ** 2
            lz_flux_sr = tau * b_tot * np.sqrt(
                np.pi / np.arccos(np.cos(long) * np.cos(lat)) /
                (np.sin(lat) ** 2
                 + (0.6 * (self.data.inst['wl_bins'] / 11e-6) ** (-0.4) * np.cos(lat)) ** 2)
            )

        lz_flux = lz_flux_sr * (np.pi * self.data.inst['hfov'] ** 2)

        lz_leak = (ap * self.data.inst['t_map']).sum(axis=(-2, -1)) / ap.sum() * lz_flux \
                  * self.data.inst['telescope_area']

        return lz_leak
