import numpy as np

from lifesim.dataio.bus import Module
from lifesim.dataio.catalog import Catalog
from lifesim.modules.util import black_body


def get_localzodi_leakage(lz_model: str,
                          lat_s: float,
                          telescope_area: float,
                          image_size: int,
                          t_map: np.ndarray,
                          radius_map: np.ndarray,
                          wl_bins: np.ndarray,
                          wl_width: np.ndarray,
                          hfov: np.ndarray):
    # TODO Implement longitude dependence of localzodi

    if lz_model == 'glasse':
        temp = 270
        epsilon = 4.30e-8

    elif lz_model == 'darwinsim':
        radius_sun_au = 0.00465047  # in AU
        tau = 4e-8
        temp_eff = 265
        temp_sun = 5777
        a = 0.22

    else:
        raise ValueError('Specified model does not exist')
    
    long = 4 / 4 * np.pi
    lat = lat_s

    ap = np.where(radius_map <= image_size / 2, 1, 0)

    if lz_model == "glasse":
        lz_flux_sr = epsilon * black_body(mode='wavelength',
                                          bins=wl_bins,
                                          width=wl_width,
                                          temp=temp)

    elif lz_model == "darwinsim":
        b_tot = black_body(mode='wavelength',
                           bins=wl_bins,
                           width=wl_width,
                           temp=temp_eff) + a \
               * black_body(mode='wavelength',
                            bins=wl_bins,
                            width=wl_width,
                            temp=temp_sun) \
               * (radius_sun_au / 1.5) ** 2
        lz_flux_sr = tau * b_tot * np.sqrt(
            np.pi / np.arccos(np.cos(long) * np.cos(lat)) /
            (np.sin(lat) ** 2 + (0.6 * (wl_bins / 11e-6) ** (-0.4) * np.cos(lat)) ** 2)
        )

    lz_flux = lz_flux_sr * (np.pi * hfov ** 2)

    lz_leak = (ap * t_map).sum(axis=(-2, -1)) / ap.sum(
    ) * lz_flux * telescope_area
    return lz_leak


class PhotonNoiseLocalzodi(Module):
    def __init__(self,
                 name: str):
        super().__init__(name=name)
        self.f_type = 'photon_noise'
        self.noise = None

    def run(self):
        self.noise = get_localzodi_leakage(lz_model=self.data['lz_model'],
                                           lat_s=self.data['lat_s'],
                                           telescope_area=self.data['telescope_area'],
                                           image_size=self.data['image_size'],
                                           t_map=self.data['t_map'],
                                           radius_map=self.data['radius_map'],
                                           wl_bins=self.data['wl_bins'],
                                           wl_width=self.data['wl_width'],
                                           hfov=self.data['hfov'])
