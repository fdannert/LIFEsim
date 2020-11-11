import numpy as np

from lifesim.dataio.bus import Module
from lifesim.dataio.catalog import Catalog
from lifesim.instrument.transmission import fast_transmission
from lifesim.modules.util import black_body


def get_stellar_leakage(radius_s: float,
                        distance_s: float,
                        temp_s: float,
                        bl: float,
                        telescope_area: float,
                        wl_bins: np.ndarray,
                        wl_width: np.ndarray,
                        image_size: int = 50,
                        map_selection: str = 'tm3'):
    if map_selection not in ['tm1', 'tm2', 'tm3', 'tm4']:
        raise ValueError('Nonexistent transmission map')

    Rs_au = 0.00465047 * radius_s
    Rs_mas = Rs_au / distance_s * 1000
    Rs_mas = float(Rs_mas)

    # TODO Instead of recalculating the transmission map for the stelar radius here, one could try
    #   to reuse the inner part of the transmission map already calculated in the get_snr function
    #   of the instrument class
    tm_star = fast_transmission(wl=wl_bins,
                                hfov_mas=Rs_mas,
                                image_size=image_size,
                                bl=bl,
                                map_selection=[map_selection])[int(map_selection[-1]) - 1]

    x_map = np.tile(np.array(range(0, image_size)), (image_size, 1))
    y_map = x_map.T
    r_square_map = (x_map - (image_size - 1) / 2) ** 2 + (y_map - (image_size - 1) / 2) ** 2
    star_px = np.where(r_square_map < (image_size / 2) ** 2, 1, 0)

    sl = (star_px * tm_star).sum(axis=(-2, -1)) / star_px.sum(
    ) * black_body(bins=wl_bins,
                   width=wl_width,
                   temp=temp_s,
                   radius=radius_s,
                   distance=distance_s,
                   mode='star') * telescope_area
    return sl


class PhotonNoiseStar(Module):
    def __init__(self,
                 name: str):
        super().__init__(name=name)
        self.f_type = 'photon_noise'
        self.noise = None

    def run(self):
        self.noise = get_stellar_leakage(radius_s=self.data['radius_s'],
                                         distance_s=self.data['distance_s'],
                                         temp_s=self.data['temp_s'],
                                         bl=self.data['bl'],
                                         telescope_area=self.data['telescope_area'],
                                         wl_bins=self.data['wl_bins'],
                                         wl_width=self.data['wl_width'])
