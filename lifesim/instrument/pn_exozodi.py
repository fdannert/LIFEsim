import numpy as np

from lifesim.dataio.bus import Module
from lifesim.modules import constants
from lifesim.modules.util import black_body


def get_exozodi_leakage(image_size: int,
                        l_sun: float,
                        distance_s: float,
                        mas_pix: float,
                        telescope_area: float,
                        radius_map: np.ndarray,
                        wl_bins: np.ndarray,
                        wl_bin_edges: np.ndarray):

    alpha = 0.34
    r_in = 0.034422617777777775 * np.sqrt(l_sun)
    r_0 = np.sqrt(l_sun)
    distance_s_au = distance_s * 648000 / np.pi  # 360 * 3600 / (2 * pi)
    sigma_zero = 7.12e-8  # Sigma_{m,0} from Kennedy 2015 (doi:10.1088/0067-0049/216/2/23)

    mas_pix = np.array([mas_pix])
    if mas_pix.shape[-1] > 1:
        mas_pix = np.reshape(mas_pix, (mas_pix.shape[-1], 1, 1))

    au_pix = mas_pix / 1e3 * distance_s
    area_m = au_pix ** 2
    r_au = radius_map * au_pix
    r_cond = ((r_au >= r_in)
              & (r_au <= image_size / 2 * au_pix))

    temp_map = np.where(r_cond,
                        278.3 * (l_sun ** 0.25) / np.sqrt(r_au), 0)

    sigma = np.where(r_cond,
                     area_m * sigma_zero *
                     (r_au / r_0) ** (-alpha), 0)

    freq_bins = constants.c / np.array([wl_bins])
    if freq_bins.shape[-1] > 1:
        freq_bins = np.reshape(freq_bins, (freq_bins.shape[-1], 1, 1))

    freq_bin_edges = constants.c / np.array(wl_bin_edges)
    freq_widths = []
    for i in range(freq_bin_edges.shape[0]-1):
        freq_widths.append(freq_bin_edges[i+1]-freq_bin_edges[i])
    freq_widths = np.array(freq_widths)

    f_nu_disk = black_body(bins=freq_bins,
                           width=freq_widths,
                           temp=temp_map,
                           mode='frequency') * area_m / (distance_s_au**2) * sigma * telescope_area

    return f_nu_disk


class PhotonNoiseExozodi(Module):
    def __init__(self,
                 name: str):
        super().__init__(name=name)
        self.f_type = 'photon_noise'
        self.noise = None

    def run(self):
        mask = self.data['c'].data.nstar == self.data['nstar']
        self.noise = get_exozodi_leakage(image_size=self.data['image_size'],
                                         l_sun=self.data['c'].data.l_sun[mask].to_numpy()[0],
                                         distance_s=self.data['c'].data.distance_s[mask].to_numpy()[0],
                                         mas_pix=self.data['mas_pix'],
                                         telescope_area=self.data['telescope_area'],
                                         radius_map=self.data['radius_map'],
                                         wl_bins=self.data['wl_bins'],
                                         wl_bin_edges=self.data['wl_bin_edges']) \
                     * self.data['c'].data['z'][mask].to_numpy()[0]
