import numpy as np

from lifesim.dataio.bus import Module


def fast_transmission(wl, hfov_mas, image_size, bl, map_selection):
    wl = np.array([wl])  # wavelength in m
    if wl.shape[-1] > 1:
        wl = np.reshape(wl, (wl.shape[-1], 1, 1))

    hfov_mas = np.array([hfov_mas])  # wavelength in m
    if hfov_mas.shape[-1] > 1:
        hfov_mas = np.reshape(hfov_mas, (hfov_mas.shape[-1], 1, 1))

    hfov = hfov_mas / (1000 * 3600 * 180) * np.pi  # hfov in radians

    angle = np.linspace(-1, 1, image_size)  # generare 1D array that spans field of view
    alpha = np.tile(angle, (image_size, 1))  # angle matrix in x-direction ("alpha")
    beta = alpha.T  # angle matrix in y-direction ("beta")
    alpha = alpha * hfov
    beta = beta * hfov
    L = bl / 2  # smaller distance of telecscopes from center line

    tm1, tm2, tm3, tm4, tm_chop = None, None, None, None, None

    if 'tm1' in map_selection:
        tm1 = np.cos(2 * np.pi * L * alpha / wl) ** 2 * np.cos(
            12 * np.pi * L * beta / wl - np.pi / 4) ** 2  # transmission map of mode 1

    if 'tm2' in map_selection:
        tm2 = np.cos(2 * np.pi * L * alpha / wl) ** 2 * np.cos(
            12 * np.pi * L * beta / wl + np.pi / 4) ** 2  # transmission map of mode 2

    if 'tm3' in map_selection:
        tm3 = np.sin(2 * np.pi * L * alpha / wl) ** 2 * np.cos(
            12 * np.pi * L * beta / wl - np.pi / 4) ** 2  # transmission map of mode 3

    if 'tm4' in map_selection:
        tm4 = np.sin(2 * np.pi * L * alpha / wl) ** 2 * np.cos(
            12 * np.pi * L * beta / wl + np.pi / 4) ** 2  # transmission map of mode 4

    if 'tm4' in map_selection:
        if (tm3 is None) or (tm4 is None):
            raise ValueError('Third and fourth transmission maps are required to caluclate the'
                             ' copped transmission map')

        # difference of transmission maps 3 and 4 = "chopped transmission"
        tm_chop = tm3 - tm4

    return tm1, tm2, tm3, tm4, tm_chop


class TransmissionMap(Module):
    def __init__(self,
                 name: str):
        super().__init__(name=name)

        self.f_type = 'transmission'

        self.tm1, self.tm2, self.tm3, self.tm4, self.tm_chop = None, None, None, None, None

        # data needed ['wl', 'hfov_mas', 'image_size', 'bl', 'map_selection']

    def run(self):
        self.tm1, self.tm2, self.tm3, self.tm4, self.tm_chop = \
            fast_transmission(wl=self.data['wl'],
                              hfov_mas=self.data['hfov_mas'],
                              image_size=self.data['image_size'],
                              bl=self.data['bl'],
                              map_selection=self.data['map_selection'])


