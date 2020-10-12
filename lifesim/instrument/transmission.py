import numpy as np

from lifesim.dataio.bus import Module


class TransmissionMap(Module):
    def __init__(self,
                 name: str):
        super().__init__(name=name)

        self.f_type = 'transmission'

        self.tm1, self.tm2, self.tm3, self.tm4, self.tm_chop = None, None, None, None, None

        # data needed ['wl', 'hfov_mas', 'image_size', 'bl', 'map_selection']

    def fast_transmission(self):
        wl = np.array([self.data['wl']])  # wavelength in m
        if wl.shape[-1] > 1:
            wl = np.reshape(wl, (wl.shape[-1], 1, 1))

        hfov_mas = np.array([self.data['hfov_mas']])  # wavelength in m
        if hfov_mas.shape[-1] > 1:
            hfov_mas = np.reshape(hfov_mas, (hfov_mas.shape[-1], 1, 1))

        hfov = hfov_mas / (1000 * 3600 * 180) * np.pi  # hfov in radians

        angle = np.linspace(-1, 1, self.data['image_size'])  # generare 1D array that spans field of view
        alpha = np.tile(angle, (self.data['image_size'], 1))  # angle matrix in x-direction ("alpha")
        beta = alpha.T  # angle matrix in y-direction ("beta")
        alpha = alpha * hfov
        beta = beta * hfov
        L = self.data['bl'] / 2  # smaller distance of telecscopes from center line

        if 'tm1' in self.data['map_selection']:
            self.tm1 = np.cos(2 * np.pi * L * alpha / wl) ** 2 * np.cos(
                12 * np.pi * L * beta / wl - np.pi / 4) ** 2  # transmission map of mode 1

        if 'tm2' in self.data['map_selection']:
            self.tm2 = np.cos(2 * np.pi * L * alpha / wl) ** 2 * np.cos(
                12 * np.pi * L * beta / wl + np.pi / 4) ** 2  # transmission map of mode 2

        if 'tm3' in self.data['map_selection']:
            self.tm3 = np.sin(2 * np.pi * L * alpha / wl) ** 2 * np.cos(
                12 * np.pi * L * beta / wl - np.pi / 4) ** 2  # transmission map of mode 3

        if 'tm4' in self.data['map_selection']:
            self.tm4 = np.sin(2 * np.pi * L * alpha / wl) ** 2 * np.cos(
                12 * np.pi * L * beta / wl + np.pi / 4) ** 2  # transmission map of mode 4

        if 'tm4' in self.data['map_selection']:
            if (self.tm3 is None) or (self.tm4 is None):
                raise ValueError('Third and fourth transmission maps are required to caluclate the'
                                 ' copped transmission map')

            # difference of transmission maps 3 and 4 = "chopped transmission"
            self.tm_chop = self.tm3 - self.tm4

    def run(self):
        self.fast_transmission()
