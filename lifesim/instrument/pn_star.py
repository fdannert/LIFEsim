import numpy as np

from lifesim.dataio.bus import Module
from lifesim.instrument.transmission import fast_transmission
from lifesim.modules.util import black_body


def get_stellar_leakage(radius_s: float,
                        distance_s: float,
                        temp_s: float,
                        bl: float,
                        telescope_area: float,
                        wl_bins: np.ndarray,
                        wl_bin_widths: np.ndarray,
                        ratio: float,
                        image_size: int = 50,
                        map_selection: str = 'tm3'):
    """
    Simulates the amount of photon noise originating from the star of the observed system leaking
    into the LIFE array measurement

    Parameters
    ----------
    image_size : int
        Number of pixels on one axis of a square detector (dimensionless). I.e. for a 512x512
        detector this value is 512
    telescope_area : float
        Area of all array apertures combined in [m^2]
    wl_bins : np.ndarray
        Central values of the spectral bins in the wavelength regime in [m]
    wl_bin_widths : np.ndarray
        Widths of the spectral wavelength bins in [m]
    distance_s : float
        Distance between the observed star and the LIFE array in [pc]
    temp_s : float
        Temperature of the observed star in [K]
    radius_s : float
        Radius of the observed star in [sun radii]
    bl : float
        Length of the shorter, nulling baseline in [m]
    map_selection : str
        Select from which arm of the array the transmission map for the calculation of the leakage
        is taken
    ratio : float
        Ratio between the nulling and the imaging baseline. E.g. if the imaging baseline is twice
        as long as the nulling baseline, the ratio will be 2

    Returns
    -------
    sl_leak
        Stellar light leakage in [s^-1] per wavelength bin

    Raises
    ______

    ValueError
        If the specified transmission map does not exits
    """

    # check if the specified map exists
    if map_selection not in ['tm1', 'tm2', 'tm3', 'tm4']:
        raise ValueError('Nonexistent transmission map')

    # convert units
    Rs_au = 0.00465047 * radius_s
    Rs_as = Rs_au / distance_s
    Rs_mas = float(Rs_as)
    Rs_rad = Rs_mas / (3600. * 180.) * np.pi

    # TODO Instead of recalculating the transmission map for the stellar radius here, one could try
    #   to reuse the inner part of the transmission map already calculated in the get_snr function
    #   of the instrument class
    # TODO: why are we not reusing the maps calculated in the instrument class
    # get a transmission map at the stellar radius
    tm_star = fast_transmission(wl_bins=wl_bins,
                                hfov=Rs_rad,
                                image_size=image_size,
                                bl=bl,
                                map_selection=[map_selection],
                                ratio=ratio)[int(map_selection[-1]) - 1]

    x_map = np.tile(np.array(range(0, image_size)), (image_size, 1))
    y_map = x_map.T
    r_square_map = (x_map - (image_size - 1) / 2) ** 2 + (y_map - (image_size - 1) / 2) ** 2
    star_px = np.where(r_square_map < (image_size / 2) ** 2, 1, 0)

    # get the stellar leakage
    sl_leak = (star_px * tm_star).sum(axis=(-2, -1)) / star_px.sum(
    ) * black_body(bins=wl_bins,
                   width=wl_bin_widths,
                   temp=temp_s,
                   radius=radius_s,
                   distance=distance_s,
                   mode='star') * telescope_area
    return sl_leak


class PhotonNoiseStar(Module):
    """
    'Plugin' module for calculating the photon noise created by the stellar light in the observed
    system. The module has the function type 'photon_noise'

    Attributes
    ----------
    noise
        Stellar light leakage in [s^-1] per wavelength bin
    """
    def __init__(self,
                 name: str):
        """
        Parameters
        ----------
        name
            Name of the PhotonNoiseStar plugin
        """
        super().__init__(name=name)
        self.f_type = 'photon_noise'
        self.noise = None

    def run(self):
        """
        The run method executes the 'plugin' module
        """
        self.noise = get_stellar_leakage(radius_s=self.data['radius_s'],
                                         distance_s=self.data['distance_s'],
                                         temp_s=self.data['temp_s'],
                                         bl=self.data['bl'],
                                         telescope_area=self.data['telescope_area'],
                                         wl_bins=self.data['wl_bins'],
                                         wl_bin_widths=self.data['wl_width'],
                                         ratio=self.data['ratio'])
