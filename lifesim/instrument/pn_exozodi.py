import numpy as np

from lifesim.dataio.bus import Module
from lifesim.modules import constants
from lifesim.modules.util import black_body


def get_exozodi_leakage(image_size: int,
                        l_sun: float,
                        distance_s: float,
                        mas_pix: np.ndarray,
                        rad_pix: np.ndarray,
                        z: float,
                        telescope_area: float,
                        radius_map: np.ndarray,
                        wl_bins: np.ndarray,
                        wl_bin_edges: np.ndarray,
                        t_map: np.ndarray,
                        hfov: np.ndarray):
    """
    Simulates the amount of photon noise originating from the exozodi of the observed system
    leaking into the LIFE array measurement

    Parameters
    ----------
    image_size : int
        Number of pixels on one axis of a square detector (dimensionless). I.e. for a 512x512
        detector this value is 512
    l_sun : float
        Luminosity of the observed star in [solar luminosities]
    distance_s : float
        Distance between the observed star and the LIFE array in [pc]
    mas_pix : np.ndarray
        Contains the size of each pixel projected to the sky in [milliarcseconds]
    z : float
        Zodi level in the observed system in [zodis]
    telescope_area : float
        Area of all array apertures combined in [m^2]
    radius_map : np.ndarray
        Contains the distance of a pixel from the center of the detector in [pix]
    wl_bins : np.ndarray
        Central values of the spectral bins in the wavelength regime in [m]
    wl_bin_edges : np.ndarray
        Edges of the spectral wavelength bins in [m]. For N bins, this array will contain N+1 edges
    t_map : np.ndarray
        Transmission map of the TM3 mode of the array created by the
        lifesim.TransmissionMap module

    Returns
    -------
    ez_leak
        Exozodi leakage in [s^-1] per wavelength bin
    """

    # calculate the parameters required by Kennedy2015
    alpha = 0.34
    r_in = 0.034422617777777775 * np.sqrt(l_sun)
    r_0 = np.sqrt(l_sun)
    distance_s_au = distance_s * 648000 / np.pi  # 360 * 3600 / (2 * pi)
    sigma_zero = 7.12e-8  # Sigma_{m,0} from Kennedy 2015 (doi:10.1088/0067-0049/216/2/23)

    # reshape the mas per pixel array for calculation (to (n, 1, 1))
    mas_pix = np.array([mas_pix])
    if mas_pix.shape[-1] > 1:
        mas_pix = np.reshape(mas_pix, (mas_pix.shape[-1], 1, 1))
    rad_pix = np.array([rad_pix])
    if rad_pix.shape[-1] > 1:
        rad_pix = np.reshape(rad_pix, (rad_pix.shape[-1], 1, 1))

    au_pix = mas_pix / 1e3 * distance_s
    area_m = au_pix ** 2
    r_au = radius_map * au_pix
    area_pix_at_aperture = telescope_area / image_size

    # identify all pixels where the radius is larges than the inner radius by kennedy2015
    r_cond = ((r_au >= r_in)
              & (r_au <= image_size / 2 * au_pix))

    # calculate the temperature at all pixel positions according to Kennedy2015 Eq. 2
    temp_map = np.where(r_cond,
                        278.3 * (l_sun ** 0.25) / np.sqrt(r_au), 0)

    # calculate the Sigma (Eq. 3) in Kennedy2015 and set everything inside the inner radius to 0
    sigma = np.where(r_cond,
                     # area_m *
                     sigma_zero * z *
                     (r_au / r_0) ** (-alpha), 0)

    # convert the wavelength bins to frequency bins and reshape like above (this approach is okay
    # because when calculating in densities it is admissible to change the units of the limits)
    freq_bins = constants.c / np.array([wl_bins])
    if freq_bins.shape[-1] > 1:
        freq_bins = np.reshape(freq_bins, (freq_bins.shape[-1], 1, 1))

    # same conversion for the bin edges
    freq_bin_edges = constants.c / np.array(wl_bin_edges)
    freq_widths = []
    for i in range(freq_bin_edges.shape[0]-1):
        freq_widths.append(freq_bin_edges[i]-freq_bin_edges[i+1])
    freq_widths = np.array(freq_widths)

    # TODO remove
    hfov = np.array(hfov)
    if hfov.shape[-1] > 1:
        hfov = np.reshape(hfov, (hfov.shape[-1], 1, 1))


    # get the black body radiation emitted by the interexoplanetary dust
    f_nu_disk = black_body(bins=freq_bins,
                           width=freq_widths,
                           temp=temp_map,
                           mode='frequency') \
                * sigma * np.pi * np.sin(hfov)**2 * telescope_area / 512**2
                # * sigma * np.pi * np.sin(rad_pix) ** 2 * area_pix_at_aperture


    ap = np.where(radius_map <= image_size/2, 1, 0)
    # add the transmission map
    ez_leak = (f_nu_disk * t_map * ap).sum(axis=(-2, -1))

    return ez_leak


class PhotonNoiseExozodi(Module):
    """
    'Plugin' module for calculating the photon noise created by the exozodiacal light of the
    observed system. The module has the function type 'photon_noise'

    Attributes
    ----------
    noise
        Exozodi leakage in [s^-1] per wavelength bin
    """
    def __init__(self,
                 name: str):
        """
        Parameters
        ----------
        name
            Name of the PhotonNoiseExozodi plugin
        """
        super().__init__(name=name)
        self.f_type = 'photon_noise'
        self.noise = None

    def run(self):
        """
        The run method executes the 'plugin' module
        """
        self.noise = get_exozodi_leakage(image_size=self.data['image_size'],
                                         l_sun=self.data['l_sun'],
                                         distance_s=self.data['distance_s'],
                                         mas_pix=self.data['mas_pix'],
                                         rad_pix=self.data['rad_pix'],
                                         z=self.data['z'],
                                         telescope_area=self.data['telescope_area'],
                                         radius_map=self.data['radius_map'],
                                         wl_bins=self.data['wl_bins'],
                                         wl_bin_edges=self.data['wl_bin_edges'],
                                         t_map=self.data['t_map'],
                                         hfov=self.data['hfov'])
