import numpy as np
from typing import Union

from lifesim.dataio.bus import Module


def fast_transmission(wl_bins: np.ndarray,
                      hfov: np.ndarray,
                      image_size: Union[int, type(None)],
                      bl: float,
                      map_selection: list,
                      ratio: float,
                      direct_mode: bool = False,
                      d_alpha: np.ndarray = None,
                      d_beta: np.ndarray = None):
    """
    Return the transmission map the LIFE array

    Parameters
    ----------
    image_size : int
        Number of pixels on one axis of a square detector (dimensionless). I.e. for a 512x512
        detector this value is 512
    bl : float
        Length of the shorter, nulling baseline in [m]
    wl_bins : np.ndarray
        Central values of the spectral bins in the wavelength regime in [m]
    hfov : np.ndarray
        Contains the half field of view of the observatory in [rad] for each of the spectral bins
    map_selection : list
        A list of strings specifying the array transmission mode from which the transmission map is
        generated. The possible options are 'tm1', 'tm2', 'tm3', 'tm4'
    ratio : float
        Ratio between the nulling and the imaging baseline. E.g. if the imaging baseline is twice
        as long as the nulling baseline, the ratio will be 2
    direct_mode: bool
        If direct mode is set to True, the variables wl_bins, d_alpha and d_beta are injected into
        the transmission calculation directly. The output is not necessarily in the the form of a
        map. The inputs hfov and image size will be ignored in this mode
    d_alpha: np.ndarray
        The x-positions of the points to be evaluated measured from the central viewing axis in
        [rad]
    d_beta: np.ndarray
        The y-positions of the points to be evaluated measured from the central viewing axis in
        [rad]

    Returns
    -------
    tm1
        Transmission map from first mode
    tm2
        Transmission map from second mode
    tm3
        Transmission map from third mode
    tm4
        Transmission map from fourth mode
    tm_chop
        The chopped transmission map calculated by subtracting tm4 from tm3
    """

    if direct_mode:
        alpha = d_alpha
        beta = d_beta
    else:
        # reshape the wl_bins and hfov arrays for calculation (to (n, 1, 1))
        wl_bins = np.array([wl_bins])  # wavelength in m
        if wl_bins.shape[-1] > 1:
            wl_bins = np.reshape(wl_bins, (wl_bins.shape[-1], 1, 1))

        hfov = np.array([hfov])  # wavelength in m
        if hfov.shape[-1] > 1:
            hfov = np.reshape(hfov, (hfov.shape[-1], 1, 1))

        # generare 1D array that spans field of view
        angle = np.linspace(-1, 1, image_size)

        # angle matrix in x-direction ("alpha")
        alpha = np.tile(angle, (image_size, 1))

        # angle matrix in y-direction ("beta")
        beta = alpha.T

        # convert angle matrices to fov units
        alpha = alpha * hfov
        beta = beta * hfov

    # smaller distance of apertures from center line
    L = bl / 2

    tm1, tm2, tm3, tm4, tm_chop = None, None, None, None, None

    if 'tm1' in map_selection:
        tm1 = np.cos(2 * np.pi * L * alpha / wl_bins) ** 2 * np.cos(
            2 * ratio * np.pi * L * beta / wl_bins - np.pi / 4) ** 2  # transmission map of mode 1

    if 'tm2' in map_selection:
        tm2 = np.cos(2 * np.pi * L * alpha / wl_bins) ** 2 * np.cos(
            2 * ratio * np.pi * L * beta / wl_bins + np.pi / 4) ** 2  # transmission map of mode 2

    if 'tm3' in map_selection:
        tm3 = np.sin(2 * np.pi * L * alpha / wl_bins) ** 2 * np.cos(
            2 * ratio * np.pi * L * beta / wl_bins - np.pi / 4) ** 2  # transmission map of mode 3

    if 'tm4' in map_selection:
        tm4 = np.sin(2 * np.pi * L * alpha / wl_bins) ** 2 * np.cos(
            2 * ratio * np.pi * L * beta / wl_bins + np.pi / 4) ** 2  # transmission map of mode 4

    # difference of transmission maps 3 and 4 = "chopped transmission"
    if 'tm_chop' in map_selection:

        # if tm3 and tm4 exist, calculate the chopped transmission directly from the difference
        if (tm3 is not None) and (tm4 is not None):
            tm_chop = tm3 - tm4

        # if they don't exist, calculate the chopped transmission from formula
        else:
            tm_chop = np.sin(2 * np.pi * L * alpha / wl_bins) ** 2 * np.sin( # chopped transm. map
            4 * ratio * np.pi * L * beta / wl_bins)

    return tm1, tm2, tm3, tm4, tm_chop


def transm_curve(bl: float,
                 wl_bins: np.ndarray,
                 angsep: float,
                 ratio: float,
                 phi_n: int = 360):
    """
    Calculates the radial transmission curve of the LIFE array

    Parameters
    ----------
    bl : float
        Length of the shorter, nulling baseline in [m]
    wl_bins : np.ndarray
        Central values of the spectral bins in the wavelength regime in [m]
    angsep : float
        Angular separation between the observed star and the observed exoplanet in [arcsec]
    phi_n : int
        Number of rotation steps used in integration
    ratio : float
        Ratio between the nulling and the imaging baseline. E.g. if the imaging baseline is twice
        as long as the nulling baseline, the ratio will be 2

    Returns
    -------
    transm_curve_chop
        Radial transmission curve corresponding to the chopped transmission map
    transm_curve_tm4
        Radial transmission curve corresponding to the transmission map of the 4th mode 'tm4'
    """

    # reshape the wl_bins arrays for calculation (to (n, 1))
    # TODO make this a method and apply it where needed
    wl_bins = np.array([wl_bins])  # wavelength in m
    if wl_bins.shape[-1] > 1:
        wl_bins = np.reshape(wl_bins, (wl_bins.shape[-1], 1))

    # convert angular separation to radians
    angsep_rad = angsep / (3600 * 180) * np.pi

    # create 1D array with azimuthal coordinates
    phi_lin = np.linspace(0, 2 * np.pi, phi_n, endpoint=False)

    # retrieve the transmission curves
    (_, _, _,
     transm_curve_tm4,
     transm_curve_chop) = fast_transmission(wl_bins=wl_bins,
                                            hfov=None,
                                            image_size=None,
                                            bl=bl,
                                            map_selection=['tm4', 'tm_chop'],
                                            ratio=ratio,
                                            direct_mode=True,
                                            d_alpha=angsep_rad * np.cos(phi_lin),
                                            d_beta=angsep_rad * np.sin(phi_lin))

    return transm_curve_chop, transm_curve_tm4


def transm_eff(bl: float,
               wl_bins: np.ndarray,
               angsep: np.ndarray,
               ratio: float):
    """
    Integrates over transmission curves to get the transmission efficiency for signal and noise

    Parameters
    ----------
    bl : float
        Length of the shorter, nulling baseline in [m]
    wl_bins : np.ndarray
        Central values of the spectral bins in the wavelength regime in [m]
    angsep : float
        Angular separation between the observed star and the observed exoplanet in [arcsec]
    ratio : float
        Ratio between the nulling and the imaging baseline. E.g. if the imaging baseline is twice
        as long as the nulling baseline, the ratio will be 2

    Returns
    -------
    transm_eff
        Transmission efficiency per spectral bin for the exoplanet signal
    transm_noise
        Transmission efficiency per spectral bin for the photon noise received from the exoplanet
        signal
    """
    tc_chop, tc_tm4 = transm_curve(bl=bl,
                                   wl_bins=wl_bins,
                                   angsep=angsep,
                                   ratio=ratio)

    # integrate over angles to get transmission efficiency
    transm_eff = np.sqrt((tc_chop ** 2).mean(axis=-1))
    transm_noise = np.sqrt((tc_tm4 ** 2).mean(axis=-1))
    return transm_eff, transm_noise


class TransmissionMap(Module):
    """
    'Plugin' module for calculating transmission trough the LIFE array. The module can be
    operated in 'map' mode, where it returns the transmission maps, or in 'efficiency' mode,
    where it returns the transmission efficiency. The module has the function type
    'transmission'

    Attributes
    ----------
    tm1
        Transmission map from first mode
    tm2
        Transmission map from second mode
    tm3
        Transmission map from third mode
    tm4
        Transmission map from fourth mode
    tm_chop
        The chopped transmission map calculated by subtracting tm4 from tm3
    transm_eff
        Transmission efficiency per spectral bin for the exoplanet signal
    transm_noise
        Transmission efficiency per spectral bin for the photon noise received from the
        exoplanet signal
    """

    def __init__(self,
                 name: str):
        """
        Parameters
        ----------
        name
            Name of the TransmissionMap plugin
        """

        super().__init__(name=name)

        self.f_type = 'transmission'

        self.tm1, self.tm2, self.tm3, self.tm4, self.tm_chop = None, None, None, None, None
        self.transm_eff, self.transm_noise = None, None

    def run(self,
            args: dict):
        """
        The run method executes the 'plugin' module

        Parameters
        ----------
        args : dict
            A dictionary of arguments passed to the run method. Here, the only accepted argument is
            of key 'mode'. If it is set to 'map' (``args['mode'] == 'map'``), the TransmissionMap
            'plugin' module will generate transmission maps. If it is set to 'efficiency'
            (``args['mode'] == 'efficiency'``), the TransmissionMap 'plugin' module will generate
            transmission efficiencies

        Raises
        ------
        ValueError
            If the mode is not recognized, i.e. is not ``'map'`` or ``'efficiency'``
        """

        # select between the modes
        if args['mode'] == 'map':
            self.tm1, self.tm2, self.tm3, self.tm4, self.tm_chop = \
                fast_transmission(wl_bins=self.data['wl_bins'],
                                  hfov=self.data['hfov'],
                                  image_size=self.data['image_size'],
                                  bl=self.data['bl'],
                                  map_selection=self.data['map_selection'],
                                  ratio=self.data['ratio'])
        elif args['mode'] == 'efficiency':
            self.transm_eff, self.transm_noise = transm_eff(bl=self.data['bl'],
                                                            wl_bins=self.data['wl_bins'],
                                                            angsep=self.data['ang_sep_as'],
                                                            ratio=self.data['ratio'])
        else:
            raise ValueError('Mode not recognized')


