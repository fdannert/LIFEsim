import numpy as np
from typing import Union

from lifesim.util import constants


def planck_law(x: np.ndarray,
               temp: Union[float, np.ndarray],
               mode: str):
    """
    Calculates the photon flux emitted from a black body according to Planck's law in the
    wavelength or frequency regime

    Parameters
    ----------
    x : np.ndarray
        The frequency of wavelength at which the photon fluxes are calculated in [Hz] or [m]
    temp : Union[float, np.ndarray]
        The temperature of the black body
    mode : str
        If ``x`` is given in [Hz], set ``mode = 'frequency'. If ``x`` is given in [m], set
        ``mode = 'wavelength'

    Raises
    ------
    ValueError
        If the mode is not recognized

    Returns
    -------
    fgamma : np.ndarray
        The photon flux at the respective wavelengths or frequencies
    """

    # select the correct mode
    if mode == 'wavelength':

        # account for the temperature being zero at some pixels
        with np.errstate(divide='ignore'):

            # the Planck law divided by the photon energy to obtain the photon flux
            fgamma = 2 * constants.c / (x**4) / \
               (np.exp(constants.h * constants.c / x / constants.k / temp) - 1)
    elif mode == 'frequency':

        # account for the temperature being zero at some pixels
        with np.errstate(divide='ignore'):

            # the Planck law divided by the photon energy to obtain the photon flux
            fgamma = np.where(temp == 0,
                              0,
                              2 * x**2 / (constants.c**2) /
                              (np.exp(constants.h * x / constants.k / temp))-1.)
    else:
        raise ValueError('Mode not recognised')

    return fgamma


def black_body(mode: str,
               bins: np.ndarray,
               width: np.ndarray,
               temp: Union[float, np.ndarray],
               radius: float = None,
               distance: float = None):
    """
    Calculates the black body photon flux in wavelength or frequency as well as for planetary or
    stellar sources

    Parameters
    ----------
    mode : str
        Defines the mode of the ``black_body`` function.
            - ``mode = 'wavelength'`` : Clean photon flux black body spectrum over wavelength is
              returned. Parameters used are ``bins``, ``width`` and ``temp``
            - ``mode = 'frequency'`` : Clean photon flux black body spectrum over frequency is
              returned. Parameters used are ``bins``, ``width`` and ``temp``
            - ``mode = 'star'`` : Photon flux black body spectrum received from a star of specified
              radius from the specified distance. All parameters are used. In this mode, the
              parameter ``bins`` needs to be in wavelength
            - ``mode = 'planet'`` : Photon flux black body spectrum received from a planet of
              specified radius from the specified distance. All parameters are used. In this mode,
              the parameter ``bins`` needs to be in wavelength
    bins : np.ndarray
        The wavelength or frequency bins at which the black body is evaluated in [m] or [Hz]
        respectively
    width : np.ndarray
        The width of the wavelength or frequency bins to integrate over the black body spectrum in
        [m] or [Hz] respectively
    temp : Union[float, np.ndarray]
        The temperature of the black body
    radius : float
        The radius of the spherical black body object. For ``mode = 'star'`` in [sun_radii], for
        ``mode = 'planet'`` in [earth_radii]
    distance : float
        The distance between the instrument and the observed object in [pc]

    Raises
    ------
    ValueError
        If the mode is not recognized

    Returns
    -------
    fgamma : np.ndarray
        The photon flux at the respective wavelengths or frequencies
    """

    if mode == 'star':
        fgamma = planck_law(x=bins,
                            temp=temp,
                            mode='wavelength') * width \
                 * np.pi * ((radius * constants.radius_sun) / (distance * constants.m_per_pc)) ** 2
    elif mode == 'planet':
        fgamma = planck_law(x=bins,
                            temp=temp,
                            mode='wavelength') * width \
                 * np.pi * ((radius * constants.radius_earth) / (distance * constants.m_per_pc)) ** 2
    elif mode == 'wavelength':
        fgamma = planck_law(x=bins,
                            temp=temp,
                            mode='wavelength') * width
    elif mode == 'frequency':
        # TODO remove hardcoded np.newaxis solution. The redim is needed for the PhotonNoiseExozodi
        #   class
        fgamma = planck_law(x=bins,
                            temp=temp,
                            mode='frequency') * width[:, np.newaxis, np.newaxis]
    else:
        raise ValueError('Mode not recognised')

    return fgamma


def import_spectrum(pathtofile: str,
                    wl_bin_edges: np.ndarray,
                    radius_p: float,
                    distance_s: float,
                    radius_spec: float,
                    distance_spec: float,
                    clean: bool = False):
    spec = np.loadtxt(pathtofile).T
    spec[0] *= 1e-6  # per micron to per m
    # spec[1] /= 3600.  # hours to seconds
    spec[1] *= 1e6  # per micron to per m

    bins = np.digitize(spec[0], wl_bin_edges)
    bins_mean = [spec[1][bins == i].mean() for i in range(1, len(wl_bin_edges))]
    bins_mean = np.array(bins_mean)

    fgamma = bins_mean * ((radius_p / radius_spec) / (distance_s / distance_spec)) ** 2

    if clean:
        return spec
    else:
        return fgamma
