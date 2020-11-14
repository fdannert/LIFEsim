import numpy as np
import scipy.integrate as integrate

from lifesim.modules import constants


def planck_law(x, temp, mode):
    if mode == 'wavelength':
        fgamma = 2 * constants.c / (x**4) / \
           (np.exp(constants.h * constants.c / x / constants.k / temp) - 1)
    elif mode == 'frequency':
        # TODO temp == 0 is unused but still calculated. It thus still throws a runtime warning
        fgamma = np.where(temp == 0,
                          0,
                          2 * x**2 / (constants.c**2) /
                          (np.exp(constants.h * x / constants.k / temp)))
    else:
        raise ValueError('Mode not recognised')
    return fgamma


# TODO The black body number flux function was changed to an integral over the planck function.
#   Compare with non-integral function (commented out) to see performance difference
def black_body(mode: str,
               bins: np.ndarray,
               width: np.ndarray,
               temp: float,
               radius: float = None,
               distance: float = None):
    # fgamma = []
    # for i in range(wl_edges.shape[0] - 1):
    #     f = integrate.quad(lambda l: planck_law(l, temp), wl_edges[i], wl_edges[i + 1])[0]
    #     fgamma.append(f)
    # fgamma = np.array(fgamma)
    # fgamma = planck_law(wl=wl_bins,
    #                     temp=temp) * wl_width
    # k1 = 2 * constants.c
    # k2 = (constants.h * constants.c) / constants.k
    # fact1 = k1 / (wl ** 4)
    # fact2 = k2 / (temp * wl)
    # fgamma = np.array(fact1 / (np.exp(fact2) - 1.0)) * 1e-6

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
                    clean: bool = False):
    spec = np.loadtxt(pathtofile).T
    spec[0] *= 1e-6  # per micron to per m
    spec[1:] /= 3600.  # hours to seconds

    bins = np.digitize(spec[0], wl_bin_edges)
    bins_mean = [spec[1][bins == i].mean() for i in range(1, len(wl_bin_edges))]
    bins_mean = np.array(bins_mean)

    fgamma = bins_mean * np.pi * (
                (radius_p * constants.radius_earth) / (distance_s * constants.m_per_pc)) ** 2

    if clean:
        return spec
    else:
        return fgamma
