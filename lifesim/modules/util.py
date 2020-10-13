import numpy as np
import scipy.integrate as integrate

from lifesim.modules import constants


def planck_law(wl, temp):
    return 2 * constants.c / (wl**4) / \
           (np.exp(constants.h * constants.c / wl / constants.k / temp) - 1)


# TODO The black body number flux function was changed to an integral over the planck function.
#   Compare with non-integral function (commented out) to see performance difference
def black_body(mode: str,
               wl_edges: np.ndarray,
               temp: float,
               radius: float = None,
               distance: float = None):
    fgamma = []
    for i in range(wl_edges.shape[0] - 1):
        f = integrate.quad(lambda l: planck_law(l, temp), wl_edges[i], wl_edges[i + 1])[0]
        fgamma.append(f)
    fgamma = np.array(fgamma)

    # k1 = 2 * constants.c
    # k2 = (constants.h * constants.c) / constants.k
    # fact1 = k1 / (wl ** 4)
    # fact2 = k2 / (temp * wl)
    # fgamma = np.array(fact1 / (np.exp(fact2) - 1.0)) * 1e-6

    if mode == 'star':
        fgamma = fgamma \
                 * np.pi * ((radius * constants.R_sun) / (distance * constants.m_per_pc)) ** 2
    elif mode == 'clean':
        pass
    else:
        raise ValueError('Mode not recognised')

    return fgamma
