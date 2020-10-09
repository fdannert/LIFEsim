import warnings

import numpy as np

from lifesim.archive.catalog import Catalog

def compute_habitable_zone(catalog: Catalog,
                           model: str = 'MS'):

    if 'stars' not in catalog.masks.keys():
        raise KeyError('Stars mask pre-calculation missing')

    if not ((model == 'MS') or (model == 'POST-MS')):
        warnings.warn(str(model) + ' is an unknown model. Using model MS')
        model = 'MS'

    s_in = np.zeros_like(catalog.data.nstar, dtype=float)
    s_out = np.zeros_like(catalog.data.nstar, dtype=float)
    l_sun = np.zeros_like(catalog.data.nstar, dtype=float)

    hz_in = np.zeros_like(catalog.data.nstar, dtype=float)
    hz_out = np.zeros_like(catalog.data.nstar, dtype=float)
    hz_center = np.zeros_like(catalog.data.nstar, dtype=float)

    for _, n in enumerate(np.where(catalog.masks['stars'])[0]):
        if (model == 'MS'):
            s0_in, s0_out = 1.7665, 0.3240
            a_in, a_out = 1.3351E-4, 5.3221E-5
            b_in, b_out = 3.1515E-9, 1.4288E-9
            c_in, c_out = -3.3488E-12, -1.1049E-12
        elif (model == 'POST-MS'):
            s0_in, s0_out = 1.1066, 0.3240
            a_in, a_out = 1.2181E-4, 5.3221E-5
            b_in, b_out = 1.5340E-8, 1.4288E-9
            c_in, c_out = -1.5018E-12, -1.1049E-12

        t_star = catalog.data.temp_s[n] - 5780

        # in units of incoming flux, normalized to stellar flux on earth_0
        s_in[catalog.data.nstar == catalog.data.nstar[n]] = \
            s0_in + a_in * t_star + b_in * t_star ** 2 + c_in * t_star ** 3
        s_out[catalog.data.nstar == catalog.data.nstar[n]] = \
            s0_out + a_out * t_star + b_out * t_star ** 2 + c_out * t_star ** 3
        l_sun[catalog.data.nstar == catalog.data.nstar[n]] = \
            catalog.data.radius_s[n] ** 2 \
            * (catalog.data.temp_s[n] / 5780) ** 4  # luminosity in L_sun

        hz_in[catalog.data.nstar == catalog.data.nstar[n]] = \
            np.sqrt(l_sun[n] / s_in[n])  # HZ inner boundery in AU
        hz_out[catalog.data.nstar == catalog.data.nstar[n]] = \
            np.sqrt(l_sun[n] / s_out[n])  # HZ outer boundery in AU
        hz_center[catalog.data.nstar == catalog.data.nstar[n]] = \
            (hz_in[n] + hz_out[n]) / 2  # HZ center  in AU

    catalog.data['s_in'] = s_in
    catalog.data['s_out'] = s_out
    catalog.data['l_sun'] = l_sun
    catalog.data['hz_in'] = hz_in
    catalog.data['hz_out'] = hz_out
    catalog.data['hz_center'] = hz_center

