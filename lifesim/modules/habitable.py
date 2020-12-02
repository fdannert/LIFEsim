import warnings

import numpy as np

# from lifesim.dataio.catalog import Catalog


def single_habitable_zone(model: str,
                          temp_s: float,
                          radius_s: float):

    # TODO which paper is this model based on?
    if model == 'MS':
        s0_in, s0_out = 1.7665, 0.3240
        a_in, a_out = 1.3351E-4, 5.3221E-5
        b_in, b_out = 3.1515E-9, 1.4288E-9
        c_in, c_out = -3.3488E-12, -1.1049E-12
    elif model == 'POST-MS':
        s0_in, s0_out = 1.1066, 0.3240
        a_in, a_out = 1.2181E-4, 5.3221E-5
        b_in, b_out = 1.5340E-8, 1.4288E-9
        c_in, c_out = -1.5018E-12, -1.1049E-12
    else:
        raise ValueError('Unknown model')

    t_star = temp_s - 5780

    # in units of incoming flux, normalized to stellar flux on earth_0
    s_in = s0_in + a_in * t_star + b_in * t_star ** 2 + c_in * t_star ** 3
    s_out = s0_out + a_out * t_star + b_out * t_star ** 2 + c_out * t_star ** 3
    l_sun = radius_s ** 2 * (temp_s / 5780) ** 4  # luminosity in L_sun

    hz_in = np.sqrt(l_sun / s_in)  # HZ inner boundery in AU
    hz_out = np.sqrt(l_sun / s_out)  # HZ outer boundery in AU
    hz_center = (hz_in + hz_out) / 2  # HZ center  in AU

    return s_in, s_out, l_sun, hz_in, hz_out, hz_center

# TODO remove frankenstein fix of catalog object
def compute_habitable_zone(catalog: object,
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

        s_in[catalog.data.nstar == catalog.data.nstar[n]], \
            s_out[catalog.data.nstar == catalog.data.nstar[n]],\
            l_sun[catalog.data.nstar == catalog.data.nstar[n]],\
            hz_in[catalog.data.nstar == catalog.data.nstar[n]],\
            hz_out[catalog.data.nstar == catalog.data.nstar[n]],\
            hz_center[catalog.data.nstar == catalog.data.nstar[n]] \
            = single_habitable_zone(model=model,
                                    temp_s=catalog.data.temp_s[n],
                                    radius_s=catalog.data.radius_s[n])

    catalog.data['s_in'] = s_in
    catalog.data['s_out'] = s_out
    catalog.data['l_sun'] = l_sun
    catalog.data['hz_in'] = hz_in
    catalog.data['hz_out'] = hz_out
    catalog.data['hz_center'] = hz_center
    catalog.data['habitable'] = np.logical_and.reduce((
        (catalog.data['semimajor_p'] > catalog.data['hz_in']).to_numpy(),
        (catalog.data['semimajor_p'] < catalog.data['hz_out']).to_numpy(),
        (catalog.data['radius_p'].ge(0.5)).to_numpy(),
        (catalog.data['radius_p'].le(1.5)).to_numpy()))

