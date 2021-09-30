import numpy as np


def single_habitable_zone(model: str,
                          temp_s: float,
                          radius_s: float):
    """
    Calculates the location of the habitable zone according to Kaltenegger+2017.

    Parameters
    ----------
    model : str
        Specifies the model for the habitable zone. The options are 'MS' and 'PRE-MS' according to
        Kaltenegger+2017 Table 1.
    temp_s : float
        Temperature of the star in [K].
    radius_s : float
        Radius of the star in [solar radii].

    Returns
    -------
    s_in
        Inner edge of the habitable zone in [earth insolation].
    s_out
        Outer edge of the habitable zone in [earth insolation].
    l_sun
        Luminosity of the star in [stellar luminosities].
    hz_in
        Inner edge of the habitable zone in [AU].
    hz_out
        Outer edge of the habitable zone in [AU].
    hz_center
        Center of the habitable zone in [AU].

    Raises
    ------
    ValueError
        If the specified model does not exits.
    """
    # model based on Kaltenegger+2017
    if model == 'MS':
        s0_in, s0_out = 1.7665, 0.3240
        a_in, a_out = 1.3351E-4, 5.3221E-5
        b_in, b_out = 3.1515E-9, 1.4288E-9
        c_in, c_out = -3.3488E-12, -1.1049E-12
        d_in, d_out = 0, 0
    elif model == 'POST-MS':
        s0_in, s0_out = 1.1066, 0.3240
        a_in, a_out = 1.2181E-4, 5.3221E-5
        b_in, b_out = 1.5340E-8, 1.4288E-9
        c_in, c_out = -1.5018E-12, -1.1049E-12
        d_in, d_out = 0, 0
    elif model == 'Kopparapu-Optimistic':
        s0_in, s0_out = 1.776, 0.32
        a_in, a_out = 2.136e-4, 5.547e-5
        b_in, b_out = 2.533e-8, 1.526e-9
        c_in, c_out = -1.332e-11, -2.874e-12
        d_in, d_out = -3.097e-15, -5.011e-16
    elif model == 'Kopparapu-Conservative':
        s0_in, s0_out = 1.107, 0.356
        a_in, a_out = 1.332e-4, 6.171e-5
        b_in, b_out = 1.58e-8, 1.698e-9
        c_in, c_out = -8.308e-12, -3.198e-12
        d_in, d_out = -1.931e-15, -5.575e-16
    else:
        raise ValueError('Unknown model')

    t_star = temp_s - 5780

    # in units of incoming flux, normalized to stellar flux on earth_0
    s_in = s0_in + a_in * t_star + b_in * t_star ** 2 + c_in * t_star ** 3 + d_in * t_star ** 4
    s_out = s0_out + a_out * t_star + b_out * t_star ** 2 + c_out * t_star ** 3 + d_out * t_star ** 4
    l_sun = radius_s ** 2 * (temp_s / 5780) ** 4  # luminosity in L_sun

    hz_in = np.sqrt(l_sun / s_in)  # HZ inner boundery in AU
    hz_out = np.sqrt(l_sun / s_out)  # HZ outer boundery in AU
    hz_center = (hz_in + hz_out) / 2  # HZ center  in AU

    return s_in, s_out, l_sun, hz_in, hz_out, hz_center
