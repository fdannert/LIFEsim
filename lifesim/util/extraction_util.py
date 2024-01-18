import numpy as np

def convert_to_minus_pi_to_pi(angle):
    return (((angle + np.pi) % (2*np.pi)) - np.pi)


def convert_to_zero_to_2pi(angle):
    if (angle<0):
        angle += 2*np.pi
    return angle