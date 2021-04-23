#%%
import lifesim as ls
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from numpy import cos, sin, pi

def true_angular_separation(theta, large_omega, inclination, L, a):
    """
    this function calculates the true angular separation between the exoplanet and the hoststar as seen from earth
    using orbital parameters which can be found in the universe catalog.
    This assumes circular orbits
    Parameters
    ----------
    theta : float
        The true anomaly of the exoplanet in radians
    large_omega : float
        The angle between the ascending node of the orbit and the line of sight vector of the observer centered on the host star in radians
    inclination : float
        The angle between the plane of reference and the orbital plane in radians
    L : float
        The distance of the star system to ours in pc
    a : float
        The length of the semi-major axis of the orbit in AU
    
    Returns
    --------
    angular separation : float
        Returns the observed angular separation in arc seconds
    """

    return (a / L) * np.sqrt(np.cos(large_omega + theta)**2 + np.cos(inclination)**2 * np.sin(large_omega + theta) ** 2)


# %%
"""
def v(theta, o, i):
    return np.array([np.cos(o)*np.cos(i)*np.sin(theta)+np.sin(o)*np.cos(theta),
            np.sin(o)*np.cos(i)*np.sin(theta) - np.cos(o)*np.cos(theta),
            -np.sin(i)*np.sin(theta)])
"""
def v(theta, o, i):
    start = np.array([cos(theta), 0, sin(theta)])
    #Rotation around x axis
    Ri = np.array([[1, 0, 0],
                    [0, cos(i), -sin(i)],
                    [0, sin(i), cos(i)]])
    """
    #Rotation around y axis
    Ro = np.array([[cos(o), 0, -sin(o)],
                    [0, 1, 0],
                    [sin(o), 0, cos(o)]])
    """
    #Rotation around z axis
    Ro = np.array([[cos(o), -sin(o), 0],
                    [sin(o), cos(o), 0],
                    [0, 0, 1]])
    return Ro.dot(Ri.dot(start))
# %%
def showOrbit(o, i):
    thetas = np.linspace(0, np.pi * 2, 360)
    vs = v(thetas, o, i)
    coords = np.array([vs[0],vs[2]])
    plt.plot(coords[0], coords[1])
    plt.ylim([-1.1,1.1])
    plt.xlim([-1.1,1.1])
    plt.show()
    plt.figure()
    distances = np.sqrt(coords[0]**2 + coords[1]**2)
    distances_smart = true_angular_separation(thetas, o, i, 1, 1)
    plt.plot(thetas,distances, label = "me", marker= "x")
    plt.plot(thetas,distances_smart,label="smart")
    plt.legend()
    plt.show()

# %%
