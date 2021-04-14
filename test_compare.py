#%%
import lifesim as ls
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from numpy import cos, sin, pi

"""
Test and compare the result of timedependent simulation vs original time independent simulation 
"""
#%%
def setup_bus():
    bus = ls.Bus()
    bus.data.options.set_scenario('baseline')
    bus.data.import_catalog("C:/Users/Stoephu/Projects/SemesterProject2021/LIFEsim/output/suitable.hdf5")
    return bus

#%%
bus_o = setup_bus() #original sim
bus_t = setup_bus() #new sim

#%%
# ---------- MEASUREMENT SIMULATION ---------- copied from lifesim_demo.py by fdannert

# create modules and add to bus
inst = ls.Instrument(name='inst')
bus_t.add_module(inst)

transm = ls.TransmissionMap(name='transm')
bus_t.add_module(transm)

exo = ls.PhotonNoiseExozodi(name='exo')
bus_t.add_module(exo)
local = ls.PhotonNoiseLocalzodi(name='local')
bus_t.add_module(local)
star = ls.PhotonNoiseStar(name='star')
bus_t.add_module(star)

# connect all modules
bus_t.connect(('inst', 'transm'))
bus_t.connect(('inst', 'exo'))
bus_t.connect(('inst', 'local'))
bus_t.connect(('inst', 'star'))

bus_t.connect(('star', 'transm'))

#%%
inst.get_snr_test()

#%%
bus_o.data.catalog["snr_1h"]
#%%
bus_t.data.catalog["snr_1h"]
# %%
snrs = []
true_anomalies = np.arange(0, 2 * np.pi, np.pi / 360)
for anomaly in true_anomalies:
    bus_t.data.catalog["true_anomaly"] = anomaly
    inst.get_snr_t()
    snrs.append(bus_t.data.catalog[["true_anomaly", "snr_1h"]].copy())
snrs
# %%

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
