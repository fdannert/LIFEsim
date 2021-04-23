#%%
import importlib
import lifesim as ls
importlib.reload(ls)
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from numpy import cos, sin, pi

"""
Test and compare the result of timedependent simulation vs original time independent simulation
"""
# %%

def setup_bus():
    bus = ls.Bus()
    bus.data.options.set_scenario('baseline')
    bus.data.import_catalog("C:/Users/Stoephu/Projects/SemesterProject2021/LIFEsim/output/suitable.hdf5")
    return bus

# %%
bus_o = setup_bus()  # original sim

bus_t = setup_bus()  # new sim

# %%
# ---------- MEASUREMENT SIMULATION ---------- copied from lifesim_demo.py by fdannert

# create modules and add to bus
inst_o = ls.Instrument(name='inst')
bus_o.add_module(inst_o)

transm_o = ls.TransmissionMap(name='transm')
bus_o.add_module(transm_o)

exo_o = ls.PhotonNoiseExozodi(name='exo')
bus_o.add_module(exo_o)
local_o = ls.PhotonNoiseLocalzodi(name='local')
bus_o.add_module(local_o)
star_o = ls.PhotonNoiseStar(name='star')
bus_o.add_module(star_o)

# connect all modules
bus_o.connect(('inst', 'transm'))
bus_o.connect(('inst', 'exo'))
bus_o.connect(('inst', 'local'))
bus_o.connect(('inst', 'star'))

bus_o.connect(('star', 'transm'))

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

# %%
inst_o.get_snr()

# %%
bus_o.data.catalog["snr_1h"]

# %%


def test_anomalies(n_angles):
    snrs = []
    true_anomalies = np.arange(0, 2 * np.pi, np.pi * 2 / n_angles)
    for anomaly in true_anomalies:
        bus_t.data.catalog["theta_p"] = anomaly
        inst.get_snr_t()
        snrs.append(bus_t.data.catalog[["theta_p", "snr_1h", "inc_p"]].copy())
    snrs
# %%


def plot_transmission_curves(rotation_time=1, rotation=1, rotation_steps=360):
    # TODO support more rotation times
    transmission_curves = []
    for i in enumerate(bus_t.data.catalog):
        tr_chop_one_wv = inst.get_transmission_curve(i)[0][0, 0]
        transmission_curves.append(tr_chop_one_wv)
        plt.plot(np.linspace(0, 2*pi, 360), tr_chop_one_wv, label=str(i))
    plt.legend()
    plt.show()
