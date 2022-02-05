# %%
import lifesim as ls
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
%matplotlib


# %%
"""
Search for suitable Planets:
1 Edge on short period
3 Face on: 1 Counter rotation, 1 "stand still", 1 Co rotation 5h/12h
# TODO stats how common those planets

"""
print("load")
bus = ls.Bus()
bus.data.options.set_scenario('baseline')
bus.data.import_catalog(
    input_path='./output/TestSet.hdf5')
db = bus.data.catalog
suitable = pd.DataFrame(data=None, columns=db.columns)
detected_n = db[(db["snr_1h"] > 10) & db["detected"] == True]
# detected_n = db.copy()
max_period = detected_n["p_orb"].max()
min_period = detected_n["p_orb"].min()
m_stars = detected_n[detected_n["stype"] == 4]
face_on_a = detected_n[(detected_n["inc_p"] <= 0.2)]
face_on_c = detected_n[detected_n["inc_p"] >= (np.pi - 0.2)]
edge_on = detected_n[(detected_n["inc_p"] >= (np.pi/2 - 0.02))
                     & (detected_n["inc_p"] >= (np.pi/2 + 0.02))]
face_max_period = face_on_c["p_orb"].max()


# %%
# real edge on 12h
suitable = suitable.append(edge_on[edge_on["p_orb"] <= (
    min_period + 0.1)].sort_values("snr_1h").iloc[[0]])
# real face on clockwise rotation 12h
suitable = suitable.append(face_on_c[face_on_c["p_orb"] <= (
    min_period + 0.1)].sort_values("snr_1h").iloc[[0]])
# artificial stand still (period >> 12h) same separation
artificial = face_on_c[face_on_c["p_orb"] <= (
    min_period + 0.1)].sort_values("snr_1h").iloc[[0]]
artificial["p_orb"] = 1000
suitable = suitable.append(artificial)
# real face on anticlockwise rotation 12h
suitable = suitable.append(face_on_a[face_on_a["p_orb"] <= (
    min_period + 0.1)].sort_values("snr_1h").iloc[[0]])
# real face on clockwise period >> 12h => bigger separation than the others
suitable = suitable.append(
    face_on_c[face_on_c["p_orb"] >= (90)].sort_values("snr_1h").iloc[[0]])

# edge_on in hz
suitable = suitable.append(edge_on[edge_on["habitable"]].sort_values("maxangsep").iloc[[-1]])


# %%
# TODO Create artificial planets to test limits.


def export():
    bus.data.catalog = suitable
    bus.data.export_catalog(
        "./output/suitable.hdf5")


# %%
def pop():
    print("plot")
    max_period = detected_n["p_orb"].max()
    min_period = detected_n["p_orb"].min()
    detected = db[db["detected"] == True]
    #bins_period = np.arange(min_period, max_period, 1/24)
    bins_period = np.arange(min_period, max_period, 1)
    plt.hist(db["p_orb"], bins_period, label="total population")
    plt.hist(detected["p_orb"], bins_period, label="detected", color="green")
    plt.hist(detected[detected["stype"] == 4]["p_orb"],
             bins_period, label="detected M", color="red")
    plt.xlim(-1, 100)
    plt.legend()
    plt.show()
    ms = detected[detected["stype"] == 4]
    print(len(ms[ms["p_orb"] <= 50]),"/",len(ms),"=",len(ms[ms["p_orb"] <= 50])/len(ms))

# %%


def inc():
    bins_inc = np.arange(0, 190, 5)
    plt.hist(db["inc_p"]/np.pi * 180, bins_inc)
    plt.xlabel("Degree")
    plt.ylabel("N")
    plt.show()


# np.cos(dtheta*t)*np.cos(dphi*t) = 1/2 (cos((dtheta - dphi)*t) + cos((dtheta + dphi)*t))
"""

"""
# 'radius_p' radius/radius_earth
# 'p_orb' Days
# 'mass_p', M/Mearth
# 'ecc_p',
# 'inc_p', rad 0:face on pi/2: edge on
# 'large_omega_p', rad
# 'small_omega_p', rad
# 'theta_p', rad
# 'albedo_bond',
# 'albedo_geom_vis',
# 'albedo_geom_mir',
# 'z',
# 'semimajor_p',
# 'sep_p',
# 'angsep',
# 'maxangsep',
# 'flux_p',
# 'fp',
# 'temp_p',
# 'radius_s',
# 'mass_s',
# 'temp_s',
# 'distance_s',
# 'ra',
# 'dec',
# 'nuniverse',
# 'nstar',
# 'stype',
# 'id' 
# "hz" in au

# %%
