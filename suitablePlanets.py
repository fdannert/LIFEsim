#%%
import lifesim as ls
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

#%%
def min_separation(large_omega_p, inc_p):
    """
    Returns value from 0 to 1, 
    where 0 is the planet travels exactly over the middle of the star(edge_on),
    where 1 is never comes closer to the star in the projected separation(face_on).
    tries to give a measure of edge_on'ness
    """
    separation = 1
    separation = np.abs(np.sin(large_omega_p)*np.sin(inc_p))
    return separation

#%%

print("load")
bus = ls.Bus()
bus.data.options.set_scenario('baseline')
bus.data.import_catalog(input_path='C:/Users/Stoephu/Projects/SemesterProject2021/LIFEsim/output/TestSet.hdf5')
db = bus.data.catalog
suitable = pd.DataFrame(data=None, columns=db.columns)
detected_n = db[db["detected"]]
max_period = detected_n["p_orb"].max()
min_period = detected_n["p_orb"].min()
m_stars = detected_n[detected_n["stype"] == 4]
edge_on = detected_n[(min_separation(detected_n["large_omega_p"],detected_n["inc_p"]) <= 0.01)]
face_on = detected_n[(min_separation(detected_n["large_omega_p"],detected_n["inc_p"]) >= 0.99)]



#%%
suitable = suitable.append(edge_on[edge_on["p_orb"] <= (min_period + 0.1)].sort_values("snr_1h").iloc[[5]])
suitable = suitable.append(face_on[face_on["p_orb"] <= (min_period + 0.1)].sort_values("snr_1h").iloc[[9]])
#%%
bus.data.catalog = suitable
bus.data.export_catalog("C:/Users/Stoephu/Projects/SemesterProject2021/LIFEsim/output/suitable.hdf5")

#%%
bus2 = ls.Bus()
bus2.data.options.set_scenario('baseline')
bus2.data.import_catalog("C:/Users/Stoephu/Projects/SemesterProject2021/LIFEsim/output/suitable.hdf5")
bus.data.catalog

#%%
print("plot")
#bins_period = np.arange(min_period, max_period, 1/24)
bins_period = np.arange(min_period, 1, 1/24)
plt.hist(db["p_orb"], bins_period)
plt.show()

bins_inc = np.arange(0, 180, 1)
plt.hist(db["inc_p"]/np.pi * 180, bins_inc)
plt.xlabel("Degree")
plt.ylabel("N")
plt.show()


"""
'radius_p' radius/radius_earth
'p_orb' Days
'mass_p', M/Mearth
'ecc_p',
'inc_p', rad 0:face on pi/2: edge on
'large_omega_p', rad
'small_omega_p', rad
'theta_p', rad
'albedo_bond',
'albedo_geom_vis',
'albedo_geom_mir',
'z',
'semimajor_p',
'sep_p',
'angsep',
'maxangsep',
'flux_p',
'fp',
'temp_p',
'radius_s',
'mass_s',
'temp_s',
'distance_s',
'ra',
'dec',
'nuniverse',
'nstar',
'stype',
'id'
"""

# %%
