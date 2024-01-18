import numpy as np

h = 6.62607e-34
k = 1.380649e-23
c = 2.99792e+8
m_per_pc = 3.086e+16
m_per_au = 1.49597870691e+11

radius_sun = 6.947e+8
temp_sun = 5778.

radius_earth = 6.371e+6
temp_earth = 276.

#some conversions and constants (signal extraction)
au2par = 4.8481368111358*10**-6
par2au = 1/au2par
rad2arcsec = 180/np.pi*3600
arcsec2rad = 1/rad2arcsec

