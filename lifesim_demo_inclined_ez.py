import astropy.units as u
import numpy as np

import lifesim

# ---------- Set-Up ----------

# parameters for the simulation
temp_p = 285.  # planet temperature in K
radius_p = 1.  # planet radius in Earth radii
distance_s = 10.  # distance to target system in pc
temp_s = 5778.
radius_s = 1.
lat_s = 0.79
z = 3.
angsep = 0.1
integration_time = 200000


# create bus
bus = lifesim.Bus()

# setting the options
bus.data.options.set_scenario('baseline')

# set options manually
#bus.data.options.set_manual(rotation_steps=15)

# ---------- Creating the Instrument ----------

# create modules and add to bus
instrument = lifesim.Instrument(name='inst')
bus.add_module(instrument)

transm = lifesim.TransmissionMap(name='transm')
bus.add_module(transm)

exoi = lifesim.PhotonNoiseExozodiInclined(name='exoi')
bus.add_module(exoi)

exo = lifesim.PhotonNoiseExozodiInclined(name='exo')
bus.add_module(exo)

local = lifesim.PhotonNoiseLocalzodi(name='local')
bus.add_module(local)
star = lifesim.PhotonNoiseStar(name='star')
bus.add_module(star)

# connect all modules
bus.connect(('inst', 'transm'))
bus.connect(('inst', 'exoi'))
bus.connect(('inst', 'local'))
bus.connect(('inst', 'star'))

bus.connect(('star', 'transm'))

# ---------- Running the Simulation ----------

# the ff. lines create the black body spectrum for the exoplanet
# do not change any parameters here
wl_edge = 1.
wl_bins = []
wl_bin_widths = []
wl_bin_edges = [wl_edge]
R = 200
wl_max = 30

while wl_edge < wl_max:

    # set the wavelength bin width according to the spectral resolution
    wl_bin_width = wl_edge / R / \
                   (1 - 1 / R / 2)

    # make the last bin shorter when it hits the wavelength limit
    if wl_edge + wl_bin_width > wl_max:
        wl_bin_width = wl_max - wl_edge

    # calculate the center and edges of the bins
    wl_center = wl_edge + wl_bin_width / 2
    wl_edge += wl_bin_width

    wl_bins.append(wl_center)
    wl_bin_widths.append(wl_bin_width)
    wl_bin_edges.append(wl_edge)

# convert everything to [m]
wl_bins = np.array(wl_bins) * 1e-6  # in m
wl_bin_widths = np.array(wl_bin_widths) * 1e-6  # in m

fgamma = lifesim.util.radiation.black_body(mode='planet',
                                           bins=wl_bins,
                                           width=wl_bin_widths,
                                           temp=temp_p,
                                           radius=radius_p,
                                           distance=distance_s
                                           ) \
         / wl_bin_widths \
         * u.photon / u.second / (u.meter ** 3)

flux_planet_spectrum = [wl_bins * u.meter, fgamma]

# --- End of Black Body calculation for planet ---

# run simulation. This function assigns every planet an SNR for 1 hour of integration time. Since
# we are currently only simulating photon noise, the SNR will scale with the integration time as
# sqrt(t)

(spectrum,
 flux_planet,
 noise) = instrument.get_spectrum_inclined(temp_s=temp_s,  # for description, view get_spectrum documentation
                                           radius_s=radius_s,
                                           distance_s=distance_s,
                                           lat_s=lat_s,
                                           z=z,
                                           angsep=angsep,
                                           flux_planet_spectrum=flux_planet_spectrum,
                                           integration_time=integration_time,
                                           safe_mode=False,
                                           inclination_ez=45./180.*np.pi,
                                           ascending_node_ez=45./180.*np.pi)

snri = np.sqrt(np.sum(spectrum[1] ** 2))

bus.disconnect(('inst', 'exoi'))
bus.connect(('inst', 'exo'))

(spectrum,
 flux_planet,
 noise) = instrument.get_spectrum_inclined(temp_s=temp_s,  # for description, view get_spectrum documentation
                                           radius_s=radius_s,
                                           distance_s=distance_s,
                                           lat_s=lat_s,
                                           z=z,
                                           angsep=angsep,
                                           flux_planet_spectrum=flux_planet_spectrum,
                                           integration_time=integration_time,
                                           safe_mode=False)

snr = np.sqrt(np.sum(spectrum[1] ** 2))

print(snr)
print(snri)

# overall snr given by root squared sum SNR = sqrt(sum(snr_n^2))

a=1
