import pandas as pd
import numpy as np
from astropy import units as u

from lifesim.util.radiation import black_body

def get_integration_time(temp_p,
                         radius_p,
                         angsep,
                         distance_s,
                         temp_s,
                         radius_s,
                         z,
                         lat_s,
                         spec_res,
                         wl_optimized,
                         target_snr,
                         opt_wl,
                         bus,
                         instrument,
                         use_blackbody,
                         return_reference=False,
                         return_time=False,):
    """
    Function to get the integration time for a given target SNR
    """

    # ---------- Creating the planet ----------

    if use_blackbody:
        fgamma = (black_body(mode='planet',
                             bins=bus.data.inst['wl_bins'],
                             width=bus.data.inst['wl_bin_widths'],
                             temp=temp_p,
                             radius=radius_p,
                             distance=distance_s
                             )
                  / bus.data.inst['wl_bin_widths']
                  * u.photon / u.second / (u.meter ** 3))

        flux_planet_spectrum = [bus.data.inst['wl_bins'] * u.meter, fgamma]

    else:
        data = pd.read_csv(
            '/Users/fdannert/Documents/projects/InLIFEsim/working/nice_requirements/input_data/Earth_PRTunits_10pc.txt',
            header=None, sep='\s+')

        lam_PRT = data[0].values * u.micron
        f_PRT = data[1].values * u.erg / u.cm ** 2 / u.s / u.Hz
        f_lifesim = f_PRT.to(u.photon / u.m ** 2 / u.s / u.micron, \
                             equivalencies=u.spectral_density(lam_PRT))

        f_lifesim = f_lifesim.to(u.photon / u.s / u.meter ** 3)
        lam_lifesim = lam_PRT.to(u.meter)

        # scale planet flux to distance
        f_lifesim *= (10 / distance_s) ** 2

        flux_planet_spectrum = [lam_lifesim, f_lifesim]

    bus.data.options.set_manual(wl_optimal=opt_wl)
    bus.data.options.set_manual(spec_res=spec_res)
    # bus.modules['inst'].apply_options()

    res_in = instrument.get_spectrum(temp_s=temp_s,
                                          radius_s=radius_s,
                                          distance_s=distance_s,
                                          lat_s=lat_s,
                                          z=z,
                                          angsep=angsep,
                                          flux_planet_spectrum=flux_planet_spectrum,
                                          integration_time=24*60*60,
                                          exposure_time=60*2,
                                          n_rot=1)

    snr_fundamental = res_in['snr'] * np.sqrt(1/24)

    print('In 1h of integration time:')
    print('Bulk SNR = ' + str(np.round(np.sqrt(np.sum(snr_fundamental ** 2)), 2)))
    if wl_optimized == 'bulk':
        integration_time_new = 60 * 60 * (target_snr / np.sqrt(np.sum(snr_fundamental ** 2))) ** 2
    else:
        # find appropriate wavelength bin from wl_optimized
        wl_id = np.argmin(np.abs(res_in.index.to_numpy() - wl_optimized))
        print(
            'SNR @ ' + str(np.round(res_in.index.to_numpy()[wl_id] * 1e6, 2)) + 'µm = ' + str(np.round(snr_fundamental[wl_id], 2)))
        integration_time_new = 60 * 60 * (target_snr / snr_fundamental[wl_id]) ** 2

    n_rot = int(integration_time_new/24/60/60)  # Convert to integer (truncates towards zero)
    if n_rot % 2 == 0:
        n_rot -= 1  # If even, decrement by 1 to make it odd
    n_rot = np.max((n_rot, 1))

    res_in = instrument.get_spectrum(temp_s=temp_s,
                                          radius_s=radius_s,
                                          distance_s=distance_s,
                                          lat_s=lat_s,
                                          z=z,
                                          angsep=angsep,
                                          flux_planet_spectrum=flux_planet_spectrum,
                                          integration_time=integration_time_new,
                                          exposure_time=60*2,
                                          n_rot=n_rot)

    snr_fundamental = res_in['snr']

    print('In ' + str(np.round(integration_time_new / (24 * 60 * 60), 2)) + 'd of integration time:')
    print('Bulk SNR = ' + str(np.round(np.sqrt(np.sum(snr_fundamental ** 2)), 2)))
    if not wl_optimized == 'bulk':
        print(
            'SNR @ ' + str(np.round(res_in.index.to_numpy()[wl_id] * 1e6, 2)) + 'µm = ' + str(np.round(snr_fundamental[wl_id], 2)))

    if return_time:
        return integration_time_new