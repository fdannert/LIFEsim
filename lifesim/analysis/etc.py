from typing import Union

import yaml
import astropy.units as u
import numpy as np

import lifesim
from lifesim.util.radiation import black_body

def etc(
        instrument_config_file:str,
        sources_config_file:str,
        target_snr: float,
        wl_optimized: Union[str, float] = 'bulk',
):
    bus = lifesim.Bus()

    bus.build_from_config(instrument_config_file)

    instrument = lifesim.Instrument(name='inst')
    bus.add_module(instrument)

    if bus.data.options.array['num_apertures'] == 2:
        transm = lifesim.TransmissionMapSBW(name='transm')
    elif bus.data.options.array['num_apertures'] == 4:
        transm = lifesim.TransmissionMap(name='transm')
    bus.add_module(transm)

    exo = lifesim.PhotonNoiseExozodi(name='exo')
    bus.add_module(exo)
    local = lifesim.PhotonNoiseLocalzodi(name='local')
    bus.add_module(local)
    star = lifesim.PhotonNoiseStar(name='star')
    bus.add_module(star)

    # connect all modules
    bus.connect(('inst', 'transm'))
    bus.connect(('inst', 'exo'))
    bus.connect(('inst', 'local'))
    bus.connect(('inst', 'star'))

    bus.connect(('star', 'transm'))

    instrument.apply_options()

    # parse the system from the sources config file (yaml)
    with open(sources_config_file) as file:
        sources = yaml.load(file, Loader=yaml.FullLoader)

    # ---------- Creating the planet ----------

    fgamma = (black_body(mode='planet',
                         bins=bus.data.inst['wl_bins'],
                         width=bus.data.inst['wl_bin_widths'],
                         temp=sources['planet']['temperature'],
                         radius=sources['planet']['radius'],
                         distance=sources['star']['distance'],
                         )
              / bus.data.inst['wl_bin_widths']
              * u.photon / u.second / (u.meter ** 3))

    flux_planet_spectrum = [bus.data.inst['wl_bins'] * u.meter, fgamma]
    #
    # bus.modules['inst'].adjust_bl_to_hz(hz_center=sources['planet']['separation'],
    #                                     distance_s=sources['star']['distance'],)
    bus.modules['inst'].apply_options()

    snr, _, _ = instrument.get_spectrum(temp_s=sources['star']['temperature'],
                                        radius_s=sources['star']['radius'],
                                        distance_s=sources['star']['distance'],
                                        lat_s=sources['star']['latitude'],
                                        z=sources['exozodi']['z'],
                                        angsep=sources['planet']['separation'] / sources['star']['distance'],
                                        flux_planet_spectrum=flux_planet_spectrum,
                                        integration_time=60 * 60,
                                        safe_mode=False)

    snr_fundamental = snr[1]

    print('In 1h of integration time:')
    print('Bulk SNR = ' + str(np.round(np.sqrt(np.sum(snr_fundamental ** 2)), 2)))
    if wl_optimized == 'bulk':
        integration_time_new = 60 * 60 * (target_snr / np.sqrt(np.sum(snr_fundamental ** 2))) ** 2
    else:
        # find appropriate wavelength bin from wl_optimized
        wl_id = np.argmin(np.abs(snr[0] - wl_optimized))
        print('SNR @ ' + str(np.round(snr[0][wl_id] * 1e6, 2)) + 'µm = ' + str(np.round(snr_fundamental[wl_id], 2)))
        integration_time_new = 60 * 60 * (target_snr / snr_fundamental[wl_id]) ** 2

    snr, _, _ = instrument.get_spectrum(temp_s=sources['star']['temperature'],
                                        radius_s=sources['star']['radius'],
                                        distance_s=sources['star']['distance'],
                                        lat_s=sources['star']['latitude'],
                                        z=sources['exozodi']['z'],
                                        angsep=sources['planet']['separation'] / sources['star']['distance'],
                                        flux_planet_spectrum=flux_planet_spectrum,
                                        integration_time=integration_time_new,
                                        safe_mode=False)

    snr_fundamental = snr[1]

    print('In ' + str(np.round(integration_time_new / (24 * 60 * 60), 2)) + 'd of integration time:')
    print('Bulk SNR = ' + str(np.round(np.sqrt(np.sum(snr_fundamental ** 2)), 2)))
    if not wl_optimized == 'bulk':
        print('SNR @ ' + str(np.round(snr[0][wl_id] * 1e6, 2)) + 'µm = ' + str(np.round(snr_fundamental[wl_id], 2)))

    print('Nulling baseline used: ' + str(np.round(bus.data.inst['bl'], 1)) + ' m')

    return integration_time_new, bus