import os
import warnings
from typing import Union
from importlib.resources import files
from io import BytesIO
import requests

import yaml
import astropy.units as u
import numpy as np
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
import planets
import pandas as pd
import whereistheplanet
from astropy.io.votable import parse_single_table
from astropy.constants import h, c

import lifesim
from lifesim.core.core import add_numpy_representers
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

class SourceConfig(object):

    def __init__(self,
                 sources_config_file: str):
        self.sources_config_file = sources_config_file
        if os.path.isfile(sources_config_file):
            # parse the system from the sources config file (yaml)
            with open(sources_config_file) as file:
                self.sources = yaml.load(file, Loader=yaml.FullLoader)

            warnings.warn('You have opened an existing config file which will be overwritten.')

        else:
            self.sources = {}

    def save(self):
        add_numpy_representers()
        with open(self.sources_config_file, 'w') as file:
            yaml.dump(self.sources, file)

    def add_star(self, star_name: str):

        print(f"--- Searching for: {star_name} ---")

        if star_name.casefold() == 'sol':
            print('Selecting 10 pc Sun Twin')

            self.sources['star'] = {'temperature': 5778,
                                    'radius': 1.,
                                    'distance': 10.,
                                    'latitude': 0.78,
                                    'l_sun': 1.,
                                    'name': 'Sol'}
        else:

            # 1. Resolve Name to Coordinates
            try:
                coord = SkyCoord.from_name(star_name)
            except Exception:
                raise ValueError("Error: Name not recognized by SIMBAD/NED.")

            custom_simbad = Simbad()
            custom_simbad.add_votable_fields('sp_type')

            print(f"Resolving {star_name} in SIMBAD...")
            simbad_result = custom_simbad.query_object(star_name)

            # Extract Spectral Type and Coordinates
            spec_type = simbad_result['sp_type'][0]

            # 2. Configure Vizier to search the TESS Input Catalog (TIC)
            # This catalog is specifically aggregated to avoid "empty" values for bright stars
            v = Vizier(
                catalog="IV/38",
                columns=['TIC', 'Teff', 'Rad', 'Dist', 'SpType', 'Vmag', 'ELAT', '_r'],
            )
            # This is the correct way to sort by distance from the center coordinates
            v.ucd = "pos.angDistance"
            v.ROW_LIMIT = 500

            # 3. Search with a wider radius (30 arcseconds)
            result = v.query_region(coord, radius=30 * u.arcsec)

            if not result:
                return f"No matches found in TIC for {star_name} within 30 arcsec."

            # 4. Sort by V-magnitude to get the brightest (most likely) star
            table = result[0]
            table.sort('_r')

            # Extract the top match
            best_match = table[0]

            print(f"Found: TIC {best_match['TIC']}")
            print(f"Spectral Type: {spec_type}")
            print(f"Temperature:   {best_match['Teff']} K")
            print(f"Radius:        {best_match['Rad']} Solar Radii")
            print(f"Distance:      {best_match['Dist']} pc")

            lum_s = best_match['Rad'] ** 2 * (best_match['Teff'] / 5780) ** 4

            print('Luminosity:     {:.2f} L_sun'.format(lum_s))
            print('')

            self.sources['star'] = {'temperature': best_match['Teff'],
                                    'radius': best_match['Rad'],
                                    'distance': best_match['Dist'],
                                    'latitude': np.deg2rad(best_match['ELAT']),
                                    'l_sun': lum_s,
                                    'name': star_name}

    def planet_solar_system(self,
                            planet_name: str,):
        # NASA Official Effective Temperatures (K)
        # Source: https://nssdc.gsfc.nasa.gov/planetary/factsheet/
        eff_temps = {
            'mercury': 437,
            'venus': 232,
            'earth': 255,
            'mars': 209,
            'jupiter': 88,
            'saturn': 95,
            'uranus': 58,
            'neptune': 55.5,
        }

        geom_albedos = {
            'mercury': 0.142,
            'venus': 0.689,
            'earth': 0.24,
            'mars': 0.17,
            'jupiter': 0.538,
            'saturn': 0.499,
            'uranus': 0.488,
            'neptune': 0.442,
        }

        print('--- Resolving: {} ---'.format(planet_name))
        # Mapping for flexibility: handles full names and common abbreviations
        lookup = {
            'mercury': planets.Mercury, 'm': planets.Mercury,
            'venus': planets.Venus, 'v': planets.Venus,
            'earth': planets.Earth, 'e': planets.Earth,
            'mars': planets.Mars, 'ma': planets.Mars,
            'jupiter': planets.Jupiter, 'j': planets.Jupiter,
            'saturn': planets.Saturn, 's': planets.Saturn,
            'uranus': planets.Uranus, 'u': planets.Uranus,
            'neptune': planets.Neptune, 'n': planets.Neptune,
        }

        # Normalize input: lowercase and strip whitespace
        key = str(planet_name).strip().lower()

        if key in lookup:
            p = lookup[key]
            print('Name: {}'.format(key.capitalize()))
            print('Effective Temperature: {} K'.format(eff_temps[key]))
            print('Radius: {} Earth radii'.format(p.R / planets.Earth.R))
            print('Semi-Major Axis: {} AU'.format(p.rAU))
            print('Geometric Albedo: {}'.format(geom_albedos[key]))
        else:
            return f"Error: '{planet_name}' not found. Are you sure that's a planet?"

        # if self.lum_s is None:
        #     print('Warning: Stellar luminosity not set. Cannot calculate scaled semimajor axis.')
        try:
            # Calculate scaled semimajor axis (a/R_star)
            sma = p.rAU * np.sqrt(self.sources['star']['l_sun'])
            print('Scaled Semi-Major Axis: {:.2f} AU'.format(sma))
        except:
            sma = p.rAU
            warnings.warn('Stellar luminosity not set. Will use solar reference.')

        self.sources['planet'] = {'radius': p.R / planets.Earth.R,
                                'temperature': eff_temps[key],
                                'sma': sma,
                                'geom_albedo': geom_albedos[key],
                                'name': key}

    def from_lband(self,
                   planet_name: str,
                   date: str):
        catalog = pd.read_csv(str(files("lifesim.analysis") / "etc_data" / "reliable_photometry.csv"))
        planet = catalog.loc[catalog['planet_name'] == planet_name]
        if len(planet) == 0:
            raise ValueError('Could not find planet named {}'.format(planet_name))
        star_name = ' '.join(planet['planet_name'].values[0].split(' ')[:-1])

        self.add_star(star_name=star_name)

        # retrieve the angular separation of the planet
        if planet['witp_name'].values[0] == 'None':
            print('No orbit data available, taking last known angular separation.')
            angsep = planet['angsep_arcsep'].values[0]
        else:
            print('Retrieving angular separation using WhereIsThePlanet.')
            _, _, sep_args, _ = whereistheplanet.predict_planet(planet['witp_name'].values[0], date)
            angsep = sep_args[0] * 1e-3
        print('Using angular separation: {}'.format(angsep))

        # retrieve the photon flux of the planet
        ph_flux = get_flux_filter(filter_id=planet['filter_id'].values[0], magnitude=planet['l_band_mag'].values[0])
        print('L-band magnitude of {} resulting in photon flux of {}'.format(planet['l_band_mag'].values[0], ph_flux))

        self.sources['planet'] = {'name': planet_name,
                                  'angsep': angsep,
                                  'sma': angsep * self.sources['star']['distance'],
                                  'ph_flux': ph_flux,  # in ph s-1 m-2 µm-1
                                  }
        a=1


# Get filter metadata (zero point, central wavelength, etc.)
def get_svo_filter_info(filter_id):
    url = "http://svo2.cab.inta-csic.es/theory/fps/fps.php"
    params = {"ID": filter_id, "FORMAT": "votable&VERB=2"}
    r = requests.get(url, params=params)
    vot = parse_single_table(BytesIO(r.content))
    # Metadata is in VOTable PARAM fields
    params_dict = {p.name: p.value for p in vot.params}
    return params_dict

def jy_to_photons(flux_jy, wavelength):
    """
    Convert flux from Jy to photons / m² / µm.

    Parameters
    ----------
    flux_jy : float or astropy.units.Quantity
        Flux in Jansky. If float, Jy are assumed.
    wavelength : float or astropy.units.Quantity
        Reference wavelength. If float, µm are assumed.

    Returns
    -------
    astropy.units.Quantity
        Flux in photons / m² / µm
    """
    if not isinstance(flux_jy, u.Quantity):
        flux_jy = flux_jy * u.Jy
    if not isinstance(wavelength, u.Quantity):
        wavelength = wavelength * u.micron

    # Step 1: Jy → W / m² / Hz
    flux = flux_jy.to(u.W / u.m**2 / u.Hz)

    # Step 2: to per wavelenth instead of per frequency
    flux = flux.to(u.J / u.m**2 / u.um / u.s, equivalencies=u.spectral_density(wavelength))

    # # Step 3: to photons
    flux = (flux *u.photon / (h * c / wavelength)).to(u.photon / u.s / u.m**2 / u.micron)
    #
    # # Step 3: W / m² / µm → photons / m² / µm  (divide by photon energy E = hc/λ)
    # photon_flux = (flux_wUm / (h * c / wavelength)).to(u.photon / u.m**2 / u.um)

    return flux

def get_flux_filter(filter_id, magnitude):
    '''
    Converts filter magnitude into photon flux with units ph s-1 m-2 µm-1
    :param filter_id:
    :param magnitude:
    :return:
    '''
    info = get_svo_filter_info(filter_id)
    wavelength = (float(info["WavelengthCen"]) * u.AA).to(u.micron)
    zp_jy = float(info["ZeroPoint"]) * u.Unit(info["ZeroPointUnit"])

    flux_ph = jy_to_photons(zp_jy, wavelength) * 10**(-magnitude/2.5)
    return flux_ph
