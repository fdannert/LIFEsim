# Top of run_tse.py, BEFORE any other imports.
# astroquery.gaia.core runs `Gaia = GaiaClass()` at import time, which
# calls get_status_messages() -> a blocking HTTP GET to ESA's TAP server.
# ESA's archive is currently hanging on that endpoint, so we stub the GET.
from astroquery.utils.tap.conn.tapconn import TapConn

def _noop_get(self, *args, **kwargs):
    class _Resp:
        status = 500
        def __iter__(self): return iter(())
    return _Resp()

TapConn.execute_tapget = _noop_get

import os
import sys
import warnings
from typing import Union
from importlib.resources import files
from io import BytesIO
import requests
import contextlib

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
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich import box as rich_box
from rich.text import Text

import lifesim
from lifesim.core.core import add_numpy_representers
from lifesim.util.radiation import black_body

def etc(
        instrument_config_file: str,
        sources_config_file: str,
        target_snr: float,
        wl_optimized: Union[str, float] = 'bulk',
        verbose: bool = True,
        additional_options: dict = None,
        return_noise: bool = False
):
    console = Console()

    def _print(*args, **kwargs):
        if verbose:
            console.print(*args, **kwargs)

    def _status(message: str):
        _print(f"  [dim]›[/dim] {message}")

    _print(Panel(
        f"[bold white]Exposure Time Calculator[/bold white]  [dim]target SNR {target_snr}[/dim]",
        border_style="bright_blue", expand=False
    ))

    # ---------- Build bus ----------
    _status("Building instrument from config...")
    bus = lifesim.Bus()
    bus.build_from_config(instrument_config_file)

    if additional_options is not None:
        bus.data.options.set_manual(**additional_options)

    instrument = lifesim.Instrument(name='inst')
    bus.add_module(instrument)

    _status("Configuring transmission map...")
    if bus.data.options.array['num_apertures'] == 2:
        transm = lifesim.TransmissionMapSBW(name='transm')
    elif bus.data.options.array['num_apertures'] == 4:
        transm = lifesim.TransmissionMap(name='transm')
    bus.add_module(transm)

    _status("Adding noise modules (exozodi, localzodi, star)...")
    exo   = lifesim.PhotonNoiseExozodi(name='exo')
    local = lifesim.PhotonNoiseLocalzodi(name='local')
    star  = lifesim.PhotonNoiseStar(name='star')
    bus.add_module(exo)
    bus.add_module(local)
    bus.add_module(star)

    if ((bus.data.options.thermal['ota_temperature'] != 0.)
            or (bus.data.options.thermal['instrument_temperature'] != 0.)
            or (bus.data.options.thermal['detector_temperature'] != 0.)):
        _status("Adding therman noise modules...")
        mirror = lifesim.PhotonNoiseThermal(name='mirror')
        bus.add_module(mirror)
        bus.connect(('inst', 'mirror'))

    bus.connect(('inst', 'transm'))
    bus.connect(('inst', 'exo'))
    bus.connect(('inst', 'local'))
    bus.connect(('inst', 'star'))
    bus.connect(('star', 'transm'))

    instrument.apply_options()

    # ---------- Load sources ----------
    if type(sources_config_file) == str:
        _status(f"Loading source config: [dim]{sources_config_file}[/dim]")
        with open(sources_config_file) as file:
            sources = yaml.load(file, Loader=yaml.FullLoader)
    else:
        sources = sources_config_file

    # ---------- Planet flux ----------
    if 'ph_flux' in sources['planet']:
        _status("Using measured L-band photon flux for planet spectrum...")
        flux_planet_spectrum = [
            bus.data.inst['wl_bins'] * u.meter,
            np.ones_like(bus.data.inst['wl_bins']) * sources['planet']['ph_flux'] * 1e6 * u.photon / u.second / (u.meter ** 3)
        ]
    else:
        _status("Computing planet blackbody spectrum...")
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

    bus.modules['inst'].apply_options()

    # ---------- 1h integration ----------
    _status("Running 1h integration to estimate SNR...")
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
    bulk_snr_1h = np.round(np.sqrt(np.sum(snr_fundamental ** 2)), 2)
    _status(f"Bulk SNR in 1h = [bold]{bulk_snr_1h}[/bold]")

    # ---------- Required integration time ----------
    if wl_optimized == 'bulk':
        integration_time_new = 60 * 60 * (target_snr / np.sqrt(np.sum(snr_fundamental ** 2))) ** 2
        _status(f"Optimising for bulk SNR → target integration time: "
                f"[bold]{np.round(integration_time_new / (24 * 60 * 60), 2)}d[/bold]")
    else:
        wl_id = np.argmin(np.abs(snr[0] - wl_optimized))
        wl_actual = np.round(snr[0][wl_id] * 1e6, 2)
        snr_at_wl  = np.round(snr_fundamental[wl_id], 2)
        integration_time_new = 60 * 60 * (target_snr / snr_fundamental[wl_id]) ** 2
        _status(f"SNR @ {wl_actual}µm in 1h = [bold]{snr_at_wl}[/bold]  →  "
                f"target integration time: [bold]{np.round(integration_time_new / (24 * 60 * 60), 2)}d[/bold]")

    # ---------- Final integration ----------
    _status(f"Running final integration ({np.round(integration_time_new / (24 * 60 * 60), 2)}d)...")
    snr, _, noise = instrument.get_spectrum(temp_s=sources['star']['temperature'],
                                        radius_s=sources['star']['radius'],
                                        distance_s=sources['star']['distance'],
                                        lat_s=sources['star']['latitude'],
                                        z=sources['exozodi']['z'],
                                        angsep=sources['planet']['separation'] / sources['star']['distance'],
                                        flux_planet_spectrum=flux_planet_spectrum,
                                        integration_time=integration_time_new,
                                        safe_mode=True)

    snr_fundamental = snr[1]

    # ---------- Results table ----------
    result_table = Table(
        box=rich_box.ROUNDED,
        border_style="bright_blue",
        show_header=True,
        header_style="bold bright_white",
        title="[bold]Results[/bold]"
    )
    result_table.add_column("Parameter", min_width=28)
    result_table.add_column("Value")

    int_time_seconds = integration_time_new
    if int_time_seconds >= 365.25 * 24 * 60 * 60:
        int_time_val = np.round(int_time_seconds / (365.25 * 24 * 60 * 60), 2)
        int_time_unit = "yrs"
    elif int_time_seconds >= 24 * 60 * 60:
        int_time_val = np.round(int_time_seconds / (24 * 60 * 60), 2)
        int_time_unit = "d"
    elif int_time_seconds >= 60 * 60:
        int_time_val = np.round(int_time_seconds / (60 * 60), 2)
        int_time_unit = "h"
    elif int_time_seconds >= 60:
        int_time_val = np.round(int_time_seconds / 60, 2)
        int_time_unit = "min"
    else:
        int_time_val = np.round(int_time_seconds, 2)
        int_time_unit = "s"

    result_table.add_row("Integration time", f"{int_time_val} {int_time_unit}")

    result_table.add_row("Bulk SNR",
                         f"{np.round(np.sqrt(np.sum(snr_fundamental ** 2)), 2)}")
    if not wl_optimized == 'bulk':
        result_table.add_row(f"SNR @ {wl_actual} µm",
                             f"{np.round(snr_fundamental[wl_id], 2)}")
    result_table.add_row("Nulling baseline", f"{np.round(bus.data.inst['bl'], 1)} m")

    _print(result_table)

    if return_noise:
        return integration_time_new, bus, noise
    else:
        return integration_time_new, bus

@contextlib.contextmanager
def _suppress_output():
    """Suppresses all stdout and stderr within the context."""
    with open(os.devnull, 'w') as devnull:
        old_stdout, old_stderr = sys.stdout, sys.stderr
        sys.stdout, sys.stderr = devnull, devnull
        try:
            yield
        finally:
            sys.stdout, sys.stderr = old_stdout, old_stderr


class SourceConfig(object):

    def __init__(self,
                 sources_config_file: str,
                 verbose: bool = True):
        self._console = Console()
        self._verbose = verbose
        self.sources_config_file = sources_config_file
        if os.path.isfile(sources_config_file):
            with open(sources_config_file) as file:
                self.sources = yaml.load(file, Loader=yaml.FullLoader)
            warnings.warn('You have opened an existing config file which will be overwritten.')
        else:
            self.sources = {}

    def _print(self, *args, **kwargs):
        if self._verbose:
            self._console.print(*args, **kwargs)

    def _status(self, message: str):
        """Print a single-line status update with a leading indicator."""
        self._print(f"  [dim]›[/dim] {message}")

    def save(self):
        add_numpy_representers()
        with open(self.sources_config_file, 'w') as file:
            yaml.dump(self.sources, file)

        # --- Overview panel after saving ---
        star = self.sources.get('star', {})
        planet = self.sources.get('planet', {})

        overview = Table(
            box=rich_box.ROUNDED,
            border_style="bright_blue",
            show_header=True,
            # header_style="bold bright_white",
            title=f"[bold]Saved[/bold] [dim]{self.sources_config_file}[/dim]"
        )
        overview.add_column("Parameter", min_width=26)
        overview.add_column("Value")

        if star:
            overview.add_section()
            overview.add_row("[bold cyan]Star[/bold cyan]", "")
            overview.add_row("  Name",        str(star.get('name', '—')))
            overview.add_row("  Temperature", f"{star.get('temperature', '—')} K")
            overview.add_row("  Radius",      f"{star.get('radius', '—')} R☉")
            overview.add_row("  Distance",    f"{star.get('distance', '—')} pc")
            overview.add_row("  Luminosity",  f"{star.get('l_sun', '—')} L☉")

        if planet:
            overview.add_section()
            overview.add_row("[bold green]Planet[/bold green]", "")
            overview.add_row("  Name",                str(planet.get('name', '—')))
            if 'temperature' in planet:
                overview.add_row("  Temperature",     f"{planet.get('temperature')} K")
            if 'radius' in planet:
                overview.add_row("  Radius",          f"{planet.get('radius')} R⊕")
            overview.add_row("  Semi-major axis",     f"{planet.get('sma', '—')} AU")
            if 'geom_albedo' in planet:
                overview.add_row("  Geometric albedo",f"{planet.get('geom_albedo')}")
            if 'angsep' in planet:
                overview.add_row("  Angular sep.",    f"{planet.get('angsep')} arcsec")
            if 'ph_flux' in planet:
                overview.add_row("  Photon flux",     f"{planet.get('ph_flux')}")

        self._print(overview)

    def add_star(self, star_name: str):

        self._print(Panel(
            f"[bold white]Star:[/bold white] [cyan]{star_name}[/cyan]",
            border_style="bright_blue", expand=False
        ))

        if star_name.casefold() == 'sol':
            self._status("Using 10 pc Sun Twin")
            self.sources['star'] = {'temperature': 5778,
                                    'radius': 1.,
                                    'distance': 10.,
                                    'latitude': 0.78,
                                    'l_sun': 1.,
                                    'name': 'Sol',}
        else:
            try:
                self._status("Resolving coordinates via SIMBAD/NED...")
                coord = SkyCoord.from_name(star_name)
            except Exception:
                raise ValueError("Error: Name not recognized by SIMBAD/NED.")

            custom_simbad = Simbad()
            custom_simbad.add_votable_fields('sp_type')

            self._status("Querying SIMBAD for spectral type...")
            simbad_result = custom_simbad.query_object(star_name)
            spec_type = simbad_result['sp_type'][0]

            self._status("Querying TESS Input Catalog (TIC) via Vizier...")
            v = Vizier(
                catalog="IV/38",
                columns=['TIC', 'Teff', 'Rad', 'Dist', 'SpType', 'Vmag', 'ELAT', 'ELON', '_r'],
            )
            v.ucd = "pos.angDistance"
            v.ROW_LIMIT = 500
            result = v.query_region(coord, radius=30 * u.arcsec)

            if not result:
                return f"No matches found in TIC for {star_name} within 30 arcsec."

            table = result[0]
            table.sort('_r')
            best_match = table[0]

            self._status(f"Best match: TIC {best_match['TIC']}  |  SpT {spec_type}  |  "
                         f"Teff {best_match['Teff']:.0f} K  |  "
                         f"d {best_match['Dist']:.1f} pc")

            lum_s = best_match['Rad'] ** 2 * (best_match['Teff'] / 5780) ** 4

            self.sources['star'] = {'temperature': best_match['Teff'],
                                    'radius': best_match['Rad'],
                                    'distance': best_match['Dist'],
                                    'latitude': np.deg2rad(best_match['ELAT']),
                                    'longitude': np.deg2rad(best_match['ELON']),
                                    'l_sun': lum_s,
                                    'name': star_name}

    def planet_solar_system(self,
                            planet_name: str):
        eff_temps = {
            'mercury': 437, 'venus': 232, 'earth': 255, 'mars': 209,
            'jupiter': 88, 'saturn': 95, 'uranus': 58, 'neptune': 55.5,
        }
        geom_albedos = {
            'mercury': 0.142, 'venus': 0.689, 'earth': 0.24, 'mars': 0.17,
            'jupiter': 0.538, 'saturn': 0.499, 'uranus': 0.488, 'neptune': 0.442,
        }

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

        self._print(Panel(
            f"[bold white]Planet:[/bold white] [cyan]{planet_name}[/cyan]",
            border_style="green", expand=False
        ))

        key = str(planet_name).strip().lower()

        if key not in lookup:
            return f"Error: '{planet_name}' not found. Are you sure that's a planet?"

        self._status(f"Looking up ephemeris data for {key.capitalize()}...")
        p = lookup[key]

        try:
            self._status("Computing scaled semi-major axis from stellar luminosity...")
            sma = p.rAU * np.sqrt(self.sources['star']['l_sun'])
        except:
            sma = p.rAU
            warnings.warn('Stellar luminosity not set. Will use solar reference.')

        self._status(f"Teff {eff_temps[key]} K  |  R {p.R / planets.Earth.R:.3f} R⊕  |  "
                     f"a {p.rAU} AU  |  albedo {geom_albedos[key]}")

        self.sources['planet'] = {'radius': p.R / planets.Earth.R,
                                  'temperature': eff_temps[key],
                                  'separation': sma,
                                  'geom_albedo': geom_albedos[key],
                                  'name': key}

    def from_lband(self,
                   planet_name: str,
                   date: str):

        self._print(Panel(
            f"[bold white]L-band photometry:[/bold white] [cyan]{planet_name}[/cyan]",
            border_style="yellow", expand=False
        ))

        self._status("Loading photometry catalog...")
        catalog = pd.read_csv(str(files("lifesim.analysis") / "etc_data" / "reliable_photometry.csv"))
        planet = catalog.loc[catalog['planet_name'] == planet_name]
        if len(planet) == 0:
            raise ValueError('Could not find planet named {}'.format(planet_name))

        star_name = ' '.join(planet['planet_name'].values[0].split(' ')[:-1])
        self.add_star(star_name=star_name)

        if planet['witp_name'].values[0] == 'noorbit':
            self._status("No orbit data available — using last known angular separation.")
            angsep = planet['angsep_arcsec'].values[0]
        else:
            self._status("Retrieving angular separation via WhereIsThePlanet...")
            with _suppress_output():
                _, _, sep_args, _ = whereistheplanet.predict_planet(
                    planet['witp_name'].values[0], date
                )
            angsep = sep_args[0] * 1e-3

        self._status(f"Angular separation: {angsep}")

        self._status("Converting L-band magnitude to photon flux...")
        ph_flux = get_flux_filter(
            filter_id=planet['filter_id'].values[0],
            magnitude=planet['l_band_mag'].values[0]
        )
        self._status(f"L-band mag {planet['l_band_mag'].values[0]}  →  flux {ph_flux}")

        self.sources['planet'] = {'name': planet_name,
                                  'angsep': angsep,
                                  'separation': angsep * self.sources['star']['distance'],
                                  'ph_flux': ph_flux.value,
                                  }

    def list_lband(self):
        catalog = pd.read_csv(str(files("lifesim.analysis") / "etc_data" / "reliable_photometry.csv"))
        return catalog['planet_name'].tolist()

    def add_exozodi(self,
                    z: float):
        self.sources['exozodi'] = {'z': z,}
        self._status("Adding exozodi...")


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
