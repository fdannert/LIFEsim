import numpy as np
from tqdm import tqdm

from lifesim.modules.options import Options
from lifesim.dataio.catalog import Catalog
from lifesim.dataio.bus import PrimaryModule
from lifesim.modules.util import black_body, import_spectrum
from lifesim.modules.habitable import single_habitable_zone


class Instrument(PrimaryModule):
    """
    The Instrument class represents the central module for simulating the LIFE array. It connects
    to 'plugin' modules which calculate signal and noise terms and distributes tasks and data
    between them. The instrument class features two socket types:
        a)  For calculation of the instrument transmission map a single socket of f_type
            'transmission'
        b)  For simulation of the photon noise sources a number (set in the options class) of
            sockets of f_type 'photon_noise'

    Attributes
    ----------
    options : lifesim.Options
        The options class containing all setting for the array and computations
    bl : float
        Length of the shorter, nulling baseline in [m]
    telescope_area : float
        Area of all array apertures combined in [m^2]
    eff_tot : float
        Total efficiency of the telescope as ratio of generated counts over incoming photons
        (dimensionless)
    wl_bins : np.ndarray
        Central values of the spectral bins in the wavelength regime in [m]
    wl_bin_widths : np.ndarray
        Widths of the spectral wavelength bins in [m]
    wl_bin_edges : np.ndarray
        Edges of the spectral wavelength bins in [m]. For N bins, this array will contain N+1 edges
    hfov : np.ndarray
        Contains the half field of view of the observatory in [rad] for each of the spectral bins
    hfov_mas : np.ndarray
        Contains the half field of view of the observatory in [milliarcseconds] for each of the
        spectral bins
    rad_pix : np.ndarray
        Contains the size of each pixel projected to the sky in [rad]
    mas_pix : np.ndarray
        Contains the size of each pixel projected to the sky in [milliarcseconds]
    apertures : np.ndarray
        Positions of the collector spacecraft relative to the beam combiner in [m]
    x_map : np.ndarray
        A map used for speeding up calculations. Contains the indices of the pixels on the detector
        in x-direction (dimensionless)
    y_map : np.ndarray
        A map used for speeding up calculations. Contains the indices of the pixels on the detector
        in x-direction (dimensionless)
    r_square_map : np.ndarray
        A map used for speeding up calculations. Contains the square of the distance of a pixel
        from the center of the detector in [pix]
    r_map : np.ndarray
        A map used for speeding up calculations. Contains the distance of a pixel
        from the center of the detector in [pix]
    """

    def __init__(self,
                 name: str,
                 options: Options):
        """
        Parameters
        ----------
        name : str
            Name of the instrument module
        options : lifesim.Options
            The options class containing all setting for the array and computations
        """

        super().__init__(name=name)

        # initialize all instrument parameters
        self.options = options
        self.bl, self.telescope_area, self.eff_tot = None, None, None
        self.wl_bins, self.wl_bin_widths, self.wl_bin_edges = None, None, None
        self.hfov, self.hfov_mas, self.rad_pix, self.mas_pix = None, None, None, None
        self.apertures = None
        self.x_map, self.y_map, self.r_square_map, self.r_map = None, None, None, None

        # create sockets
        self.add_socket(name='transmission_generator',
                        f_type='transmission',
                        data={})
        for i in range(self.options.other['n_plugins']):
            self.add_socket(name='p_noise_source_' + str(i),
                            f_type='photon_noise',
                            data={})

        # set options
        self.apply_options(options)

    def apply_options(self,
                      options: Options):
        """
        Applies the options given to the instrument module and recalculates all necessary values

        Parameters
        ----------
        options : lifesim.Options
            The options class containing all setting for the array and computations
        """
        self.options = options

        # Get array parameters from options for faster calculation
        self.bl = self.options.array['baseline']

        self.telescope_area = np.pi * (self.options.array['diameter'] / 2.) ** 2 * 4.
        self.eff_tot = self.options.array['quantum_eff'] * self.options.array['throughput']

        self.wl_bins, self.wl_bin_widths, self.wl_bin_edges = self.get_wl_bins_const_spec_res()

        # fov = wl / D -> hfov=wl/(2*D)
        # TODO remove the double usage of mas and rad, stick to only one
        self.hfov = self.wl_bins / (2. * self.options.array['diameter'])

        self.hfov_mas = self.hfov * (3600000. * 180.) / np.pi
        self.rad_pix = (2 * self.hfov) / self.options.other['image_size']  # Radians per pixel
        self.mas_pix = (2 * self.hfov_mas) / self.options.other['image_size']  # mas per pixel

        # apertures defines the telescope positions (and *relative* radius)
        self.apertures = np.array([[-self.bl / 2, -self.options.array['ratio'] * self.bl / 2., 1.],
                                   [self.bl / 2, -self.options.array['ratio'] * self.bl / 2., 1.],
                                   [self.bl / 2, self.options.array['ratio'] * self.bl / 2., 1.],
                                   [-self.bl / 2, self.options.array['ratio'] * self.bl / 2., 1.]])

        # coordinate maps for faster calculations
        self.x_map = np.tile(np.array(range(0, self.options.other['image_size'])),
                             (self.options.other['image_size'], 1))
        self.y_map = self.x_map.T
        self.r_square_map = ((self.x_map - (self.options.other['image_size'] - 1) / 2) ** 2
                             + (self.y_map - (self.options.other['image_size'] - 1) / 2) ** 2)
        self.r_map = np.sqrt(self.r_square_map)

        # push the new options and parameters to the sockets
        self.update_socket(name='transmission_generator',
                           data={'wl_bins': self.wl_bins,
                                 'hfov': self.hfov,
                                 'image_size': self.options.other['image_size'],
                                 'bl': self.bl,
                                 'map_selection': 'tm3',
                                 'ratio': self.options.array['ratio']})
        for i in range(self.options.other['n_plugins']):
            self.update_socket(name='p_noise_source_' + str(i),
                               data={'lz_model': self.options.models['localzodi'],
                                     'image_size': self.options.other['image_size'],
                                     'radius_map': self.r_map,
                                     'wl_bins': self.wl_bins,
                                     'wl_width': self.wl_bin_widths,
                                     'wl_bin_edges': self.wl_bin_edges,
                                     'hfov': self.hfov,
                                     'telescope_area': self.telescope_area,
                                     'mas_pix': self.mas_pix,
                                     'rad_pix': self.rad_pix,
                                     'ratio': self.options.array['ratio']})

    def get_wl_bins_const_spec_res(self):
        """
        Create the wavelength bins for the given spectral resolution and wavelength limits
        """
        wl_edge = self.options.array['wl_min']
        wl_bins = []
        wl_bin_widths = []
        wl_bin_edges = [wl_edge]

        while wl_edge < self.options.array['wl_max']:

            # set the wavelength bin width according to the spectral resolution
            wl_bin_width = wl_edge / self.options.array['spec_res'] / \
                           (1 - 1 / self.options.array['spec_res'] / 2)

            # make the last bin shorter when it hits the wavelength limit
            if wl_edge + wl_bin_width > self.options.array['wl_max']:
                wl_bin_width = self.options.array['wl_max'] - wl_edge

            # calculate the center and edges of the bins
            wl_center = wl_edge + wl_bin_width / 2
            wl_edge += wl_bin_width

            wl_bins.append(wl_center)
            wl_bin_widths.append(wl_bin_width)
            wl_bin_edges.append(wl_edge)

        # convert everything to [m]
        wl_bins = np.array(wl_bins) * 1e-6  # in m
        wl_bin_widths = np.array(wl_bin_widths) * 1e-6  # in m
        wl_bin_edges = np.array(wl_bin_edges) * 1e-6  # in m

        return wl_bins, wl_bin_widths, wl_bin_edges

    # TODO does not take the inclination into account!
    def adjust_bl_to_hz(self,
                        hz_center: float,
                        distance_s: float):
        """
        Adjusts the baseline of the array to be optimal for observations in the habitable zone of
        the target star for the selected optimal wavelength

        Parameters
        ----------
        hz_center : float
            Separation of the center of the habitable zone in [AU]
        distance_s : float
            Distance between the observed star and the LIFE array in [pc]
        """

        # convert the habitable zone to radians
        hz_center_rad = hz_center / distance_s / (3600 * 180) * np.pi  # in rad

        # put first transmission peak of optimal wl on center of HZ
        self.bl = 0.589645 / hz_center_rad * self.options.other['wl_optimal']*10**(-6)

        # make sure that the baseline does not exeed the set baseline limits
        self.bl = np.maximum(self.bl, self.options.array['bl_min'])
        self.bl = np.minimum(self.bl, self.options.array['bl_max'])

        # update the position of the apertures
        self.apertures = np.array([[-self.bl / 2, -6 * self.bl / 2., 1.],
                                   [self.bl / 2, -6 * self.bl / 2., 1.],
                                   [self.bl / 2, 6 * self.bl / 2., 1.],
                                   [-self.bl / 2, 6 * self.bl / 2., 1.]])

    # TODO Re-add functionality for calculating the SNR without certain noise term
    def get_snr(self,
                c: Catalog,
                time: int):
        """
        Calculates the signal-to-noise ration for all planets within the catalog if the are
        observed by the LIFE array for the given observing time

        Parameters
        ----------
        c : lifesim.Catalog
            The catalog for which the SNRs will be calculated and added to
        time : float
            The observation time per planet in [s]
        """

        snr_tot = np.zeros_like(c.data.nstar, dtype=float)
        c.data['snr_1h'] = np.zeros_like(c.data.nstar, dtype=float)

        # Push catalog to photon noise sockets
        for j in range(self.options.other['n_plugins']):
            self.update_socket(name='p_noise_source_'+str(j),
                               data={'c': c})

        # iterate over all stars
        for i, n in enumerate(tqdm(np.where(c.masks['stars'])[0])):
            nstar = c.data.nstar[n]
            index_s = np.argwhere((c.data.nstar.to_numpy() == nstar))[0]

            # adjust baseline of array and give new baseline to transmission generator plugin
            self.adjust_bl_to_hz(hz_center=float(c.data.hz_center[index_s]),
                                 distance_s=float(c.data.distance_s[index_s]))
            self.update_socket(name='transmission_generator',
                               data={'bl': self.bl})

            # get transmission map
            self.run_socket(name='transmission_generator',
                            mode='map')
            tm3 = self.sockets['transmission_generator'].tm3

            # update and run the photon noise plugins
            for j in range(self.options.other['n_plugins']):
                self.update_socket(name='p_noise_source_'+str(j),
                                   data={'bl': self.bl,
                                         't_map': tm3,
                                         'radius_s': float(c.data.radius_s[index_s]),
                                         'distance_s': float(c.data.distance_s[index_s]),
                                         'temp_s': float(c.data.temp_s[index_s]),
                                         'lat_s': float(c.data.lat[index_s]),
                                         'l_sun': float(c.data.l_sun[index_s]),
                                         'z': float(c.data.z[index_s])})

                # in most cases, more sockets are initialized than plugins are needed
                self.run_socket(name='p_noise_source_'+str(j))

            # get the new transmission maps
            for _, n_p in enumerate(np.argwhere((c.data['nstar'].to_numpy() == nstar))):
                self.update_socket(name='transmission_generator',
                                   data={'ang_sep_as': c.data['angsep'].to_numpy()[n_p]})
                self.run_socket(name='transmission_generator',
                                mode='efficiency')

                # calculate the photon flux originating from the planet
                flux_planet_thermal = black_body(mode='planet',
                                                 bins=self.wl_bins,
                                                 width=self.wl_bin_widths,
                                                 temp=c.data['temp_p'].to_numpy()[n_p],
                                                 radius=c.data['radius_p'].to_numpy()[n_p],
                                                 distance=c.data['distance_s'].to_numpy()[n_p])
                flux_planet = self.sockets['transmission_generator'].transm_eff * \
                              flux_planet_thermal * time * self.eff_tot
                noise_planet = self.sockets['transmission_generator'].transm_noise * \
                               flux_planet_thermal * time * self.eff_tot

                # calculate the noise from the background sources
                noise_bg = 0
                for p in range(self.options.other['n_plugins']):
                    if self.sockets['p_noise_source_'+str(p)] is not None:
                        noise_bg += self.sockets['p_noise_source_'+str(p)].noise
                noise_bg = noise_bg * time * self.eff_tot

                # Add up the noise and caluclate the SNR
                noise = noise_bg + noise_planet
                snr_tot[n_p] = np.sqrt((flux_planet ** 2 / noise).sum())
            c.data['snr_1h'] = snr_tot

    # TODO: fix units in documentation
    def get_spectrum(self,
                     temp_s: float,  # in K
                     radius_s: float,  # in R_sun
                     distance_s: float,  # in pc
                     lat_s: float,  # in radians
                     z: float,  # in zodis
                     angsep: float,  # in arcsec
                     flux_planet_spectrum: list,  # in ph m-3 s-1 over m
                     integration_time: float):
        """
        Simulates the spectrum of an exoplanet as seen by LIFE array, i.e. with the respective
        noise terms

        Parameters
        ----------
        pathtofile : str
            Path to the .txt file holding the input spectrum. The input spectrum need to be in the
            units of a photon flux density per micron, i.e. [s^-1 m^-2 microns^-1]
        temp_s : float
            Temperature of the observed star in [K]
        radius_s : float
            Radius of the observed star in [sun radii]
        distance_s : float
            Distance between the observed star and the LIFE array in [pc]
        lat_s : float
            Ecliptic latitude of the observed star in [rad]
        z : float
            Zodi level in the observed system in [zodis]
        angsep : float
            Angular separation between the observed star and the observed exoplanet in [arcsec]
        radius_p : float
            Radius of the observed exoplanet in [earth radii]
        radius_spec : float
            Radius of the observed exoplanet that was assumed in the creation of the input spectrum
            in [earth radii]
        distance_spec : float
            Distance between the observed star and the observer that was assumed in the creation of
            the input spectrum in [pc]
        integration_time : float
            Time that the LIFE array spends for observing the observed planet in [s]

        Returns
        -------
        Tuple[wl_bins, snr_spec]
            Some text
        flux_planet
            Some text
        noise
            Some text
        """

        # calculate the habitable zone of the specified star
        s_in, s_out, l_sun, \
            hz_in, hz_out, \
            hz_center = single_habitable_zone(model=self.options.models['habitable'],
                                              temp_s=temp_s,
                                              radius_s=radius_s)

        # adjust the baseline of the array to the habitable zone
        self.adjust_bl_to_hz(hz_center=hz_center,
                             distance_s=distance_s)

        # get transmission map
        self.update_socket(name='transmission_generator',
                           data={'bl': self.bl,
                                 'ang_sep_as': angsep})
        self.run_socket(name='transmission_generator',
                        mode='map')
        self.run_socket(name='transmission_generator',
                        mode='efficiency')
        tm3 = self.sockets['transmission_generator'].tm3

        # update and run the photon noise plugins
        for j in range(self.options.other['n_plugins']):
            self.update_socket(name='p_noise_source_' + str(j),
                               data={'bl': self.bl,
                                     't_map': tm3,
                                     'radius_s': radius_s,
                                     'distance_s': distance_s,
                                     'temp_s': temp_s,
                                     'lat_s': lat_s,
                                     'l_sun': l_sun,
                                     'z': z})

            self.run_socket(name='p_noise_source_' + str(j))

        bins = np.digitize(flux_planet_spectrum[0].value, self.wl_bin_edges)
        bins_mean = [flux_planet_spectrum[1].value[bins == i].mean()
                     for i in range(1, len(self.wl_bin_edges))]
        bins_mean = np.array(bins_mean)
        flux_planet_spectrum = bins_mean

        # calculate the signal and photon noise flux received from the planet
        flux_planet = flux_planet_spectrum \
                          * self.sockets['transmission_generator'].transm_eff \
                          * integration_time \
                          * self.eff_tot \
                          * self.telescope_area \
                          * self.wl_bin_widths
        noise_planet = flux_planet_spectrum \
                           * self.sockets['transmission_generator'].transm_noise \
                           * integration_time \
                           * self.eff_tot \
                           * self.telescope_area \
                           * self.wl_bin_widths

        # calculate the noise from the background sources
        noise_bg = 0
        for p in range(self.options.other['n_plugins']):
            if self.sockets['p_noise_source_' + str(p)] is not None:
                noise_bg += self.sockets['p_noise_source_' + str(p)].noise
        noise_bg = noise_bg * integration_time * self.eff_tot

        # Add up the noise and caluclate the SNR
        noise = (noise_bg + noise_planet) * 2
        snr_spec = np.sqrt((flux_planet ** 2 / noise))

        return ([self.wl_bins, snr_spec],
                flux_planet_spectrum,
                noise)



