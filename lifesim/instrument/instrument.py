import time

import numpy as np
from tqdm import tqdm

from lifesim.modules.options import Options
from lifesim.archive.catalog import Catalog
from lifesim.dataio.bus import PrimaryModule
from lifesim.modules.util import black_body, import_spectrum
from lifesim.modules.habitable import single_habitable_zone


class Instrument(PrimaryModule):
    def __init__(self,
                 name: str,
                 options: Options):
        super().__init__(name=name)

        # initialize all instrument parameters
        self.options = options
        self.bl, self.telescope_area, self.eff_tot = None, None, None
        self.wl_bins, self.wl_bin_widths, self.wl_bin_edges = None, None, None
        self.hfov, self.hfov_mas, self.rpp, self.mas_pix = None, None, None, None
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
        self.options = options

        # Get array parameters from options for faster calculation
        self.bl = self.options.array['baseline']

        self.telescope_area = np.pi * (self.options.array['diameter'] / 2.) ** 2 * 4.
        self.eff_tot = self.options.array['quantum_eff'] * self.options.array['throughput']

        self.wl_bins, self.wl_bin_widths, self.wl_bin_edges = self.get_wl_bins_const_spec_res()

        # fov = wl / D -> hfov=wl/(2*D)
        self.hfov = self.wl_bins / (2. * self.options.array['diameter'])

        self.hfov_mas = self.hfov * (3600000. * 180.) / np.pi
        self.rpp = (2 * self.hfov) / self.options.other['image_size']  # Radians per pixel
        self.mas_pix = (2 * self.hfov_mas) / self.options.other['image_size']  # mas per pixel

        # apertures defines the telescope positions (and *relative* radius)
        self.apertures = np.array([[-self.bl / 2, -6 * self.bl / 2., 1.],
                                   [self.bl / 2, -6 * self.bl / 2., 1.],
                                   [self.bl / 2, 6 * self.bl / 2., 1.],
                                   [-self.bl / 2, 6 * self.bl / 2., 1.]])

        # self.planet_offset = [int(self.options.other['image_size'] / 2), int(self.options.other['image_size'] * 3 / 4)]
        # self.t_step = self.t_tot / self.rot_steps

        # coordinate maps for faster calculations
        self.x_map = np.tile(np.array(range(0, self.options.other['image_size'])),
                             (self.options.other['image_size'], 1))
        self.y_map = self.x_map.T
        self.r_square_map = ((self.x_map - (self.options.other['image_size'] - 1) / 2) ** 2
                             + (self.y_map - (self.options.other['image_size'] - 1) / 2) ** 2)
        self.r_map = np.sqrt(self.r_square_map)

        # push the new options and parameters to the sockets
        self.update_socket(name='transmission_generator',
                           data={'wl': self.wl_bins,
                                 'hfov_mas': self.hfov_mas,
                                 'image_size': self.options.other['image_size'],
                                 'bl': self.bl,
                                 'map_selection': 'tm3'})
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
                                     'mas_pix': self.mas_pix})

    def get_wl_bins_const_spec_res(self):
        wl_edge = self.options.array['wl_min']
        wl_bins = []
        wl_bin_widths = []
        wl_bin_edges = [wl_edge]

        while wl_edge < self.options.array['wl_max']:
            wl_bin_width = wl_edge / self.options.array['spec_res'] / \
                           (1 - 1 / self.options.array['spec_res'] / 2)

            if wl_edge + wl_bin_width > self.options.array['wl_max']:
                wl_bin_width = self.options.array['wl_max'] - wl_edge

            wl_center = wl_edge + wl_bin_width / 2
            wl_edge += wl_bin_width

            wl_bins.append(wl_center)
            wl_bin_widths.append(wl_bin_width)
            wl_bin_edges.append(wl_edge)

        wl_bins = np.array(wl_bins) * 1e-6  # in m
        wl_bin_widths = np.array(wl_bin_widths) * 1e-6  # in m
        wl_bin_edges = np.array(wl_bin_edges) * 1e-6  # in m

        return wl_bins, wl_bin_widths, wl_bin_edges

    def adjust_bl_to_hz(self,
                        hz_center: float,
                        distance_s: float):
        hz_center_rad = hz_center / distance_s / (3600 * 180) * np.pi  # in rad

        # put first transmission peak of optimal wl on center of HZ
        self.bl = 0.589645 / hz_center_rad * self.options.other['wl_optimal']*10**(-6)

        self.bl = np.maximum(self.bl, self.options.array['bl_min'])
        self.bl = np.minimum(self.bl, self.options.array['bl_max'])

        self.apertures = np.array([[-self.bl / 2, -6 * self.bl / 2., 1.],
                                   [self.bl / 2, -6 * self.bl / 2., 1.],
                                   [self.bl / 2, 6 * self.bl / 2., 1.],
                                   [-self.bl / 2, 6 * self.bl / 2., 1.]])

    # TODO Re-add functionality for calculating the SNR without certain noise term
    def get_snr(self,
                c: Catalog,
                time: int):

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
                # in most cases, more sockets are initialized than plugins are needed. Running
                # empty sockets currently throws exceptions
                # TODO change socket to catch running an empty socket. Throw and catch specific
                #   exception
                try:
                    self.run_socket(name='p_noise_source_'+str(j))
                except AttributeError:
                    pass

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

    def get_spectrum(self,
                     pathtofile: str,
                     temp_s: float,  # in K
                     radius_s: float,  # in R_sun
                     distance_s: float,  # in pc
                     lat_s: float,  # in radians
                     z: float,  # in zodis
                     angsep: float,  # in arcsec
                     radius_p: float,  # in R_earth
                     radius_spec: float,  # in R_earth radius of the planet in the reference spec
                     distance_spec: float,  # in pc distance at which the spectrum is simulated
                     integration_time: float
                     ):
        s_in, s_out, l_sun, \
            hz_in, hz_out, \
            hz_center = single_habitable_zone(model=self.options.models['habitable'],
                                              temp_s=temp_s,
                                              radius_s=radius_s)
        self.adjust_bl_to_hz(hz_center=hz_center,
                             distance_s=distance_s)

        # TODO: REMOVE REMOVE REMOVE
        self.bl = 14.54118465607078

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

        flux_planet_spectrum = import_spectrum(pathtofile=pathtofile,
                                               wl_bin_edges=self.wl_bin_edges,
                                               radius_p=radius_p,
                                               distance_s=distance_s,
                                               radius_spec=radius_spec,
                                               distance_spec=distance_spec)
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
        noise = noise_bg + noise_planet
        snr_spec = np.sqrt((flux_planet ** 2 / noise))

        return [self.wl_bins, snr_spec]



