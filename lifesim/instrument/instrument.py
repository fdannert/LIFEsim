import time

import numpy as np

from lifesim.modules.options import Options
from lifesim.archive.catalog import Catalog
from lifesim.dataio.bus import PrimaryModule


class Instrument(PrimaryModule):
    def __init__(self,
                 name: str,
                 options: Options):
        super().__init__(name=name)
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

        self.add_socket(name='transmission_generator',
                        f_type='transmission',
                        data={'wl': self.wl_bins,
                              'hfov_mas': self.hfov_mas,
                              'image_size': self.options.other['image_size'],
                              'bl': self.bl,
                              'map_selection': 'tm3'})

        self.add_socket(name='p_noise_source_1',
                        f_type='photon_noise',
                        data={'lz_model': self.options.models['localzodi'],
                              'image_size': self.options.other['image_size'],
                              'radius_map': self.r_map,
                              'wl_bins': self.wl_bins,
                              'wl_width': self.wl_bin_widths,
                              'wl_bin_edges': self.wl_bin_edges,
                              'hfov': self.hfov,
                              'telescope_area': self.telescope_area,
                              'mas_pix': self.mas_pix})

        # List of data for photon noise plugin:
        #   nstar, catalog, bl, wl_bins, wl_width
        #   lz_model, lat, image_size, transmission map, radius_map, wl_bins, wl_width, hfov

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
                        nstar: int,
                        c: Catalog):

        HZcenter_rad = c.data.hz_center[c.data.nstar == nstar].to_numpy()[0] \
                       / c.data.distance_s[c.data.nstar == nstar].to_numpy()[0] \
                       / (3600 * 180) * np.pi  # in rad

        # put first transmission peak of optimal wl on center of HZ
        self.bl = 0.589645 / HZcenter_rad * self.options.other['wl_optimal']

        self.bl = np.maximum(self.bl, self.options.array['bl_min'])
        self.bl = np.minimum(self.bl, self.options.array['bl_max'])

        self.apertures = np.array([[-self.bl / 2, -6 * self.bl / 2., 1.],
                                   [self.bl / 2, -6 * self.bl / 2., 1.],
                                   [self.bl / 2, 6 * self.bl / 2., 1.],
                                   [-self.bl / 2, 6 * self.bl / 2., 1.]])

    # TODO Re-add functionality for calculating the SNR without certain noise term
    def get_snr(self,
                c: Catalog):

        c.data['snr_1h'] = np.zeros_like(c.data.nstar, dtype=float)
        self.update_socket(name='p_noise_source_1',
                           data={'c': c})

        for i, n in enumerate(np.where(c.masks['stars'])[0]):
            nstar = c.data.nstar[n]
            self.adjust_bl_to_hz(nstar=nstar,
                                 c=c)
            t = time.time()
            self.update_socket(name='transmission_generator',
                               data={'bl': self.bl})
            self.run_socket(name='transmission_generator')
            tm3 = self.sockets['transmission_generator'].tm3
            print(time.time()-t)
            self.update_socket(name='p_noise_source_1',
                               data={'nstar': nstar,
                                     'bl': self.bl,
                                     't_map': tm3})
            self.run_socket(name='p_noise_source_1')
            if i == 20:
                break

            a=1

    # def get_SNR1h(self, n_max=np.inf, multiprocessing=False):
    #
    #     self.planets.SNR1h = np.zeros(self.n_p)
    #     self.planets.SNR1h_sl = np.zeros(self.n_p)
    #     self.planets.SNR1h_lz = np.zeros(self.n_p)
    #     self.planets.SNR1h_ez = np.zeros(self.n_p)
    #     self.planets.is_observed = np.zeros(self.n_p)
    #     self.planets.t_obs = np.zeros(self.n_p)
    #
    #     inds_s, inds_p, t_obs_s = self.get_obs_times(tmax_d=np.inf, tslew_h=10,
    #                                                  dmax_Mstars=10, mode="const1h")
    #
    #     for i in tqdm(range(min(n_max, len(inds_s)))):
    #         Nstar = inds_s[i]
    #         ind_p = inds_p[i]
    #         self.options.t_tot = t_obs_s[i]
    #
    #         self.star = Star.init_from_planets(self.planets, ind_p)
    #         self.lz = Localzodi(model=self.options.lz_model)
    #         self.ez1z = Disk(star=self.star,
    #                          options=self.options,
    #                          maspp=self.inst.maspp,
    #                          maps=True,
    #                          norm1z=True)
    #
    #         self.inst.adjust_bl_to_HZ(self.star, self.options)
    #         self.inst.add_transmission_maps(self.options, map_selection="tm3_only")
    #
    #         # get photon noise of astro sources. Factor 2 is due to symmetry of TMs 3&4
    #         N_sl = 2 * self.inst.get_stellar_leakage(self.star, self.options)
    #         N_lz = 2 * self.inst.get_localzodi_leakage(self.lz, self.star, self.options)
    #         N_ez1z = 2 * self.inst.get_exozodi_leakage(self.ez1z, self.options)
    #
    #         planet_inds_s = np.argwhere((self.planets.Nstar == Nstar))
    #
    #         mult_factor = (self.options.t_tot * self.inst.telescope_area *
    #                        self.inst.wl_bin_widths * self.inst.eff_tot)
    #
    #         for i_p in planet_inds_s:
    #             self.planets.is_observed[i_p] = True
    #             self.planets.t_obs[i_p] = t_obs_s[i]
    #
    #             planet = Planet.init_from_planets(self.planets, i_p)
    #
    #             trans_eff_p = tm_gen.transm_eff(self.inst.bl,
    #                                             self.inst.wl_bins,
    #                                             planet.ang_sep)
    #
    #             Fp = planet.fgamma(wl=self.inst.wl_bins)
    #             S_p = trans_eff_p * Fp * mult_factor
    #             N_bg = (N_sl + N_lz + planet.zodis * N_ez1z) * mult_factor
    #             N_p = S_p
    #             N_comb = N_bg + N_p
    #             snr_tot = np.sqrt((S_p ** 2 / N_comb).sum())
    #
    #             self.planets.SNR1h[i_p] = snr_tot
    #             self.planets.SNR1h_sl[i_p] = np.sqrt((S_p ** 2 / (N_sl * mult_factor)).sum())
    #             self.planets.SNR1h_lz[i_p] = np.sqrt((S_p ** 2 / (N_lz * mult_factor)).sum())
    #             self.planets.SNR1h_ez[i_p] = np.sqrt(
    #                 (S_p ** 2 / (N_ez1z * planet.zodis * mult_factor)).sum())

