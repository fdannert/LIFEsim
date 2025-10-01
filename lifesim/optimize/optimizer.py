import numpy as np
from tqdm import tqdm

from lifesim.core.modules import OptimizationModule

class Optimizer(OptimizationModule):

    def __init__(self,
                 name: str):
        super().__init__(name=name)

    def find_phase(self,
                   recalc: bool = False):
        self.run_socket(s_name='instrument',
                        method='apply_options')
        self.data.catalog['snr_new'] = 0.
        theta_p = np.linspace(start=0,
                              stop=2 * np.pi,
                              num=self.data.options.optimization['N_pf'])

        if recalc:
            self.data.catalog['snr_phase'] = None

            for n_p in tqdm(range(self.data.catalog.shape[0])):
                snr_phase = np.zeros_like(theta_p)

                self.run_socket(s_name='instrument',
                                method='adjust_bl_to_hz',
                                hz_center=float(self.data.catalog.hz_center.iloc[n_p]),
                                distance_s=float(self.data.catalog.distance_s.iloc[n_p]))

                if n_p == 3200:
                    break
                for i, theta in enumerate(theta_p):
                    self.data.catalog.angsep.iat[n_p] = (self.data.catalog.semimajor_p.iloc[n_p]
                                                         / self.data.catalog.distance_s.iloc[n_p]
                                                         * np.sqrt(
                        np.cos(self.data.catalog.small_omega_p.iloc[n_p] + theta) ** 2
                        + np.cos(self.data.catalog.inc_p.iloc[n_p]) ** 2
                        * np.sin(self.data.catalog.small_omega_p.iloc[n_p] + theta) ** 2))

                    transm_eff, transm_noise = self.run_socket(s_name='transmission',
                                                               method='transmission_efficiency',
                                                               index=n_p)
                    noise = (self.data.catalog.noise_astro.iloc[n_p][0]
                             + self.data.catalog.planet_flux_use.iloc[n_p][0] * transm_noise * 2)
                    flux_planet = self.data.catalog.planet_flux_use.iloc[n_p][0] * transm_eff
                    snr_phase[i] = np.sqrt((flux_planet ** 2 / noise).sum())

                self.data.catalog.snr_phase.iat[n_p] = [snr_phase]

        for n_p in tqdm(range(self.data.catalog.shape[0])):
            i = np.argmax(self.data.catalog.snr_phase.iloc[n_p][0])
            self.data.catalog.theta_p.iat[n_p] = theta_p[i]
            self.data.catalog.angsep.iat[n_p] = (self.data.catalog.semimajor_p.iloc[n_p]
                                                 / self.data.catalog.distance_s.iloc[n_p]
                                                 * np.sqrt(
                        np.cos(self.data.catalog.small_omega_p.iloc[n_p]
                               + self.data.catalog.theta_p.iloc[n_p]) ** 2
                        + np.cos(self.data.catalog.inc_p.iloc[n_p]) ** 2
                        * np.sin(self.data.catalog.small_omega_p.iloc[n_p]
                                 + self.data.catalog.theta_p.iloc[n_p]) ** 2))
            self.data.catalog.snr_new.iat[n_p] = self.data.catalog.snr_phase.iloc[n_p][0][i]

    def ahgs(self):
        # sum of detected planets per stype
        self.data.optm['sum_detected'] = np.zeros(5)

        self.data.optm['num_universe'] = np.unique(self.data.catalog.nuniverse).shape[0]

        # initialize optimization limits

        # handle the legacy mode
        if self.data.options.optimization['limit_mode'] == 'legacy':
            self.data.options.optimization['experiments'] = {}
            if self.data.options.optimization['habitable']:
                rmin = 0.5
                rmax = 1.5
            else:
                rmin = 0.
                rmax = np.inf
            for stype, limit in self.data.options.optimization['limit'].items():
                self.data.options.optimization['experiments'][stype] = {
                    'radius_p_min': rmin,
                    'radius_p_max': rmax,
                    'temp_s_min': self.data.catalog[self.data.catalog.stype==stype].temp_s.min(),
                    'temp_s_max': self.data.catalog[self.data.catalog.stype==stype].temp_s.max(),
                    'in_HZ': self.data.options.optimization['habitable'],
                    'sample_size': self.data.options.optimization['limit'][stype],
                }

        if self.data.options.optimization['experiments'] is None:
            self.data.catalog['is_interesting'] = True
            self.data.optm['hit_limit'] = None
            self.data.optm['exp_detected'] = None
        else:
            # assign targets to experiments
            self.data.catalog['is_interesting'] = False
            self.data.optm['hit_limit'] = {}
            self.data.optm['exp_detected'] = {}
            for exp in self.data.options.optimization['experiments'].keys():
                mask_exp = ((self.data.catalog.radius_p
                             >= self.data.options.optimization['experiments'][exp]['radius_p_min'])
                            & (self.data.catalog.radius_p
                               <= self.data.options.optimization['experiments'][exp]['radius_p_max'])
                            & (self.data.catalog.temp_s
                               >= self.data.options.optimization['experiments'][exp]['temp_s_min'])
                            & (self.data.catalog.temp_s
                               <= self.data.options.optimization['experiments'][exp]['temp_s_max']))

                if self.data.options.optimization['experiments'][exp]['in_HZ']:
                    mask_exp = (mask_exp
                                & (self.data.catalog['habitable']))

                self.data.catalog['exp_' + exp] = mask_exp

                self.data.catalog['is_interesting'] = np.logical_or(mask_exp, self.data.catalog['is_interesting'])

                self.data.optm['hit_limit'][exp] = False
                self.data.optm['exp_detected'][exp] = 0

        self.data.optm['tot_time'] = 0  # in sec

        # add new columns to catalog
        self.data.catalog['detected'] = False
        self.data.catalog['snr_current'] = 0.
        self.data.catalog['int_time'] = 0.
        self.data.catalog['t_slew'] = -self.data.options.array['t_slew']
        self.data.catalog['t_detected'] = 0.

        self.run_socket(s_name='slope',
                        method='distribute_time')
        print('')
