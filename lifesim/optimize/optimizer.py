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
        self.data.optm['hit_limit'] = np.zeros(5)
        self.data.optm['sum_detected'] = np.zeros(5)
        self.data.optm['num_universe'] = self.data.catalog.nuniverse.max() + 1
        self.data.optm['hit_limit'] = ((self.data.optm['sum_detected']
                                       / (self.data.optm['num_universe']))
                                       >= self.data.options.optimization['limit'][1][:])
        self.data.optm['tot_time'] = 0

        self.data.catalog['detected'] = False
        self.data.catalog['snr_current'] = 0.
        self.data.catalog['int_time'] = 0.
        self.data.catalog['t_slew'] = -self.data.options.array['t_slew']
        self.run_socket(s_name='slope',
                        method='distribute_time')


