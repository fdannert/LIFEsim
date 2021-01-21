import numpy as np
from tqdm import tqdm

from lifesim.core.modules import OptimizationModule

class Optimizer(OptimizationModule):

    def __init__(self,
                 name: str):
        super().__init__(name=name)

    def find_phase(self):
        self.run_socket(s_name='instrument',
                        method='apply_options')
        theta_p = np.linspace(start=0,
                              stop=2 * np.pi,
                              num=self.data.options.optimization['N_pf'])
        self.data.catalog['snr_phase'] = None

        for n_p in tqdm(range(self.data.catalog.shape[0])):
            snr_phase = np.zeros_like(theta_p)
            if n_p == 10:
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

