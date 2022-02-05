import numpy as np
from typing import Union

from lifesim.instrument.transmission import TransmissionMap


class OrbitalTransmissionMap(TransmissionMap):
    """
    Module for calculating nulling interferometer transmission maps for orbiting exoplanets.
    """

    def projected_vector(self, theta, inc, rad):
        """
        Calculates the 3d rotated star to planet vector from inclincation and true anomaly.
        It returns the rotated vector projected on to the image plane.

        inc 0 => Face on counter-clockwise, co-rotating to the shift of the observatory
        inc pi/2 => Edge ON
        inc pi => Face on clockwise, counter-rotating to the shift of the observatory
        Returns
        X and Z value of the exoplanet after rotation scaled with maxangsep in radiens
        """
        start = np.array([np.cos(theta), 0*theta, np.sin(theta)])
        Ri = np.array([[1, 0, 0],
                       [0, np.cos(inc), -np.sin(inc)],
                       [0, np.sin(inc), np.cos(inc)]])
        vs = rad*Ri.dot(start)
        return vs[0], vs[2]

    def get_transmission_curve(self, index, time_dependent=True):
        """
        Calculates the radial transmission curve of the LIFE array for a orbiting exoplanet.

        Parameters
        ----------
        bl : float
            Length of the shorter, nulling baseline in [m]
        wl_bins : np.ndarray
            Central values of the spectral bins in the wavelength regime in [m]
        angsep : float
            Angular separation between the observed star and the observed exoplanet in [arcsec]
        phi_n : int
            Number of rotation steps used in integration
        ratio : float
            Ratio between the nulling and the imaging baseline. E.g. if the imaging baseline is twice
            as long as the nulling baseline, the ratio will be 2

        Returns
        -------
        transm_curve_chop
            Radial transmission curve corresponding to the chopped transmission map
        transm_curve_tm4
            Radial transmission curve corresponding to the transmission map of the 4th mode 'tm4'
        """

        angsep = self.data.catalog["maxangsep"].iloc[index]
        theta = self.data.catalog["theta_p"].iloc[index]
        inclination = self.data.catalog["inc_p"].iloc[index]
        n_rotations = self.data.inst["rotations"]
        rotation_steps = self.data.inst["rotation_steps"]
        angsep_rad = angsep / (3600 * 180) * np.pi
        if time_dependent:
            true_anom = self.data.catalog.iloc[index]["theta_p"]

            # period of planet orbit in days
            p_orb = self.data.catalog.iloc[index]["p_orb"]

            # Instrument rotation period 1h/ 5h / 20h / (12h) in seconds
            rotation_period = self.data.inst["rotation_period"]

            # calculates the time length of one integration step in seconds
            rotation_step_time = rotation_period / rotation_steps

            percentage_orbit_step = rotation_step_time / (p_orb * 24 * 60 * 60)

            # delta_theta: change per time_step in true_anomaly
            delta_theta = percentage_orbit_step * 2 * np.pi

            # movement of planet inside the observation period
            thetas = np.linspace(true_anom, true_anom + delta_theta * (rotation_steps * n_rotations),
                                 n_rotations * rotation_steps, endpoint=False, dtype="float")

            # calculating positions of planet and using numpy broadcasting get transmission map values
            xs, zs = self.projected_vector(thetas, inclination, angsep_rad)
            return self.transmission_curve_t(xs, zs, n_rotations, rotation_steps)

        else:
            x, z = self.projected_vector(theta, inclination, angsep_rad)
            return self.transmission_curve_t(x, z, n_rotations, rotation_steps)

    def transmission_efficiency(self,
                                index: Union[int, None],
                                time_dependent=True):
        """
        Integrates over transmission curves to get the transmission efficiency for signal and noise

        Parameters
        ----------
        bl : float
            Length of the shorter, nulling baseline in [m]
        wl_bins : np.ndarray
            Central values of the spectral bins in the wavelength regime in [m]
        angsep : float
            Angular separation between the observed star and the observed exoplanet in [arcsec]
        ratio : float
            Ratio between the nulling and the imaging baseline. E.g. if the imaging baseline is
            twice as long as the nulling baseline, the ratio will be 2

        Returns
        -------
        transm_eff
            Transmission efficiency per spectral bin for the exoplanet signal
        transm_noise
            Transmission efficiency per spectral bin for the photon noise received from the exoplanet
            signal

        """

        tc_chop, tc_tm4 = self.get_transmission_curve(index, time_dependent)

        # integrate over angles to get transmission efficiency
        transm_eff = np.sqrt((tc_chop ** 2).mean(axis=(-2, -1)))
        transm_noise = np.sqrt((tc_tm4 ** 2).mean(axis=(-2, -1)))
        return transm_eff, transm_noise
