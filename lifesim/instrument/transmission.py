import numpy as np
from typing import Union

from lifesim.core.modules import TransmissionModule


class TransmissionMap(TransmissionModule):
    def __init__(self,
                 name: str):
        super().__init__(name)

    def transmission_map(self,
                         map_selection: list,
                         direct_mode: bool = False,
                         d_alpha: np.ndarray = None,
                         d_beta: np.ndarray = None,
                         hfov: np.ndarray = None,
                         image_size: int = None):
        """
        Return the transmission map the LIFE array

        Parameters
        ----------
        image_size : int
            Number of pixels on one axis of a square detector (dimensionless). I.e. for a 512x512
            detector this value is 512
        bl : float
            Length of the shorter, nulling baseline in [m]
        wl_bins : np.ndarray
            Central values of the spectral bins in the wavelength regime in [m]
        hfov : np.ndarray
            Contains the half field of view of the observatory in [rad] for each of the spectral
            bins
        map_selection : list
            A list of strings specifying the array transmission mode from which the transmission
            map is generated. The possible options are 'tm1', 'tm2', 'tm3', 'tm4'
        ratio : float
            Ratio between the nulling and the imaging baseline. E.g. if the imaging baseline is
            twice as long as the nulling baseline, the ratio will be 2
        direct_mode: bool
            If direct mode is set to True, the variables wl_bins, d_alpha and d_beta are injected
            into the transmission calculation directly. The output is not necessarily in the the
            form of a map. The inputs hfov and image size will be ignored in this mode
        d_alpha: np.ndarray
            The x-positions of the points to be evaluated measured from the central viewing axis in
            [rad]
        d_beta: np.ndarray
            The y-positions of the points to be evaluated measured from the central viewing axis in
            [rad]

        Returns
        -------
        tm1
            Transmission map from first mode
        tm2
            Transmission map from second mode
        tm3
            Transmission map from third mode
        tm4
            Transmission map from fourth mode
        tm_chop
            The chopped transmission map calculated by subtracting tm4 from tm3
        """

        if hfov is None:
            hfov = self.data.inst['hfov']
        if image_size is None:
            image_size = self.data.options.other['image_size']

        # reshape the wl_bins and hfov arrays for calculation (to (n, 1, 1))
        wl_bins = np.array([self.data.inst['wl_bins']])  # wavelength in m
        if wl_bins.shape[-1] > 1:
            wl_bins = np.reshape(wl_bins, (wl_bins.shape[-1], 1, 1))

        if direct_mode:
            alpha = d_alpha
            beta = d_beta
        else:
            hfov = np.array([hfov])  # wavelength in m
            if hfov.shape[-1] > 1:
                hfov = np.reshape(hfov, (hfov.shape[-1], 1, 1))

            # generare 1D array that spans field of view
            angle = np.linspace(-1, 1, image_size)

            # angle matrix in x-direction ("alpha")
            alpha = np.tile(angle, (image_size, 1))

            # angle matrix in y-direction ("beta")
            beta = alpha.T

            # convert angle matrices to fov units
            alpha = alpha * hfov
            beta = beta * hfov

        # smaller distance of apertures from center line
        L = self.data.inst['bl'] / 2

        tm1, tm2, tm3, tm4, tm_chop = None, None, None, None, None

        # transmission map of mode 1
        if 'tm1' in map_selection:
            tm1 = np.cos(2 * np.pi * L * alpha / wl_bins) ** 2 * np.cos(
                2 * self.data.options.array['ratio'] * np.pi * L * beta / wl_bins - np.pi / 4) ** 2

        # transmission map of mode 2
        if 'tm2' in map_selection:
            tm2 = np.cos(2 * np.pi * L * alpha / wl_bins) ** 2 * np.cos(
                2 * self.data.options.array['ratio'] * np.pi * L * beta / wl_bins + np.pi / 4) ** 2

        # transmission map of mode 3
        if 'tm3' in map_selection:
            tm3 = np.sin(2 * np.pi * L * alpha / wl_bins) ** 2 * np.cos(
                2 * self.data.options.array['ratio'] * np.pi * L * beta / wl_bins - np.pi / 4) ** 2

        # transmission map of mode 4
        if 'tm4' in map_selection:
            tm4 = np.sin(2 * np.pi * L * alpha / wl_bins) ** 2 * np.cos(
                2 * self.data.options.array['ratio'] * np.pi * L * beta / wl_bins + np.pi / 4) ** 2

        # difference of transmission maps 3 and 4 = "chopped transmission"
        if 'tm_chop' in map_selection:

            # if tm3 and tm4 exist, calculate the chopped transmission directly from the difference
            if (tm3 is not None) and (tm4 is not None):
                tm_chop = tm3 - tm4

            # if they don't exist, calculate the chopped transmission from formula
            else:
                tm_chop = np.sin(2 * np.pi * L * alpha / wl_bins) ** 2 * np.sin(
                    # chopped transm. map
                    4 * self.data.options.array['ratio'] * np.pi * L * beta / wl_bins)

        return tm1, tm2, tm3, tm4, tm_chop

    def transmission_efficiency(self,
                                index: Union[int, type(None)]):
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

        if index is None:
            angsep = self.data.single['angsep']
        else:
            angsep = self.data.catalog.angsep.iloc[index]

        tc_chop, tc_tm4 = self.transmission_curve(angsep=angsep)

        # integrate over angles to get transmission efficiency
        transm_eff = np.sqrt((tc_chop ** 2).mean(axis=(-2, -1)))
        transm_noise = np.sqrt((tc_tm4 ** 2).mean(axis=(-2, -1)))
        return transm_eff, transm_noise

    def transmission_curve(self,
                           angsep: float,
                           phi_n: int = 360):
        """
        Calculates the radial transmission curve of the LIFE array

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

        # convert angular separation to radians
        angsep_rad = angsep / (3600 * 180) * np.pi

        # create 1D array with azimuthal coordinates
        phi_lin = np.linspace(0, 2 * np.pi, phi_n, endpoint=False)

        # retrieve the transmission curves
        (_, _, _,
         transm_curve_tm4,
         transm_curve_chop) = self.transmission_map(map_selection=['tm4', 'tm_chop'],
                                                    direct_mode=True,
                                                    d_alpha=angsep_rad *
                                                    np.cos(phi_lin),
                                                    d_beta=angsep_rad * np.sin(phi_lin))

        return transm_curve_chop, transm_curve_tm4

    def projected_vector(self, theta, inc, rad):
        """
        Calculating the star to planet vector and projecting it to the 2d image.
        The planet orbits in XZ axis and the inclination rotates this orbit around the X-axis
        inc 0 > Face on Counterclockwise
        inc pi/2 > Edge ON
        inc pi > Face on clockwise
        Returns
        X and Z value after rotation scaled with maxangsep in radiens
        """
        start = np.array([np.cos(theta), 0*theta, np.sin(theta)])
        Ri = np.array([[1, 0, 0],
                       [0, np.cos(inc), -np.sin(inc)],
                       [0, np.sin(inc), np.cos(inc)]])
        vs = rad*Ri.dot(start)
        return vs[0], vs[2]

    def get_transmission_curve(self, index, time_dependent=True):
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

            # calculate positions of planet and using numpy broadcasting position in transmission map
            xs, zs = self.projected_vector(thetas, inclination, angsep_rad)
            return self.transmission_curve_t(xs, zs, n_rotations, rotation_steps)

        else:
            # same as above just 1 value instead of numpy array
            x, z = self.projected_vector(theta, inclination, angsep_rad)
            return self.transmission_curve_t(x, z, n_rotations, rotation_steps)

    def transmission_efficiency_t(self,
                                  index: Union[int, type(None)],
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
