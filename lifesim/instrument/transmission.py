import numpy as np
from typing import Union

from lifesim.core.modules import TransmissionModule


class TransmissionMap(TransmissionModule):
    """
    Module for calculating nulling interferometer transmission maps.
    """
    def __init__(self,
                 name: str):
        """
        Parameters
        ----------
        name : str
            Name of the module.
        """

        super().__init__(name)

    def transmission_map(self,
                         map_selection: list,
                         direct_mode: bool = False,
                         d_alpha: np.ndarray = None,
                         d_beta: np.ndarray = None,
                         hfov: np.ndarray = None,
                         image_size: int = None):
        """
        Return the transmission map of a double-Bracewell configuration for the LIFE array.

        Parameters
        ----------
        map_selection : list
            A list of strings specifying the array transmission mode from which the transmission
            map is generated. The possible options are 'tm1', 'tm2', 'tm3', 'tm4'.
        direct_mode: bool
            If direct mode is set to True, the variables wl_bins, d_alpha and d_beta are injected
            into the transmission calculation directly. The output is not necessarily in the the
            form of a map. The inputs hfov and image size will be ignored in this mode.
        d_alpha: np.ndarray
            The x-positions of the points to be evaluated measured from the central viewing axis in
            [rad].
        d_beta: np.ndarray
            The y-positions of the points to be evaluated measured from the central viewing axis in
            [rad].
        hfov : np.ndarray
            Contains the half field of view of the observatory in [rad] for each of the spectral
            bins. If no value is given, `data.inst['hfov']` is used.
        image_size : int
            Number of pixels on one axis of a square detector (dimensionless). I.e. for a 512x512
            detector this value is 512. If no value is given, `data.options.other['image_size']` is
            used.

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

        Notes
        -----
        The following additional parameters are required for the calculation of the transmission
        maps and should be specified either in `data.catalog` or in `data.single`.

        data.inst['wl_bins'] : np.ndarray
            Central values of the spectral bins in the wavelength regime in [m].
        data.inst['bl'] : float
            Length of the shorter, nulling baseline in [m].
        data.options.array['ratio'] : float
            Ratio between the nulling and the imaging baseline. E.g. if the imaging baseline is
            twice as long as the nulling baseline, the ratio will be 2.
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
                # chopped transm. map
                tm_chop = np.sin(2 * np.pi * L * alpha / wl_bins) ** 2 * np.sin(
                    4 * self.data.options.array['ratio'] * np.pi * L * beta / wl_bins)

        return tm1, tm2, tm3, tm4, tm_chop

    def transmission_efficiency(self,
                                index: Union[int, type(None)]):
        """
        Integrates over transmission curves to get the transmission efficiency for signal and
        noise.

        Parameters
        ----------
        index: Union[int, type(None)]
            Specifies the planet for which to calculate the transmission efficiency. If an integer
            n is given, the noise will be calculated for the n-th row in the `data.catalog`. If
            `None` is given, the noise is caluculated for the parameters located in `data.single`.

        Returns
        -------
        transm_eff
            Transmission efficiency per spectral bin for the exoplanet signal
        transm_noise
            Transmission efficiency per spectral bin for the photon noise received from the
            exoplanet signal

        Notes
        -----
        The following additional parameters are required for the calculation of the transmission
        efficiency and should be specified either in `data.catalog` or in `data.single`.

        data.single['angsep']
            Angular separation between the observed star and the observed exoplanet in [arcsec].
        data.inst['wl_bins'] : np.ndarray
            Central values of the spectral bins in the wavelength regime in [m].
        data.inst['bl'] : float
            Length of the shorter, nulling baseline in [m].
        data.options.array['ratio'] : float
            Ratio between the nulling and the imaging baseline. E.g. if the imaging baseline is
            twice as long as the nulling baseline, the ratio will be 2.
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
        angsep : float
            Angular separation between the observed star and the observed exoplanet in [arcsec].
        phi_n : int
            Number of rotation steps used in integration.

        Returns
        -------
        transm_curve_chop
            Radial transmission curve corresponding to the chopped transmission map.
        transm_curve_tm4
            Radial transmission curve corresponding to the transmission map of the 4th mode 'tm4'.
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
                                                    d_alpha=angsep_rad * np.cos(phi_lin),
                                                    d_beta=angsep_rad * np.sin(phi_lin))

        return transm_curve_chop, transm_curve_tm4
