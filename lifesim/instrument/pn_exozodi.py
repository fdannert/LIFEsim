import numpy as np
from typing import Union

from lifesim.core.modules import PhotonNoiseModule, PhotonNoiseRotationModule
from lifesim.util import constants
from lifesim.util.radiation import black_body


class PhotonNoiseExozodi(PhotonNoiseModule):
    """
    This class simulates the noise contribution of an exozodi disk to the interferometric
    measurement of LIFE.
    """

    def __init__(self,
                 name: str):
        super().__init__(name=name)
        """
        Parameters
        ----------
        name : str
            Name of the module.
        """

    def noise(self,
              index: Union[int, type(None)]):
        """
        Simulates the amount of photon noise originating from the exozodi of the observed system
        leaking into the LIFE array measurement.

        Parameters
        ----------
        index: Union[int, type(None)]
            Specifies the planet for which to calculate the noise contribution. If an integer n is
            given, the noise will be calculated for the n-th row in the `data.catalog`. If `None`
            is given, the noise is caluculated for the parameters located in `data.single`.

        Returns
        -------
        ez_leak
            Exozodi leakage in [photon s-1] per wavelength bin.

        Notes
        -----
        All of the following parameters are needed for the calculation of the exozodi noise
        contribution and should be specified either in `data.catalog` or in `data.single`.

        l_sun : float
            Luminosity of the observed star in [solar luminosities].
        distance_s : float
            Distance between the observed star and the LIFE array in [pc].
        z : float
            Zodi level in the observed system in [zodis].
        mas_pix : np.ndarray
            Contains the size of each pixel projected to the sky in [milliarcseconds].
        rad_pix : np.ndarray
            Contains the size of each pixel projected to the sky in [radians].
        data.inst['radius_map'] : np.ndarray
            Contains the distance of a pixel from the center of the detector in [pix].
        data.options.other['image_size']
            Number of pixels on one axis of a square detector (dimensionless). I.e. for a 512x512
            detector this value is 512.
        wl_bins : np.ndarray
            Central values of the spectral bins in the wavelength regime in [m].
        wl_widths : np.ndarray
            Widths of the spectral wavelength bins in [m].
        data.inst['telescope_area'] : float
            Area of all array apertures combined in [m^2].
        data.inst['t_map'] : np.ndarray
            Transmission map of the TM3 mode of the array created by the
            lifesim.TransmissionMap module.
        """

        # read from catalog or single data depending on index specification
        if index is None:
            l_sun = self.data.single['l_sun']
            distance_s = self.data.single['distance_s']
            z = self.data.single['z']
        else:
            l_sun = self.data.catalog.l_sun.iloc[index]
            distance_s = self.data.catalog.distance_s.iloc[index]
            z = self.data.catalog.z.iloc[index]

        # calculate the parameters required by Kennedy2015
        alpha = 0.34
        r_in = 0.034422617777777775 * np.sqrt(l_sun)
        r_0 = np.sqrt(l_sun)
        sigma_zero = 7.11889e-8  # Sigma_{m,0} from Kennedy+2015 (doi:10.1088/0067-0049/216/2/23)

        # reshape the mas per pixel array for calculation (to (n, 1, 1))
        mas_pix = np.array([self.data.inst['mas_pix']])
        if mas_pix.shape[-1] > 1:
            mas_pix = np.reshape(mas_pix, (mas_pix.shape[-1], 1, 1))
        rad_pix = np.array([self.data.inst['rad_pix']])
        if rad_pix.shape[-1] > 1:
            rad_pix = np.reshape(rad_pix, (rad_pix.shape[-1], 1, 1))

        au_pix = mas_pix / 1e3 * distance_s

        # the radius as measured from the central star for every pixel in [AU]
        r_au = self.data.inst['radius_map'] * au_pix

        # identify all pixels where the radius is larges than the inner radius by Kennedy+2015
        r_cond = ((r_au >= r_in)
                  & (r_au <= self.data.options.other['image_size'] / 2 * au_pix))

        # calculate the temperature at all pixel positions according to Kennedy2015 Eq. 2
        temp_map = np.where(r_cond,
                            278.3 * (l_sun ** 0.25) / np.sqrt(r_au), 0)

        # calculate the Sigma (Eq. 3) in Kennedy2015 and set everything inside the inner radius to 0
        sigma = np.where(r_cond,
                         sigma_zero * z *
                         (r_au / r_0) ** (-alpha), 0)

        wl_bins = np.array([self.data.inst['wl_bins']])
        if wl_bins.shape[-1] > 1:
            wl_bins = np.reshape(wl_bins, (wl_bins.shape[-1], 1, 1))

        wl_bin_widths = np.array([self.data.inst['wl_bin_widths']])
        if wl_bin_widths.shape[-1] > 1:
            wl_bin_widths = np.reshape(wl_bin_widths, (wl_bin_widths.shape[-1], 1, 1))

        # get the black body radiation emitted by the interexoplanetary dust
        f_nu_disk = black_body(bins=wl_bins,
                               width=wl_bin_widths,
                               temp=temp_map,
                               mode='wavelength') \
                    * sigma * rad_pix ** 2 * self.data.inst['telescope_area']

        ap = np.where(self.data.inst['radius_map']
                      <= self.data.options.other['image_size'] / 2, 1, 0)
        # add the transmission map
        ez_leak = (f_nu_disk * self.data.inst['t_map'] * ap).sum(axis=(-2, -1))

        return ez_leak


class PhotonNoiseExozodiInclined(PhotonNoiseRotationModule):
    """
    This class simulates the noise contribution of an exozodi disk to the interferometric
    measurement of LIFE.
    """

    def __init__(self,
                 name: str):
        super().__init__(name=name)
        """
        Parameters
        ----------
        name : str
            Name of the module.
        """

    def noise(self,
              index: Union[int, type(None)]):
        """
        Simulates the amount of photon noise originating from the exozodi of the observed system
        leaking into the LIFE array measurement.

        Parameters
        ----------
        index: Union[int, type(None)]
            Specifies the planet for which to calculate the noise contribution. If an integer n is
            given, the noise will be calculated for the n-th row in the `data.catalog`. If `None`
            is given, the noise is caluculated for the parameters located in `data.single`.

        Returns
        -------
        ez_leak
            Exozodi leakage in [photon s-1] per wavelength bin.

        Notes
        -----
        All of the following parameters are needed for the calculation of the exozodi noise
        contribution and should be specified either in `data.catalog` or in `data.single`.

        l_sun : float
            Luminosity of the observed star in [solar luminosities].
        distance_s : float
            Distance between the observed star and the LIFE array in [pc].
        z : float
            Zodi level in the observed system in [zodis].
        mas_pix : np.ndarray
            Contains the size of each pixel projected to the sky in [milliarcseconds].
        rad_pix : np.ndarray
            Contains the size of each pixel projected to the sky in [radians].
        data.inst['radius_map'] : np.ndarray
            Contains the distance of a pixel from the center of the detector in [pix].
        data.options.other['image_size']
            Number of pixels on one axis of a square detector (dimensionless). I.e. for a 512x512
            detector this value is 512.
        wl_bins : np.ndarray
            Central values of the spectral bins in the wavelength regime in [m].
        wl_widths : np.ndarray
            Widths of the spectral wavelength bins in [m].
        data.inst['telescope_area'] : float
            Area of all array apertures combined in [m^2].
        data.inst['t_map'] : np.ndarray
            Transmission map of the TM3 mode of the array created by the
            lifesim.TransmissionMap module.
        """

        if index is None:
            l_sun = self.data.single['l_sun']
            distance_s = self.data.single['distance_s']
            z = self.data.single['z']
            inclination_ez = self.data.single['inclination_ez']
            ascending_node_ez = self.data.single['ascending_node_ez']
        else:
            l_sun = self.data.catalog.l_sun.iloc[index]
            distance_s = self.data.catalog.distance_s.iloc[index]
            z = self.data.catalog.z.iloc[index]
            inclination_ez = self.data.catalog.inc_p.iloc[index]
            ascending_node_ez = self.data.catalog.large_omega_p[index]

        # calculate the parameters required by Kennedy2015
        alpha = 0.34
        r_in = 0.034422617777777775 * np.sqrt(l_sun)
        r_0 = np.sqrt(l_sun)
        sigma_zero = 7.11889e-8  # Sigma_{m,0} from Kennedy 2015 (doi:10.1088/0067-0049/216/2/23)

        # reshape the mas per pixel array for calculation (to (n, 1, 1))
        mas_pix = np.array([self.data.inst['mas_pix']])
        if mas_pix.shape[-1] > 1:
            mas_pix = np.reshape(mas_pix, (mas_pix.shape[-1], 1, 1))
        rad_pix = np.array([self.data.inst['rad_pix']])
        if rad_pix.shape[-1] > 1:
            rad_pix = np.reshape(rad_pix, (rad_pix.shape[-1], 1, 1))

        au_pix = mas_pix / 1e3 * distance_s
        r_au = self.data.inst['radius_map'] * au_pix

        # @Felix: Inputs von mir sind zwischen den beiden '###...' Zeilen
        ###################
        # obtain dimensions of r_au (array of depth x length^2 matrices for different wavelengths)
        length = len(r_au[1, :, 1])
        depth = len(r_au[:, 1, 1])

        # fill maps array with positions -> array entry [:,i,j] contains array [i,j] to indicate its position
        maps = np.moveaxis(np.mgrid[:length, :length], 0, -1)
        maps = np.moveaxis(np.repeat(maps[:, :, np.newaxis], depth, axis=2), 2, 0)

        # shift the maps such that the center of each map is at [0,0,0] and
        # each of the length^2 array entries contains its 3-d position on the center-origin unit coordinate system
        # we consider a plane with z-coordinate 0 s.t. each 2-d position vector gets an additional zero-coordinate

        shifted_maps = np.zeros((depth, length, length, 3))

        shifted_maps[:, :, :, 0] = maps[:, :, :, 0] - (length - 1) / 2
        shifted_maps[:, :, :, 1] = - maps[:, :, :, 1] + (length - 1) / 2

        radii = np.sqrt((shifted_maps ** 2).sum(axis=3)) * au_pix

        # exo_incl_deg = self.data.single['inclination']
        # exo_asc_deg_input = 0
        # exo_asc_deg = -(90 + exo_asc_deg_input)
        # Comment Felix: In the end, the inclination and ascending node should be read from the
        # catalog, i.e. from self.data.catalog.inc_p and self.data.catalog.large_omega_p. Exactly
        # like it is done in the very beginning of this function after the 'if index is None'
        # statement

        # incl = exo_incl_deg * np.pi / 180
        # asc = exo_asc_deg * np.pi / 180

        rotation_steps = self.data.inst['rotation_steps']
        ez_leak = np.zeros((self.data.inst['wl_bins'].shape[0], rotation_steps.shape[0]))

        for i in range(rotation_steps.shape[0]):
            print(i)
            incl = -(np.pi/2 + inclination_ez)
            asc = ascending_node_ez + rotation_steps[i]

            # define rotation matrix with euler angles \omega = 0, \Omega and i
            R_incl = np.array([[1, 0, 0],
                               [0, np.cos(incl), np.sin(incl)],
                               [0, -np.sin(incl), np.cos(incl)]])

            R_asc = np.array([[np.cos(asc), np.sin(asc), 0],
                              [-np.sin(asc), np.cos(asc), 0],
                              [0, 0, 1]])

            R = np.matmul(R_incl, R_asc)
            R_proj = np.array(((1, 0, 0),
                               (0, 1, 0),
                               (R[2, 0] / R[2, 2], R[2, 1] / R[2, 2], 0)))

            # transform origin-centered true-to-scale 'sky-plane' maps into skewed origin-centered 'disk-system' maps
            rotated_maps = np.matmul(R_proj, shifted_maps[:, :, :, :, np.newaxis])[:, :, :, :, 0]

            # calculate distance from origin to rotated pixels in 3d-space in order to be able to
            # calculate temperatures of pixels at their actual position in space and not projected position
            rotated_radii = np.sqrt((rotated_maps ** 2).sum(axis=3)) * au_pix

            # project skewed maps back onto 2-d sky-plane
            # projected_maps = rotated_maps[:, :, :, :2]

            a = np.array((1, 0, 0))
            b = np.array((0, 1, 0))
            area_frac = np.linalg.norm(np.cross(np.matmul(R_proj, a), np.matmul(R_proj, b)))

            ###################

            # identify all pixels where the radius is larges than the inner radius by kennedy2015
            r_cond = ((rotated_radii >= r_in)
                      & (radii <= self.data.options.other['image_size'] / 2 * au_pix))

            # calculate the temperature at all pixel positions according to Kennedy2015 Eq. 2
            temp_map = np.where(r_cond,
                                278.3 * (l_sun ** 0.25) / np.sqrt(rotated_radii), 0)

            # calculate the Sigma (Eq. 3) in Kennedy2015 and set everything inside the inner radius to
            # 0, include area fraction here, as this is the area density
            sigma = np.where(r_cond,
                             sigma_zero * z * area_frac *
                             (rotated_radii / r_0) ** (-alpha), 0)

            wl_bins = np.array([self.data.inst['wl_bins']])
            if wl_bins.shape[-1] > 1:
                wl_bins = np.reshape(wl_bins, (wl_bins.shape[-1], 1, 1))

            wl_bin_widths = np.array([self.data.inst['wl_bin_widths']])
            if wl_bin_widths.shape[-1] > 1:
                wl_bin_widths = np.reshape(wl_bin_widths, (wl_bin_widths.shape[-1], 1, 1))

            # get the black body radiation emitted by the interexoplanetary dust
            # does the 1/area make sense?
            f_nu_disk = (black_body(bins=wl_bins,
                                    width=wl_bin_widths,
                                    temp=temp_map,
                                    mode='wavelength')
                         * sigma * rad_pix ** 2 * self.data.inst['telescope_area'])

            # plt.imshow(f_nu_disk[0,0,:,:])
            ap = np.where(self.data.inst['radius_map']
                          <= self.data.options.other['image_size'] / 2, 1, 0)
            # add the transmission map
            ez_leak[:, i] = (f_nu_disk * self.data.inst['t_map'] * ap).sum(axis=(-2, -1))

        return ez_leak