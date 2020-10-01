import warnings

import numpy as np

# todo Documentation
class Catalog(object):
    """Gets and prints the spreadsheet's header columns

    Parameters
    ----------

    Returns
    -------
    """

    def __init__(self):
        """Gets and prints the spreadsheet's header columns

        Parameters
        ----------

        Returns
        -------
        """

        self.nuniverse = None
        self.radius_p = None  # Rearth
        self.p_orb = None  # d
        self.mass_p = None  # Mearth
        self.ecc_p = None
        self.inc_p = None  # rad
        self.large_omega_p = None  # rad
        self.small_omega_p = None  # rad
        self.theta_p = None  # rad
        self.albedo_bond = None
        self.albedo_geom_vis = None
        self.albedo_geom_mir = None
        self.z = None
        self.semimajor_p = None  # au
        self.sep_p = None  # au
        self.angsep = None  # arcsec
        self.maxangsep = None  # arcsec
        self.flux_p = None  # Searth
        self.fp = None
        self.temp_p = None  # K
        self.nstar = None
        self.radius_s = None  # Rsun
        self.mass_s = None  # Msun
        self.temp_s = None  # K
        self.distance_s = None  # pc
        self.stype = None
        self.ra = None  # deg
        self.dec = None  # deg

    def remove_distance(self,
                        stype: str,
                        dist: float,
                        mode: str):
        """Removes planets around stellar type above or below a certain distance from the sample

        Parameters
        ----------
        stype : str
            Planets around stars of the specified stellar type are removed. Possible options are
            'A', 'F', 'G', 'K', 'M'
        dist : float
            Specifies the distance over or under which the planets are removed in pc
        mode : str
            'larger':   planets around stars with distances larger than the specified distance are
                        removed
            'smaller':  planets around stars with distances smaller than the specified distance are
                        removed

        Returns
        -------
        """

        if not np.isin(stype, np.array(('A', 'F', 'G', 'K', 'M'))):
            raise ValueError('Stellar type not recognised')
        if mode == 'larger':
            mask = np.logical_and(self.stype == stype,
                                  self.distance_s >= dist)
        elif mode == 'smaller':
            mask = np.logical_and(self.stype == stype,
                                  self.distance_s <= dist)
        else:
            warnings.warn('Mode ' + mode + ' not available. Using mode larger.')
            mask = np.logical_and(self.stype == stype,
                                  self.distance_s >= dist)

        mask = np.invert(mask)

        self.nuniverse = self.nuniverse[mask]
        self.radius_p = self.radius_p[mask]
        self.p_orb = self.p_orb[mask]
        self.mass_p = self.mass_p[mask]
        self.ecc_p = self.ecc_p[mask]
        self.inc_p = self.inc_p[mask]
        self.large_omega_p = self.large_omega_p[mask]
        self.small_omega_p = self.small_omega_p[mask]
        self.theta_p = self.theta_p[mask]
        self.albedo_bond = self.albedo_bond[mask]
        self.albedo_geom_vis = self.albedo_geom_vis[mask]
        self.albedo_geom_mir = self.albedo_geom_mir[mask]
        self.z = self.z[mask]
        self.semimajor_p = self.semimajor_p[mask]
        self.sep_p = self.sep_p[mask]
        self.angsep = self.angsep[mask]
        self.maxangsep = self.maxangsep[mask]
        self.flux_p = self.flux_p[mask]
        self.fp = self.fp[mask]
        self.temp_p = self.temp_p[mask]
        self.nstar = self.nstar[mask]
        self.radius_s = self.radius_s[mask]
        self.mass_s = self.mass_s[mask]
        self.temp_s = self.temp_s[mask]
        self.distance_s = self.distance_s[mask]
        self.stype = self.stype[mask]
        self.ra = self.ra[mask]
        self.dec = self.dec[mask]