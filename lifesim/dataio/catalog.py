import warnings
import sys

import numpy as np
import pandas as pd
from astropy.io import fits

from lifesim.modules.habitable import compute_habitable_zone

# todo Documentation
class Catalog(object):
    """Gets and prints the spreadsheet's header columns

    Parameters
    ----------

    Returns
    -------
    """

    def __init__(self,
                 input_path: str):
        """Gets and prints the spreadsheet's header columns

        Parameters
        ----------

        Returns
        -------
        """

        self.data = pd.DataFrame(columns=['radius_p',
                                          'p_orb',
                                          'mass_p',
                                          'ecc_p',
                                          'inc_p',
                                          'large_omega_p',
                                          'small_omega_p',
                                          'theta_p',
                                          'albedo_bond',
                                          'albedo_geom_vis',
                                          'albedo_geom_mir',
                                          'z',
                                          'semimajor_p',
                                          'sep_p',
                                          'angsep',
                                          'maxangsep',
                                          'flux_p',
                                          'fp',
                                          'temp_p',
                                          'radius_s',
                                          'mass_s',
                                          'temp_s',
                                          'distance_s',
                                          'ra',
                                          'dec',
                                          'nuniverse',
                                          'nstar',
                                          'stype',
                                          'id'])

        self.stype_key = {'A': 0,
                          'F': 1,
                          'G': 2,
                          'K': 3,
                          'M': 4}

        self.masks = {}
        self.read_ppop(input_path=input_path)
        self.stype = None
        self.get_stype()
        self.create_mask()
        compute_habitable_zone(catalog=self,
                               model='MS')

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
            mask = np.logical_and(self.data.stype == stype,
                                  self.data.distance_s >= dist)
        elif mode == 'smaller':
            mask = np.logical_and(self.data.stype == stype,
                                  self.data.distance_s <= dist)
        else:
            warnings.warn('Mode ' + mode + ' not available. Using mode larger.')
            mask = np.logical_and(self.data.stype == stype,
                                  self.data.distance_s >= dist)

        self.data = self.data.drop(np.where(mask)[0])
        self.data.index = np.arange(0, self.data.shape[0], 1)

    def get_stype(self):
        self.stype = np.zeros_like(self.data.nstar, dtype=str)
        for _, k in enumerate(self.stype_key):
            self.stype[self.data.stype == self.stype_key[k]] = k

    def read_ppop(self,
                  input_path: str):
        """Read the contents of the P-Pop output .txt file to a catalog

        Parameters
        ----------
        input_path : str
            path to the P-Pop output file in a .txt format

        Returns
        -------
        """

        if input_path[-4:] == '.txt':
            table = open(input_path, 'r')
            table_lines = table.readlines()

            nuniverse = []
            radius_p = []  # Rearth
            p_orb = []  # d
            mass_p = []  # Mearth
            ecc_p = []
            inc_p = []  # rad
            large_omega_p = []  # rad
            small_omega_p = []  # rad
            theta_p = []  # rad
            albedo_bond = []
            albedo_geom_vis = []
            albedo_geom_mir = []
            z = []
            semimajor_p = []  # au
            sep_p = []  # au
            angsep = []  # arcsec
            maxangsep = []  # arcsec
            flux_p = []  # Searth
            fp = []
            temp_p = []  # K
            nstar = []
            radius_s = []  # Rsun
            mass_s = []  # Msun
            temp_s = []  # K
            distance_s = []  # pc
            stype = []
            ra = []  # deg
            dec = []  # deg

            nlines = len(table_lines)

            # The second line (i = 1) contains the column names of the new
            # Ppop while the first line (i = 0) contains the column names
            # of the old Ppop.

            tempLine = table_lines[1].split('\t')
            col_nuniverse = np.where(np.array(tempLine) == 'Nuniverse')[0][0]
            col_radius_p = np.where(np.array(tempLine) == 'Rp')[0][0]
            col_p_orb = np.where(np.array(tempLine) == 'Porb')[0][0]
            col_mass_p = np.where(np.array(tempLine) == 'Mp')[0][0]
            col_ecc_p = np.where(np.array(tempLine) == 'ep')[0][0]
            col_inc_p = np.where(np.array(tempLine) == 'ip')[0][0]
            col_large_omega_p = np.where(np.array(tempLine) == 'Omegap')[0][0]
            col_small_omega_p = np.where(np.array(tempLine) == 'omegap')[0][0]
            col_theta_p = np.where(np.array(tempLine) == 'thetap')[0][0]
            col_albedo_bond = np.where(np.array(tempLine) == 'Abond')[0][0]
            col_albedo_geom_vis = np.where(np.array(tempLine) == 'AgeomVIS')[0][0]
            col_albedo_geom_mir = np.where(np.array(tempLine) == 'AgeomMIR')[0][0]
            col_z = np.where(np.array(tempLine) == 'z')[0][0]
            col_semimajor_p = np.where(np.array(tempLine) == 'ap')[0][0]
            col_sep_p = np.where(np.array(tempLine) == 'rp')[0][0]
            col_angsep = np.where(np.array(tempLine) == 'AngSep')[0][0]
            col_maxangsep = np.where(np.array(tempLine) == 'maxAngSep')[0][0]
            col_flux_p = np.where(np.array(tempLine) == 'Fp')[0][0]
            col_fp = np.where(np.array(tempLine) == 'fp')[0][0]
            col_temp_p = np.where(np.array(tempLine) == 'Tp')[0][0]
            col_nstar = np.where(np.array(tempLine) == 'Nstar')[0][0]
            col_radius_s = np.where(np.array(tempLine) == 'Rs')[0][0]
            col_mass_s = np.where(np.array(tempLine) == 'Ms')[0][0]
            col_temp_s = np.where(np.array(tempLine) == 'Ts')[0][0]
            col_distance_s = np.where(np.array(tempLine) == 'Ds')[0][0]
            col_stype = np.where(np.array(tempLine) == 'Stype')[0][0]
            col_ra = np.where(np.array(tempLine) == 'RA')[0][0]
            col_dec = np.where(np.array(tempLine) == 'Dec')[0][0]

            for i, line in enumerate(table_lines[2:]):

                if (i % 10000) == 0:
                    sys.stdout.write('\rProcessed line %.0f of %.0f' % (i, nlines))
                    sys.stdout.flush()

                tempLine = line.split('\t')
                nuniverse += [int(tempLine[col_nuniverse])]
                radius_p += [float(tempLine[col_radius_p])]  # Rearth
                p_orb += [float(tempLine[col_p_orb])]  # d
                mass_p += [float(tempLine[col_mass_p])]  # Mearth
                ecc_p += [float(tempLine[col_ecc_p])]
                inc_p += [float(tempLine[col_inc_p])]  # rad
                large_omega_p += [float(tempLine[col_large_omega_p])]  # rad
                small_omega_p += [float(tempLine[col_small_omega_p])]  # rad
                theta_p += [float(tempLine[col_theta_p])]  # rad
                albedo_bond += [float(tempLine[col_albedo_bond])]
                albedo_geom_vis += [float(tempLine[col_albedo_geom_vis])]
                albedo_geom_mir += [float(tempLine[col_albedo_geom_mir])]
                z += [float(tempLine[col_z])]
                semimajor_p += [float(tempLine[col_semimajor_p])]  # au
                sep_p += [float(tempLine[col_sep_p])]  # au
                angsep += [float(tempLine[col_angsep])]  # arcsec
                maxangsep += [float(tempLine[col_maxangsep])]  # arcsec
                flux_p += [float(tempLine[col_flux_p])]  # Searth
                fp += [float(tempLine[col_fp])]
                temp_p += [float(tempLine[col_temp_p])]  # K
                nstar += [int(tempLine[col_nstar])]
                radius_s += [float(tempLine[col_radius_s])]  # Rsun
                mass_s += [float(tempLine[col_mass_s])]  # Msun
                temp_s += [float(tempLine[col_temp_s])]  # K
                distance_s += [float(tempLine[col_distance_s])]  # pc
                stype += [str(tempLine[col_stype])]
                ra += [float(tempLine[col_ra])]  # deg
                dec += [float(tempLine[col_dec])]  # deg

            sys.stdout.write('\rProcessed line %.0f of %.0f' % (nlines, nlines))
            sys.stdout.flush()
            print('')

            self.data['nuniverse'] = np.array(nuniverse).astype(int)
            self.data['radius_p'] = np.array(radius_p).astype(float)
            self.data['p_orb'] = np.array(p_orb)
            self.data['mass_p'] = np.array(mass_p)
            self.data['ecc_p'] = np.array(ecc_p)
            self.data['inc_p'] = np.array(inc_p)
            self.data['large_omega_p'] = np.array(large_omega_p)
            self.data['small_omega_p'] = np.array(small_omega_p)
            self.data['theta_p'] = np.array(theta_p)
            self.data['albedo_bond'] = np.array(albedo_bond)
            self.data['albedo_geom_vis'] = np.array(albedo_geom_vis)
            self.data['albedo_geom_mir'] = np.array(albedo_geom_mir)
            self.data['z'] = np.array(z)
            self.data['semimajor_p'] = np.array(semimajor_p)
            self.data['sep_p'] = np.array(sep_p)
            self.data['angsep'] = np.array(angsep)
            self.data['maxangsep'] = np.array(maxangsep)
            self.data['flux_p'] = np.array(flux_p)
            self.data['fp'] = np.array(fp)
            self.data['temp_p'] = np.array(temp_p)
            self.data['nstar'] = np.array(nstar)
            self.data['radius_s'] = np.array(radius_s)
            self.data['mass_s'] = np.array(mass_s)
            self.data['temp_s'] = np.array(temp_s)
            self.data['distance_s'] = np.array(distance_s)
            self.data['stype'] = np.array(stype)
            self.data['ra'] = np.array(ra)
            self.data['dec'] = np.array(dec)
            self.data['id'] = np.arange(0, len(dec), 1)

        elif input_path[-5:] == '.fits':
            hdu = fits.open(input_path)
            stype_int = np.zeros_like(hdu[1].data.Nstar.astype(int), dtype=int)
            for _, k in enumerate(self.stype_key.keys()):
                stype_int[hdu[1].data.Stype.astype(str) == k] = self.stype_key[k]

            self.data = pd.DataFrame({'nuniverse': hdu[1].data.Nuniverse.astype(int),
                                      'nstar': hdu[1].data.Nstar.astype(int),
                                      'stype': stype_int,
                                      'id': np.arange(0, hdu[1].data.Dec.shape[0], 1).astype(int),
                                      'radius_p': hdu[1].data.Rp.astype(float),
                                      'p_orb': hdu[1].data.Porb.astype(float),
                                      'mass_p': hdu[1].data.Mp.astype(float),
                                      'ecc_p': hdu[1].data.ep.astype(float),
                                      'inc_p': hdu[1].data.ip.astype(float),
                                      'large_omega_p': hdu[1].data.Omegap.astype(float),
                                      'small_omega_p': hdu[1].data.omegap.astype(float),
                                      'theta_p': hdu[1].data.thetap.astype(float),
                                      'albedo_bond': hdu[1].data.Abond.astype(float),
                                      'albedo_geom_vis': hdu[1].data.AgeomVIS.astype(float),
                                      'albedo_geom_mir': hdu[1].data.AgeomMIR.astype(float),
                                      'z': hdu[1].data.z.astype(float),
                                      'semimajor_p': hdu[1].data.ap.astype(float),
                                      'sep_p': hdu[1].data.rp.astype(float),
                                      'angsep': hdu[1].data.AngSep.astype(float),
                                      'maxangsep': hdu[1].data.maxAngSep.astype(float),
                                      'flux_p': hdu[1].data.Fp.astype(float),
                                      'fp': hdu[1].data.fp.astype(float),
                                      'temp_p': hdu[1].data.Tp.astype(float),
                                      'radius_s': hdu[1].data.Rs.astype(float),
                                      'mass_s': hdu[1].data.Ms.astype(float),
                                      'temp_s': hdu[1].data.Ts.astype(float),
                                      'distance_s': hdu[1].data.Ds.astype(float),
                                      'ra': hdu[1].data.RA.astype(float),
                                      'dec': hdu[1].data.Dec.astype(float)})
            hdu.close()

    def create_mask(self):
        _, temp = np.unique(self.data.nstar, return_index=True)
        self.masks['stars'] = np.zeros_like(self.data.nstar, dtype=bool)
        self.masks['stars'][temp] = True

    def safe_add(self,
                 name: str,
                 data: np.ndarray):
        if data.ndim != 1:
            raise ValueError('Only one-dimensional data can be added')
        if data.shape[0] != self.data.nstar.shape[0]:
            raise ValueError('Creation of ragged database is supressed')
        if name in self.data.keys():
            raise ValueError('Data can not be overwritten in safe mode')
        self.data[name] = data
