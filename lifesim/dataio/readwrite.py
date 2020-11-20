import sys

import numpy as np
import pandas as pd
from astropy.io import fits

from lifesim.dataio.catalog import Catalog


# TODO remove this file (it is not used)
def read_ppop(tc: Catalog,
              input_path: str):
    """Read the contents of the P-Pop output .txt file to a catalog

    Parameters
    ----------
    tc : object
        catalog to which the planet population should be written
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

        tc.data['nuniverse'] = np.array(nuniverse).astype(int)
        tc.data['radius_p'] = np.array(radius_p).astype(float)
        tc.data['p_orb'] = np.array(p_orb)
        tc.data['mass_p'] = np.array(mass_p)
        tc.data['ecc_p'] = np.array(ecc_p)
        tc.data['inc_p'] = np.array(inc_p)
        tc.data['large_omega_p'] = np.array(large_omega_p)
        tc.data['small_omega_p'] = np.array(small_omega_p)
        tc.data['theta_p'] = np.array(theta_p)
        tc.data['albedo_bond'] = np.array(albedo_bond)
        tc.data['albedo_geom_vis'] = np.array(albedo_geom_vis)
        tc.data['albedo_geom_mir'] = np.array(albedo_geom_mir)
        tc.data['z'] = np.array(z)
        tc.data['semimajor_p'] = np.array(semimajor_p)
        tc.data['sep_p'] = np.array(sep_p)
        tc.data['angsep'] = np.array(angsep)
        tc.data['maxangsep'] = np.array(maxangsep)
        tc.data['flux_p'] = np.array(flux_p)
        tc.data['fp'] = np.array(fp)
        tc.data['temp_p'] = np.array(temp_p)
        tc.data['nstar'] = np.array(nstar)
        tc.data['radius_s'] = np.array(radius_s)
        tc.data['mass_s'] = np.array(mass_s)
        tc.data['temp_s'] = np.array(temp_s)
        tc.data['distance_s'] = np.array(distance_s)
        tc.data['stype'] = np.array(stype)
        tc.data['ra'] = np.array(ra)
        tc.data['dec'] = np.array(dec)
        tc.data['id'] = np.arange(0, len(dec), 1)

    elif input_path[-5:] == '.fits':
        hdu = fits.open(input_path)
        stype_int = np.zeros_like(hdu[1].data.Nstar.astype(int), dtype=int)
        for _, k in enumerate(tc.stype_key.keys()):
            stype_int[hdu[1].data.Stype.astype(str) == k] = tc.stype_key[k]

        tc.data = pd.DataFrame({'nuniverse': hdu[1].data.Nuniverse.astype(int),
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

        tc.parameters = pd.DataFrame({'nuniverse': hdu[1].data.Nuniverse.astype(int),
                                      'nstar': hdu[1].data.Nstar.astype(int),
                                      'stype': hdu[1].data.Stype.astype(str),
                                      'id': np.arange(0, hdu[1].data.Dec.shape[0], 1).astype(int)})

        hdu.close()
