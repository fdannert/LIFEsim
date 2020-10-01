import sys

import numpy as np


def read_ppop(tc: object,
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

    table = open(input_path, 'r')
    table_lines = table.readlines()

    tc.nuniverse = []
    tc.radius_p = []  # Rearth
    tc.p_orb = []  # d
    tc.mass_p = []  # Mearth
    tc.ecc_p = []
    tc.inc_p = []  # rad
    tc.large_omega_p = []  # rad
    tc.small_omega_p = []  # rad
    tc.theta_p = []  # rad
    tc.albedo_bond = []
    tc.albedo_geom_vis = []
    tc.albedo_geom_mir = []
    tc.z = []
    tc.semimajor_p = []  # au
    tc.sep_p = []  # au
    tc.angsep = []  # arcsec
    tc.maxangsep = []  # arcsec
    tc.flux_p = []  # Searth
    tc.fp = []
    tc.temp_p = []  # K
    tc.nstar = []
    tc.radius_s = []  # Rsun
    tc.mass_s = []  # Msun
    tc.temp_s = []  # K
    tc.distance_s = []  # pc
    tc.stype = []
    tc.ra = []  # deg
    tc.dec = []  # deg

    nlines = len(table_lines)

    # The second line (i = 1) contains the column names of the new
    # P-pop while the first line (i = 0) contains the column names
    # of the old P-pop.

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
        tc.nuniverse += [int(tempLine[col_nuniverse])]
        tc.radius_p += [float(tempLine[col_radius_p])]  # Rearth
        tc.p_orb += [float(tempLine[col_p_orb])]  # d
        tc.mass_p += [float(tempLine[col_mass_p])]  # Mearth
        tc.ecc_p += [float(tempLine[col_ecc_p])]
        tc.inc_p += [float(tempLine[col_inc_p])]  # rad
        tc.large_omega_p += [float(tempLine[col_large_omega_p])]  # rad
        tc.small_omega_p += [float(tempLine[col_small_omega_p])]  # rad
        tc.theta_p += [float(tempLine[col_theta_p])]  # rad
        tc.albedo_bond += [float(tempLine[col_albedo_bond])]
        tc.albedo_geom_vis += [float(tempLine[col_albedo_geom_vis])]
        tc.albedo_geom_mir += [float(tempLine[col_albedo_geom_mir])]
        tc.z += [float(tempLine[col_z])]
        tc.semimajor_p += [float(tempLine[col_semimajor_p])]  # au
        tc.sep_p += [float(tempLine[col_sep_p])]  # au
        tc.angsep += [float(tempLine[col_angsep])]  # arcsec
        tc.maxangsep += [float(tempLine[col_maxangsep])]  # arcsec
        tc.flux_p += [float(tempLine[col_flux_p])]  # Searth
        tc.fp += [float(tempLine[col_fp])]
        tc.temp_p += [float(tempLine[col_temp_p])]  # K
        tc.nstar += [int(tempLine[col_nstar])]
        tc.radius_s += [float(tempLine[col_radius_s])]  # Rsun
        tc.mass_s += [float(tempLine[col_mass_s])]  # Msun
        tc.temp_s += [float(tempLine[col_temp_s])]  # K
        tc.distance_s += [float(tempLine[col_distance_s])]  # pc
        tc.stype += [str(tempLine[col_stype])]
        tc.ra += [float(tempLine[col_ra])]  # deg
        tc.dec += [float(tempLine[col_dec])]  # deg

    sys.stdout.write('\rProcessed line %.0f of %.0f' % (nlines, nlines))
    sys.stdout.flush()
    print('')

    tc.nuniverse = np.array(tc.nuniverse)
    tc.radius_p = np.array(tc.radius_p)
    tc.p_orb = np.array(tc.p_orb)
    tc.mass_p = np.array(tc.mass_p)
    tc.ecc_p = np.array(tc.ecc_p)
    tc.inc_p = np.array(tc.inc_p)
    tc.large_omega_p = np.array(tc.large_omega_p)
    tc.small_omega_p = np.array(tc.small_omega_p)
    tc.theta_p = np.array(tc.theta_p)
    tc.albedo_bond = np.array(tc.albedo_bond)
    tc.albedo_geom_vis = np.array(tc.albedo_geom_vis)
    tc.albedo_geom_mir = np.array(tc.albedo_geom_mir)
    tc.z = np.array(tc.z)
    tc.semimajor_p = np.array(tc.semimajor_p)
    tc.sep_p = np.array(tc.sep_p)
    tc.angsep = np.array(tc.angsep)
    tc.maxangsep = np.array(tc.maxangsep)
    tc.flux_p = np.array(tc.flux_p)
    tc.fp = np.array(tc.fp)
    tc.temp_p = np.array(tc.temp_p)
    tc.nstar = np.array(tc.nstar)
    tc.radius_s = np.array(tc.radius_s)
    tc.mass_s = np.array(tc.mass_s)
    tc.temp_s = np.array(tc.temp_s)
    tc.distance_s = np.array(tc.distance_s)
    tc.stype = np.array(tc.stype)
    tc.ra = np.array(tc.ra)
    tc.dec = np.array(tc.dec)
