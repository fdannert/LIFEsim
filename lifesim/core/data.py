import sys
import warnings
import time

import numpy as np
import xarray as xr
import pandas as pd
from astropy.io import fits
from tqdm import tqdm
import pickle

from astropy.coordinates import SkyCoord, BarycentricMeanEcliptic

from lifesim.util.options import Options
from lifesim.util.habitable import single_habitable_zone


# TODO: automatically add data storage for all
class Data(object):
    """
    The data class is the central storage class for catalogs, options, parameters and data. Any
    data used in simulations should be stored in this class. Via the bus, access to the data class
    is given to all modules.

    Attributes
    ----------
    inst : dict
        Data used for simulation of the instrument.
    catalog : pd.DataFrame
        Catalog containing all exoplanets in the sample.
    single : dict
        Data used for the spectral simulation of single exoplanets.
    other : dict
        Data storage for any other pertinent data.
    options : Options
        Location of the Options class. All options and free parameters used in a LIFEsim simulation
        must be stored here.
    """
    def __init__(self):
        self.inst = {}
        self.catalog = None
        self.noise_catalog = None
        self.noise_catalog_pivot = None
        self.single = {}
        self.other = {}
        self.options = Options()
        self.optm = {}

    def catalog_delete(self):
        self.catalog = None

    def catalog_from_ppop(self,
                          input_path: str,
                          overwrite: bool = False):
        """
        Read the contents of the P-Pop output file (in .txt or .fits format) to a catalog. Note that reading catalogs
        in .fits format is significantly faster.

        Parameters
        ----------
        input_path : str
            Path to the P-Pop output file in a .txt or .fits format.
        overwrite : bool
            If set to true, existing catalogs can overwritten.

        Raises
        ------
        ValueError
            If the data class already has an initialized catalog and overwrite is set to False.
        """

        # make sure that no catalog exists
        if (self.catalog is not None) and (not overwrite):
            raise ValueError('A catalog has already been imported. Delete the old catalog or set '
                             'overwrite=True')

        print('Loading catalog from P-Pop...')

        self.options.other['database_path'] = input_path

        # initialize catalog
        self.catalog = pd.DataFrame(columns=['radius_p',
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
                                             'id',
                                             'name_s'])

        # check the format of the input file
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
            name_s = []

            nlines = len(table_lines)

            # The second line (i = 1) contains the column names of the new
            # Ppop while the first line (i = 0) contains the column names
            # of the old Ppop.

            tempLine = table_lines[1].split('\t')

            get_name = ('name' in tempLine)

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

            if get_name:
                col_name_s = np.where(np.array(tempLine) == 'name')[0][0]

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

                if get_name:
                    name_s += [str(tempLine[col_name_s])]

            sys.stdout.write('\rProcessed line %.0f of %.0f' % (nlines, nlines))
            sys.stdout.flush()
            print('')

            # remove the newline string from the star names
            if get_name:
                name_s = [name.replace('\n', '') for name in name_s]
            else:
                name_s = ['None'] * len(dec)

            # save the data to the pandas DataFrame
            self.catalog['nuniverse'] = np.array(nuniverse).astype(int)
            self.catalog['radius_p'] = np.array(radius_p).astype(float)
            self.catalog['p_orb'] = np.array(p_orb)
            self.catalog['mass_p'] = np.array(mass_p)
            self.catalog['ecc_p'] = np.array(ecc_p)
            self.catalog['inc_p'] = np.array(inc_p)
            self.catalog['large_omega_p'] = np.array(large_omega_p)
            self.catalog['small_omega_p'] = np.array(small_omega_p)
            self.catalog['theta_p'] = np.array(theta_p)
            self.catalog['albedo_bond'] = np.array(albedo_bond)
            self.catalog['albedo_geom_vis'] = np.array(albedo_geom_vis)
            self.catalog['albedo_geom_mir'] = np.array(albedo_geom_mir)
            self.catalog['z'] = np.array(z)
            self.catalog['semimajor_p'] = np.array(semimajor_p)
            self.catalog['sep_p'] = np.array(sep_p)
            self.catalog['angsep'] = np.array(angsep)
            self.catalog['maxangsep'] = np.array(maxangsep)
            self.catalog['flux_p'] = np.array(flux_p)
            self.catalog['fp'] = np.array(fp)
            self.catalog['temp_p'] = np.array(temp_p)
            self.catalog['nstar'] = np.array(nstar)
            self.catalog['radius_s'] = np.array(radius_s)
            self.catalog['mass_s'] = np.array(mass_s)
            self.catalog['temp_s'] = np.array(temp_s)
            self.catalog['distance_s'] = np.array(distance_s)
            self.catalog['stype'] = pd.Series(stype, dtype=pd.StringDtype())
            self.catalog['ra'] = np.array(ra)
            self.catalog['dec'] = np.array(dec)
            self.catalog['id'] = np.arange(0, len(dec), 1)
            self.catalog['name_s'] = pd.Series(name_s, dtype=pd.StringDtype())

        # check the format of the input file
        elif input_path[-5:] == '.fits':
            hdu = fits.open(input_path)

            get_name = ('name' in hdu[1].columns.names)

            # save the data to the pandas DataFrame, make sure it is saved in the correct type
            # some errors were produced here by not respecting the endianess of the data (numpy
            # usually works with little endian)
            if get_name:
                self.catalog = pd.DataFrame({'nuniverse': hdu[1].data.Nuniverse.astype(int),
                                             'nstar': hdu[1].data.Nstar.astype(int),
                                             'stype': pd.Series(hdu[1].data.Stype, dtype=pd.StringDtype()),
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
                                             'dec': hdu[1].data.Dec.astype(float),
                                             'lat': hdu[1].data.lat.astype(float),
                                             'lon': hdu[1].data.lon.astype(float),
                                             'name_s': pd.Series(hdu[1].data.name, dtype=pd.StringDtype())})
            else:
                self.catalog = pd.DataFrame({'nuniverse': hdu[1].data.Nuniverse.astype(int),
                                             'nstar': hdu[1].data.Nstar.astype(int),
                                             'stype': pd.Series(hdu[1].data.Stype, dtype=pd.StringDtype()),
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
                                             'dec': hdu[1].data.Dec.astype(float),
                                             'lat': hdu[1].data.lat.astype(float),
                                             'lon': hdu[1].data.lon.astype(float),
                                             'name_s': pd.Series(['None']*hdu[1].data.shape[0], dtype=pd.StringDtype())})
            hdu.close()

        # create mask returning only unique stars
        _, temp = np.unique(self.catalog.nstar, return_index=True)
        star_mask = np.zeros_like(self.catalog.nstar, dtype=bool)
        star_mask[temp] = True

        # TODO: why is this commented out? AFAIK P-Pop uses equitorial coordinates
        # transform from equitorial to ecliptic coordinates
        coord = SkyCoord(self.catalog.ra, self.catalog.dec, frame='icrs', unit='deg')
        coord_ec = coord.transform_to(BarycentricMeanEcliptic())
        self.catalog['lon'] = np.array(coord_ec.lon.radian)
        self.catalog['lat'] = np.array(coord_ec.lat.radian)

        # add the inner/ outer edges and centers of the habitable zone
        s_in = np.zeros_like(self.catalog.nstar, dtype=float)
        s_out = np.zeros_like(self.catalog.nstar, dtype=float)
        l_sun = np.zeros_like(self.catalog.nstar, dtype=float)

        hz_in = np.zeros_like(self.catalog.nstar, dtype=float)
        hz_out = np.zeros_like(self.catalog.nstar, dtype=float)
        hz_center = np.zeros_like(self.catalog.nstar, dtype=float)

        for _, n in enumerate(np.where(star_mask)[0]):
            s_in[self.catalog.nstar == self.catalog.nstar[n]], \
            s_out[self.catalog.nstar == self.catalog.nstar[n]], \
            l_sun[self.catalog.nstar == self.catalog.nstar[n]], \
            hz_in[self.catalog.nstar == self.catalog.nstar[n]], \
            hz_out[self.catalog.nstar == self.catalog.nstar[n]], \
            hz_center[self.catalog.nstar == self.catalog.nstar[n]] \
                = single_habitable_zone(model=self.options.models['habitable'],
                                        temp_s=self.catalog.temp_s[n],
                                        radius_s=self.catalog.radius_s[n])

        self.catalog['s_in'] = s_in
        self.catalog['s_out'] = s_out
        self.catalog['l_sun'] = l_sun
        self.catalog['hz_in'] = hz_in
        self.catalog['hz_out'] = hz_out
        self.catalog['hz_center'] = hz_center
        self.catalog['habitable'] = np.logical_and.reduce((
            (self.catalog['semimajor_p'] > self.catalog['hz_in']).to_numpy(),
            (self.catalog['semimajor_p'] < self.catalog['hz_out']).to_numpy(),
            (self.catalog['radius_p'].ge(0.5)).to_numpy(),
            (self.catalog['radius_p'].le(1.5)).to_numpy()))

        print('[Done]')

    # TODO: Definition of stype here is wrong. It should be an int not a string.
    #   Think about this a bit more. It is important to keep the ints in the DataFrame, but that
    #   decreases usability. Maybe an option is to use intermediate masks for the stellar types.
    def catalog_remove_distance(self,
                                stype: str,
                                dist: float,
                                mode: str):
        """
        Removes planets around stellar type above or below a certain distance from the sample. A
        warning is raise if the mode is not recognized.

        Parameters
        ----------
        stype : str
            Planets around stars of the specified stellar type are removed. Possible options are
            'A, 'F', 'G', 'K' and 'M'.
        dist : float
            Specifies the distance over or under which the planets are removed in pc.
        mode : str
            'larger':   planets around stars with distances larger than the specified distance are
                        removed.
            'smaller':  planets around stars with distances smaller than the specified distance are
                        removed.
        """

        # check if stellar type is valid
        if not np.isin(stype, np.array(('A', 'F', 'G', 'K', 'M'))):
            raise ValueError('Stellar type not recognised')

        # create masks selecting closer or more far away planets
        if mode == 'larger':
            mask = np.logical_and(self.catalog.stype == stype,
                                  self.catalog.distance_s >= dist)
        elif mode == 'smaller':
            mask = np.logical_and(self.catalog.stype == stype,
                                  self.catalog.distance_s <= dist)
        else:
            warnings.warn('Mode ' + mode + ' not available. Using mode larger.')
            mask = np.logical_and(self.catalog.stype == stype,
                                  self.catalog.distance_s >= dist)

        # drop the planets to remove from the data
        self.catalog = self.catalog.drop(np.where(mask)[0])
        self.catalog.index = np.arange(0, self.catalog.shape[0], 1)

    def catalog_safe_add(self,
                         name: str,
                         data: np.ndarray):
        """
        Add data to the catalog while making sure that none of the catalogs vital properties are
        disrespected and that no data is deleted.

        Parameters
        ----------
        name : str
            Name key of the data in the pandas DataFrame.
        data : np.ndarray
            The data added to the catalog.

        Raises
        ------
        ValueError
            If the data is not one-dimensional.
            If the array does not have the same length as the other data arrays in the catalog.
            If the catalog already contains data under the given name.
        """

        if data.ndim != 1:
            raise ValueError('Only one-dimensional data can be added')
        if data.shape[0] != self.catalog.nstar.shape[0]:
            raise ValueError('Creation of ragged database is supressed')
        if name in self.catalog.keys():
            raise ValueError('Data can not be overwritten in safe mode')
        self.catalog[name] = data

    def export_catalog(self):
        """
        Save the catalog to an file in the hdf-format.

        Parameters
        ----------
        output_path : str
            path to the new file.

        Raises
        ------
        ValueError
            If not catalog exists in this data class.
        """
        if self.catalog is None:
            raise ValueError('No catalog found')

        self.str_to_obj(reverse=False)
        self.catalog.to_hdf(path_or_buf=self.options.other['output_path']
                                        + self.options.other['output_filename'] + '_catalog.hdf5',
                            key='catalog',
                            mode='w')
        self.str_to_obj(reverse=True)

        print('Main Catalog Stored')

        print('Exporting Noise Catalog...')
        if self.noise_catalog is not None:
            self.noise_catalog.to_netcdf(path=self.options.other['output_path']
                                              + self.options.other['output_filename'] + '_noise.nc',
                                         mode='w',
                                         engine='h5netcdf')
        print('[Done]')
        # if (self.options.other['pickle_mode'] == True) and (self.noise_catalog is not None):
        #     file = open(self.options.other['output_path']
        #                 + self.options.other['output_filename'] + '_noise_pickle.pickle', 'wb')
        #     pickle.dump(self.noise_catalog, file)
        #     file.close()
        # elif self.options.other['large_file']:
        #     if self.noise_catalog_pivot is not None:
        #         self.store_pivot_noise_catalog()
        #         print('Pivoted Noise Catalog Stored')
        #     elif self.noise_catalog is not None:
        #         self.pivot_noise_catalog(to_wavelength=True)
        #         self.store_pivot_noise_catalog()
        #         print('Pivoted Noise Catalog Stored')
        #     else:
        #         print('No Noise Catalog Stored')
        # else:
        #     if self.noise_catalog is not None:
        #         self.store_noise_catalog()
        #         print('Standard Noise Catalog Stored')
        #     elif self.noise_catalog_pivot is not None:
        #         self.pivot_noise_catalog(to_wavelength=False)
        #         self.store_noise_catalog()
        #         print('Standard Noise Catalog Stored')
        #     else:
        #         print('No Noise Catalog Stored')


        # if self.noise_catalog is not None:
        #     print('Exporting Noise Catalog...')
        #     if self.options.other['large_file']:
        #         self.pivot_noise_catalog(to_wavelength=True)
        #         store = pd.HDFStore(self.options.other['output_path']
        #                             + self.options.other['output_filename'] + '_noise_large.hdf5')
        #         for k in self.noise_catalog_pivot.keys():
        #             store.put(key='id_' + k.replace('.', ''), value=self.noise_catalog_pivot[k].astype(float))
        #         store.put(key='wl_keys', value=pd.Series(
        #             ['id_' + k.replace('.', '') for k in self.noise_catalog_pivot.keys()]
        #         ))
        #         store.close()
        #     else:
        #         store = pd.HDFStore(self.options.other['output_path']
        #                             + self.options.other['output_filename'] + '_noise.hdf5')
        #         for k in self.noise_catalog.keys():
        #             store.put(key='id_' + k, value=self.noise_catalog[k].astype(float))
        #         store.close()
        #     print('[Done]')
        #     self.noise_catalog.to_hdf(path_or_buf=self.options.other['output_path']
        #                                           + self.options.other['output_filename']
        #                                           + '_noise.hdf5', key='noise_catalog', mode='w')

    # def store_pivot_noise_catalog(self):
    #     store = pd.HDFStore(self.options.other['output_path']
    #                         + self.options.other['output_filename'] + '_noise_large.hdf5')
    #     for k in self.noise_catalog_pivot.keys():
    #         store.put(key='id_' + k.replace('.', ''), value=self.noise_catalog_pivot[k].astype(float))
    #     store.put(key='wl_keys', value=pd.Series(
    #         ['id_' + k.replace('.', '') for k in self.noise_catalog_pivot.keys()]
    #     ))
    #     store.close()
    #
    # def store_noise_catalog(self):
    #     store = pd.HDFStore(self.options.other['output_path']
    #                         + self.options.other['output_filename'] + '_noise.hdf5')
    #     for k in self.noise_catalog.keys():
    #         store.put(key='id_' + k, value=self.noise_catalog[k].astype(float))
    #     store.close()

    def import_catalog(self,
                       input_path: str,
                       overwrite: bool = False,
                       noise_catalog: bool = False):
        """
        Import catalog from external file of hdf-format.

        Parameters
        ----------
        input_path : str
            path to the P-Pop output file in a .txt or .fits format.
        overwrite : bool
            if set to true, existing catalogs can overwritten.

        Raises
        ------
        ValueError
            If the data class already has an initialized catalog and overwrite is set to False.
        """

        print('Importing Catalog...')
        if (self.catalog is not None) and (not overwrite):
            raise ValueError('Can not overwrite existing catalog')

        self.options.other['database_path'] = input_path

        print('Beginning Import...')
        t0 = time.time()
        self.catalog = pd.read_hdf(path_or_buf=input_path,
                                   key='catalog')
        print('Import completed (Time: ' + str(time.time()-t0) + '), changing string object types...')
        t0 = time.time()
        self.str_to_obj(reverse=True)

        print('[Done] (Time: ' + str(time.time()-t0) + ')')

        if noise_catalog:
            print('Importing Noise Catalog...')
            with xr.open_dataarray(input_path[:-5] + '_noise.nc',
                                   engine='h5netcdf') as file:
                self.noise_catalog = file
            print('[Done]')
            # if self.options.other['pickle_mode']:
            #     file = open(input_path[:-5] + '_noise_pickle.pickle', 'rb')
            #     self.noise_catalog = pickle.load(file)
            #     file.close()
            # elif self.options.other['large_file']:
            #     store = pd.HDFStore(input_path[:-5] + '_noise_large.hdf5')
            #     self.noise_catalog_pivot = {}
            #     wl_keys = store.get('wl_keys')
            #     for wl_key in tqdm(wl_keys):
            #         self.noise_catalog_pivot[wl_key[3:-1] + '.' + wl_key[-1]] = store.get(wl_key)
            #     store.close()
            # else:
            #     store = pd.HDFStore(input_path[:-5] + '_noise.hdf5')
            #     self.noise_catalog = {}
            #     for id in tqdm(self.catalog.id):
            #         self.noise_catalog[str(id)] = store.get('id_' + str(id))
            #     store.close()

    # def pivot_noise_catalog(self,
    #                         to_wavelength: bool):
    #     print('Pivoting Noise Catalog...')
    #     if to_wavelength:
    #         self.noise_catalog_pivot = {}
    #         idx = list(self.noise_catalog.keys())
    #         wl_ids = self.noise_catalog[idx[0]].index.values
    #         columns = self.noise_catalog[idx[0]].columns.values
    #         for wl_id in tqdm(wl_ids):
    #             pd_table = pd.DataFrame(columns=columns, index=idx)
    #             for id in idx:
    #                 pd_table.loc[id] = self.noise_catalog[id].loc[wl_id]
    #             self.noise_catalog_pivot[wl_id] = pd_table
    #         self.noise_catalog = None
    #
    #     else:
    #         self.noise_catalog = {}
    #         wl_ids = list(self.noise_catalog_pivot.keys())
    #         idx = self.noise_catalog_pivot[wl_ids[0]].index.values
    #         columns = self.noise_catalog_pivot[wl_ids[0]].columns.values
    #         for id in tqdm(idx):
    #             pd_table = pd.DataFrame(columns=columns, index=wl_ids)
    #             for wl_id in wl_ids:
    #                 pd_table.loc[wl_id] = self.noise_catalog_pivot[wl_id].loc[id]
    #             self.noise_catalog[id] = pd_table
    #         self.noise_catalog_pivot = None
    #     print('')
    #     print('[Done]')


    def noise_catalog_from_catalog(self):
        pass
        # self.noise_catalog = pd.DataFrame(columns=['signal',  # planet signal
        #                                            'noise',  # overall noise contribution
        #                                            'wl',  # wavelength bin
        #                                            'pn_sgl',  # stellar geometric leakage
        #                                            'pn_ez',  # exozodi leakage
        #                                            'pn_lz',  # localzodi leakage
        #                                            'pn_dc',  # dark current
        #                                            'pn_tbd',  # thermal background detector
        #                                            'pn_tbpm',  # thermal background primary mirror
        #                                            'pn_pa',  # polarization angle
        #                                            'pn_snfl',  # stellar null floor leakage
        #                                            'pn_ag_cld',  # agnostic cold instrumental photon noise
        #                                            'pn_ag_ht',  # agnostic hot instrumental photon noise
        #                                            'pn_ag_wht',  # agnostic white instrumental photon noise
        #                                            'pn',  # photon noise
        #                                            'sn_fo_a',  # first order amplitude
        #                                            'sn_fo_phi',  # first order phase
        #                                            'sn_fo_x',  # first order x position
        #                                            'sn_fo_y',  # first order y position
        #                                            'sn_fo',  # systematic noise first order
        #                                            'sn_so_aa',  # second order amplitude-amplitude term
        #                                            'sn_so_phiphi',  # second order phase-phase term
        #                                            'sn_so_aphi',  # amplitude phase cross term
        #                                            'sn_so_polpol',  # second order polarization-polarization term
        #                                            'sn_so',  # systematic noise second order
        #                                            'sn',  # systematic noise
        #                                            'fundamental',  # fundamental noise (astrophysical)
        #                                            'instrumental',  # instrumental noise
        #                                            'snr'  # signal to noise ratio
        #                                            ],
        #                                   index=self.catalog.id)

        # self.noise_catalog =

    def str_to_obj(self,
                   reverse: bool):
        """
        Converts all string type columns in the catalog between type 'object' (needed for saving to hdf5) and type
        'pandas.StringDtype' (needed for fast computation).

        Parameters
        ----------
        reverse : bool
            if reveres is set true, the type will be converted 'object' -> 'pandas.StringDtype'
        """
        if not reverse:
            for key in list(set(self.catalog.keys()[np.where(self.catalog.dtypes == 'string')])
                            & {'name_s', 'stype'}):
                self.catalog[key] = self.catalog[key].astype(object)
        else:
            for key in list(set(self.catalog.keys()[np.where(self.catalog.dtypes == 'object')])
                            & {'name_s', 'stype'}):
                self.catalog[key] = self.catalog[key].astype(pd.StringDtype())
