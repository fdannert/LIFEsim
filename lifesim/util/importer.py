from astropy import units as u
from astropy import constants as c
import numpy as np
from typing import Union


# TODO: Documentation
class SpectrumImporter(object):

    def __init__(self):
        self.pathtotxt = None
        self.y_string = None
        self.x_string = None
        self.radius_p_spectrum = None  # in earth radii
        self.radius_p_target = None  # in earth radii
        self.integration_time = None  # in s
        self.distance_s_spectrum = None  # in pc
        self.distance_s_target = None  # in pc

        self.x_data = None
        self.y_data = None

        self.x_raw = None
        self.y_raw = None

        # order: update_x/update_y -> import_spectrum -> to_wavelength -> to_counts

    def to_wavelength(self):
        if self.x_data.unit.decompose().bases == [u.s]:
            target = self.y_data.unit / u.m * u.Hz

            self.x_data = self.x_data.to(u.m, equivalencies=u.spectral())
            self.y_data = self.y_data.to(target, equivalencies=u.spectral_density(self.x_data)).decompose()

    def to_counts(self):
        if u.kg in self.y_data.unit.decompose().bases:
            if u.m in self.x_data.unit.bases:
                self.y_data = self.y_data / (c.c * c.h / self.x_data)
            else:
                self.y_data = self.y_data / (c.h * self.x_data)
            self.y_data = self.y_data.decompose()

    def to_per_second(self):
        if self.integration_time != 0.:
            self.y_data = (self.y_data / (self.integration_time * u.s)).decompose()

    def adjust_to_target(self):
        if self.distance_s_spectrum == 0:
            self.y_data = self.y_data / (self.distance_s_target * u.pc)**2 / np.pi
        elif self.distance_s_spectrum is None:
            pass
        else:
            self.y_data = self.y_data * (self.distance_s_spectrum / self.distance_s_target)**2
        if self.radius_p_spectrum == 0:
            self.y_data = self.y_data * np.pi * (self.radius_p_target * c.R_earth)**2
        elif self.radius_p_spectrum is None:
            pass
        else:
            self.y_data = self.y_data * (self.radius_p_target / self.radius_p_spectrum)**2
        self.y_data = self.y_data.decompose()

    def remove_sr(self):
        if u.rad in self.y_data.unit.decompose().bases:
            if (self.y_data.unit.powers[np.where(np.equal(np.array(self.y_data.unit.decompose().bases), u.rad))[0][0]]
                    == -2):
                self.y_data = self.y_data * u.sr * np.pi
            else:
                raise ValueError('Given units cannot be converted to units required by LIFEsim')

    # TODO: Implement that the user can select which column is x and which is y
    # TODO: Check whether the vector is [N,2] or [2,N] and transpose accordingly
    # TODO: Ignore any non-numeric values during import
    def import_spectrum(self):
        spec = np.loadtxt(self.pathtotxt).T
        self.x_data = (spec[0] * u.Unit(self.x_string)).decompose()
        self.y_data = (spec[1] * u.Unit(self.y_string)).decompose()

        self.x_raw = spec[0] * u.Unit(self.x_string)
        self.y_raw = spec[1] * u.Unit(self.y_string)

    def update_x(self,
                 x_string: str):
        self.x_string = x_string
        if (u.Unit(x_string).decompose().bases != [u.s]) and (u.Unit(x_string).decompose().bases != [u.m]):
            # TODO: Instead of error raising prob. better to set a flag and warn user in GUI
            raise ValueError('Unit not accepted')

    # TODO: Test if the input is valid?
    def update_y(self,
                 y_string: str):
        self.y_string = y_string

    def update_pathtotext(self,
                          pathtotext: str):
        self.pathtotxt = pathtotext

    def do_import(self,
                  pathtotext: str,
                  x_string: str,
                  y_string: str,
                  radius_p_spectrum: Union[float, type(None)],
                  radius_p_target: float,
                  distance_s_spectrum: Union[float, type(None)],
                  distance_s_target: float,
                  integration_time: float):
        self.radius_p_spectrum = radius_p_spectrum
        self.radius_p_target = radius_p_target
        self.distance_s_target = distance_s_target
        self.distance_s_spectrum = distance_s_spectrum
        self.integration_time = integration_time

        self.update_pathtotext(pathtotext)
        self.update_x(x_string)
        self.update_y(y_string)
        self.import_spectrum()
        self.to_wavelength()
        self.to_counts()
        self.to_per_second()
        self.remove_sr()
        self.adjust_to_target()

        if (self.y_data.unit != (u.ph / u.m ** 3 / u.s)) or (self.x_data.unit != u.m):
            raise ValueError('Given units cannot be converted to units required by LIFEsim')

# all settings
# in_string = 'erg cm-2 s-1 Hz-1'
# x_string = 'Hz'
# radius_p = None
# observation_time = None
# distance_s = None
#
#
# target = u.photon / u.m / u.m**2 / u.s
# base = u.Unit(in_string)
# x_unit = u.Unit(x_string)


if __name__ == '__main__':
    # res = np.loadtxt('/home/felix/Documents/MA/LIFEsim/docs/_static/example_spectrum.txt')
    # x = res[:, 0]
    # y = res[:, 1]
    #
    # x_unit = x * u.micron
    # y_unit = y * u.photon / u.meter**2 / u.micron / u.s
    #
    # target = u.photon / u.meter**2 / u.Hz / u.s
    #
    # x_new = x_unit.to(u.Hz, equivalencies=u.spectral())
    # y_new = y_unit.to(target, equivalencies=u.spectral_density(x_new))
    #
    # int_new = np.zeros(y_new.shape[0]) * u.photon / u.meter**2 / u.s
    # int_old = np.zeros(y_new.shape[0]) * u.photon / u.meter**2 / u.s
    # for i in range(y_new.shape[0]-1):
    #     int_new[i] = y_new[i] * abs(x_new[i] - x_new[i+1])
    #     int_old[i] = y_unit[i] * abs(x_unit[i] - x_unit[i + 1])
    # a=1

    si = SpectrumImporter()
    si.update_x('mHz')
    si.update_y('erg cm-2 s-1 Hz-1')
    si.update_pathtotext(pathtotext='/home/felix/Documents/MA/LIFEsim/docs/_static/example_spectrum.txt')
    si.import_spectrum()
    si.to_wavelength()
    si.to_counts()

    a=1