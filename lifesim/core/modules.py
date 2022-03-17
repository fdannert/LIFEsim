import abc
import numpy as np
from typing import Union

from lifesim.core.core import Module


class InstrumentModule(Module):
    """
    Module to handle the central simulation of on-sky sources and the nulling-interferometric
    combination of their emissions.
    """

    def __init__(self,
                 name: str):
        super().__init__(name=name)

        # the instrument module needs control over the module responsible for creating the
        # transmission maps. Only one such module should exist per instrument or simulation.
        self.add_socket(s_name='transmission',
                        s_type=TransmissionModule,
                        s_number=1)

        # the photon noise from the individual sources in calculated in separate modules. The
        # instrument module needs access to their data and control over when they execute
        # TODO: not optimal to hardcode the s_number like this. Maybe find a way to draw it from
        #   options
        self.add_socket(s_name='photon_noise_star',
                        s_type=PhotonNoiseStarModule,
                        s_number=5)
        self.add_socket(s_name='photon_noise_universe',
                        s_type=PhotonNoiseUniverseModule,
                        s_number=5)

    @abc.abstractmethod
    def get_snr(self):
        """
        The get_snr function should take the existing P-Pop catalog and add the signal-to-noise
        ratio after one hour of integration to this catalog.
        """

        pass

    @abc.abstractmethod
    def get_spectrum(self,
                     temp_s: float,
                     radius_s: float,
                     distance_s: float,
                     lat_s: float,
                     z: float,
                     angsep: float,
                     flux_planet_spectrum: list,
                     integration_time: float):
        """
        The get_spectrum function should calculate the signal-to-noise ratio of a single planet
        after one hour of integration time.

        Parameters
        ----------
        temp_s : float
            Temperature of the observed star in [Kelvin]
        radius_s : float
            Radius of the observed star in [sun radii]
        distance_s : float
            Distance between the observed star and the LIFE array in [pc]
        lat_s : float
            Ecliptic latitude of the observed star in [rad]
        z : float
            Zodi-level in the observed system in [zodis]
        angsep : float
            Angular separation between the observed star and the observed exoplanet in [arcsec]
        flux_planet_spectrum : np.ndarray
            Array containing the spectrum of the planet. `flux_planet_spectrum[0]` should contain
            wavelength values in [m] and `flux_planet_spectrum[1]` should contain flux values in
            [ph m-3 s-1]
        integration_time : float
            Time that the LIFE array spends for observing the observed planet in [s]
        """

        pass


class PhotonNoiseStarModule(Module):
    """
    Module for simulating astrophysical sources and their photon shot noise contribution specific
    to a single star to the interferometric measurement.
    """
    @abc.abstractmethod
    def noise(self,
              index: Union[int, type(None)]):
        """
        Calculates the photon shot noise contribution.

        Parameters
        ----------
        index : Union[int, type(None)]
            If an integer is given, the photon noise of the planet corresponding to the respective
            interger row position in the catalog is given. If `None` is given, the photon noise is
            calculated for the parameters found in `bus.data.single`.
        """
        pass


class PhotonNoiseUniverseModule(Module):
    """
    Module for simulating astrophysical sources and their photon shot noise contribution specific
    to a single universe to the interferometric measurement.
    """
    @abc.abstractmethod
    def noise(self,
              index: Union[int, type(None)]):
        """
        Calculates the photon shot noise contribution.

        Parameters
        ----------
        index : Union[int, type(None)]
            If an integer is given, the photon noise of the planet corresponding to the respective
            interger row position in the catalog is given. If `None` is given, the photon noise is
            calculated for the parameters found in `bus.data.single`.
        """
        pass


class TransmissionModule(Module):
    """
    Module for calculating transmission maps of the nulling-interferometer.
    """
    @abc.abstractmethod
    def transmission_map(self,
                         map_selection: list):
        """
        Calculates the transmission map.

        Parameters
        ----------
        map_selection : list
            The list specifies which of the transmission maps should be calculated. In the current
            implementation, a double Bracewell interferometer is simulated. Therefore, one can
            choose to return the transmission maps of the two single Bracewell arms (TM1-4) or the
            chopped transmission map (TM Chop).
        """
        pass

    @abc.abstractmethod
    def transmission_efficiency(self,
                                index: Union[int, type(None)]):
        """
        Calculates the transmission efficiency on a single specified angular distance from the
        central star.

        Parameters
        ----------
        index : Union[int, type(None)]
            If an integer is given, the transmission efficiency of the planet corresponding to the
            respective interger row position in the catalog is given. If `None` is given, the
            transmission efficiency is calculated for the parameters found in `bus.data.single`.
        """
        pass

class OptimizationModule(Module):
    """
    Module for distributing the observing time during the search phase.
    """
    def __init__(self,
                 name: str):
        super().__init__(name=name)
        # needs control over the instrument module
        self.add_socket(s_name='instrument',
                        s_type=InstrumentModule,
                        s_number=1)
        # needs control over the transmission map
        self.add_socket(s_name='transmission',
                        s_type=TransmissionModule,
                        s_number=1)
        # needs control over the slope module
        self.add_socket(s_name='slope',
                        s_type=SlopeModule,
                        s_number=1)

    @abc.abstractmethod
    def find_phase(self):
        """
        This function calculates the optimal phase for observation. This is need for the simulation
        of multi-epoch observations.
        """
        pass


class SlopeModule(Module):
    """
    The slope module can find optimal observing time distributions using the aHGS algorithm in
    connection with the algorithm developed by Stark+.
    """
    @abc.abstractmethod
    def distribute_time(self):
        """
        Run the aHGS algorithm.
        """
        pass
