import abc
import numpy as np
from typing import Union

from lifesim.core.core import Module


class InstrumentModule(Module):

    def __init__(self,
                 name: str):
        super().__init__(name=name)
        self.add_socket(s_name='transmission',
                        s_type=TransmissionModule,
                        s_number=1)
        # TODO: not optimal to hardcode the s_number like this. Maybe find a way to draw it from
        #   options
        self.add_socket(s_name='photon_noise',
                        s_type=PhotonNoiseModule,
                        s_number=5)

    @abc.abstractmethod
    def get_snr(self):
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
        pass


class PhotonNoiseModule(Module):

    @abc.abstractmethod
    def noise(self,
              index: Union[int, type(None)]):
        pass


class TransmissionModule(Module):

    @abc.abstractmethod
    def transmission_map(self,
                         map_selection: list):
        pass

    @abc.abstractmethod
    def transmission_efficiency(self,
                                index: Union[int, type(None)]):
        pass

class OptimizationModule(Module):

    def __init__(self,
                 name: str):
        super().__init__(name=name)
        self.add_socket(s_name='instrument',
                        s_type=InstrumentModule,
                        s_number=1)
        self.add_socket(s_name='transmission',
                        s_type=TransmissionModule,
                        s_number=1)
        self.add_socket(s_name='slope',
                        s_type=SlopeModule,
                        s_number=1)

    @abc.abstractmethod
    def find_phase(self):
        pass


class SlopeModule(Module):

    @abc.abstractmethod
    def distribute_time(self):
        pass
