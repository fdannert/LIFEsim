import abc
from lifesim.dataio.core import Module


class InstrumentModule(Module):

    def __init__(self,
                 name: str):
        super().__init__(name=name)
        self.add_socket(s_name='transmission',
                        s_type=type(TransmissionModule),
                        s_number=1)
        # TODO: not optimal to hardcode the s_number like this. Maybe find a way to draw it from
        #   options
        self.add_socket(s_name='photon_noise',
                        s_type=type(PhotonNoiseModule),
                        s_number=5)

    @abc.abstractmethod
    def get_snr(self):
        pass

    @abc.abstractmethod
    def get_spectrum(self):
        pass


class PhotonNoiseModule(Module):

    @property
    @abc.abstractmethod
    def noise(self):
        pass


class TransmissionModule(Module):

    @property
    @abc.abstractmethod
    def transmission_map(self):
        pass

    @property
    @abc.abstractmethod
    def transmission_efficiency(self):
        pass
