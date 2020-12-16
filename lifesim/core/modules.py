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
              image_size: int,
              l_sun: float,
              distance_s: float,
              mas_pix: np.ndarray,
              rad_pix: np.ndarray,
              z: float,
              telescope_area: float,
              radius_map: np.ndarray,
              wl_bins: np.ndarray,
              wl_bin_edges: np.ndarray,
              wl_bin_widths: np.ndarray,
              t_map: np.ndarray,
              hfov: np.ndarray,
              lz_model: str,
              lat_s: float,
              radius_s: float,
              temp_s: float,
              bl: float,
              ratio: float,):
        pass


class TransmissionModule(Module):

    @abc.abstractmethod
    def transmission_map(self,
                         wl_bins: np.ndarray,
                         hfov: np.ndarray,
                         image_size: Union[int, type(None)],
                         bl: float,
                         map_selection: list,
                         ratio: float,
                         direct_mode: bool = False,
                         d_alpha: np.ndarray = None,
                         d_beta: np.ndarray = None):
        pass

    @abc.abstractmethod
    def transmission_efficiency(self,
                                bl: float,
                                wl_bins: np.ndarray,
                                angsep: np.ndarray,
                                ratio: float):
        pass
