from __future__ import annotations
import abc
from abc import ABC

from lifesim.core.data import Data


class Module(ABC):
    def __init__(self,
                 name: str):
        self.name = name
        self.data = None
        self.sockets = {}
        self.catalog = None

    def add_socket(self,
                   s_name: str,
                   s_type: type,
                   s_number: int = 1):

        # check if socket already exists
        self.socket_exists(s_name=s_name,
                           should_exist=False)

        # create socket with s_name as key
        self.sockets[s_name] = {}

        # add the maximum number of connected modules and the expected types to the socket dict
        self.sockets[s_name]['s_type'] = s_type
        self.sockets[s_name]['s_number'] = s_number

        # add the socket data
        self.sockets[s_name]['data'] = {}

        # create list of modules
        self.sockets[s_name]['modules'] = [None] * s_number

    def run_socket(self,
                   method: str,
                   s_name: str,
                   **kwargs):

        self.socket_exists(s_name=s_name,
                           should_exist=True)

        if self.sockets[s_name]['s_number'] == 1:
            return getattr(self.sockets[s_name]['modules'][0], method)(**kwargs)
        else:
            return [getattr(module, method)(**kwargs)
                    for module in filter(None, self.sockets[s_name]['modules'])]

    def update_socket(self,
                      s_name: str,
                      s_data: dict):

        # check if the socket exists
        self.socket_exists(s_name=s_name,
                           should_exist=True)

        # set the new/updated data to the socket data
        self.sockets[s_name]['data'].update(s_data)

    def connect_module(self,
                       module: Module):

        s_name = None
        for key in self.sockets.keys():
            # if type(module) == self.sockets[key]['s_type']:
            if isinstance(module, self.sockets[key]['s_type']):
                s_name = key
                break
        if s_name is None:
            raise ValueError('Cannot find socket that accepts module of type ' + str(type(module)))

        for i in range(self.sockets[s_name]['s_number']):
            if self.sockets[s_name]['modules'][i] is None:
                self.sockets[s_name]['modules'][i] = module
                break
            if i == self.sockets[s_name]['s_number'] - 1:
                raise ValueError('All available socket spots for this module type are already '
                                 'occupied')

    def disconnect_module(self,
                          module: Module):

        # search the sockets for the given module
        s_index = None
        for key in self.sockets.keys():
            if module in self.sockets[key]['modules']:

                # if the module is found, remove it from the socket
                s_index = self.sockets[key]['modules'].index(module)
                self.sockets[key]['modules'][s_index] = None
                break

        # if the module is not found raise an error
        if s_index is None:
            raise ValueError('The module ' + module.name + ' is not connected to the module '
                             + self.name)

    def socket_exists(self,
                      s_name: str,
                      should_exist: bool):
        if should_exist and (s_name not in self.sockets.keys()):
            raise ValueError('The socket ' + s_name + ' does not exist')
        elif (not should_exist) and (s_name in self.sockets.keys()):
            raise ValueError('A socket with the name ' + s_name + ' already exists')


class Bus(object):
    def __init__(self):
        self.modules = {}
        self.data = Data()

    @classmethod
    def init_from_config(cls):
        pass

    def save_to_config(self):
        pass

    def connect(self):
        pass

    def disconnect(self):
        pass

    def add_module(self,
                   module: Module):
        """
        Adds a given module to the bus

        Parameters
        ----------
        module : lifesim.Module
            Module that will be added to the bus

        Raises
        ------
        ValueError
            If the name of the module is already in use
        """

        # add module to the list of modules
        if module.name in self.modules.keys():
            raise ValueError('The module name ' + module.name + ' is already in use in this bus')
        self.modules[module.name] = module

        # point module to the data storage class
        module.data = self.data
