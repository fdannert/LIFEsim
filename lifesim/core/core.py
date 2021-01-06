from __future__ import annotations
import abc
from abc import ABC
from typing import Tuple
import warnings

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
        self.connections = []

    @classmethod
    def init_from_config(cls):
        pass

    def save_to_config(self):
        pass

    def connect(self,
                module_names: Tuple[str, str]):
        """
        Connects a 'plugin' to a 'primary' module

        Parameters
        ----------
        module_names : Tuple[str, str]
            A tuple containing the names of the 'primary' and 'plugin' modules which will be
            connected. Must contain one 'primary' and one 'plugin' module

        Raises
        ------
        ValueError
            If at least one of the specified modules does not exists. If two modules of same module
            type are given
        """

        # check if modules exist
        if module_names[0] not in self.modules.keys():
            raise ValueError('The module ' + module_names[0] + ' does not exist')
        elif module_names[1] not in self.modules.keys():
            raise ValueError('The module ' + module_names[1] + ' does not exist')

        i = 0
        try:
            self.modules[module_names[0]].connect_module(self.modules[module_names[1]])
            self.connections.append({module_names[0], module_names[1]})
        except ValueError:
            i += 1
        try:
            self.modules[module_names[1]].connect_module(self.modules[module_names[0]])
            self.connections.append({module_names[1], module_names[0]})
        except ValueError:
            i += 1
        if i == 2:
            raise ValueError('The modules '
                             + module_names[0] + ' and '
                             + module_names[1] + ' could not be connected')

    def disconnect(self,
                   module_names: Tuple[str, str]):
        if set(module_names) not in self.connections:
            warnings.warn('The modules '
                          + module_names[0] + ' and '
                          + module_names[1] + ' are not connected')

        i = 0
        try:
            self.modules[module_names[0]].disconnect_module(self.modules[module_names[1]])
            self.connections.remove({module_names[0], module_names[1]})
        except ValueError:
            i += 1
        try:
            self.modules[module_names[1]].disconnect_module(self.modules[module_names[0]])
            self.connections.remove({module_names[1], module_names[0]})
        except ValueError:
            i += 1
        if i == 2:
            warnings.warn('The modules '
                          + module_names[0] + ' and '
                          + module_names[1] + ' could not be disconnected')

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

    # def connect(self,
    #             module_names: Tuple[str, str]):
    #     """
    #     Connects a 'plugin' to a 'primary' module
    #
    #     Parameters
    #     ----------
    #     module_names : Tuple[str, str]
    #         A tuple containing the names of the 'primary' and 'plugin' modules which will be
    #         connected. Must contain one 'primary' and one 'plugin' module
    #
    #     Raises
    #     ------
    #     ValueError
    #         If at least one of the specified modules does not exists. If two modules of same module
    #         type are given
    #     """
    #
    #     # check if modules exist
    #     if (module_names[0] not in self.modules) or (module_names[1] not in self.modules):
    #         raise ValueError('At least one of the specified modules does not exist')
    #
    #     types = [self.modules[m].m_type for m in module_names]
    #
    #     # check if modules are of different type
    #     if not (set(types) == {'primary', 'plugin'}):
    #         raise ValueError('Can only connect a plugin to a master module')
    #
    #     # identify which module is primary and which is plugin
    #     name_primary = module_names[types.index('primary')]
    #     name_plugin = module_names[(types.index('primary') + 1) % 2]
    #
    #     # connect the plugin to the primary module
    #     self.modules[name_primary].connect_plugin(self.modules[name_plugin])
    #
    #     # create a cable for the conntection
    #     self.cables.append(Cable(name_primary=name_primary,
    #                              name_plugin=name_plugin))
    #
    # def disconnect(self,
    #                module_names: Tuple[str, str]):
    #     """
    #     Disconnects a 'plugin' and a 'primary' module
    #
    #     Parameters
    #     ----------
    #     module_names : Tuple[str, str]
    #         A tuple containing the names of the 'primary' and 'plugin' modules which will be
    #         disconnected. Must contain one 'primary' and one 'plugin' module. Raises a warning if
    #         the specified modules are not connected.
    #
    #     Raises
    #     ------
    #     ValueError
    #         If at least one of the specified modules does not exists.
    #     """
    #
    #     # check if the modules exist
    #     if (module_names[0] not in self.modules) or (module_names[1] not in self.modules):
    #         raise ValueError('At least one of the specified modules does not exist')
    #
    #     # identify which module is primary and which is plugin
    #     types = [self.modules[m].m_type for m in module_names]
    #     name_primary = module_names[types.index('primary')]
    #     name_plugin = module_names[(types.index('primary') + 1) % 2]
    #
    #     # find the appropriate cable
    #     cable_found = False
    #     for cab in self.cables:
    #         if (cab.name_primary == name_primary) \
    #                 and (cab.name_plugin == name_plugin):
    #             cable_found = True
    #             # remove cable and disconnect plugin
    #             self.cables.remove(cab)
    #             self.modules[name_primary].disconnect_plugin(self.modules[name_plugin])
    #
    #     # give a warning if the connection does not exist
    #     if not cable_found:
    #         warnings.warn('The specified connection does not exist or was already disconnected')
