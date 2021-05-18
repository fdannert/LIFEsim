from __future__ import annotations
import abc
from abc import ABC
from typing import Tuple
import warnings

from lifesim.core.data import Data


class Module(ABC):
    """
    A `Module` is the most basic component of LIFEsim. A single `Module` should be used to
    encapsulate a physical, (partially) self-contained model. Communication between two or more
    Modules is done through sockets. See the LIFEsim ULM diagram or the inheritance and
    implementation in :file:`lifesim.core.modules.py` and :file:`lifesim.instrument.*`

    Attributes
    ----------
    name : str
        Name of the module. The name is used as an identifier for different instances of modules
        and should therefore be unique to a single instance of the `Bus` class.
    data : dict
        Deprecated. Store data in the respective location in the `Bus` class (`bus.data`)
    sockets : dict
        The sockets dictionary stores the properties of the sockets under the key equal to the
        socket name

    Examples
    --------
    When any module is created within LIFEsim, it should inherit the ``Module`` class and define
    stub methods as follows

    >>> from lifesim.core.core import Module
    >>> class NewModule(Module):
    ...     def __init__(self,
    ...                  name: str):
    ...         ().__init__(name=name)
    ...     @abc.abstractmethod
    ...     def new_method(self):
    ...         pass
    """

    def __init__(self,
                 name: str):
        self.name = name
        self.data = None
        self.sockets = {}

    def add_socket(self,
                   s_name: str,
                   s_type: type,
                   s_number: int = 1):
        """
        Add a socket to the module.

        Parameters
        ----------
        s_name : str
            The name of the socket. This name is used as the key in the `sockets` dictionary. Must
            be unique within a single module.
        s_type : type
            The socket will only accept connections with modules of the given type.
        s_number : int
            Specifies the number of module instances that can be connected to this socket.

        Raises
        ------
        ValueError
            If the a socket with the specified socket name already exists.

        Examples
        --------
        In this example, a newly coded module needs an interface to the TransmissionModule. A
        suitable socket can be added as follows

        >>> import abc
        >>> from lifesim.core.core import Module
        >>> from lifesim.core.modules import TransmissionModule
        >>> class NewModule(Module):
        ...     def __init__(self,
        ...                  name: str):
        ...         ().__init__(name=name)
        ...
        ...         self.add_socket(s_name='transm',
        ...                         s_type=TransmissionModule,
        ...                         s_number=1)
        ...     @abc.abstractmethod
        ...     def new_method(self):
        ...         pass
        """

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
        """
        A generic command for execution of a specific method in all modules connected to the
        socket. This can be used to specify workflows and hierachy between different modules.

        Parameters
        ----------
        s_name : str
            Specifies the name of the socket that should send the run command.
        method : str
            Specifies the method that should be called in the connected module(s).
        kwargs :
            Can be used to pass any argument to the specified method within the connected
            module(s).

        Returns
        -------
        Passes trough any returns from the specified methods within the connected module(s). If
        multiple modules are connected, the return variables are encapsulated in a list.

        Raises
        ------
        ValueError
            If there exist no sockets with the specified socket names.

        Example
        -------
        Continuing the above example in ``add_socket``, any module added to the socket can be run
        by using

        >>> class New(NewModule):
        ...         def __init__(self,
        ...                      name: str):
        ...             super().__init__(name=name)
        ...
        ...         def new_method(self):
        ...             _, _, result, _, _ = self.run_socket(s_name='transm',
        ...                                                  method='transmission_map',
        ...                                                  map_selection='tm3')
        """

        self.socket_exists(s_name=s_name,
                           should_exist=True)

        if self.sockets[s_name]['s_number'] == 1:
            return getattr(self.sockets[s_name]['modules'][0], method)(**kwargs)
        else:
            # TODO: Attach the name of the module to the returned variables in the list s.t. the
            #       source of the individual returned variables can be identified
            return [getattr(module, method)(**kwargs)
                    for module in filter(None, self.sockets[s_name]['modules'])]

    def update_socket(self,
                      s_name: str,
                      s_data: dict):
        """
        Deprecated. All data should be handled through the bus. Pushes changed data to all modules
        connected to the bus.

        Parameters
        ----------
        s_name : str
            Specifies the name of the socket that should push data.
        s_data : dict
            Dictionary containing the updated data.

        Raises
        ------
        ValueError
            If there exist no sockets with the specified socket names.
        """

        # check if the socket exists
        self.socket_exists(s_name=s_name,
                           should_exist=True)

        # set the new/updated data to the socket data
        self.sockets[s_name]['data'].update(s_data)

    def connect_module(self,
                       module: Module):
        """
        Connect the specified module to a suitable and available (the number of module that can
        connect to a single socket can be specified in the creation of the socket) socket. If
        multiple sockets are suitable and available, the module is connected to only the first such
        socket (in the order as they appear in the socket dictionary).

        Parameters
        ----------
        module : Module
            The module that should be connected to a socket.

        Raises
        ------
        ValueError
            If no socket accepting modules of the same type as the specified module exists.
        ValueError
            If none of the suitable sockets have an spot available.
        """

        # TODO: The described functionality is not quite reached since only the first socket of
        #  acceptable type is used for any connection. In the case that the first socket of correct
        #  type does not have free spots but later sockets of correct type do, this will lead to a
        #  non-connection of the module.
        s_name = None
        # cycle through all sockets
        for key in self.sockets.keys():
            # if the type of the specified module is the same as the type of the socket safe the
            # name of the socket an break
            if isinstance(module, self.sockets[key]['s_type']):
                s_name = key
                break
        if s_name is None:
            # raise error if no sockets of the required type exist
            raise ValueError('Cannot find socket that accepts module of type ' + str(type(module)))

        # try to find an empty spot for the module in the previously identified socket
        for i in range(self.sockets[s_name]['s_number']):
            if self.sockets[s_name]['modules'][i] is None:
                self.sockets[s_name]['modules'][i] = module
                break
            if i == self.sockets[s_name]['s_number'] - 1:
                raise ValueError('All available socket spots for this module type are already '
                                 'occupied')

    def disconnect_module(self,
                          module: Module):
        """
        Disconnect the specified module.

        Parameters
        ----------
        module : Module
            The module that should be disconnected.

        Raises
        ------
        ValueError
            If the module is not connected.
        """

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
        """
        Checks for the existence of a socket of the specified name.

        Parameters
        ----------
        s_name : str
            Name of the socket.
        should_exist : bool
            If set to `true`, the exception is raised if the socket of given name does not exist.
            If set to `false`, the exception is raised if the socket of given name does exist.

        Raises
        ------
        ValueError
            Depending of `should_exist`, the error is raised if a socket of given name does or does
            not exist.
        """
        if should_exist and (s_name not in self.sockets.keys()):
            raise ValueError('The socket ' + s_name + ' does not exist')
        elif (not should_exist) and (s_name in self.sockets.keys()):
            raise ValueError('A socket with the name ' + s_name + ' already exists')


class Bus(object):
    """
    The `Bus` is responsible for handling and interfacing with modules, for creating connections
    between modules, for data storage and for option handling. Every simulation run with LIFEsim
    should be based of a bus.

    Attributes
    ----------
    modules : dict
        A dictionary containing all modules used in a simulation.
    data : Data
        The data and options of a simulation are stored in here in a `Data` class.
    connections : list
        A list keeping track over all connections made between modules via sockets.
    """
    def __init__(self):
        self.modules = {}
        self.data = Data()
        self.connections = []

    # TODO: Implement the use of configuration files.
    @classmethod
    def init_from_config(cls):
        pass

    def save_to_config(self):
        pass

    def connect(self,
                module_names: Tuple[str, str]):
        """
        Establish a connection between two specified modules using available sockets

        Parameters
        ----------
        module_names : Tuple[str, str]
            A tuple containing the names of the two modules which will be connected.

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

        # try to connect the modules in both directions
        i = 0
        try:
            self.modules[module_names[0]].connect_module(self.modules[module_names[1]])

            # note the establishment of the connection in the connection list of the bus
            self.connections.append({module_names[0], module_names[1]})
        except ValueError:
            i += 1
        try:
            self.modules[module_names[1]].connect_module(self.modules[module_names[0]])

            # note the establishment of the connection in the connection list of the bus
            self.connections.append({module_names[1], module_names[0]})
        except ValueError:
            i += 1

        # if a connection is not possible in either direction, raise an error that the connection
        # was unsuccessful
        if i == 2:
            raise ValueError('The modules '
                             + module_names[0] + ' and '
                             + module_names[1] + ' could not be connected')

    def disconnect(self,
                   module_names: Tuple[str, str]):
        """
        Disconnect two previously connected modules.

        Parameters
        ----------
        module_names : Tuple[str, str]
            A tuple containing the names of the two modules which will be disconnected.

        Raises
        ------
        Warning
            If there exists no connection between the two specified modules.
        Warning
            If the connected modules cannot be disconnected.
        """

        # check if the modules are connected
        if set(module_names) not in self.connections:
            warnings.warn('The modules '
                          + module_names[0] + ' and '
                          + module_names[1] + ' are not connected')

        # try to disconnect the two modules both ways
        i = 0
        try:
            self.modules[module_names[0]].disconnect_module(self.modules[module_names[1]])

            # if the disconnect is successful, remove it from the connections list
            self.connections.remove({module_names[0], module_names[1]})
        except ValueError:
            i += 1
        try:
            self.modules[module_names[1]].disconnect_module(self.modules[module_names[0]])

            # if the disconnect is successful, remove it from the connections list
            self.connections.remove({module_names[1], module_names[0]})
        except ValueError:
            i += 1
        if i == 2:
            # raise a warning if the modules cannot be disconnected
            warnings.warn('The modules '
                          + module_names[0] + ' and '
                          + module_names[1] + ' could not be disconnected')

    def add_module(self,
                   module: Module):
        """
        Adds a given module to the bus

        Parameters
        ----------
        module : Module
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
