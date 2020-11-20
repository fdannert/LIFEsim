from typing import Tuple
import warnings


class Cable(object):
    """
    Stores connection information between modules in the bus

    Attributes
    ----------
    name_primary : str
        Name of the 'primary' module to which the cable is connected
    name_plugin : str
        Name of the 'plugin' module to which the cable is connected
    """
    def __init__(self,
                 name_primary: str,
                 name_plugin: str):
        """
        Parameters
        ----------
        name_primary : str
            Name of the 'primary' module to which the cable is connected
        name_plugin : str
            Name of the 'plugin' module to which the cable is connected
        """
        self.name_primary = name_primary
        self.name_plugin = name_plugin


class Module(object):
    """
    Modules represent computational cores on the bus. Individual modules should be used to
    perform calculations regarding a single model or theme. Here, modules are initialized as
    'plugins' which should be linked to 'primary' modules.

    Attributes
    ----------
    name : str
        The unique name of the module
    m_type : str
        The module type differentiating between a 'primary' or 'plugin' module
    f_type : str
        The function type specifies which tasks, given by another primary module, a plugin
        module can perform.
    data : object
        A common collection place for the modules data
    """
    def __init__(self,
                 name: str):
        """
        Parameters
        ----------
        name : str
        The unique name of the module
        """
        self.name = name
        self.m_type = 'plugin'
        self.f_type = None  # function type of plugin
        self.data = None


class PrimaryModule(Module):
    """
    Modules represent computational cores on the bus. Individual modules should be used to
    perform calculations regarding a single model or theme. 'Primary' modules are used to handle
    tasks performed by connected 'plugin' modules.

    Attributes
    ----------
    name : str
        The unique name of the module
    m_type : str
        The module type differentiating between a 'primary' or 'plugin' module
    sockets : dict
        The sockets provide points to which 'plugin' modules can connect
    socket_f_type : dict
        The function type of the sockets specify what tasks the connected 'plugin' module should be
         able to perform
    socket_data : dict
        A repository used to transfer data between the 'primary' and 'plugin' modules
    """
    def __init__(self,
                 name: str):
        """
        Parameters
        ----------
        name : str
            The unique name of the module
        """
        super().__init__(name=name)
        self.m_type = 'primary'
        self.sockets = {}
        self.socket_f_type = {}
        self.socket_data = {}

    def connect_plugin(self,
                       plugin: Module):
        """
        Tries to connect a given 'plugin' module to the primary module. The 'plugin' module is
        connected to an empty socket with a matching function type.

        Parameters
        ----------
        plugin : lifesim.Module
            'Plugin' module which is connected to the 'primary' module

        Raises
        ------
        ValueError
            If there is no socket with a function type matching the function type of the given
            'plugin' module
        AttributeError
            If all available sockets of matching function type are already occupied
        """

        # check if primary does accept the f_type of the plugin
        if plugin.f_type not in self.socket_f_type.keys():
            raise ValueError('Primary module does not accept plugin of this type')

        # try a connection with all sockets of the correct f_type
        connected = False
        for _, name in enumerate(self.socket_f_type[plugin.f_type]):
            if self.sockets[name] is None:
                self.sockets[name] = plugin
                plugin.data = self.socket_data[name]
                connected = True
                break

        # raise error if no connection could be made
        if not connected:
            raise AttributeError('All available sockets are already occupied')

    def disconnect_plugin(self,
                          plugin: Module):
        """
        Tries to disconnect a given 'plugin' module from the primary module.

        Parameters
        ----------
        plugin : lifesim.Module
            'Plugin' module which should be disconnected from the 'primary' module

        Raises
        ------
        ValueError
            If no disconnection can be done because there is not fitting 'plugin' module connected
            to the 'primary' module
        """

        # check if specified plugin is connected to the primary
        if plugin not in self.sockets.values():
            raise ValueError('Plugin ' + str(plugin.name) + ' not found')

        # disconnect the plugin and clear the plugin data
        for socket in self.sockets.items():
            if socket[1] is plugin:
                plugin.data = None
                self.sockets[socket[0]] = None

    def add_socket(self,
                   name: str,
                   f_type: str,
                   data: dict):
        """
        Adds a socket of a given function type to the 'primary' module

        Parameters
        ----------
        name : str
            Name of the socket
        f_type : str
            Function type of the socket. The function type will define which tasks the 'plugin'
            module connected to this socket is expected to perform. This means that sockets of a
            ceratin function type will only accept a connection by a 'plugin' module with the same
            function type
        data : dict
            Data which needs to be transferred to and from the 'plugin' module

        Raises
        ------
        ValueError
            If a socket of the same name already exists for this 'primary' module
        """

        # check if socket already exists
        if name in self.sockets.keys():
            raise ValueError('This socket already exists')

        # initialize socket
        self.sockets[name] = None

        # add socket f_type to the list of f_types
        if f_type in self.socket_f_type.keys():
            self.socket_f_type[f_type].append(name)
        else:
            self.socket_f_type[f_type] = [name]

        # store data to the socket
        self.socket_data[name] = data

    def run_socket(self,
                   name: str,
                   **kwargs):
        """
        Executes the connected 'plugin' module. Produces a warning if no 'plugin' module is
        connected to the specified socket

        Parameters
        ----------
        name : str
            Name of the socket which is executed
        kwargs
            The keyword arguments will be forwarded to the 'plugin' module

        Raises
        ------
        ValueError
            If the specified socket does not exist
        """

        # check if socket exists
        if name not in self.sockets.keys():
            raise ValueError('This socket does not exist')

        # check if a plugin is connected to the socket
        elif self.sockets[name] is None:
            warnings.warn('No module run since no plugin connected to socket ' + name)

        # run with or without kwargs
        else:
            if len(kwargs) == 0:
                self.sockets[name].run()
            else:
                self.sockets[name].run(kwargs)

    def update_socket(self,
                      name: str,
                      data: dict):
        """
        Updates or adds data given to the specified socket

        Parameters
        ----------
        name : str
            Name of the socket which is executed
        data : dict
            A dictionary containing the data that will be updated or added

        Raises
        ------
        ValueError
            If the specified socket does not exist
        """

        # check if the socket exists
        if name not in self.sockets.keys():
            raise ValueError('This socket does not exist')

        # set the new/updated data to the socket data
        self.socket_data[name].update(data)


class Bus(object):
    """
    The bus is the underlying structure to which both 'primary' and 'plugin' modules are connected.
    It also stores the connections between the individual 'primary' and 'plugin' modules

    Attributes
    ----------
    modules : dict
        A dictionary containing all connected modules with the respective module names
    cables : list
        A list containing all cables between modules
    """
    def __init__(self):
        """

        """
        self.modules = {}
        self.cables = []

    def add_module(self,
                   module: Module):
        """
        Adds a given 'primary' or 'plugin' module to the bus

        Parameters
        ----------
        module: lifesim.Module
            Module that will be added to the bus
        """
        self.modules[module.name] = module

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
        if (module_names[0] not in self.modules) or (module_names[1] not in self.modules):
            raise ValueError('At least one of the specified modules does not exist')

        types = [self.modules[m].m_type for m in module_names]

        # check if modules are of different type
        if not (set(types) == {'primary', 'plugin'}):
            raise ValueError('Can only connect a plugin to a master module')

        # identify which module is primary and which is plugin
        name_primary = module_names[types.index('primary')]
        name_plugin = module_names[(types.index('primary') + 1) % 2]

        # connect the plugin to the primary module
        self.modules[name_primary].connect_plugin(self.modules[name_plugin])

        # create a cable for the conntection
        self.cables.append(Cable(name_primary=name_primary,
                                 name_plugin=name_plugin))

    def disconnect(self,
                   module_names: Tuple[str, str]):
        """
        Disconnects a 'plugin' and a 'primary' module

        Parameters
        ----------
        module_names : Tuple[str, str]
            A tuple containing the names of the 'primary' and 'plugin' modules which will be
            disconnected. Must contain one 'primary' and one 'plugin' module. Raises a warning if
            the specified modules are not connected.

        Raises
        ------
        ValueError
            If at least one of the specified modules does not exists.
        """

        # check if the modules exist
        if (module_names[0] not in self.modules) or (module_names[1] not in self.modules):
            raise ValueError('At least one of the specified modules does not exist')

        # identify which module is primary and which is plugin
        types = [self.modules[m].m_type for m in module_names]
        name_primary = module_names[types.index('primary')]
        name_plugin = module_names[(types.index('primary') + 1) % 2]

        # find the appropriate cable
        cable_found = False
        for cab in self.cables:
            if (cab.name_primary == name_primary) \
                    and (cab.name_plugin == name_plugin):
                cable_found = True
                # remove cable and disconnect plugin
                self.cables.remove(cab)
                self.modules[name_primary].disconnect_plugin(self.modules[name_plugin])

        # give a warning if the connection does not exist
        if not cable_found:
            warnings.warn('The specified connection does not exist or was already disconnected')
