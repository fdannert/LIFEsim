# todo Documentation

from typing import Tuple
import warnings


class Cable(object):
    """Store connection information between modules in the bus

    """
    def __init__(self,
                 name_primary: str,
                 name_plugin: str):
        self.name_primary = name_primary
        self.name_plugin = name_plugin


class Module(object):
    def __init__(self,
                 name: str):
        self.name = name
        self.m_type = 'plugin'
        self.f_type = None  # functionality type of plugin
        self.data = None


class PrimaryModule(Module):
    def __init__(self,
                 name: str):
        super().__init__(name=name)
        self.m_type = 'primary'
        self.sockets = {}
        self.socket_f_type = {}
        self.socket_data = {}

    def connect_plugin(self,
                       plugin: Module):
        if plugin.f_type not in self.socket_f_type.keys():
            raise ValueError('Primary module does not accept plugin of this type')
        connected = False
        for _, name in enumerate(self.socket_f_type[plugin.f_type]):
            if self.sockets[name] is None:
                self.sockets[name] = plugin
                plugin.data = self.socket_data[name]
                connected = True
                break
        if not connected:
            raise AttributeError('All available sockets are already occupied')

    def disconnect_plugin(self,
                          plugin: Module):
        if plugin not in self.sockets.values():
            raise ValueError('Plugin ' + str(plugin.name) + ' not found')
        for socket in self.sockets.items():
            if socket[1] is plugin:
                plugin.data = None
                self.sockets[socket[0]] = None

    def add_socket(self,
                   name: str,
                   f_type: str,
                   data: dict):
        if name in self.sockets.keys():
            raise ValueError('This socket already exists')
        self.sockets[name] = None
        if f_type in self.socket_f_type.keys():
            self.socket_f_type[f_type].append(name)
        else:
            self.socket_f_type[f_type] = [name]
        self.socket_data[name] = data

    def run_socket(self,
                   name: str,
                   **kwargs):
        if name not in self.sockets.keys():
            raise ValueError('This socket does not exist')
        elif self.sockets[name] is None:
            warnings.warn('No module run since no plugin connected to socket ' + name)
        else:
            if len(kwargs) == 0:
                self.sockets[name].run()
            else:
                self.sockets[name].run(kwargs)

    def update_socket(self,
                      name: str,
                      data: dict):
        if name not in self.sockets.keys():
            raise ValueError('This socket does not exist')
        self.socket_data[name].update(data)


class Bus(object):
    def __init__(self):
        self.modules = {}
        self.cables = []

    def add_module(self,
                   module: Module):
        self.modules[module.name] = module

    def connect(self,
                module_names: Tuple[str, str]):
        if (module_names[0] not in self.modules) or (module_names[1] not in self.modules):
            raise ValueError('At least one of the specified modules does not exist')

        types = [self.modules[m].m_type for m in module_names]

        if not (set(types) == {'primary', 'plugin'}):
            raise ValueError('Can only connect a plugin to a master module')

        name_primary = module_names[types.index('primary')]
        name_plugin = module_names[(types.index('primary') + 1) % 2]

        # connect the plugin to the primary module
        self.modules[name_primary].connect_plugin(self.modules[name_plugin])

        self.cables.append(Cable(name_primary=name_primary,
                                 name_plugin=name_plugin))

    def disconnect(self,
                   module_names: Tuple[str, str]):
        if (module_names[0] not in self.modules) or (module_names[1] not in self.modules):
            raise ValueError('At least one of the specified modules does not exist')

        types = [self.modules[m].m_type for m in module_names]
        name_primary = module_names[types.index('primary')]
        name_plugin = module_names[(types.index('primary') + 1) % 2]

        cable_found = False
        for cab in self.cables:
            if (cab.name_primary == name_primary) \
                    and (cab.name_plugin == name_plugin):
                cable_found = True
                self.cables.remove(cab)
                self.modules[name_primary].disconnect_plugin(self.modules[name_plugin])

        if not cable_found:
            warnings.warn('The specified connection does not exist or was already disconnected')





