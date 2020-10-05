# todo comment

from typing import Tuple


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
        self.plugs = {}
        self.m_type = 'plugin'


class PrimaryModule(Module):
    def __init__(self,
                 name: str):
        super().__init__(name=name)
        self.m_type = 'primary'

    def connect_plugin(self,
                      plugin: Module):
        self.plugs[plugin.name] = plugin


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



if __name__ == '__main__':
    b = Bus()
    pm = PrimaryModule('Prim')
    plugin = Module('Plug')
    b.add_module(pm)
    b.add_module(plugin)
    b.connect(('Prim', 'Plug'))






