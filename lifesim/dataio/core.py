from abc import ABC


class Module(ABC):
    def __init__(self):
        self.name = None
        self.data = {}
        self.sockets = {}
        self.catalog = None

    def add_socket(self):
        pass

    def run_socket(self):
        pass

    def update_socket(self):
        pass

    def connect_module(self):
        pass

    def disconnect_module(self):
        pass


class Bus(object):
    def __init__(self):
        self.modules = {}

    @classmethod
    def init_from_config(cls):
        pass

    def save_to_config(self):
        pass

    def connect(self):
        pass

    def disconnect(self):
        pass
