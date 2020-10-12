import warnings

import numpy as np


class Options(object):
    def __init__(self):
        self.diameter = None
        self.quantum_eff = None
        self.throughput = None
        self.wl_limits = None
        self.spec_res = None
        self.baseline = None
        self.image_size = None
        self.wl_optimal = None

        self.array = {'diameter': 0.,
                      'quantum_eff': 0.,
                      'throughput': 0.,
                      'wl_min': 0.,
                      'wl_max': 0.,
                      'spec_res': 0,
                      'baseline': 0.,
                      'bl_min': 0.,
                      'bl_max': 0.}

        self.other = {'image_size': 0,
                      'wl_optimal': 0.}

    def set_scenario(self,
                     case: str):

        self.array['quantum_eff'] = 0.7
        self.array['throughput'] = 0.05
        self.array['spec_res'] = 20.
        self.array['baseline'] = 20.
        self.array['bl_min'] = 10.
        self.array['bl_max'] = 100.

        self.other['image_size'] = 512
        self.other['wl_optimal'] = 15

        if case == 'baseline':
            self.array['diameter'] = 2.
            self.array['wl_min'] = 4.
            self.array['wl_max'] = 18.5

        elif case == 'pessimistic':
            self.array['diameter'] = 1.
            self.array['wl_min'] = 6.
            self.array['wl_max'] = 17.

        elif case == 'optimistic':
            self.array['diameter'] = 3.5
            self.array['wl_min'] = 3.
            self.array['wl_max'] = 20.

        else:
            warnings.warn('Option case not recognised, no options set')

    def set_manual(self,
                   diameter: float,
                   quantum_eff: float,
                   throughput: float,
                   wl_limits: np.ndarray,
                   spec_res: float,
                   baseline: float,
                   image_size: int):
        self.array['diameter'] = diameter
        self.array['quantum_eff'] = quantum_eff
        self.array['throughput'] = throughput
        self.array['wl_limits'] = wl_limits
        self.array['spec_res'] = spec_res
        self.array['baseline'] = baseline
        self.other['image_size'] = image_size
