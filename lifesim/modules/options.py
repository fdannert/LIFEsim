import warnings

import numpy as np


class Options(object):
    def __init__(self):
        self.array = {'diameter': 0.,
                      'quantum_eff': 0.,
                      'throughput': 0.,
                      'wl_min': 0.,
                      'wl_max': 0.,
                      'spec_res': 0,
                      'baseline': 0.,
                      'bl_min': 0.,
                      'bl_max': 0.,
                      'ratio': 0.}

        self.other = {'image_size': 0,
                      'wl_optimal': 0.,
                      'n_plugins': 0}

        self.models = {'localzodi': '',
                       'habitable': ''}

    def set_scenario(self,
                     case: str):

        self.array['quantum_eff'] = 0.7
        self.array['throughput'] = 0.05
        self.array['spec_res'] = 20.
        self.array['baseline'] = 20.
        self.array['bl_min'] = 10.
        self.array['bl_max'] = 100.
        self.array['ratio'] = 6.

        self.other['image_size'] = 512
        self.other['wl_optimal'] = 15
        self.other['n_plugins'] = 5

        self.models['localzodi'] = 'darwinsim'
        self.models['habitable'] = 'MS'

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

    def set_manual(self, **kwargs):
        for i, key in enumerate(kwargs.keys()):
            option_set = False
            for sub_dict in [self.array, self.other, self.models]:
                if key in sub_dict:
                    sub_dict[key] = kwargs[key]
                    option_set = True
                    break
            if not option_set:
                raise ValueError(str(key) + ' is an unknown option')

