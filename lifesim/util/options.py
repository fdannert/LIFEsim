import warnings

import numpy as np


# TODO: Rename habitable variable
class Options(object):
    """
    The Options class contains all options and settings needed for the simulation of the LIFE
    array. Any settings are exclusively found in the Options class and must only be set through
    changing attributes of the Options class itself.

    Attributes
    ----------
    array : dict
        Options concerning the LIFE array itself. They are
            - ``'diameter'`` : Diameter of single aperture in [m].
            - ``'quantum_eff'`` : Quantum efficiency of the detector in (# of electrons/# of
              photons).
            - ``'throughput'`` : Optical throughput of the whole instrument in [%].
            - ``'wl_min'`` : Minimum wavelength of the spectrometer in [microns].
            - ``'wl_max'`` : Maximum wavelength of the spectrometer in [microns].
            - ``'spec_res'`` : Spectral resolution of the spectrometer (dimensionless).
            - ``'baseline'`` : Length of the shorter nulling baseline in [m].
            - ``'bl_min'`` : Minimum allowed length of the shorter nulling baseline in [m].
            - ``'bl_max'`` : Maximum allowed length of the shorter nulling baseline in [m].
            - ``'ratio'`` : Ratio between the nulling and the imaging baseline. E.g. if the imaging
              baseline is twice as long as the nulling baseline, the ratio will be 2.
            - ``'t_slew'`` : Slew time required for the array to shift from one target system to
              the next in [s].
            - ``'t_efficiency'`` : Time efficiency of the observation accounting for overheads.
              E.g. if the time efficiency is 0.8, 80% of the on-target observation time will be
              actually spend integrating photons.
    other : dict
        Options concerning simulation parameters. They are
            - ``'image_size'`` : Number of pixels (in one axis) which will be simulated.
              Corresponds to the pixel resolution of the detector. I.e. if image size is 512, the
              detector will be simulated with 512^2 pixels.
            - ``'wl_optimal'`` : The wavelength to which the baseline is optimized in [micron].
            - ``'n_plugins'`` : Number of sockets the instrument class will feature.
    models : dict
        Options concerning different models used in the simulation. They are
            - ``'localzodi'`` : Model for the localzodi, possible options are ``'glasse'`` and
              ``'darwinsim'``
            - ``'habitable'`` : Model used for calculating the habitable zone, possible options are
              ``'MS'`` and ``'POST_MS'``
    optimization : dict
        Options concerning the methods used to optimally distribute the observing time.
            - ``'N_pf'`` : Number of sampling locations per orbit.
            - ``'snr_target'`` : Planets with a larger signal-to-noise ratio than `'snr_target'`
              are counted as detections.
            - ``'limit'`` : Limits the number of wanted detections for the different stellar host
              star types.
            - ``'habitable'`` : If true, the integration time is optimized towards planets residing
              in the habitable zone.
            - ``'t_search'`` : Duration of the search phase in [s].
    """
    def __init__(self):
        """
        """
        self.array = {'diameter': 0.,
                      'quantum_eff': 0.,
                      'throughput': 0.,
                      'wl_min': 0.,
                      'wl_max': 0.,
                      'spec_res': 0,
                      'baseline': 0.,
                      'bl_min': 0.,
                      'bl_max': 0.,
                      'ratio': 0.,
                      't_slew': 0.,
                      't_efficiency': 0.}

        self.other = {'image_size': 0,
                      'wl_optimal': 0.,
                      'n_plugins': 0}

        self.models = {'localzodi': '',
                       'habitable': ''}

        self.optimization = {'N_pf': 0.,
                             'snr_target': 0.,
                             'limit': None,
                             'habitable': False,
                             't_search': 0.}

    def set_scenario(self,
                     case: str):
        """
        Sets the options according the chosen scenario

        Parameters
        ----------
        case : str
            The scenario after which to set the options. Possible scenarios are ``'optimistic'``,
            ``'baseline'`` and ``'pessimistic'``
        """

        self.array['quantum_eff'] = 0.7
        self.array['throughput'] = 0.05
        self.array['spec_res'] = 20.
        self.array['baseline'] = 20.
        self.array['bl_min'] = 10.
        self.array['bl_max'] = 100.
        self.array['ratio'] = 6.
        self.array['t_slew'] = 10. * 60. * 60.
        self.array['t_efficiency'] = 0.8

        self.other['image_size'] = 256  # TODO: or 512?
        self.other['wl_optimal'] = 15
        self.other['n_plugins'] = 5

        self.models['localzodi'] = 'darwinsim'
        self.models['habitable'] = 'MS'

        self.optimization['N_pf'] = 25
        self.optimization['snr_target'] = 7
        self.optimization['limit'] = np.array(((0, 1, 2, 3, 4),
                                               (np.inf, np.inf, np.inf, np.inf, np.inf)))
        self.optimization['habitable'] = True
        self.optimization['t_search'] = 2.5 * 365. * 24. * 60. * 60.

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
        """
        For manually setting one or more options

        Parameters
        ----------
        kwargs : dict
            Dictionary containing the options with the option name as key and the option itself as
            value. E.g. ``kwarg = {'diameter': 2., 'n_plugins': 3, 'localzodi': 'glasse'}``

        Raises
        ------
        ValueError
            If attempting to set an unknown option
        """

        # cycle through all keys
        for i, key in enumerate(kwargs.keys()):
            option_set = False

            # check if the key exists in any of the options dictionaries
            for sub_dict in [self.array, self.other, self.models, self.optimization]:
                if key in sub_dict:

                    # set the option
                    sub_dict[key] = kwargs[key]
                    option_set = True
                    break

            # raise error if no option was set
            if not option_set:
                raise ValueError(str(key) + ' is an unknown option')

