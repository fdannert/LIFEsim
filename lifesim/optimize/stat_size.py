

import lifesim as ls


def run_stat_size(stat_size: float):
    bus = ls.Bus()

    # setting the options
    bus.data.options.set_scenario('baseline')
    bus.data.options.set_manual(stat_size=stat_size,
                                multi_scaler=None)


    # loading and preparing the catalog
    bus.data.import_catalog(input_path='/home/felix/Documents/MA/Outputs/LIFEsim_development/'
                                        'SNR_ahgs.hdf5')
    bus.data.catalog_remove_distance(stype=3, dist=5, mode='larger')
    bus.data.stars_from_catalog()
    bus.data.universe_from_catalog(nuniverse=10)


    # add all modules
    bld = ls.BlindSimulationModule(name='blind')
    bus.add_module(bld)
    zt = ls.PredictionIntervalModule(name='ZTest')
    bus.add_module(zt)
    obs = ls.SimpleObservationModule(name='obs')
    bus.add_module(obs)
    inst = ls.Instrument(name='inst')
    bus.add_module(inst)
    transm = ls.TransmissionMap(name='transm')
    bus.add_module(transm)
    exo = ls.PhotonNoiseExozodi(name='exo')
    bus.add_module(exo)
    local = ls.PhotonNoiseLocalzodi(name='local')
    bus.add_module(local)
    star = ls.PhotonNoiseStar(name='star')
    bus.add_module(star)
    int = ls.SimpleInterestModule(name='inter')
    bus.add_module(int)


    # connect all modules
    bus.connect(('inst', 'transm'))
    bus.connect(('inst', 'exo'))
    bus.connect(('inst', 'local'))
    bus.connect(('inst', 'star'))
    bus.connect(('star', 'transm'))
    bus.connect(('inter', 'ZTest'))
    bus.connect(('inter', 'blind'))
    bus.connect(('blind', 'inst'))
    bus.connect(('blind', 'obs'))
    bus.connect(('obs', 'inst'))


    # run simulation
    bld.run_simulation()


def run_grid(input):
    stat_size = input[0]
    multi_scaler = input[1]
    localzodi_scaler = input[2]

    bus = ls.Bus()

    # setting the options
    bus.data.options.set_scenario('baseline')
    bus.data.options.set_manual(stat_size=stat_size,
                                multi_scaler=multi_scaler,
                                localzodi_scaler=localzodi_scaler)

    # loading and preparing the catalog
    bus.data.import_catalog(input_path='/home/felix/Documents/MA/Outputs/LIFEsim_development/'
                                       'SNR_ahgs.hdf5')
    bus.data.catalog_remove_distance(stype=3, dist=5, mode='larger')
    bus.data.stars_from_catalog()
    bus.data.universe_from_catalog(nuniverse=10)

    # add all modules
    bld = ls.BlindSimulationModule(name='blind')
    bus.add_module(bld)
    zt = ls.PredictionIntervalModule(name='ZTest')
    bus.add_module(zt)
    obs = ls.SimpleObservationModule(name='obs')
    bus.add_module(obs)
    inst = ls.Instrument(name='inst')
    bus.add_module(inst)
    transm = ls.TransmissionMap(name='transm')
    bus.add_module(transm)
    exo = ls.PhotonNoiseExozodi(name='exo')
    bus.add_module(exo)
    local = ls.PhotonNoiseLocalzodi(name='local')
    bus.add_module(local)
    star = ls.PhotonNoiseStar(name='star')
    bus.add_module(star)
    int = ls.SimpleInterestModule(name='inter')
    bus.add_module(int)

    # connect all modules
    bus.connect(('inst', 'transm'))
    bus.connect(('inst', 'exo'))
    bus.connect(('inst', 'local'))
    bus.connect(('inst', 'star'))
    bus.connect(('star', 'transm'))
    bus.connect(('inter', 'ZTest'))
    bus.connect(('inter', 'blind'))
    bus.connect(('blind', 'inst'))
    bus.connect(('blind', 'obs'))
    bus.connect(('obs', 'inst'))

    # run simulation
    bld.run_simulation()

    bus.data.export_catalog(output_path='/home/felix/Documents/MA/Outputs/LIFEsim_development/'
                                        'Grid_Test_ss_' + str(stat_size).replace('.', '_')
                                        + '_ls_' + str(localzodi_scaler).replace('.', '_')
                                        + '_ms_' + str(multi_scaler).replace('.', '_')
                                        + '.hdf5')
