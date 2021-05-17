import lifesim as ls

# ---------- SET-UP ----------

# create bus
bus = ls.Bus()

# setting the options
bus.data.options.set_scenario('baseline')

# set options manually
# bus.data.options.set_manual(diameter=4.)
bus.data.options.set_manual(bl_max=200.)

# loading and preparing the catalog
bus.data.catalog_from_ppop(input_path='/home/felix/Documents/MA/lifeOS/Data/baselineSample.fits')
bus.data.catalog_remove_distance(stype=0, mode='larger', dist=0.)  # remove all A stars
bus.data.catalog_remove_distance(stype=4, mode='larger', dist=10.)  # remove M stars > 10pc to speed up calculation

# ---------- MEASUREMENT SIMULATION ----------

# create modules and add to bus
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

# connect all modules
bus.connect(('inst', 'transm'))
bus.connect(('inst', 'exo'))
bus.connect(('inst', 'local'))
bus.connect(('inst', 'star'))

bus.connect(('star', 'transm'))


# run simulation. This function assigns every planet an SNR for 1 hour of integration time. Since
# we are currently only simulating photon noise, the SNR will scale with the integration time as sqrt(t)
inst.get_snr()

# ---------- DISTRIBUTE THE OBSERVATION TIME ----------
# After every planet is given an SNR, we want to distribute the time available in the search phase
# such that we maximize the number of detections.

# optimizing the result
opt = ls.Optimizer(name='opt')
bus.add_module(opt)
ahgs = ls.AhgsModule(name='ahgs')
bus.add_module(ahgs)

bus.connect(('transm', 'opt'))
bus.connect(('inst', 'opt'))
bus.connect(('opt', 'ahgs'))

opt.ahgs()

# save the catalog to a file
bus.data.export_catalog(output_path='/home/felix/Documents/MA/Outputs/Thesis/'
                                    'Baseline_AHGS_habitable_maxbl200.hdf5')


# ---------- READ SAVED CATALOG ----------
# import a previously saved catalog
# bus = ls.Bus()
# bus.data.options.set_scenario('baseline')
# bus.data.import_catalog(input_path='/home/felix/Documents/MA/Outputs/Thesis/Baseline_AHGS_habitable.hdf5')
