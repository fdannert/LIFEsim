# `ScienceYield` Wrapper

The purpose of the science yield wrapper is to collect all functions necessary to run a yield with LIFEsim.
While detection yields with single catalogs or instrument setups can be quickly run using 
`instrument.snr()` and `opt.ahgs()`, this wrapper is aimed at
- yields for varying instrument setups (e.g. aperture size)
- yields including the characterization yield

Central tools to achieve this are
- a table of all synthetic planets, to which LIFEsim adds the `snr_1h` column. We call this _snr table_.
- the optimizer, which based on the `snr_1h` column creates an optimized observing sequence

## Current State
So far, the ScienceYield wrapper has grown organically. 
Here, we provide an overview of the current state.

---

## `ScienceYield`

### Constructor

#### `__init__(self, ...)`

**Description:**  
Sets up the working environment, mainly paths to file locations.

**Arguments:**

| Argument        | Type | Description                                  |
|-----------------|------|----------------------------------------------|
| `config_path`   | str  | Path to the config file in .yaml format      |
| `catalog_path`  | str  | Path to the input catalog                    |
| `output_path`   | str  | Location where the snr catalogs are saved    |
| `n_cpu`         | int  | Number of CPUs to be used in multiprocessing |
| `cat_from_ppop` | bool | True if input catalog is in .txt format      |

**Example:**
```python
import lifesim
ywrap = lifesim.ScienceYield(config_path='/path/CONFIG_251104.yaml',
                 catalog_path='/path/Bryson2021_hab2low.txt',
                 output_path='/path/hab2low_500/',
                 n_cpu=64)
```

#### `_compute_snrs(self, ...)`

**Description:**  
Calculates on individual SNR table for one setting of the instrument. 
Saves a catalog with `snr_1h` and a config file called `output_filename_catalog.hdf5` and `output_filename.yaml` to `output_path`.
This function is mainly meant to be an internal helper function, it is usually not required to use it.

**Arguments:**

| Argument          | Type | Description                                                                     |
|-------------------|------|---------------------------------------------------------------------------------|
| `output_path`     | str  | Path to where this individual run is saved                                      |
| `output_filename` | str  | Filename of catalog and config table                                            |
| `run_maxsep`      | bool | If true, all planets will be placed in their maximum separation                 |
| `diameter`        | int  | Mirror diameter for this run. If set to None, diameter from config file is used |

**Files created:**

| Name / Location                           | Content and Purpose            |
|-------------------------------------------|--------------------------------|
| `output_path/output_filename_catalog.hdf5`| snr table for the single run   |
| `output_path/output_filename.yaml`        |config file for the single run |

**Example:**
(continued from above)
```python
ywrap.compute_snrs(output_path='path/',
                   output_filename='single_test',
                   run_maxsep=False,
                   diameter=None)
```

#### `run_aperture_sweep_snr(self, ...)`

**Description:**
One of the key instrument parameters is the sensitivity.
For LIFE, we usually fix the photon-conversion-efficiency, and the aperture is varied.
This function creates snr tables for an array of mirror diameters. 

**Arguments:**

| Argument           | Type             | Description                               |
|--------------------|------------------|-------------------------------------------|
| `mirror_diameters` | list/np.ndarray  | Array of all mirror sizes to be evaluated |
| `run_name`         | str              | Name of the run (for file creation)       |

**Files created:**

For every diameter in `mirror_diameters`, the following files will be created (here we take 2.5m as example)

| Name / Location                                                         | Content and Purpose                                                                     |
|-------------------------------------------------------------------------|-----------------------------------------------------------------------------------------|
| `self.output_path/run_name/diam_2_5/`                                   | new parent directory is created                                                         |
| `self.output_path/run_name/diam_2_5/sweep_diam_2_5.yaml`                | _config file_ for run with the actual planet position for the detection phase           |
| `self.output_path/run_name/diam_2_5/sweep_diam_2_5_catalog.hdf5`        | &uarr; and its _snr catalog_                                                            |
| `self.output_path/run_name/diam_2_5/sweep_diam_maxsep_2_5.yaml`         | _config file_ for run with the maximum planet separation for the characterization phase |
| `self.output_path/run_name/diam_2_5/sweep_diam_maxsep_2_5_catalog.hdf5` | &uarr; and its _snr catalog_                                                            |

**Example:**
(continued from above)
```python
ywrap.run_aperture_sweep_snr(mirror_diameters=np.arange(2.0, 5.1, 0.25),
                             run_name='aperture_size_2_to_5')
```

#### `run_optimizer_sweep(self, ...)`

**Description:**
Runs the optimizer on all snr tables in a given source folder.

**Arguments:**

| Argument      | Type | Description                                                           |
|---------------|------|-----------------------------------------------------------------------|
| `run_name`    | str  | Name of the run (for file creation)                                   |
| `source_name` | str  | Name of the folder in self.output_path that contains the _snr tables_ |

**Files created:**

| Name / Location                                                | Content and Purpose                                            |
|----------------------------------------------------------------|----------------------------------------------------------------|
| `self.output_path/run_name/`                                   | new parent directory is created, if it does not already exists |
| `self.output_path/run_name/subdir/`                            | New directory for every directory in the source folder         |
| `self.output_path/run_name/subdir/sweep_diam_2_5.yaml`         | _config file_ for run of optimizer                             |
| `self.output_path/run_name/subdir/sweep_diam_2_5_catalog.hdf5` | &uarr; and its _snr catalog_ with optimized sequence           |

**Example:**
(continued from above)
```python
ywrap.run_optimizer_sweep(run_name='opt_aperture_size_2_to_5',
                          source_name='aperture_size_2_to_5')
```

---

## Helper Functions

#### `compute_yields_mp()`

**Description:**
Runs the optimizer on all snr tables in a given source folder.

**Arguments:**

| Argument          | Type     | Description                                                                   |
|-------------------|----------|-------------------------------------------------------------------------------|
| `output_filename` | str      | Name of output files                                                          |
| `output_path`     | str      | Path to which catalog and config file will be saved                           |
| `catalog_path`    | str      | Path to the input catalog file                                                |
| `config_path`     | str      | Path to the config file                                                       |
| `uni_sel`         | int/None | Number of universes to select from all universes. If `None`, all are selected |
| `return_yields`   | bool     | Set true to return experiment yield of this run                               |

**Files created:**

| Name / Location                           | Content and Purpose              |
|-------------------------------------------|----------------------------------|
| `output_path/output_filename_catalog.hdf5`| snr table for the single run     |
| `output_path/output_filename.yaml`        | config file for the single run   |


# What happens in the characterization ahgs



