# hpsOpenMM

-------------------------------------
### A coarse-grained simulation engine empowered by openMM

OpenMM codebase for IDP which uses HPS-Urry/HPS-KR force field



`hps:` hydropathy scale. Currently, there are two models are supported:

1) `hps_urry:` Hydropathy according to Urry scale (Recommended).
2) `hps_kr:`  Kapcha-Rossy scale (default).
   This model has parameters for nucleic acids and post-translational modification residues.
3) other models can be easily implemented by defining them in `hps/parameters/model_parameters.py`

The package is ready for studying various problems such as, conformation dynamics of single chain, LLPS ...


-------------------------------------

## Requirements:

- **OpenMM 7.7** (select cuda version that compatible with your nvidia driver)
- **Parmed**

## How to use hpsOpenMM:

#### Linux:
- Create conda environment: `conda create -n hpsOpenMM python=3.10`
- activate `hpsOpenMM` env : `conda activate hpsOpenMM`
- Install openMM 7.7: `conda install openMM=7.7 conda install openmm=7.7 cudatoolkit=10.2`
- `(Note that cudatoolkit=10.2 for plgrid, 11.6 for masli6)`
- Download folder and place in target location, for example: `PATH_TO_CODE/hpsOpenMM/`
- Add folder in Python path (in `.bashrc` file): `export PYTHONPATH=$PYTHONPATH:PATH_TO_CODE/hpsOpenMM/`
- The standard example can be found at `example/standard_example`
- Run simulation: 
  - goto example folder, e.g examples/standard_example/: 
  - edit simulation config file: `md.ini`
  - execute command: `python run_simulation.py -f md.ini`
#### Windows:

- No idea (no time to test) !!!

#### MacOS:

- No money to test !!!

## Notes:

- *Note that in cluster, when submit job, the environment may not load `.bashrc`, so need to
  load conda environment in job file:*
  `source PATH_TO_ANNACONDA/anaconda3/etc/profile.d/conda.sh`
- Activate environment (e.g py310): `conda activate py310`

-------------------------------------

## Bugs

- If you encounter any bugs, please report the issue to Quyen Vu (`vuqv.phys@gmail.com`)

## Acknowledgments

We would like to thank:

- OpenMM team for developing an excellent and open source engine.
- SBMOpenMM's developer team for developing an excellent model and open sourcing the software.
- Sugita's group for making a clear tutorial to introduce many features of OpenMM
- OpenMM Github community for many useful discussions

## Do you want to cite this work?

Vu, Quyen. (2022). hpsOpenMM (1.2). Zenodo. https://doi.org/10.5281/zenodo.6976690

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6976690.svg)](https://doi.org/10.5281/zenodo.6976690)