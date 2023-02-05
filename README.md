# hpsOpenMM

-------------------------------------
### A coarse-grained simulation engine empowered by openMM

OpenMM codebase for IDPs simulation using hydropathy scale.



`hps:` hydropathy scale. Currently, there are four models are supported:

1) `hps_urry:` Hydropathy according to Urry scale (default, Recommended).
2) `hps_kr:`  Kapcha-Rossy scale.
   This model has parameters for nucleic acids and post-translational modification residues.
3) `hps_ss:` `hps_urry` with bonded potential.  
4) `mpipi`: another model that using Wang-Frenkel short range potential instead of LJ 12-6
4) other models can be easily implemented by defining them in `hps/parameters/model_parameters.py`

Model summary
-------------
| Model      | component support                     | Implemented             | Tested          |
|------------|---------------------------------------|-------------------------|-----------------|
| `hps_kr`   | protein, RNA, phosphorylation protein | protein, p-protein, RNA | protein |
| `hps_urry` | protein, DNA                          | protein                 | protein |
| `hps_ss`   | protein                               | protein                 | protein |
| `mpipi`    | protein, RNA                          | protein                 | protein |

The package is ready for studying various problems such as, conformation dynamics of single chain, LLPS ...

Check out documentation for more details: [here](https://qvv5013.github.io/docs-hpsOpenMM/)

-------------------------------------

## Requirements:

- **OpenMM 7.7** (select cuda version that compatible with your nvidia driver)
- **Parmed**

## How to use hpsOpenMM:

#### Linux:
The main requirements is openMM 7.7
For my specific case, I don't know why I can not install openMM 7.7 with python version < 3.10
For that reason, I will create a virtual environment with Python 3.10 and then install openMM 7.7.

- Create conda environment: `conda create -n hpsOpenMM python=3.10`
- activate `hpsOpenMM` env : `conda activate hpsOpenMM`
- Install openMM 7.7: `conda install openmm=7.7 cudatoolkit=10.2`
  (Note that cudatoolkit=10.2 for plgrid, 11.6 for masli6)
- Download folder and place in target location, for example: `PATH_TO_CODE/hpsOpenMM/`
- Add hps module in Python path so that Python know what `hps` is (in `.bashrc` file): `export PYTHONPATH=$PYTHONPATH:PATH_TO_CODE/hpsOpenMM/`
- The standard example can be found at `example/standard_example`
- Run simulation: 
  - I give you two options to perform simulations:
    * First option:
      - goto example folder, e.g examples/standard_example/: 
      - edit simulation config file: `md.ini`
      - execute command: `python run_simulation.py -f md.ini`
    * Second option:
      * add python environment created above in the beginning of `hps-simulation.py` script: by changing its first line.
      * make the `hps-simulation.py` script to be executable: `chmod +x /PATH_TO_HPS/hps-simulation.py`
      * Create an `alias` to `hps-simulation.py`. e.g: I modify my `.bashrc`: 
        `alias hps-simulation='/home/qvv5013/work3/code/hpsOpenMM/hps/hps-simulation.py '`
      * in your simulation directory, prepare control file contains simulation parameters.
      * run simulation: `hps-simulation -f md.ini`
       
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