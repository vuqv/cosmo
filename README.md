# hpsOpenMM

OpenMM codebase for IDP which uses HPS-Urry/HPS-KR force field

This is a repository contains the source code of hps-urry/kr model.

`hps:` hydropathy scale. Currently, there are two models are supported:

1) `urry:` Hydropathy according to Urry scale.
2) `kr:`  Kapcha-Rossy scale (default).
   This model has parameters for nucleic acids and post-translational modification.

-------------------------------------

## Requirements:

- **OpenMM 7.7**
- **Parmed**

## How to use hpsOpenMM:

#### Linux:

- Download folder and place in target location, for example: `PATH_TO_CODE/hpsOpenMM/`
- Add folder in Python path (in `.bashrc` file): `export PYTHONPATH=$PYTHONPATH:PATH_TO_CODE/hpsOpenMM/`
- The standard example can be found at `example/standard_example`

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
- Sugita's group to make a clear tutorial to introduce many features of OpenMM
- OpenMM Github community for many useful discussions

## Do you want to cite this work?

- We will make this source code open to everyone to use.
