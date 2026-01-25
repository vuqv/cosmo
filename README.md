# COSMO: COarse-grained Simulation of intrinsically disordered prOteins with OpenMM

COSMO is a coarse-grained simulation engine for intrinsically disordered proteins and
related biomolecules, built on OpenMM.

## Supported models

Currently, these models are supported:

1) `hps_urry`: Hydropathy according to the Urry scale (default, recommended).
2) `hps_kr`: Kapcha-Rossy scale. Includes parameters for nucleic acids and
   post-translational modifications.
3) `hps_ss`: `hps_urry` with bonded potential.
4) `mpipi`: Wang-Frenkel short-range potential instead of LJ 12-6.
5) Additional models can be added by defining them in
   `cosmo/parameters/model_parameters.py`.

### Models summary

| Model      | Components supported                  | Implemented             | Tested  |
|------------|---------------------------------------|-------------------------|---------|
| `hps_kr`   | protein, RNA, phosphorylation protein | protein, p-protein, RNA | protein |
| `hps_urry` | protein, DNA                          | protein                 | protein |
| `hps_ss`   | protein                               | protein                 | protein |
| `mpipi`    | protein, RNA                          | protein                 | protein |

COSMO can be used to study single-chain conformations, LLPS, and related phenomena.

## Documentation and tutorials

- Documentation: https://vuqv.github.io/docs-cosmo/
- Tutorial: https://vuqv.github.io/posts/openMM/cosmo_tutorial.html
- Additional notes: https://vuqv.github.io/

## Requirements

- OpenMM >= 7.7 (choose a CUDA toolkit compatible with your NVIDIA driver)
- ParmEd

Notes:
- `getStepCount()` is not available or reliable before OpenMM 7.7 and is required to
  restart simulations.
- OpenMM 8.2 is recommended for better performance.

## Installation and usage (Linux)

1) Create a conda environment:
   `conda create -n py310 python=3.10`
2) Activate it:
   `conda activate py310`
3) Install OpenMM 7.7 or later:
   `conda install -c conda-forge openmm=7.7 cudatoolkit=10.2`
   - Conda may pick a newer CUDA toolkit by default. Choose a version compatible with
     your NVIDIA driver.
4) Download this repository to a target path, for example: `PATH_TO_CODE/cosmo/`
5) Add the module to your Python path in `.bashrc`:
   `export PYTHONPATH=$PYTHONPATH:PATH_TO_CODE/cosmo/`

Remember to replace `PATH_TO_CODE` with your actual path.

## Example

The standard example is at `examples/standard_example`. You will need a control file
(for example, `md.ini`). See the parameter reference here:
https://vuqv.github.io/docs-cosmo/usage/simulation_control.html

Example `md.ini`:

```
[OPTIONS]
md_steps = 10_000 # number of steps
dt = 0.01 ; time step in ps
nstxout = 100 ; number of steps to write checkpoint = nstxout
nstlog = 100 ; number of steps to print log
nstcomm = 100 ; frequency for center of mass motion removal
; select model, available options: hps_kr, hps_urry, hps_ss or mpipi
model = mpipi

; control temperature coupling
tcoupl = yes
ref_t = 310 ; Kelvin - reference temperature
tau_t = 0.01 ; ps^-1

; pressure coupling
pcoupl = no
ref_p = 1
frequency_p = 25

; periodic boundary condition: if pcoupl is yes then pbc must be yes
pbc = yes
; if pbc=yes, then use box_dimension to specify x or [x, y, z] in nm
box_dimension = 30 ; [30, 30, 60]

; input
protein_code = ASYN
pdb_file = asyn.pdb
; output
checkpoint = asyn.chk
; use GPU/CPU
device = GPU
; if CPU is specified, then use ppn
ppn = 4
; restart simulation
restart = no
minimize = yes ; if not restart, then minimize will be loaded
```

## Running a simulation

Option 1: run from an example directory:

1) Go to `examples/standard_example/`.
2) Edit `md.ini`.
3) Run: `python run_simulation.py -f md.ini`

Option 2: use the `cosmo-simulation.py` wrapper:

1) Set the Python path at the top of `cosmo-simulation.py` to your environment:
   `$HOME/anaconda3/envs/py310/bin/python`
2) Make it executable:
   `chmod +x /PATH_TO_COSMO/cosmo-simulation.py`
3) Add an alias in `.bashrc`, for example:
   `alias cosmo-simulation='$HOME/work3/code/cosmo/cosmo/cosmo-simulation.py'`
4) In your simulation directory, prepare a control file.
5) Run: `cosmo-simulation -f md.ini`

## Windows

Not tested yet.

## macOS

Not tested yet.

## Cluster note

When submitting jobs on a cluster, `.bashrc` may not load. Source conda manually in the
job script:

`source PATH_TO_ANACONDA/anaconda3/etc/profile.d/conda.sh`

Then activate your environment, for example:

`conda activate py310`

## Bugs

If you encounter issues, please report them to Quyen Vu (`vuqv.phys@gmail.com`).
Please do not contact the model authors for COSMO-specific bugs.

## Acknowledgments

This project builds on the original work of Prof. Jeetain Mittal's group (hps family)
and Prof. Rosana Collepardo-Guevara's group (mpipi model).

This software was developed with time and financial support from the Institute of Physics,
Polish Academy of Sciences (Prof. Mai Suan Li), and Penn State University (Prof. Edward O'Brien).

## Cite this work

If you use this software in your publications, please cite the following sources.

- Vu, Q. V.; Sitarik, I.; Li, M. S.; O’Brien, E. P. Noncovalent Lasso Entanglements are
  Common in Experimentally Derived Intrinsically Disordered Protein Ensembles and
  Strongly Influenced by Protein Length and Charge. J Phys Chem B 129, 4682–4691 (2025).

- `hps` family (`hps-urry`, `hps-kr`, and `hps-ss`):
  - Dignon, G. L.; Zheng, W.; Kim, Y. C.; Best, R. B.; Mittal, J. Sequence Determinants
    of Protein Phase Behavior from a Coarse-Grained Model. PLoS Comput. Biol. 2018,
    1–23. https://doi.org/10.1101/238170.
  - Regy, R. M.; Thompson, J.; Kim, Y. C.; Mittal, J. Improved Coarse-Grained Model for
    Studying Sequence Dependent Phase Separation of Disordered Proteins. Protein Sci.
    2021, 30 (7), 1371–1379. https://doi.org/10.1002/pro.4094.
  - Rizuan, A.; Jovic, N.; Phan, T. M.; Kim, Y. C.; Mittal, J. Developing Bonded
    Potentials for a Coarse-Grained Model of Intrinsically Disordered Proteins. J. Chem.
    Inf. Model. 2022, 62 (18), 4474–4485. https://doi.org/10.1021/acs.jcim.2c00450.

- `mpipi` model:
  - Joseph, J. A.; Reinhardt, A.; Aguirre, A.; Chew, P. Y.; Russell, K. O.; Espinosa,
    J. R.; Garaizar, A.; Collepardo-Guevara, R. Physics-Driven Coarse-Grained Model for
    Biomolecular Phase Separation with near-Quantitative Accuracy. Nat. Comput. Sci.
    2021, 1 (11), 732–743. https://doi.org/10.1038/s43588-021-00155-3.