# `COSMO`: COarse-grained Simulation of intrinsically disordered prOteins with openMM 

-------------------------------------
<b>Package under refactoring</b>:
- Initially, the package was named `hpsOpenMM` because it used the hps scale to simulate IDP. 
However, we are now implementing more models, such as mpipi, which do not use the hps scale anymore. 
Therefore, a different name is needed. We decided to change the name to `cosmo`.

- If the name `hpsOpenMM` still appears in the documentation/source code, it means that we missed some instances of 
renaming.
-------------------------------------
### A coarse-grained simulation engine empowered by openMM

Currently, there are four models are supported:

1) `hps_urry:` Hydropathy according to Urry scale (default, Recommended).
2) `hps_kr:`  Kapcha-Rossy scale.
   This model has parameters for nucleic acids and post-translational modification residues.
3) `hps_ss:` `hps_urry` with bonded potential.  
4) `mpipi`: another model that using Wang-Frenkel short range potential instead of LJ 12-6
5) Other models can be easily implemented by defining them in `cosmo/parameters/model_parameters.py`

Model summary
-------------
| Model      | components support                    | Implemented             | Tested  |
|------------|---------------------------------------|-------------------------|---------|
| `hps_kr`   | protein, RNA, phosphorylation protein | protein, p-protein, RNA | protein |
| `hps_urry` | protein, DNA                          | protein                 | protein |
| `hps_ss`   | protein                               | protein                 | protein |
| `mpipi`    | protein, RNA                          | protein                 | protein |

The package is ready for studying various problems such as, conformation dynamics of single chain, LLPS ...

Checkout documentation for more details: [here](https://qvv5013.github.io/docs-hpsOpenMM/). 
A simple example can be found [here](https://qvv5013.github.io/posts/openMM/hpsOpenMM_tutorial.html).

More information can be found in my personal blog: https://qvv5013.github.io/

-------------------------------------

## Requirements:

- **OpenMM >=7.7<sup>a,b</sup>** (select cuda version that compatible with your nvidia driver)
- **Parmed**
---
<sup>a</sup>: function `getStepCount()` does not work as expected (or is not implemented in versions earlier than 7.7).
This function is necessary when restarting simulations.

<sup>b</sup>: I recommend to upgrade to openMM 8.0 for better performance.

## How to use COSMO:

#### Linux:
The main requirements is `openMM >= 7.7`. Other packages are requires as well (see `requirements.txt`)

- Create conda environment: `conda create -n py310 python=3.10`
- activate `py310` env : `conda activate py310`
- Install `openMM 7.7` (or later): `conda install -c conda-forge openmm=7.7 cudatoolkit=10.2`
  * conda will try to install the latest version of `cudatoolkit` and 
  sometime it will not work. </br> You should select version that is compatible with your nvidia-driver (if you have NVIDIA GPU)
- Download folder and place in target location, for example: </br>`PATH_TO_CODE/cosmo/`
- Add `cosmo` module in Python path so that Python know what `cosmo` is (in `.bashrc` file): 
`export PYTHONPATH=$PYTHONPATH:PATH_TO_CODE/cosmo/`

<u>REMEMBER TO CHANGE `PATH_TO_CODE` TO YOUR SPECIFIC.</u>

#### Example:
- The standard example can be found at `example/standard_example`. 
You will need a control parameter file (e.g `md.ini`). Check [here](https://qvv5013.github.io/docs-hpsOpenMM/usage/simulation_control.html) for more information. 

Example of `md.ini`:
--------------------

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
ref_t = 310 ; Kelvin- reference temperature
tau_t = 0.01 ; ps^-1

;pressure coupling
pcoupl = no
ref_p = 1
frequency_p = 25

; Periodic boundary condition: if pcoupl is yes then pbc must be yes.
pbc = yes
; if pbc=yes, then use box_dimension option to specify box_dimension = x or [x, y, z], unit of nanometer
box_dimension = 30 ; [30, 30, 60]

; input
protein_code = ASYN
pdb_file = asyn.pdb
; output
checkpoint = asyn.chk
;Use GPU/CPU
device = GPU
; If CPU is specified, then use ppn variable
ppn = 4
;Restart simulation
restart = no
minimize = yes ;if not restart, then minimize will be loaded, otherwise, minimize=False

```
-------------------------------
- Run simulation: 
  - I give you two options to perform simulations:
    * First option:
      - goto example folder, e.g `examples/standard_example/`: 
      - edit simulation config file: `md.ini`
      - execute command: `python run_simulation.py -f md.ini`
    * Second option:
      * add python environment created above in the beginning of `cosmo-simulation.py` script: 
       `/home/qvv5013/anaconda3/envs/hpsopenmm/bin/python` - point to the environment you created above
      * make the `cosmo-simulation.py` script to be executable: `chmod +x /PATH_TO_COSMO/cosmo-simulation.py`
      * Create an `alias` to `cosmo-simulation.py`. e.g: I modify my `.bashrc`: 
        `alias cosmo-simulation='/home/qvv5013/work3/code/cosmo/cosmo/cosmo-simulation.py '`
      * in your simulation directory, prepare control file contains simulation parameters.
      * run simulation: `cosmo-simulation -f md.ini`
       
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


## Cite this work

This software is based on the original work of Prof. Jeetain Mittal's group (hps family) and 
Prof. Rosana Collepardo-Guevara's group (mpipi model). We have not published any paper using this software yet. 
If you have used it in your publications, please cite the source accordingly.

* hps family:
  - (`hps-kr`) Dignon, G. L.; Zheng, W.; Kim, Y. C.; Best, R. B.; Mittal, J. Sequence Determinants of Protein Phase Behavior from a Coarse-Grained Model. PLoS Comput. Biol. 2018, 1–23. https://doi.org/10.1101/238170.
  - (`hps-urry`) Regy, R. M.; Thompson, J.; Kim, Y. C.; Mittal, J. Improved Coarse-Grained Model for Studying Sequence Dependent Phase Separation of Disordered Proteins. Protein Sci. 2021, 30 (7), 1371–1379. https://doi.org/10.1002/pro.4094.
  - (`hps-ss`) Rizuan, A.; Jovic, N.; Phan, T. M.; Kim, Y. C.; Mittal, J. Developing Bonded Potentials for a Coarse-Grained Model of Intrinsically Disordered Proteins. J. Chem. Inf. Model. 2022, 62 (18), 4474–4485. https://doi.org/10.1021/acs.jcim.2c00450.

* mpipi
  - (`mpipi`)  Joseph, J. A.; Reinhardt, A.; Aguirre, A.; Chew, P. Y.; Russell, K. O.; Espinosa, J. R.; Garaizar, A.; Collepardo-Guevara, R. Physics-Driven Coarse-Grained Model for Biomolecular Phase Separation with near-Quantitative Accuracy. Nat. Comput. Sci. 2021, 1 (11), 732–743. https://doi.org/10.1038/s43588-021-00155-3.