# COSMO: COarse-grained Simulation of intrinsically disordered prOteins with OpenMM

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21361272.svg)](https://doi.org/10.5281/zenodo.21361272)

COSMO is a coarse-grained simulation engine for intrinsically disordered proteins and
related biomolecules, built on OpenMM.

One coarse-grained (one bead per residue / nucleotide) model powers **two complementary
workflows**:

- **A · Equilibrium simulation** — study single-chain conformations, liquid–liquid phase
  separation (LLPS / condensates), and protein–RNA complexes with the HPS / `mpipi` force
  fields.
- **B · Co-translational synthesis** — grow the nascent chain N→C one residue at a time
  with codon-resolved (O'Brien) kinetics, so the IDP is sampled *as it emerges* from the
  ribosome exit tunnel (analytic-tunnel or explicit CG-ribosome variants).

Part B reuses the Part A force fields, so start with A if you are new here.

## Supported models

Currently, these models are supported:

| Model      | Description                                                                 | Components supported                  | Implemented             |
|------------|----------------------------------------------------------------------------|---------------------------------------|-------------------------|
| `hps_urry` | Hydropathy according to the Urry scale (default, recommended).             | protein, DNA                          | protein                 |
| `hps_kr`   | Kapcha-Rossy scale; includes parameters for nucleic acids and PTMs.        | protein, RNA, phosphorylation protein | protein, p-protein, RNA |
| `hps_ss`   | `hps_urry` with bonded potential.                                          | protein                               | protein                 |
| `mpipi`    | Wang-Frenkel short-range potential instead of LJ 12-6.                     | protein, RNA                          | protein, RNA            |
| *(custom)* | New model — defined in `cosmo/parameters/model_parameters.py`.             | user-defined                          | user-defined            |

Additional models can be added easily by defining them in
`cosmo/parameters/model_parameters.py`, following the existing entries above.

COSMO can be used to study single-chain conformations, LLPS, and related phenomena.

## Documentation and tutorials

- Documentation: https://vuqv.github.io/cosmo/
- Tutorials: hands-on, ready-to-run examples in [`tutorials/`](tutorials/) (start
  with [tutorial 1](tutorials/01_single_chain_quickstart/); see the
  [tutorials overview](tutorials/README.md)).
- Additional notes may be given at: https://vuqv.github.io/blog/

## Requirements

- **Python** ≥ 3.9
- **OpenMM** ≥ 7.7 (the MD engine; install from the `conda-forge` channel, and choose a
  CUDA toolkit compatible with your NVIDIA driver)
- **ParmEd** ≥ 3.4
- **NumPy** ≥ 1.22, **pandas** ≥ 1.4
- **MDAnalysis** ≥ 2.2, **mdtraj** ≥ 1.9.7 (trajectory / structure I/O)

These are the exact runtime dependencies of `import cosmo`, declared in
[`pyproject.toml`](pyproject.toml). Floors are the oldest versions known to work;
there are no upper caps — the package runs on current releases (e.g. NumPy 2.x,
OpenMM 8.x). The standalone analysis tools under `scripts/` need a few extra packages
(scipy, matplotlib, numba, ...); those live in `requirements.txt`, not in
`pyproject.toml`.

Notes:
- `getStepCount()` is not available or reliable before OpenMM 7.7 and is required to
  restart simulations.
- OpenMM 8.2 is recommended for better performance.

## Installation and usage (Linux)

1) Create and activate a fresh conda/mamba environment with a recent Python
   (e.g. 3.9+):
   `conda create -n cosmo python=3.10 && conda activate cosmo`
2) Install OpenMM 7.7 or later:
   `conda install -c conda-forge openmm=7.7 cudatoolkit=10.2`
   - Conda may pick a newer CUDA toolkit by default. Choose a version compatible with
     your NVIDIA driver.
3) Install the remaining dependencies. The full list (numpy, parmed, MDAnalysis,
   mdtraj, pandas, plus analysis-only extras used by `scripts/`) is in
   `requirements.txt`:
   `conda install --file requirements.txt`
   - If you use the pip install in step 5, COSMO's runtime dependencies are pulled
     in automatically from `pyproject.toml`; `requirements.txt` still adds the
     analysis-only packages.
4) Download this repository to a target path, for example: `PATH_TO_CODE/cosmo/`
5) Install COSMO. Two main ways:

   **(I) Add to PYTHONPATH (no install)** — add the repo root to your Python path
   in `.bashrc`:
   `export PYTHONPATH=$PYTHONPATH:PATH_TO_CODE/cosmo/`
   This exposes `import cosmo` and the module form of every tool
   (`python -m cosmo.mdrun`, etc.), but not the console commands.

   **(II) Install with pip** — from the repo root (the directory with
   `pyproject.toml`). This additionally registers the `cosmo-mdrun` and other
   console commands (see [Console commands](#console-commands)) on your CLI:
   - editable (recommended for development), reflects source edits immediately:
     `pip install -e .`
   - regular install, copies the package into `site-packages` (source edits
     require a reinstall):
     `pip install .`

Remember to replace `PATH_TO_CODE` with your actual path.

## Console commands

`pip install` registers these entry points (each also has a module form,
`python -m <module>`):

| Command           | Module                 | Purpose                                                              |
| ----------------- | ---------------------- | ------------------------------------------------------------------- |
| `cosmo-mdrun`     | `cosmo.mdrun`          | Run an IDP / condensate simulation from an `md.ini` control file.    |
| `cosmo-csp`       | `cosmo.csp.protocol`   | Co-translational synthesis on an explicit coarse-grained ribosome.   |
| `cosmo-cylinder`  | `cosmo.csp.cylinder`   | Co-translational synthesis through an analytic (cylindrical) tunnel. |
| `cosmo-csp-movie` | `cosmo.csp.movie`      | Stitch per-residue/-stage synthesis trajectories into one VMD movie. |
| `cosmo-make-mrna` | `cosmo.csp.synth_mrna` | Pre-generate a fastest/slowest synonymous-codon mRNA for a protein.  |

(`cosmo-simulation` is a back-compat alias for `cosmo-mdrun`.)

## Example

A ready-to-run example is at `tutorials/01_single_chain_quickstart`. You will need a control file
(for example, `md.ini`). See the parameter reference here:
https://vuqv.github.io/cosmo/usage/simulation_control.html

Example `md.ini`:

```
[OPTIONS]
md_steps = 10_000 # number of steps
dt = 0.01 ; time step in ps
nstxout = 100 ; steps between trajectory/checkpoint frames
nstlog = 100 ; steps between log lines
nstcomm = 100 ; frequency for center-of-mass motion removal
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
pdb_file = asyn.pdb
; output  (all files -> <output_dir>/<outname>.*, default traj/traj.*)
output_dir = traj
outname = asyn
; use GPU/CPU
device = GPU
; if CPU is specified, then use ppn
ppn = 4
; restart simulation
restart = no
minimize = yes ; if not restart, then minimize will be loaded
```

## Running a simulation

All you need is a control file (`md.ini`). The canonical runner reads it, builds
the coarse-grained model, runs the simulation, and writes the trajectory, log,
checkpoint and final structure. The following are equivalent — pick whichever
suits your install:

```bash
cosmo-mdrun -f md.ini            # console command (after `pip install -e .`)
python -m cosmo.mdrun -f md.ini  # module form, works with the PYTHONPATH install
```

To run from a tutorial directory (e.g. `tutorials/01_single_chain_quickstart/`), edit its
`md.ini` and run any of the commands above; a thin `run_simulation.py` wrapper is
also provided there, so `python run_simulation.py -f md.ini` does the same thing.

> The runner is `cosmo.mdrun.mdrun`. Control-file parsing lives in
> `cosmo.read_simulation_config` and the build/run steps in `cosmo.engine`, so you
> can also drive a custom workflow from Python (see `cosmo/csp/` — the co-translational
> synthesis subsystem — for a specialized driver built on these pieces).

## Co-translational protein synthesis

Beyond equilibrium simulation of a full-length chain, COSMO can **grow the nascent chain
N→C one residue at a time** with codon-resolved (O'Brien) kinetics, so the IDP is sampled
*as it emerges* from the ribosome exit tunnel. This subsystem lives in `cosmo/csp/` and
mirrors the sibling [`topo`](https://github.com/vuqv/topo) package's synthesis design; it
reuses the same HPS / `mpipi` force fields as the equilibrium runner. Two runners differ
only in how the ribosome is represented:

```bash
cosmo-csp -f csp.ini             # explicit rigid coarse-grained ribosome (A-/P-site anchors)
cosmo-cylinder -f cylinder.ini   # analytic cylindrical exit tunnel (no ribosome beads)
```

Each residue is timed from its codon and — for `cosmo-csp` — split into O'Brien's three
kinetic sub-stages (peptidyl-transfer / translocation / tRNA-binding). After the final
residue, up to two optional **post-synthesis phases** run in order (each `0` = skip):

- **`stall_steps`** — hold the finished chain at the PTC with the C-terminus restraint /
  tRNA tether still **on** (mimics ribosome stalling), then
- **`ejection_steps`** — release the restraint so the finished chain diffuses free.

Runs are resumable (`resume = auto`) and write per-residue `L_<L>/` folders plus optional
`stall/` and `ejection/` folders. Stitch the whole growth into one VMD-playable movie with
`cosmo-csp-movie -o <out_root>`. Ready-to-run examples:
[tutorials/07_csp_cylinder/](tutorials/07_csp_cylinder/) (analytic tunnel) and
[tutorials/08_csp_cg_ribosome/](tutorials/08_csp_cg_ribosome/) (CG ribosome); see the
[synthesis-control reference](https://vuqv.github.io/cosmo/usage/synthesis_control.html)
for every `csp.ini` / `cylinder.ini` key.

## Windows

Not tested yet.

## macOS

Not tested yet.

## Cluster note

When submitting jobs on a cluster, `.bashrc` may not load. Source conda manually in the
job script:

`source PATH_TO_ANACONDA/anaconda3/etc/profile.d/conda.sh`

Then activate your environment, for example:

`conda activate cosmo`

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
If your work is published with the help of this package, please also give it a star
on GitHub in addition to citing it — it helps others find the project.

- **Paper (preferred):** Vu, Q. V.; Sitarik, I.; Li, M. S.; O’Brien, E. P. Noncovalent
  Lasso Entanglements are Common in Experimentally Derived Intrinsically Disordered
  Protein Ensembles and Strongly Influenced by Protein Length and Charge. J Phys Chem B
  129, 4682–4691 (2025).

- **Software:** Vu, Q. (2026). *COSMO: COarse-grained Simulation of intrinsically
  disordered prOteins with OpenMM* (Version 2026.1) [Computer software]. Zenodo.
  https://doi.org/10.5281/zenodo.21361272 (concept DOI, all versions; use
  [10.5281/zenodo.21361273](https://doi.org/10.5281/zenodo.21361273) to pin 2026.1).
  GitHub's "Cite this repository" button reads
  [`CITATION.cff`](CITATION.cff); see the [citation guide](docs/citation.rst) for the
  full list of models and datasets to cite.

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