# TODO

Current work is under **Active**; lower-priority ideas under **Later / not urgent**;
finished items under **Done**.

## Active

### Add DNA support in the `hps_urry` model
The README lists `hps_urry` as supporting protein + DNA, but only protein is
implemented/tested. The HPS-scale models type residues by per-residue radius
(σ) and hydropathy (λ), so DNA needs those parameters for the nucleotide beads.

* Add the DNA nucleotide residues (DA/DT/DG/DC) to `parameters['hps_urry']` in
  `cosmo/parameters/model_parameters.py` with `mass`, `charge`, radius and
  hydropathy values (only the 20 amino acids are defined today;
  `bond_length_nucleic = 0.5` already exists). `parameters['hps_kr']` already
  carries nucleic-acid entries — use it as the template for structure/keys.
* Ensure `setRadiusPerResidueType` / `setHPSPerResidueType` (in
  `cosmo/core/system.py`) and `addAshbaughHatchForces` cover the DNA residues.
* Confirm the CA/P topology path keeps phosphate `P` beads and applies
  `bond_length_nucleic` for nucleic bonds (`getBonds` in `system.py`).
* Validate: add a DNA (or protein+DNA) case to `benchmarks/benchmark_energies.py`
  and update the README models table (`hps_urry` DNA → implemented/tested).

## Later / not urgent

Ideas to revisit when there's time; not on the critical path.

### Co-translational synthesis — remaining v2 extras (`cosmo/translation/`)
The subsystem (v1 + v2, incl. the tRNA tether and tunnel wall) is done (see **Done**
below); the post-elongation phase (ejection / stallation) is now implemented too
(`cosmo/csp/`). One refinement remains. Full checklist in
`tutorials/07_translation/TASKS.md` (phase 5b):

* **Analysis** — Rg / native-contact fraction vs. nascent length; tutorial worked
  example.

### Allow configurable `n_copies`
Let `md.ini` configure the number of chain copies in the system. Two distinct
modes, both worth supporting:

* **Non-interacting copies** — replicate the same chain as independent,
  non-interacting trajectories packed into one OpenMM `System` to better utilize
  the GPU (parallel sampling / replicas for better statistics). Copies must not
  see each other (exclusions / separated coordinates / no cross-copy forces).
* **Interacting copies** — multiple chains that *do* interact, i.e. the
  condensate / phase-separation setup standard in the IDP field. Run many chains
  together and scan concentration / temperature to locate the LLPS boundary
  (phase separation into a dense + dilute phase).

Design notes:
* Add an `n_copies` (and a mode flag, e.g. `copies_interacting = true/false`)
  option to `SimulationConfig` / `md.ini [OPTIONS]`.
* Wire into the model/system builder so copies are laid out without overlap
  (box packing) and, for the non-interacting mode, fully excluded from each other.

## Done

1. **Regression benchmark** — `benchmarks/` (`benchmark_energies.py`,
   `reference_energies.json`, `README.md`). Records per-force-group initial
   energies for `hps_urry`/`hps_kr`/`hps_ss`/`mpipi`; deterministic and exits
   non-zero on mismatch for CI.
2. **Align with the `topo` design** — config/engine/mdrun runner layer
   (`cosmo/utils/config.py`, `cosmo/engine.py`, `cosmo/mdrun/`) + `pyproject.toml`
   console entry points (`cosmo-mdrun`). Removed `cosmo.dynamics.Dynamics` and
   `cosmo-simulation.py`; example `run_simulation.py` scripts are thin wrappers
   (specialized ones keep only their custom bits). Tests in `tests/test_config.py`.
3. **Tutorials** — `tutorials/` numbered, self-contained folders (ready-to-run
   inputs + step-by-step `README.md`), wired into the docs under a Tutorials
   toctree.
4. **Furo full-width docs** — `docs/conf.py` uses `html_theme = 'furo'` with
   full-width content (`--content-width: 100%` + `docs/_static/custom.css`);
   removed the RTD/rst2pdf/GitPython machinery.
5. **Mpipi RNA support** — verified the protein implementation against the paper
   (Joseph et al. 2021, SI Tables 11–12: WF form, ε/σ matrices, μ exceptions,
   charges all correct). Added A/C/G/U beads (ids 20–23, −0.75e, standard masses)
   + aliases; extended the WF tables to 24×24 (μ=3 for RNA pairs, ν=1, R=3σ);
   raised the WF neighbour cutoff to max R_ij; made the Debye length
   model-dependent (mpipi → 0.795 nm). Validated: U–U WF == analytic Eq. (4) to
   1e-7, protein+RNA builds finite, protein benchmark unchanged. Tests in
   `tests/test_mpipi_rna.py`.
6. **Co-translational synthesis subsystem** (`cosmo/translation/`, mirrors
   `topo/translation/`) — replaces the removed `examples/growing/` prototype.
   * **Build step v1** — nascent-only elongation loop (`elongate.py`): per-length
     rebuild on the first `L` residues (sequence-based; no STRIDE/contacts),
     cold-start layout, A→P new-residue placement, C-terminus restraint to the
     P-anchor, reusing `cosmo.engine`. `cosmo-elongate` CLI.
   * **Build step v2** — rigid (mass-0) ribosome scenery (`ribosome.py`):
     ribosome↔nascent excluded volume + Yukawa, nascent-only output; dual RNA
     representation (`rna_model = topo` 3/4-bead P/R/BR, or `cosmo` 1-bead),
     `cg_ribosome.py` + `truncate_ribosome.py` to build them. `mpipi` dummy-id path.
     **v2 extras:** O'Brien tRNA tether (bond + orienting angle) and one-sided
     planar tunnel wall (forward-only extrusion), both on by default.
   * **Movie** — `make_movie.py` / `cosmo-elongate-movie` stitches the per-length
     trajectories into one VMD-playable growing-chain movie.
   * Validated on α-synuclein under both `rna_model` reps (ribosome frozen,
     C-terminus tracks the P-anchor, benchmark unchanged). Self-contained example
     in `tutorials/07_translation/`; status in its `TASKS.md`.
7. **Removed the deprecated `buildHPSModel` alias** — `buildCoarseGrainModel` is the
   sole builder name. Its only consumer (`examples/growing/`) is gone;
   `tests/test_mpipi_rna.py` and the docs/`CLAUDE.md` updated accordingly.
