# TODO

Current work is under **Active**; lower-priority ideas under **Later / not urgent**.

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
