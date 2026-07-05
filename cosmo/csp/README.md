# `cosmo.csp` — Continuous Synthesis Protocol (O'Brien) on cosmo's IDP force field

The per-codon, three-stage co-translational synthesis protocol of
`continuous_synthesis_v6.py` (Yang Jiang, Dan Nissley, Ed O'Brien), ported to cosmo.
`cosmo.csp` (the CSP runner, `cosmo-csp`) times every residue from its mRNA codon and
splits it into the three sub-stages of the elongation cycle. It is a thin outer loop
over the shared low-level MD engine `cosmo.csp.core` (`run_length`, `RunParams`), the
rigid-ribosome scenery `cosmo.csp.ribosome`, and the timing core `cosmo.csp.kinetics`.

```bash
cosmo-csp -f csp.ini
python -m cosmo.csp -f csp.ini
```

This package **mirrors the sibling `topo.csp`** package (same module layout, class /
function names, INI keys and CLI) but runs on cosmo's **sequence-based IDP force
field** (HPS / mpipi) instead of topo's structure-based Gō model. It therefore drops
topo's Gō machinery entirely: **no STRIDE, no native-contact map, no `domain.yaml`
nscales, no build-once-subset**. A length-`L` nascent model is simply
`cosmo.models.buildCoarseGrainModel` on the **first `L` residues** of the sequence (all
cosmo forces are sequence-local or pairwise-by-type, so the restriction is exact — see
[`cosmo/translation/PLAN.md`](../translation/PLAN.md) §2). topo's per-stage dt-halving
stability guard (against stiff Gō-contact divergence) is likewise omitted — cosmo's
soft HPS potentials have no such stiff wells.

> `cosmo.csp` is a **new, self-contained** package. The older `cosmo.translation`
> (single-stage constant-schedule `cosmo-elongate`) is left in place and untouched.

## What it does

For nascent length `L = L0 .. L_max`, each residue is added through three sub-stages,
each run for an **exponentially-sampled** dwell time (first-passage-time sampling):

| stage | biology | dwell-time mean (s) | C-terminus restraint target |
|-------|---------|---------------------|-----------------------------|
| 1 | peptidyl transfer / A-site delivery | `time_stage_1` | → A-target |
| 2 | translocation begins | `time_stage_2` | → A-target |
| 3 | tRNA binding / waiting | `codon_total − stage1 − stage2` | → **P-target** |

The A→P restraint switch between stages 2 and 3 reproduces translocation. Stage 3's
final structure seeds the next residue. Because the target must switch A↔P, CSP drives
the **position-restraint** path (the O'Brien tRNA tether is forced off). Dwell times
map to integration steps the O'Brien way:

```
t_sim_ns = dwell_s * 1e9 / scale_factor          # in-vivo -> in-silico
steps    = int(t_sim_ns / (dt_ps * 1e-3))        # -> MD steps
```

## Files

| file | role |
|------|------|
| `kinetics.py` | pure timing core (no OpenMM): codon tables, FPT sampling, `scale_factor`→steps, the 3-stage split. Ported near-verbatim from `topo.csp.kinetics`. |
| `core.py` | shared per-length MD engine (a library): `RunParams`, `run_length` (build / seed / restrain / minimize / run / finalize), `read_anchor`, coordinate seeding. Reused by both runners. |
| `protocol.py` | **the CSP runner** (`cosmo-csp`): per-codon, 3-stage loop + INI (`CSPConfig`, `read_csp_config`, `run_continuous_synthesis`, `csp()`). Calls `core.run_length` three times per residue. |
| `cylinder.py` | **the cylinder runner** (`cosmo-cylinder`): analytic exit tunnel (a bore in an infinite wall, `add_tunnel_cylinder`) instead of explicit beads. Same codon kinetics, **one MD segment per residue** (no A/P sub-stages). |
| `ribosome.py` | rigid CG ribosome scenery + ribosome↔nascent excluded volume / electrostatics; the tRNA tether and planar tunnel wall. |
| `cg_ribosome.py` | all-atom → CG ribosome (both `rna_model` reps: topo P/R/BR or cosmo 1-bead). |
| `truncate_ribosome.py` | crop the CG ribosome around the exit tunnel. |
| `movie.py` | stitch the per-residue/-stage trajectories into one VMD movie (`cosmo-csp-movie`; auto-detects the 3-stage vs flat layout). |
| `data/` | bundled `ecoli_trans_times_310K.txt` — the default E. coli (310 K) codon-time table (organism-universal). |
| `structures/` | truncated CG ribosomes (topo + cosmo rRNA reps). |

## Public API

```python
from cosmo.csp import run_continuous_synthesis, read_csp_config, RunParams, kinetics

cfg = read_csp_config("csp.ini")
run_continuous_synthesis(cfg.pdb_file, cfg.ribosome, L0=cfg.L0, L_max=cfg.L_max,
                         out_root=cfg.outdir, mrna=cfg.mrna,
                         codon_time_table_path=cfg.codon_time_table_path,
                         params=cfg.params)
```

## Control file (`csp.ini`)

A single `[OPTIONS]` section. Required: `pdb_file`, `ribosome`. `L0` (default `1`) and
`L_max` (default = full length) are optional. Per-codon timing also requires `mrna`
(`codon_times` is optional — defaults to the bundled E. coli 310 K table).

**Kinetic keys:** `mrna`, `codon_times` (a table path for per-codon timing, **or** a
positive number of seconds for a uniform codon time), `scale_factor`, `time_stage_1`,
`time_stage_2`, `random_seed`, `max_steps_per_stage` / `min_steps_per_stage` (test
clamps), `ejection_steps` / `dissociation_steps` (post-synthesis free runs).

**MD / ribosome keys:** `model` (nascent force field — **`hps_kr` default**), `dt`,
`ref_t`, `tau_t`, `nstout`, `device`, `ppn`, `constraints` (default `None` — flexible
bonds), `restraint_k`, `minimize`, `tunnel_wall`. The wall plane is auto-derived from the
structure; output is always nascent-only.

**PTC geometry (always optimized).** `optimal_ptc_targets` places the A-site seed /
stage-1-2 restraint target and the P-site / stage-3 target **exactly one peptide bond
apart** (the model's `bond_length_protein` — **0.380 nm** for `hps_kr`) and clear of the
ribosome excluded volume, by minimizing the soft O'Brien tRNA-bond/angle/improper
restraints + the 12-10-6 wall. Each new residue is delivered with its peptide bond at
equilibrium (stage-1 energies drop by ~50× versus seeding at the raw `AtR`/`PtR`-76 `R`
anchor beads, which sit ~0.9 nm apart). The C-terminus is held by a **position
restraint** to these points (not a bond to the tRNA beads). This is always on; there is
no knob.

**Ribosome ↔ nascent excluded volume — O'Brien 12-10-6.** The rigid ribosome interacts
with the nascent chain through the O'Brien form
`U = ε·[13(R/r)¹² − 18(R/r)¹⁰ + 4(R/r)⁶]` (ε = 0.000132 kcal/mol) with the **sum**
combining rule `R = Rmin/2ᵢ + Rmin/2ⱼ`, plus the Yukawa electrostatics. The per-bead
`Rmin/2` values are **inherited verbatim from topo** (O'Brien's CG collision radii) and
live in `cosmo.parameters.model_parameters` under **`hps_kr`** (per-AA + the rRNA
`P`/`R`/`BR` beads). The **nascent IDP↔IDP** interaction is unchanged — the `hps_kr`
Ashbaugh–Hatch potential — so `hps_kr` is the CSP force field (it carries both).

**Units:** OpenMM defaults — nm, ps, kJ/mol, K, kJ/mol/nm². In-vivo dwell times: seconds.

## Output

`<outdir>/L_<L>/stage_<1,2,3>/` per residue (each a standalone cosmo run: `traj.dcd`,
`traj.log`, `traj.psf`, `traj.chk`, `traj_final.pdb`, `traj_runinfo.log`,
`native_1_<L>.pdb`), plus a per-residue `dwell_times.dat`, and optional `ejection/` /
`dissociation/`. Stage 3's `traj_final.pdb` seeds the next residue. The cylinder runner
writes the flat `<outdir>/L_<L>/` layout (one segment per residue).

See [`../../sandbox/validate/`](../../sandbox/validate/) for a runnable proof-of-concept
on α-synuclein (both runners + the movie stitcher).
