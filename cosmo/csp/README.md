# `cosmo.csp` — Continuous Synthesis Protocol (O'Brien) on cosmo's IDP force field

The per-codon, three-stage protein synthesis protocol of
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
cosmo forces are sequence-local or pairwise-by-type, so the restriction is exact). Like
topo, each stage is wrapped
in a **per-stage dt-halving stability guard** (`STABILITY_POTE_LIMIT_KJ` /
`STABILITY_MAX_ATTEMPTS` in `core.py`): the divergence it guards is **non-native excluded
volume** (the new residue seeded at the fixed A-site landing in the stiff repulsive core of
the ribosome↔nascent 12-10-6 EV or the nascent Ashbaugh-Hatch potential — both present in
cosmo *and* topo), not a native-contact well. A diverging stage (PotE → ~1e13 kJ/mol) is
re-run at `dt/2` with `2×` steps (dwell `= n_steps · dt` preserved) until it integrates
cleanly.

> `cosmo.csp` is a **self-contained** package and the sole protein synthesis
> subsystem. It replaced the older single-stage `cosmo.translation` (`cosmo-elongate`)
> package, which has been removed.

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
| `ribosome.py` | rigid CG ribosome scenery + ribosome↔nascent excluded volume / electrostatics; the tRNA tether and planar tunnel wall. Loads a CG ribosome PDB (topo P/R/BR rep) prepared with the sibling `topo` package. |
| `movie.py` | stitch the per-residue/-stage trajectories into one VMD movie (`cosmo-csp-movie`; auto-detects the 3-stage vs flat layout). |
| `synth_mrna.py` | build a `fastest`/`slowest`/`median` synonymous-codon mRNA for a protein from its PDB + a codon dwell-time table (`cosmo-make-mrna`). |

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
`L_max` (default = full length) are optional. Per-codon timing requires both `mrna` and
a `codon_times` table path (there is no bundled default — pick one under
`assets/csp/codon_dwell_times/`). Setting `mrna = fastest`, `slowest` or `median`
auto-builds a synonymous-codon mRNA from the protein + table.

**Kinetic keys:** `mrna` (an mRNA file, **or** `fastest`/`slowest`/`median` to auto-build a
synonymous-codon mRNA — each residue's fastest/slowest/median-dwell-time codon per the
`codon_times` table, written next to the PDB), `codon_times` (a table path for per-codon
timing — required, no
bundled default — **or** a positive number of seconds for a uniform codon time),
`scale_factor`, `time_stage_1`,
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
`Rmin/2` values are **inherited verbatim from topo** (O'Brien's CG collision radii).
They are **steric** radii — force-field-independent — so they live in
`cosmo.parameters.model_parameters` as standalone tables (`OBRIEN_RMIN_2_NM` per-AA +
`OBRIEN_RNA_RMIN_2_BEADS` for the rRNA `P`/`R`/`BR` beads), **not** attached to any IDP
model, and are used for **every** nascent force field alike. Only the bead *charges* are
model-dependent (the residue's formal charge, e.g. mpipi's partial charges). The
**nascent IDP↔IDP** interaction is whatever `model` the run selects (`hps_kr` /
`hps_urry` / `mpipi`) — any of them works with the CSP excluded volume.

**Units:** OpenMM defaults — nm, ps, kJ/mol, K, kJ/mol/nm². In-vivo dwell times: seconds.

## Output

`<outdir>/L_<L>/` — **one folder per residue** (consolidated layout): a shared
`traj.psf` + `native_1_<L>.pdb` (functions of `L` only, written once), per-stage
`traj_s{1,2,3}.dcd` + `.log`, one folded `traj_runinfo.log` (a section per stage), and a
single `traj_final.pdb` (stage-3 final — it seeds the next residue and is the
resume-reload target). No per-stage `.chk`. Plus a per-residue `dwell_times.dat` and
optional `ejection/` / `dissociation/`. The cylinder runner writes the flat
`<outdir>/L_<L>/traj.*` layout (one segment per residue).

**Resume.** Re-invoking `cosmo-csp` / `cosmo-cylinder` on an interrupted `<outdir>`
continues from the last completed residue — the schedule + PTC geometry are re-read from
`dwell_times.dat` (no RNG redraw / SLSQP re-solve) and the seed is reloaded from the last
`traj_final.pdb`, tracked by `progress.log`. On by default (`resume = auto`; `yes`/`no`
or `--fresh` override). See `cosmo.csp.resume`.

See [`../../sandbox/validate/`](../../sandbox/validate/) for a runnable proof-of-concept
on α-synuclein (both runners + the movie stitcher).
