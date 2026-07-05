# Synthesis on a coarse-grained ribosome

The **Continuous Synthesis Protocol (CSP)** is cosmo's **codon-resolved, kinetic** runner
for co-translational synthesis on an **explicit coarse-grained ribosome**. It times
**every residue from its mRNA codon** and splits each elongation cycle into **three
kinetic sub-stages** — reproducing O'Brien's `continuous_synthesis_v6.py` protocol, but
on cosmo's **sequence-based IDP force field** (HPS / mpipi) instead of a structure-based
Gō model. (Codon-resolved kinetics are what make the model physically meaningful; a
fixed per-residue step count is not.) For the analytic-tunnel variant — the same codon
kinetics with the explicit ribosome replaced by a cylindrical bore — see
{doc}`cylinder_synthesis`.

- **CLI:** `cosmo-csp -f csp.ini` (or `python -m cosmo.csp -f csp.ini`)
- **Movie tool:** `cosmo-csp-movie -o <out_root> [--ribosome ribo.pdb]`
- **Proof of concept:** `sandbox/validate/` — a smoke run (`csp.ini`, L = 5 → 8) on
  α-synuclein plus the cylinder variant.
- **Architecture:** CSP is a thin outer loop. The per-length MD work — building the
  length-`L` model, seeding coordinates, restraints, running one stage — lives in the
  shared low-level engine `cosmo.csp.core` (`run_length`, `RunParams`); the
  rigid-ribosome scenery and tunnel wall live in `cosmo.csp.ribosome`; the timing lives
  in `cosmo.csp.kinetics`. CSP adds only the kinetics and the
  three-`run_length`-calls-per-residue loop.

```{note}
**This is the cosmo port of `topo.csp`.** It mirrors the sibling `topo` package's module
layout, INI keys and CLI, but the nascent chain is an **intrinsically disordered protein**
under the HPS / mpipi force field — so there is **no STRIDE, no native-contact map, no
`domain.yaml`, and no build-once-subset machinery**. A length-`L` model is simply
`cosmo.models.buildCoarseGrainModel` on the **first `L` residues** of the sequence. The
older single-stage `cosmo.translation` (`cosmo-elongate`) package is left untouched.
```

---

## Quick start

All paths in the INI are relative to the working directory; run from the example folder.
A GPU is recommended for production (the explicit-ribosome system has ~4,600 rigid beads).

```bash
cd sandbox/validate
cosmo-csp -f csp.ini              # 3-stage synthesis -> synth_out_csp/

# stitch the per-stage trajectories into one VMD movie
cosmo-csp-movie -o synth_out_csp --ribosome ribosome_trunc.pdb
vmd -e synth_out_csp/movie.tcl
```

`cosmo-csp` writes, per residue `L` and sub-stage `s`, a standalone trajectory under
`<outdir>/L_<L>/stage_<s>/`, an optional `ejection/` (and `dissociation/`) phase, and a
per-residue dwell-time log `<outdir>/dwell_times.dat`.

---

## Theory

### 1. What is being modeled: co-translational synthesis

In a living cell a protein is **synthesized vectorially, N-terminus first**, by the
ribosome, one amino acid at a time, while the growing ("nascent") chain threads out
through the ribosomal **exit tunnel** (~80 Å long, ~10–20 Å wide) and begins to fold or
collapse **co-translationally** — *as it emerges*. The kinetics of synthesis matter: how
long the ribosome dwells on each codon sets how much time each segment of the chain has to
sample conformations before the next residue is added. Rare codons (decoded slowly) act as
"translational pauses". CSP reproduces this by growing a coarse-grained protein
bead-by-bead out of a coarse-grained ribosome, **timing each residue from its mRNA codon**.

Because cosmo's chain is an **IDP** (no native fold), CSP here studies the
co-translational behaviour of *disordered* chains — extrusion, tunnel confinement,
compaction — rather than folding toward a native structure.

#### The real elongation cycle (one amino acid added)

Bacterial translation elongation repeats a three-step biochemical cycle per codon:

1. **Aminoacyl-tRNA selection / decoding.** A ternary complex delivers an aa-tRNA to the
   ribosomal **A site**; correct codon–anticodon pairing triggers accommodation. This is
   the **codon-dependent, highly variable, usually rate-limiting** step (cognate-tRNA
   abundance / codon-usage bias).
2. **Peptidyl transfer.** The peptidyl-transferase center (PTC) transfers the nascent
   peptide from the **P-site** tRNA onto the A-site aminoacyl-tRNA. The chain is now one
   residue longer and attached to the A-site tRNA. Fast (~0.3 ms).
3. **Translocation.** EF-G ratchets the ribosome forward by one codon: the tRNAs move
   **A→P**, the A site is freed. ~few ms.

CSP partitions the per-codon dwell time into these three pieces and reproduces them as
**three MD sub-stages per residue**.

### 2. The simulation model (one elongation step)

- **Nascent protein — a sequence-based IDP chain.** One bead per residue at the Cα
  position. Interactions are set by the **sequence**, not a structure: the short-range
  **Ashbaugh–Hatch** (HPS) or **Wang–Frenkel** (mpipi) pair potential plus Debye–Hückel
  electrostatics. There are **no native contacts** — the chain does not fold toward a
  target structure. Bonds are **flexible harmonic**.
- **Ribosome — rigid scenery.** The truncated CG 50S + tRNAs (~4,600 mass-0 beads) is
  fixed in space, but its **excluded-volume and electrostatic interactions with the
  nascent chain are on**. The tunnel axis is aligned with **+x** (the chain exits +x).
- **The PTC anchors.** Two ribosome beads are singled out as fixed reference points: the
  **P-anchor** (P-site tRNA residue-76 `R` bead) and the **A-anchor** (A-site tRNA
  residue-76 `R` bead) — where the peptidyl-tRNA (P) and incoming aminoacyl-tRNA (A) hold
  the chain's C-terminus.
- **C-terminus restraint — a harmonic position restraint.** The current C-terminal bead
  is restrained to one of the anchors with `U = k·|r − r₀|²`,
  `k = restraint_k = 83680 kJ/mol/nm²` (= 200 kcal/mol/Å²). **Switching the restraint
  target A→P is how translocation is reproduced.** (The `k` is a *per-particle* parameter
  so this force coexists with the tunnel wall, whose global constant is also `k`.)
- **Tunnel wall — a one-sided plane.** Because the 50S is **truncated** to a shell around
  the exit tunnel, there are no ribosome beads below the PTC. A one-sided half-harmonic
  wall `U = k·min(x − x₀, 0)²` (fixed stiffness 8368 kJ/mol/nm² = 20 kcal/mol/Å²) supplies
  the missing floor: the chain can only extrude **forward** (+x) and cannot slip below the
  synthesis point into the truncated region. **The plane `x₀` is auto-derived** from the
  ribosome structure (the lower of the two C-terminus hold planes,
  `min(A-target.x, P-target.x)`), so it can never go stale when you switch structures.
- **Thermostat.** Langevin dynamics at `ref_t = 300 K`, friction `tau_t`, timestep `dt`.

**Length-`L` model (no build-once-subset).** cosmo's forces are all sequence-local or
pairwise-by-type, so the length-`L` model is exactly
`buildCoarseGrainModel` on **residues 1..L** of the sequence — bonds, Yukawa, the
short-range HPS/mpipi term (and, for `hps_ss`, the local angle/torsion). Going `L → L+1`
just adds residue `L+1`'s terms. No STRIDE, no contact map, no matrix injection.

### 3. Ribosome ↔ nascent excluded volume: the O'Brien 12-10-6

The rigid ribosome interacts with the nascent chain through **excluded volume +
electrostatics only** — no attractive/native contacts. The excluded volume is O'Brien's
**12-10-6** form (inherited verbatim from topo):

```text
U = ε·[13(R/r)¹² − 18(R/r)¹⁰ + 4(R/r)⁶] ,   R = Rmin/2ᵢ + Rmin/2ⱼ  (sum rule)
ε = 0.000132 kcal/mol ,   cutoff 2.0 nm / switch 1.8 nm ,   {nascent}×{ribosome} only
```

The per-bead `Rmin/2` (O'Brien's structure-based CG collision radii) live in
`cosmo.parameters.model_parameters` under the **`hps_kr`** model — per-amino-acid values
plus the rRNA `P`/`R`/`BR` beads. This is why **`hps_kr` is the default CSP force field**:
it carries both the `Rmin/2` table for the ribosome wall *and* the Ashbaugh–Hatch
potential for the **nascent IDP↔IDP** interaction. Electrostatics fold into the existing
**Yukawa** force, extended over the ribosome charges (rRNA phosphate −1e, charged
residues) on `{nascent}×{nascent}` + `{nascent}×{ribosome}` (no intra-ribosome
electrostatics; the rigid ribosome's own interactions are constant and never computed).

### 4. The three stages: biology ↔ simulation

Each amino acid is added through **one elongation cycle**, split into three kinetic
sub-steps; CSP runs **one MD segment per sub-step**. For nascent length `L` each sub-stage
is a standalone short simulation (its own `L_<L>/stage_<s>/` folder); stage 3's final
structure seeds the next residue's stage 1.

| stage | real process | what the simulation does | C-terminus restrained to | mean dwell |
|-------|--------------|--------------------------|--------------------------|------------|
| **1** | Peptidyl transfer | new bead `L` **placed at the A-target**, bonded to `L−1`; minimize; run MD | **A-target** | `time_stage_1 = 0.34 ms` |
| **2** | Translocation (onset) | continue from stage 1, **still held at the A-target**; run MD (minimize skipped) | **A-target** | `time_stage_2 = 4.20 ms` |
| **3** | Translocation completes + wait | **switch the restraint A→P**, then run MD | **P-target** | remainder = (next codon total) − stage 1 − stage 2 |

```{note}
**Mechanics vs. timing.** The restraint switch (an instantaneous A→P geometric move)
happens at the **start of stage 3**, while the **duration** charged to translocation is
**stage 2**. Explicit A/P tRNA bonded geometry is not modelled, and the tRNA *tether* is
**forced off** for CSP — the switchable A↔P **position restraint** is what reproduces
translocation. The **timing** (three codon-resolved dwell times per residue) is faithful
to O'Brien; the per-stage **mechanics** are a reduced model.
```

### 5. From codon to MD steps (the kinetics)

The timing core is `cosmo.csp.kinetics` (pure Python, no OpenMM), **identical to topo's**.
For every residue it answers: *how many integration steps does each sub-stage run?*

**(a) Per-codon mean translation time.** The mRNA is split into codons; a lookup table
maps each codon to its **mean in-vivo translation time** in seconds — the codon's
intrinsic **mean first-passage time (mFPT)**, `τ(codon)`.

```{note}
**The codon-time table is organism-universal** (a property of organism + temperature, not
of the protein), so cosmo ships one: the **Fluitt *E. coli* table at 310 K** (61 sense + 3
stop codons, mean ≈ 0.068 s ≈ 15 aa/s) as `cosmo/csp/data/ecoli_trans_times_310K.txt`,
used **by default** whenever `csp.ini` gives no `codon_times` key. Set `codon_times` to a
table path only to override it. See `cosmo.csp.kinetics.default_codon_time_table_path`.
```

**(b) First-passage-time sampling.** Each stage is gated by a single rate-limiting event,
so its waiting time is **exponentially distributed**. Each sub-stage's dwell is drawn as
`t = −mean·ln(U)`, `U ∼ Uniform(0,1)`. A fixed `random_seed` makes the schedule
reproducible.

```{important}
**`time_stage_1` / `time_stage_2` (and the stage-3 remainder) are *means*, not the
per-residue dwell.** The actual time on each residue is a fresh exponential draw with that
mean, so individual residues get shorter/longer dwells and only their average equals the
configured value.
```

**(c) The three-stage split** for nascent length `L` (1-indexed). Draw three independent
`U₁, U₂, U₃`:

```text
t1(L) = −time_stage_1                                  · ln(U₁)   # peptidyl transfer
t2(L) = −time_stage_2                                  · ln(U₂)   # translocation
t3(L) = −(τ(next codon) − time_stage_1 − time_stage_2) · ln(U₃)   # wait to decode next codon
```

Stage 3 uses the *next* codon's mean (having just added residue `L`, the ribosome now
waits for `L+1`'s tRNA — which is why the codon list must extend to `L_max + 1`). If a
fast next codon makes the remainder ≤ 0 it is floored to `1e-9 s`.

**(d) In-vivo seconds → in-silico steps.** A time-compression factor maps seconds to steps:

```text
t_sim (ns) = t_s · 1e9 / scale_factor
n_steps    = t_sim (ns) / dt(ns) ,   dt(ns) = dt_ps · 1e-3
```

A **larger `scale_factor` ⇒ fewer steps per residue ⇒ a faster run**, while preserving the
*relative* timing of fast vs. slow codons. Step counts may additionally be clamped to
`[min_steps_per_stage, max_steps_per_stage]` for tractability — a clamp on **MD steps
only**; the sampled dwell **times in seconds** are recorded untouched in `dwell_times.dat`.

### 6. PTC geometry (always optimized)

The A-site seed / stage-1-2 restraint target and the P-site / stage-3 target are placed
by `optimal_ptc_targets` **exactly one peptide bond apart** (cosmo's CG bond length —
**0.380 nm** for `hps_kr`, *not* topo's 0.381 nm) and clear of the ribosome excluded
volume, by minimizing the soft O'Brien tRNA-bond/angle/improper restraints + the 12-10-6
wall over a deterministic multistart. Each new residue is delivered with its peptide bond
at equilibrium, which drops the stage-1 potential energy by ~50× versus seeding at the raw
`AtR`/`PtR`-76 `R` anchor beads (which sit ~0.9 nm apart, badly stretching the bond). The
C-terminus is held by a **position restraint** to these points — not a bond to the tRNA
beads. This is always on; there is no knob.

```{note}
cosmo **defaults to flexible harmonic bonds** (`constraints = None`) because backbone
flexibility is physically meaningful for disordered chains, but it **also supports rigid
bonds** (`constraints = AllBonds`) — every CA/P bond becomes a distance constraint pinned
at its equilibrium length, removing the fast bond-stretch mode so a larger `dt` can be
used. Constraints act only on the bonds; the **non-bonded** potentials are untouched.
So `AllBonds` does **not** by itself prevent the stiff-EV blow-up: the **non-native
excluded volume** is still stiff — the new residue is seeded at the fixed A-site target,
and a recently-added residue that has not cleared that region can land in the repulsive
core of the ribosome↔nascent 12-10-6 EV or the nascent Ashbaugh-Hatch potential (both
present in cosmo and topo). That diverges at the configured `dt` (PotE → ~1e13 kJ/mol),
which is why cosmo keeps topo's **per-stage dt-halving stability guard**: it re-runs the
stage at `dt/2` with `2×` steps (identical dwell) until it is stable. The always-on PTC
optimization reduces how often this fires (equilibrium seeding lowers the stage-1 energy
~50×), but it is a *quality* improvement, not a full substitute for the guard.
```

### 7. After the last residue: ejection (and dissociation)

Once the final residue is added, the simulation runs a **post-synthesis ejection phase**
(`ejection_steps`): the C-terminus restraint is **released** while the rigid ribosome and
tunnel wall remain, so the chain diffuses out along +x. An optional **dissociation** phase
(`dissociation_steps`) continues the free protein away from the ribosome.

---

## Configuration reference (`csp.ini`)

CSP reads a single INI control file with one `[OPTIONS]` section
(`cosmo.csp.protocol.read_csp_config`). **Units are OpenMM defaults** — nm, ps, kJ/mol, K,
kJ/mol/nm² — and **dwell times are in seconds**. Integers may use `_` separators.

```{tip}
For a compact tabular reference of every `csp.ini` option, see {doc}`synthesis_control`.
```

### Inputs & schedule

| Key | Required | Default | Meaning |
|-----|----------|---------|---------|
| `pdb_file` | **yes** | — | All-atom / CA native PDB of the target protein; the CG model is built from its first `L` residues. |
| `ribosome` | **yes** | — | Truncated CG ribosome PDB (P-/A-anchors + rigid scenery). |
| `model` | no | `hps_kr` | Nascent force field. **`hps_kr`** carries the `Rmin/2` table for the ribosome 12-10-6 wall; other models fall back to `hps_kr` for those radii only. |
| `L0` | no | `1` | Start nascent-chain length. |
| `L_max` | no | full length | Final nascent length. |
| `mrna` | cond. | — | mRNA file (one codon per residue). Required for per-codon timing (unless `codon_times` is a number). |
| `codon_times` | no | bundled E. coli 310 K | A **table path** = per-codon timing; a **positive number of seconds** = uniform codon time (no `mrna` needed); omit = bundled Fluitt table. A table filename must **not** be a bare number. |
| `outdir` | no | `synth_out` | Output root. |

There is **no `domain_def`, `stride_output_file`, or `nascent_ev_radii` key** — those are
topo's Gō-model inputs and do not apply to cosmo's sequence-based chain.

### O'Brien kinetics

| Key | Default | Meaning |
|-----|---------|---------|
| `scale_factor` | `4331293` | In-vivo-s → in-silico-ns compression (larger = fewer steps = faster). |
| `time_stage_1` | `0.00034` | Mean peptidyl-transfer dwell, **seconds**. |
| `time_stage_2` | `0.004201` | Mean translocation dwell, **seconds**. |
| `random_seed` | — | Seed for the FPT sampler. |
| `max_steps_per_stage` | — (uncapped) | **Testing only** — upper clamp on each stage's step count. |
| `min_steps_per_stage` | `1` | **Testing only** — lower clamp. |
| `ejection_steps` | `0` | Post-synthesis ejection phase (steps); `0` = skip. |
| `dissociation_steps` | `0` | Post-synthesis dissociation phase (steps); `0` = skip. |

```{warning}
**`max_steps_per_stage` / `min_steps_per_stage` are testing-only.** They clamp the MD step
count so examples finish quickly, which **breaks the physical timescale mapping**. Leave
them **unset** in production so step counts come entirely from the kinetics. The sampled
dwell **times in seconds** are always written to `dwell_times.dat` regardless.
```

### MD / ribosome mechanics (`RunParams` fields)

| Key | Default | Meaning |
|-----|---------|---------|
| `dt` | `0.01` | Timestep, ps. |
| `ref_t` | `300` | Temperature, K. |
| `tau_t` | `0.01` | Langevin friction, 1/ps. |
| `nstout` | `50` | Trajectory/log output interval (steps). |
| `device` | `CPU` | `GPU` / `CPU`. |
| `ppn` | `1` | CPU threads (CPU platform). |
| `constraints` | `None` | Bond treatment: `None` (flexible harmonic bonds, default) or `AllBonds` (rigid distance constraints — larger-timestep path). |
| `restraint_k` | `83680` | C-terminus harmonic restraint constant, kJ/mol/nm². |
| `minimize` | `yes` | Energy-minimize the seeded structure before each stage's MD. |
| `tunnel_wall` | `yes` | One-sided tunnel wall (floor below the synthesis point); plane auto-placed. |

There is **no `rigid_ribosome` key** (supplying the `ribosome` PDB *is* the signal to load
it as rigid scenery), output is **always nascent-only**, and `trna_tether` is **forced
off** by the CSP runner (CSP needs the switchable A↔P position restraint).

---

## Outputs

```text
<outdir>/
├── L_<L>/stage_<s>/        # one folder per residue L and sub-stage s ∈ {1,2,3}
│   ├── traj.dcd            # (nascent-only) trajectory for that stage
│   ├── traj_final.pdb      # last conformation (seeds the next stage/residue)
│   ├── traj.log            # energies
│   └── traj.psf, traj.chk, traj_runinfo.log, native_1_<L>.pdb
├── ejection/               # post-synthesis ejection phase (if ejection_steps > 0)
├── dissociation/           # post-synthesis free run (if dissociation_steps > 0)
└── dwell_times.dat         # per-residue dwell-time log
```

**`dwell_times.dat`** records, per residue, the codon, the three sampled dwell **times in
seconds** (`t1`/`t2`/`t3`), their nanosecond equivalents, and the integer MD step counts —
the physical schedule, independent of any step clamp.

### Console progress log

`cosmo-csp` prints one compact line per residue plus one per sub-stage, each reporting the
wall-clock time and the **total system potential energy of the last integrated step**:

```text
L=  5  uniform  dwell      0.05 s  steps   40/  40/  40
  L=  5  stage 1 peptidyl-transfer       40 steps    0.15 s  PE=  +6.7544e+01 kJ/mol
  L=  5  stage 2 translocation           40 steps    0.11 s  PE=  +5.9xxxe+01 kJ/mol
  L=  5  stage 3 tRNA-binding            40 steps    0.14 s  PE=  +5.6842e+01 kJ/mol
```

Set **`COSMO_CSP_VERBOSE=1`** to restore the full per-stage banners (build block,
minimization, run metadata). MDAnalysis' cosmetic warnings are silenced for the run.

**Movie.** Each stage writes a standalone trajectory; `cosmo-csp-movie` stitches them — in
synthesis order, padding every frame to the final length and overlaying the static
ribosome — into one VMD-playable movie (auto-detects the 3-stage vs flat layout):

```bash
cosmo-csp-movie -o <outdir> --ribosome ribosome_trunc.pdb
vmd -e <outdir>/movie.tcl
```

---

## Python API

```python
from cosmo.csp.protocol import run_continuous_synthesis, read_csp_config

# (a) drive it from an INI, exactly like the CLI:
cfg = read_csp_config("csp.ini")
run_continuous_synthesis(
    cfg.pdb_file, cfg.ribosome,
    L0=cfg.L0, L_max=cfg.L_max, out_root=cfg.outdir,
    mrna=cfg.mrna, codon_time_table_path=cfg.codon_time_table_path,
    params=cfg.params,
)

# (b) or construct parameters directly (the ribosome PDB is always rigid scenery;
#     the tunnel-wall plane is auto-derived from it):
from cosmo.csp.core import RunParams
params = RunParams(model="hps_kr",
                   scale_factor=4331293.0, random_seed=1, ejection_steps=50000)
run_continuous_synthesis("asyn.pdb", "ribosome_trunc.pdb",
                         L0=5, L_max=10, mrna="mrna.txt", params=params)
```

See the {doc}`API reference <../cosmo.csp>` for the autodocumented modules.

---

## See also

- {doc}`cylinder_synthesis` — the analytic-tunnel variant (`cosmo-cylinder`): the same
  codon kinetics with a nascent-only system and a single MD segment per residue.
- {doc}`synthesis_control` — the concise `csp.ini` control-options reference.
- {doc}`../cosmo.csp` — the API reference (`cosmo.csp.core`, `cosmo.csp.ribosome`,
  `cosmo.csp.kinetics`).
