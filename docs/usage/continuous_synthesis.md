# Synthesis on a coarse-grained ribosome

The **Continuous Synthesis Protocol (CSP)** is cosmo's **codon-resolved, kinetic** runner
for protein synthesis on an **explicit coarse-grained ribosome**. It times
**every residue from its mRNA codon** and splits each elongation cycle into **three
kinetic sub-stages** — reproducing O'Brien's `continuous_synthesis_v6.py` protocol, but
on cosmo's **sequence-based IDP force field** (HPS / mpipi) instead of a structure-based
Gō model. (Codon-resolved kinetics are what make the model physically meaningful; a
fixed per-residue step count is not.) For the analytic-tunnel variant — the same codon
kinetics with the explicit ribosome replaced by a cylindrical bore — see
{doc}`cylinder_synthesis`.

- **CLI:** `cosmo-csp -f csp.ini` (or `python -m cosmo.csp -f csp.ini`)
- **Movie tool:** `cosmo-csp-movie -o <out_root> [--ribosome ribo.pdb]`
- **Worked example:** `tutorials/08_csp_cg_ribosome/` — a tutorial-scaled run on
  α-synuclein grown on the real *E. coli* 4V9D CG ribosome (the analytic-tunnel
  variant is `tutorials/07_csp_cylinder/`). A larger production configuration lives
  in `sandbox/validate/` and `sandbox/Ecoli/`.
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
`cosmo.models.buildCoarseGrainModel` on the **first `L` residues** of the sequence.
```

---

## Quick start

All paths in the INI are relative to the working directory; run from the example folder.
A GPU is recommended for production (the explicit-ribosome system has ~4,600 rigid beads).

```bash
cd tutorials/08_csp_cg_ribosome
cosmo-csp -f csp.ini              # 3-stage synthesis -> synth_out_csp/

# stitch the per-stage trajectories into one VMD movie
cosmo-csp-movie -o synth_out_csp --ribosome 4v9d_50S_PtR_5jte_AtR_model_cg_trunc.pdb
cd synth_out_csp && vmd -e movie.tcl        # movie.tcl loads its files by basename
```

`cosmo-csp` writes, per residue `L` and sub-stage `s`, a standalone trajectory under
`<outdir>/L_<L>/` (one folder per residue, per-stage `traj_s<s>.dcd`), an optional `ejection/` phase, and a
per-residue dwell-time log `<outdir>/dwell_times.dat`.

---

## Theory

### 1. What is being modeled: protein synthesis

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
- **C-terminus restraint — position restraint (default) or O'Brien tRNA tether.** The
  current C-terminal bead is held at the A-/P-site either by a harmonic **position
  restraint** to the target *point* (`U = k·|r − r₀|²`, `k = restraint_k = 83680
  kJ/mol/nm²` = 200 kcal/mol/Å²; default) or, with `trna_tether = yes`, by the full
  **O'Brien tRNA tether** — a bond + two orienting angles + an improper to the A76 tRNA
  beads (plus a backbone orienting angle `prev–N–R` — the `hps_ss` bistable backbone
  potential, applied across the peptide–tRNA junction for **all** models), which also
  fixes the chain's *orientation*. **Switching the A-/P-site hold A→P is
  how translocation is reproduced** — the position restraint moves its target point; the
  tether re-attaches to the P-site beads in stage 3. (The `k` is a *per-particle* parameter
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

The per-bead `Rmin/2` (O'Brien's structure-based CG collision radii) are
**model-independent** steric radii: they live in the standalone `OBRIEN_RMIN_2_NM`
(per-amino-acid) and `OBRIEN_RNA_RMIN_2_BEADS` (rRNA `P`/`R`/`BR`) tables in
`cosmo.parameters.model_parameters`, **decoupled from any force field**. So the ribosome
wall is identical for **every** nascent model — CSP runs on `hps_kr`, `hps_urry` or
`mpipi` alike (`hps_kr` is merely the default); the selected model only sets the
**nascent IDP↔IDP** interaction (Ashbaugh–Hatch or Wang–Frenkel). Electrostatics fold into the existing
**Yukawa** force, extended over the ribosome charges (rRNA phosphate −1e, charged
residues) on `{nascent}×{nascent}` + `{nascent}×{ribosome}` (no intra-ribosome
electrostatics; the rigid ribosome's own interactions are constant and never computed).

### 4. The three stages: biology ↔ simulation

Each amino acid is added through **one elongation cycle**, split into three kinetic
sub-steps; CSP runs **one MD segment per sub-step**. For nascent length `L` each sub-stage
is a short simulation writing its own `traj_s<s>.dcd` in the shared `L_<L>/` folder; stage 3's final
structure seeds the next residue's stage 1.

| stage | real process | what the simulation does | C-terminus restrained to | mean dwell |
|-------|--------------|--------------------------|--------------------------|------------|
| **1** | Peptidyl transfer | new bead `L` **placed at the A-target**, bonded to `L−1`; minimize; run MD | **A-target** | `time_stage_1 = 0.34 ms` |
| **2** | Translocation (onset) | continue from stage 1, **still held at the A-target**; run MD (minimize skipped) | **A-target** | `time_stage_2 = 4.20 ms` |
| **3** | Translocation completes + wait | **switch the restraint A→P**, then run MD | **P-target** | remainder = (next codon total) − stage 1 − stage 2 |

```{note}
**Mechanics vs. timing.** The restraint switch (an instantaneous A→P geometric move)
happens at the **start of stage 3**, while the **duration** charged to translocation is
**stage 2**. With the default position restraint only the C-terminus *point* is held;
explicit A/P tRNA bonded geometry (bond + angles + improper to the A-/P-site tRNA beads,
switched A→P across the stages) **is** modelled when `trna_tether = yes`. The **timing**
(three codon-resolved dwell times per residue) is faithful to O'Brien; the per-stage
**mechanics** are a reduced model.
```

### 4a. C-terminus restraint: position restraint vs. tRNA tether

The C-terminus must be held at the PTC and translocated A→P each residue. cosmo offers two
mechanisms, chosen by the `trna_tether` key (see {doc}`synthesis_control`); both drive the
*same* A→P schedule across the three stages, but the tether additionally controls the
chain's **orientation**.

**Position restraint (`trna_tether = no`, default).** A single harmonic spring pins the
C-terminal bead to the A/P target *point* (the equilibrium-optimized A- or P-site location),
`U = k·|r − r₀|²`. Stages 1–2 target the A-point, stage 3 targets the P-point; the A→P
switch is just changing `r₀`. This is the validated default — simplest, and it needs no
particular tRNA beads.

```{figure} img/csp_restraint_position.svg
:alt: A single harmonic spring holds the C-terminus at the A- or P-site target point; sites read E, P, A left to right.
:width: 560px
:align: center

**Position restraint** (`trna_tether = no`). A single harmonic spring pins the C-terminal
bead to the A- or P-site target *point*; translocation is just moving the target `r₀` from
the A-point (stages 1–2) to the P-point (stage 3) — right→left in the E–P–A frame.
```

**O'Brien tRNA tether (`trna_tether = yes`).** Reproduces the covalent attachment of the
nascent C-terminus to the tRNA with the full bonded geometry — not a point restraint but a
set of bonded terms to the **tRNA's 3′-terminal nucleotide, residue 76 (A76)**.

*Which residue?* The truncated ribosome keeps only the acceptor-stem 3′ end of each tRNA
(a handful of residues: `AtR` 73–76, `PtR` 72–76), each nucleotide represented by a few CG
beads — `P` (phosphate), `R` (ribose), `BR1`/`BR2` (base). **The tether attaches to residue
76 only** — the 3′-terminal adenosine A76, which is exactly the CCA-3′ site where the amino
acid esterifies to the tRNA. The other tRNA residues are rigid excluded-volume scenery,
*not* part of the tether. It uses three beads *of A76*:

- `R` — the **ribose** of A76: the C-terminus **bonds** to it.
- `P` — the **phosphate** of A76: used in the orienting angle `N–R–P`.
- `BR2` — the second **base** bead of A76 (present only because A76 is a purine): used in
  `N–R–BR2` and the improper.

For the current C-terminus `N` (the newest residue) it adds:

| term | force | connectivity (beads of A76) | A-site (`AtR`) | P-site (`PtR`) |
|------|-------|-----------------------------|----------------|----------------|
| bond | harmonic | `N → R` | 0.427 nm | 0.476 nm |
| angle | harmonic | `N–R–P` | 106° | 117° |
| angle | harmonic | `N–R–BR2` | 127° | 130° |
| improper | periodic (`CustomTorsion`) | `N–R–P–BR2` | 128° | −161° |
| backbone | double-Gaussian angle | `prev–N–R` | the `hps_ss` backbone-angle potential (bistable, θ ≈ 92°/130°) across the junction; added for **all models** — a dedicated copy for the angle-free `hps_kr`/`hps_urry`/`mpipi` |

(bond/angle stiffness = 200 kcal/mol/Å²; angle/improper = 25 kcal/mol/rad².) The two
orienting angles + the improper fix the residue's **bearing in the A76 frame** — this is
what the plain point restraint cannot do.

```{figure} img/csp_restraint_tether_geom.svg
:alt: The C-terminus bonds to the ribose R of tRNA A76, with orienting angles to the phosphate P and base BR2, plus an improper dihedral.
:width: 560px
:align: center

**tRNA tether geometry** (`trna_tether = yes`). The C-terminus (`N`) **bonds** to the
**ribose (`R`)** of A76; two orienting angles (`N–R–P`, `N–R–BR2`) and the improper
(`N–R–P–BR2`) fix its bearing in the A76 frame, and a backbone angle (`prev–N–R`) keeps the
terminal segment physically oriented (the `hps_ss` double-Gaussian backbone potential, added
for all models). Values shown are the P-site (`PtR`) set;
A-site (`AtR`) values are in the table above. Only A76 is tethered — the rest of the tRNA
is rigid scenery.
```

**Per-stage tethering (tRNA-tether mode).** Each stage rebuilds the length-`L` system and
re-attaches the tether to the appropriate site — there is **no sliding spring**; the A→P
"switch" is realized by the tether referencing the P-site beads (`PtR`) instead of the
A-site beads (`AtR`):

| stage | C-terminus (residue `L`) | previous residue (`L−1`) |
|-------|--------------------------|--------------------------|
| **1** | A-site (`AtR`) | **P-site (`PtR`)** |
| **2** | A-site (`AtR`) | — free — |
| **3** | **P-site (`PtR`)** | — free — |

- **Stage 1** double-tethers *both* ends of the freshly-formed peptide bond (the new
  residue at A, the previous one at P), so the new bond starts at its equilibrium PTC
  geometry the instant the residue appears.
- **Stage 2** keeps the new residue at A; the previous residue is released (held only by
  the nascent-chain backbone + the tunnel wall).
- **Stage 3** re-tethers the new residue A→P — starting from stage 2's coordinates, the
  minimize + MD carries it from the A-site to the P-site: **this is the translocation.**
- **Cold start** (the first residue `L0`, no previous residue): tethered to the P-site in
  all three stages.

```{figure} img/csp_restraint_ap_stages.svg
:alt: tRNA sites drawn E, P, A left to right; across the three stages the new residue moves from the A-site to the P-site, right to left.
:width: 600px
:align: center

**The A→P switch across the three stages.** Each stage rebuilds the system and re-attaches
the tether (no sliding spring). Stage 1 double-tethers both ends of the new peptide bond
(new residue `N` at A, previous at P); stage 2 keeps `N` at A (previous released); stage 3
re-tethers `N` to the P-site — the minimize + MD carry it A→P (right→left), **this is the
translocation.** Sites are drawn E–P–A per convention; CSP models only the P- and A-site
tRNAs (no E-site tRNA).
```

```{note}
The tRNA tether requires a **well-formed A/P tRNA** in the ribosome PDB (segids `AtR`/`PtR`,
resid 76, beads `R`/`P`/`BR2`; the acceptor must be a purine, which carries `BR2`). A site
missing its `P` or `BR2` bead simply skips that angle/improper. The A/P target *points*
(used by both modes, and for the tunnel-wall plane) come from the always-on PTC-geometry
optimization.
```

### 5. From codon to MD steps (the kinetics)

The timing core is `cosmo.csp.kinetics` (pure Python, no OpenMM), **identical to topo's**.
For every residue it answers: *how many integration steps does each sub-stage run?*

**(a) Per-codon mean translation time.** The mRNA is split into codons; a lookup table
maps each codon to its **mean in-vivo translation time** in seconds — the codon's
intrinsic **mean first-passage time (mFPT)**, `τ(codon)`.

```{note}
**The codon-time table is an organism + temperature property** (not of the protein), so a
per-codon run needs an explicit `codon_times` table path — **there is no bundled default**.
The shared library under `assets/csp/codon_dwell_times/<organism>/` holds tables per
organism, including the **Fluitt *E. coli* table at 310 K** (Fluitt, Pienaar & Viljoen,
*Comput. Biol. Chem.* 2007; 61 sense + 3 stop codons, codon-usage-weighted mean ≈ 0.061 s
≈ 16.5 aa/s) at `assets/csp/codon_dwell_times/ecoli/ecoli_codon_dwell_times_310K.txt`. Both
`mrna` (the protein-specific sequence) and `codon_times` (the table) are mandatory for
per-codon timing; see {doc}`codon_dwell_times`.
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
`AtR`/`PtR`-76 `R` anchor beads (which sit ~0.9 nm apart, badly stretching the bond). By default the
C-terminus is held by a **position restraint** to these points; with `trna_tether = yes`
it is instead bonded to the A-/P-site tRNA beads (the full O'Brien tether). The PTC-target
optimization itself is always on; there is no knob for it.

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

### 7. After the last residue: ejection

Once the final residue is added, the simulation runs a **post-synthesis ejection phase**
(`ejection_steps`): the C-terminus restraint is **released** while the rigid ribosome and
tunnel wall remain, so the chain diffuses out along +x.

---

## Configuration

CSP reads a single INI control file with one `[OPTIONS]` section
(`cosmo.csp.protocol.read_csp_config`). **Every control key — for both `cosmo-csp` and
`cosmo-cylinder` — is documented in one place:** {doc}`synthesis_control` (grouped into
*shared*, *coarse-grained-ribosome-only*, and *cylinder-only* keys, with types and
defaults). Units are OpenMM defaults (nm, ps, kJ/mol, K, kJ/mol/nm²) and dwell times are
in seconds.

A few CSP-specific behaviors worth calling out (the physics behind the keys):

- **`pdb_file` may be all-atom or a Cα-only CG structure** — only the Cα positions and
  residue names are read (cosmo has no STRIDE / native-contact step), so there is **no**
  `domain_def` / `stride_output_file` / `nascent_ev_radii` key (those are topo's Gō-model
  inputs).
- **PTC geometry is always optimized** (A/P targets one peptide bond apart, EV-clear), so
  each new residue is delivered with its bond at equilibrium — no knob.
- **The ribosome is always rigid scenery** (supplying the `ribosome` PDB *is* the signal;
  no `rigid_ribosome` key) and output is always nascent-only. `trna_tether` defaults off
  (the switchable A↔P position restraint); set `trna_tether = yes` for the full O'Brien
  tRNA tether, which switches A→P by re-attaching to the P-site beads in stage 3.

---

## Outputs

```text
<outdir>/
├── L_<L>/                  # ONE folder per residue L (consolidated layout)
│   ├── traj.psf            # nascent topology (shared across the 3 stages; f(L) only)
│   ├── native_1_<L>.pdb    # length-L native structure (shared)
│   ├── traj_s1.dcd         # (nascent-only) trajectory, stage 1  (s2/s3 likewise)
│   ├── traj_s1.log         # energies, stage 1
│   ├── traj_runinfo.log    # folded run-info: one [run:...]/[result:...] per stage
│   └── traj_final.pdb      # stage-3 final — seeds L+1 and is the resume-reload target
├── ejection/               # post-synthesis ejection phase (if ejection_steps > 0)
├── dwell_times.dat         # per-residue dwell-time log / schedule (#PTC header)
└── progress.log            # append-only DONE/RUNNING resume status
```

Each residue's three sub-stages share **one** `L_<L>/` directory: `traj.psf` and
`native_1_<L>.pdb` depend only on `L` and are written once; trajectories stay split per
stage (`traj_s{1,2,3}.dcd`); only stage 3 writes `traj_final.pdb`. There is no per-stage
`.chk` — per-residue **resume** reloads `traj_final.pdb`. Re-invoking `cosmo-csp` on an
interrupted `<outdir>` continues from the last completed residue (`resume = auto` by
default; `--fresh` to override); see `cosmo.csp.resume`.

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
cosmo-csp-movie -o <outdir> --ribosome 4v9d_50S_PtR_5jte_AtR_model_cg_trunc.pdb
cd <outdir> && vmd -e movie.tcl        # movie.tcl loads its files by basename
```

For the cylinder-tunnel equivalent and more detail, see {doc}`synthesis_visualization`.

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
params = RunParams(model="hps_urry",   # any IDP model works; hps_kr is the default
                   scale_factor=4331293.0, random_seed=1, ejection_steps=50000)
run_continuous_synthesis("asyn.pdb", "4v9d_50S_PtR_5jte_AtR_model_cg_trunc.pdb",
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
