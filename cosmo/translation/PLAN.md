# Plan: `cosmo.translation` — co-translational synthesis runner

A new package `cosmo/translation/`, sibling of `cosmo/mdrun/`, driven by
`cosmo.read_simulation_config` + `cosmo.engine`. It mirrors the sibling project's
`topo/translation/` (`~/work/src/topo/topo/translation/`, see its `DESIGN.md` /
`README.md` / `FILES.md`): architecture follows topo; physics is cosmo's IDP
force field (HPS/mpipi), not topo's Gō model.


## Decisions (settled)

- Config key `rna_model = cosmo | topo`; **default `topo`** (3/4-bead P/R/BR).
- Ribosome source structure: **`4v9d_50S_PtR_5jte_AtR_model.pdb`** (topo's
  X-aligned 50S + P-/A-tRNA), copied into `cosmo/translation/structures/`.
- `P`/`R`/`BR` bead params kept **local** to `cosmo/translation/` (scenery-only).
- Nascent chain is **always protein**; `rna_model` governs the ribosome scenery only.
- Nascent-chain force-field default **`hps_urry`** (`hps_kr`/`hps_ss`/`mpipi` also
  selectable; `mpipi` triggers the §2.1 dummy-id path).
- First test inputs: **`tutorials/07_translation/asyn.pdb` (α-synuclein), `L0 = 5`**.
- **No coordinate shifting** — drive every length through `cfg.init_position`
  (see §4).
- **Status:** v1 (nascent-only loop) and v2 (rigid ribosome: excluded volume +
  Yukawa, both `rna_model` reps; O'Brien tRNA tether + one-sided tunnel wall) are
  implemented and validated. **Deferred:** ejection/stallation post-phase and the
  analysis layer — see `tutorials/07_translation/TASKS.md`.
- **Tunnel-wall plane is auto:** `tunnel_wall_x0` defaults to the C-terminal-AA
  addition plane (PTC) = **P-anchor x + tether bond length** (`TRNA_TETHER_BOND_NM`
  = 0.476 nm) — i.e. where each new residue is restrained/added; an explicit
  `tunnel_wall_x0` only overrides it.


## 0. What this replaces

- **`examples/growing/`** — the legacy ad-hoc prototype (`grow_nascent*.py`:
  whole-complex rebuild, single restraint at origin, planar wall) has been
  **removed**; this package replaces it with a config/engine-driven runner.
- **`cosmo/utils/crop_ribosome.py`** (`filter_pdb`: x<0 drop, cylinder crop, x≥60
  keep) — fold into `cosmo/translation/truncate_ribosome.py`.


## 1. The model

A protein is synthesized vectorially (N→C, one residue at a time); the nascent
chain emerges through the ribosomal exit tunnel. Modeled as a **sequence of
standalone simulations**, one per length `L = L0 .. N_full`, each rebuilding a
one-residue-longer system and continuing from the previous step's final
coordinates (no live `Context` resizing). The nascent chain uses cosmo's IDP force
field, so this studies co-translational behavior of disordered chains.


## 2. cosmo is sequence-based, not structure-based

The core simplification vs. topo — cosmo has **no STRIDE and no Gō native
contacts**, which deletes topo's most intricate machinery:

| | topo (Gō) | cosmo (HPS/mpipi) |
|---|---|---|
| Per-residue interactions | native 12-10-6 contacts from structure | Ashbaugh–Hatch (HPS) / Wang–Frenkel (mpipi), from sequence |
| STRIDE / native contact map | required | **none** |
| Build-once-subset `L×L` block | required (heart of `elongate.py`) | **not needed** |
| Length-`L` model | top-left `L×L` of precomputed matrices | **first `L` residues of the sequence** |

A length-`L` build is the ordinary `models.buildCoarseGrainModel` on the first `L`
residues — bonds, Yukawa, the short-range non-bonded (Ashbaugh–Hatch /
Wang–Frenkel), and (`hps_ss`) the local Gaussian angle/torsion terms. All cosmo
forces are sequence-local or pairwise-by-type, so restricting to residues `1..L` is
exact and needs no masking. There is no `precompute_contacts` and no matrix
injection. Ribosome↔nascent stays excluded-volume + electrostatics only — faithful
to topo's v2 without any Gō machinery.

### 2.1 `mpipi` exception — tabulated Wang–Frenkel is id-indexed

The HPS models use continuous per-particle σ/λ, but **`mpipi`'s WF potential is a
tabulated `CustomNonbondedForce` indexed by per-particle residue `id`** (0–19
protein, 20–23 RNA; `Discrete2DFunction` lookups). So under `mpipi` only, v2 must:

- give each appended ribosome bead an **in-range dummy `id` (e.g. 0)** in the WF
  force (OpenMM requires one `addParticle` per System particle), and restrict the
  WF force to interaction group `{nascent}×{nascent}` so dummy ids are never read;
- keep nascent ordering `0..L-1` so each nascent particle's `id` aligns with its
  residue type (already true — nascent built first, ribosome appended at `L..N-1`).

The ribosome's real σ/charge come from the separate cross-force + Yukawa (§6), not
the WF table. Under the default `hps_urry` this path is skipped entirely.


## 3. Package layout

```
cosmo/translation/
  __init__.py            # re-export run_elongation, ElongationParams,
                         #   read_elongate_config, ElongateConfig
  __main__.py            # python -m cosmo.translation -> elongate CLI
  elongate.py            # elongation runner (v1 + v2)
  ribosome.py            # rigid ribosome scenery + cross-forces (v2)
  cg_ribosome.py         # all-atom -> CG ribosome (both rna_model reps)
  truncate_ribosome.py   # crop CG ribosome around the tunnel (subsumes crop_ribosome.py)
  make_movie.py          # stitch per-length DCDs -> growing movie
  structures/            # 4v9d_50S_PtR_5jte_AtR_model.pdb + CG/truncated outputs
```

Console entry points (`pyproject.toml`): `cosmo-elongate`, `cosmo-elongate-movie`.

Module → topo source for the implementer:

| cosmo module | topo source | notes |
|---|---|---|
| `elongate.py` | `topo/translation/elongate.py` | drop `precompute_contacts` / matrix injection (§2); keep loop, anchors, placement, restraint, reporters |
| `ribosome.py` | `topo/translation/ribosome.py` | add dual-`rna_model` bead handling (§5) |
| `cg_ribosome.py` | `topo/translation/cg_ribosome.py` | add the 1-bead RNA mapping (§5) |
| `truncate_ribosome.py` | `topo/translation/truncate_ribosome.py` | subsumes `cosmo/utils/crop_ribosome.py` |
| `make_movie.py` | `topo/translation/make_movie.py` | near-verbatim |


## 4. Reuse from cosmo

The runner is an outer loop over length `L` calling cosmo's helpers each iteration:

- **`cosmo.engine`** — `build_system` / `setup_simulation` / `attach_reporters` /
  `finalize_simulation`; one build+run+finalize per length.
- **`cosmo.read_simulation_config` / `SimulationConfig`** — base run options.
  Translation-specific keys live in an `ElongateConfig` wrapper (mirror topo's
  `read_elongate_config`), keeping the shared config unbloated.
- **`cosmo.reporter.cosmoReporter` + `cosmo.runinfo`** — logging + provenance, one
  set per length.
- **`init_position` as-is coordinates (no shift).** The positive-octant shift
  applies only to the fresh `pdb_file` path; the `init_position` branch loads
  coordinates as-is (`engine.py:176-185`). The runner drives **every** length
  (including the L0 cold-start) through `cfg.init_position = seed.pdb`, which
  preserves the absolute tunnel/anchor/ribosome geometry. `shift_positions` is moot
  on this path.


## 5. RNA representation (`rna_model`)

The rigid ribosome's rRNA can be represented two ways. This is sound because
ribosome↔nascent interactions are electrostatics + excluded volume only (no
attractive contacts to RNA), so bead resolution changes only where the −1e
phosphate charges sit and the excluded-volume σ granularity — never the
interaction type.

| | `rna_model = topo` (default) | `rna_model = cosmo` |
|---|---|---|
| beads per nucleotide | 3 pyrimidine / 4 purine (`P`, `R`, `BR1`, `BR2`) | 1 (phosphate `P`) |
| charge | −1e on `P`; `R`/`BR` neutral | −1e (−0.75e under mpipi) |
| excluded-volume σ | 0.71 nm per P/R/BR (O'Brien radii) | cosmo nucleotide radius (~0.82–0.85 nm) |
| param source | local table (§ below) | cosmo `hps_kr`/mpipi RNA entries (present) |
| bead count | matches O'Brien crop (~4,576) | ~3× fewer |

Implementation:

1. **`cg_ribosome.py`** produces the chosen representation: protein → Cα; RNA →
   P/R/BR centroids (`topo`, port topo's mapping) or one phosphate `P` bead per
   nucleotide (`cosmo`). Preserve `segID` (AtR/PtR/23S/5S/L*) in both so anchors and
   selections work.
2. **`ribosome.py`** loads the truncated CG ribosome and looks up charge + σ per
   bead from the representation's table. The cross-force code is identical for both:
   one excluded-volume `CustomNonbondedForce` `{nascent}×{ribosome}` + Yukawa
   extended over ribosome charges.
3. **`P`/`R`/`BR` params** (mass, charge, radius 0.71 nm) live in a local
   ribosome-bead table in `cosmo/translation/` — scenery-only, never typed into the
   HPS/mpipi pair tables.

Anchors are representation-independent: the P-/A-anchors are the PtR/AtR residue-76
ribose position — the `R` bead under `topo`, the single bead under `cosmo`. The
anchor reader keys off `segID` + resid 76.


## 6. Build steps

### v1 — elongation loop, no ribosome forces
System = nascent chain only; the truncated ribosome supplies two fixed points:
P-anchor (PtR resid-76) and A-anchor (AtR resid-76). Loop `L = L0 .. N_full`:

1. Build length-`L` model on the first `L` residues (§2).
2. Seed coordinates: `L == L0` → extended chain along the tunnel axis (+x) from the
   P-anchor, one bond length apart, with a tiny ±perpendicular zig-zag to avoid the
   180° angle singularity; `L > L0` → residues `1..L-1` from the previous final
   structure, residue `L` at the A-anchor + buffer.
3. Harmonic restraint on residue `L` only → P-anchor (`CustomExternalForce`,
   `k ≈ 200 kcal/mol/Å² = 83680 kJ/mol/nm²`); the new residue migrates A→P.
4. Boltzmann velocities, minimize, run `n_steps`, save final → seed `L+1`.
5. Per-length output folder `<outdir>/L_<L>/` (traj/log/psf/chk/final/seed).

Bonds: flexible harmonic — cosmo's CA/P default, matching topo's required v1 setting
for free (the new bond is seeded ~1 nm off equilibrium; harmonic + minimize absorb
it).

Acceptance: runs without crashing; chain grows N→C one residue/step; new residue
migrates A→P; per-length files written.

### v2 — rigid ribosome scenery
Append the truncated CG ribosome at indices `L..N-1`, mass = 0 (frozen, coords
as-is). Cross-interactions:

- **Excluded volume:** dedicated `CustomNonbondedForce`, pure `ε·(σ_ij/r)¹²`
  (`ε = 0.000132 kcal/mol`, σ from the bead radii of the chosen rep + Cα radii),
  group `{nascent}×{ribosome}`, cutoff 2.0 nm / switch 1.8 nm.
- **Electrostatics:** extend cosmo's Yukawa over ribosome charges (rRNA phosphates
  −1e, charged residues), groups `{nascent}×{nascent}` + `{nascent}×{ribosome}`, no
  intra-ribosome electrostatics.
- **`mpipi` only:** dummy in-range `id` per ribosome bead + `{nascent}×{nascent}`
  group on the WF force (§2.1).
- **Exclusions:** OpenMM CPU requires all `CustomNonbondedForce`s to share one
  exclusion list — the ribosome–NC force copies the nascent bonded exclusions
  (verify against cosmo's existing exclusion handling).

Acceptance: stable runs, ribosome frozen (≈0 drift), nascent chain stays in/around
the tunnel and extrudes +x as it grows — under both `rna_model` settings.


## 7. Control file (`elongate.ini`)

`[OPTIONS]`, parsed by `read_elongate_config` wrapping `read_simulation_config`.
Required: `pdb_file` (nascent), `ribosome`, `L0`. Mirror topo's keys (`L_max`,
`n_steps`, `dt`, `ref_t`, `tau_t`, `nstout`, `restraint_k`, `buffer`, `minimize`,
`rigid_ribosome`, `nascent_only_output`, `device`, `ppn`, `outdir`) plus:

- `rna_model = cosmo | topo` — ribosome RNA representation (§5).
- `model = hps_urry|hps_kr|hps_ss|mpipi` — nascent-chain force field (default
  `hps_urry`).

Drop topo-only keys (`domain_def`, `stride_output_file`, Gō-specific). Defer the
v2-extra keys (`trna_tether`, `tunnel_wall*`, `ptc_offset`, `post_elongation*`)
until those features land.


## 8. Phasing

1. Scaffold `cosmo/translation/` (package, CLI, `ElongateConfig`); fold
   `crop_ribosome.py` → `truncate_ribosome.py`.
2. **v1** elongation loop on `asyn.pdb` (`L0 = 5`), per-length outputs.
3. **`cg_ribosome.py`** for both `rna_model` reps; truncate; structure inventory.
4. **v2** rigid ribosome + cross-forces, validated under both `rna_model`.
5. Later: v2 extras (tether, wall, post-elongation), `make_movie.py`, tutorial +
   docs.


## 9. References

Shared physics from topo: O'Brien et al. *JACS* 2012 (RNC force field), *JACS*
2011 (P/R/BR rRNA + truncation). cosmo nascent-chain physics: Dignon/Mittal HPS
(`hps_urry`/`hps_kr`), Joseph et al. *Nat. Comput. Sci.* 2021 (`mpipi`). Procedure
adopted; all code implemented in cosmo.
