# Plan — tutorial 10: co-translational synthesis with **codon translation kinetics**

> Status: **proposed, awaiting approval.** Port of the validated legacy
> `tutorials/growing_reference/` (hps_kr, works on legacy code) to the current
> cosmo codebase. No code written yet — this document is for review.

## 1. Motivation / what makes tutorial 10 different

Tutorials 07–09 grow the nascent chain on a **constant** schedule (the same
`n_steps` per residue). The `growing_reference` adds the scientifically
distinctive piece: **realistic per-codon translation kinetics**. The number of MD
steps spent on each residue is *variable*, set by that codon's measured median
translation time:

```
nsteps(residue k) = codon_time_ms(codon_k) * 1e9 / (dt_ps * scale_factor)
```

where `codon_time_ms` comes from `cosmo.utils.ctf_utils.codon_translation_time[organism][codon]['Median_time']`,
and `scale_factor` compresses real translation time into a tractable simulation
time. Fast-translated codons get few steps, slow ("rare") codons get many — so the
simulation reproduces **translation-rate-dependent** co-translational behaviour
(e.g. pausing at slow codons). This is the feature being ported.

## 2. The validated legacy method (to preserve)

From `tutorials/growing_reference/grow_nascent.py` (+ the `_82` / `_nascent_dcd`
variants, which only add nascent-only PSF/DCD writers):

1. **Full ribosome frozen as explicit beads.** Build the whole complex with
   `buildHPSModel(complex, frozen_indices=<all ribosome atoms>,
   except_chains=<all ribosome chains>, nb_exclusions=<all intra-ribosome pairs>)`.
   Frozen (mass-0) ribosome; ribosome–ribosome nonbonded excluded; ribosome chains
   excluded from bonding. The nascent chain interacts normally with everything.
   **Port change (§4/§6):** the reference's `nb_exclusions` was the **full O(n²)
   list of intra-ribosome pairs**; the port replaces it with **interaction groups**
   (the efficient, established mechanism — see §6).
2. **Grow by appending one CA at a time** at the origin (chain `8`) to the ribosome
   PDB → `complex_l{N}.pdb` (`append_amino_acid_to_complex`).
3. **PTC mechanics.** The newest residue is harmonically restrained at
   `(0.38, 0, 0)` nm (`k = 1e4 kJ/mol/nm²`); every earlier nascent bead is held at
   `x > 0.38 nm` by a one-sided planar wall (same `k`). New residues enter at the
   origin and are pushed out (+x).
4. **Codon-kinetics schedule** (§1): variable `nsteps` per residue.
5. Defaults: `model = hps_kr`, `organism = ecoli`, `translation_type = fast`,
   `nascent_chain_id = 8`, `scale_factor = 433129300`, `ref_t = 310 K`,
   `dt = 0.01 ps`, `pbc = no`.

Inputs (copied from the reference, validated): `nascent.pdb` (140 lines, chain 8),
`rib.pdb` (1972 atoms), `md.ini`.

## 3. What actually breaks under the current codebase (verified)

Only **two** calls are broken; everything else still exists and works:

| Legacy call | Status now | Fix |
|---|---|---|
| `cosmo.models.buildHPSModel(...)` | alias **removed** | call `buildCoarseGrainModel(...)` — `frozen_indices` / `except_chains` args **still exist and work** (verified in `cosmo/core/models.py`); see §6 re: `nb_exclusions` |
| `from distutils.util import strtobool` | removed (py3.12+) | `from cosmo.utils.config import strtobool` |

Still present and unchanged (verified): `parse_nascent_sequence`,
`generate_mRNA_sequence`, `write_pdb_with_chain_ids`,
`ctf_utils.codon_translation_time` (keys: `ecoli`, `yeast`), `crop_ribosome`,
`cosmo.system().coarseGrainingStructure()/.getAtoms()`.

→ The port is therefore **low-risk**: the scientific method needs no change, only
the two API fixes plus modernization/cleanup.

## 4. Proposed deliverables — `tutorials/10_translation_kinetics/`

1. **`grow_nascent.py`** — one consolidated module (replaces the 3 legacy variants):
   - `append_amino_acid_to_complex(...)` — kept as-is.
   - Build via `cosmo.models.buildCoarseGrainModel(..., frozen_indices,
     except_chains)` (the build-call rename) — **no `nb_exclusions`**; exclude
     ribosome–ribosome nonbonded with **interaction groups** instead (§6).
   - Restraint `(0.38,0,0)` + planar wall `x>0.38` — kept (the reference geometry).
   - **Codon-kinetics scheduler — kept verbatim** (the scientific core; Phase 5).
   - Modernizations (see §5 decisions): `strtobool` import fix; **reuse
     `cosmo.engine`** for per-iteration setup/reporters/finalize (standard
     `cosmoReporter` log + `_runinfo.log` provenance); **trajectory writes the
     nascent chain only** — nascent-only PSF + DCD (the `_82` variant /
     `NascentDCDReporter` pattern); the simulation and per-iteration final PDB still
     cover the whole complex, but the frozen ribosome is **never** written to the
     DCD/PSF; **output layout `L_<L>/traj.*`** so the shipped
     `cosmo-elongate-movie` works.
   - A small config dataclass + `read_config()` for the translation keys
     (`organism`, `scale_factor`, `translation_type`, `nascent_chain_id`,
     `rib_structure`, `nascent_structure`) layered on the standard MD keys —
     mirrors the tut-09 `cylinder.py` pattern.
2. **`md.ini`** — ported, with comments (model `hps_kr`, ecoli, fast,
   `n_steps = 1000`; **no `scale_factor`** — that arrives in Phase 5).
3. **`nascent.pdb`, `rib.pdb`** — copied from the reference (validated inputs).
   `nascent.pdb` is **CG, CA-only**; `rib.pdb` is **already cropped** — used as-is.
4. **`README.md`** — codon-kinetics vs constant schedule; the `crop_ribosome` prep
   step; run command; visualization (`cosmo-elongate-movie`).

## 5. Decisions (resolved)

1. **Per-iteration setup:** **reuse `cosmo.engine`** (standard log + provenance,
   aligns with 07–09). ✔ confirmed.
2. **Output layout:** **`L_<L>/traj.*`** so `cosmo-elongate-movie` works directly.
   ✔ confirmed.
3. **Scope:** **cosmo only** (hps_kr is a cosmo IDP force field); not topo.
4. **Method:** **keep the legacy frozen-full-ribosome + origin-PTC method**; do
   *not* re-architect onto `elongate.py`'s anchor/tunnel approach (genuinely
   different ribosome treatment; the legacy one is what was validated).
5. **Inputs are already simulation-ready:** `nascent.pdb` is the **CG (CA-only)**
   nascent chain and `rib.pdb` is **already cropped** — so the tutorial does **no**
   nascent coarse-graining and **no** `crop_ribosome` step (the latter is mentioned
   in the README only as the optional upstream prep that produced `rib.pdb`).
6. **`scale_factor` is not used in this test.** With the constant 1000-steps/residue
   schedule there is no real-time→sim-time mapping, so `scale_factor` is **omitted
   from `md.ini`** and only introduced in Phase 5 alongside the codon scheduler.

## 6. Ribosome–ribosome exclusion: interaction groups (not O(n²) `nb_exclusions`)

The reference excluded intra-ribosome nonbonded by passing **all ~1.9M pairs**
(O(n²) in the ~1972 ribosome atoms) as `nb_exclusions`. It works but is heavy to
build/store. The port replaces it with **interaction groups** — the established,
efficient mechanism already used by `cosmo/translation/ribosome.py:append_ribosome`
(lines 254–262).

Mechanism (verified): `buildCoarseGrainModel` takes no interaction-group argument,
so add them **post-build** on the model's `CustomNonbondedForce` attributes
(`cgModel.ashbaugh_HatchForce` for the HPS-scale `hps_kr`, and `cgModel.yukawaForce`):

```
nascent_idx, ribo_idx = <nascent atoms>, <frozen ribosome atoms>
for f in (cgModel.ashbaugh_HatchForce, cgModel.yukawaForce):
    f.addInteractionGroup(nascent_idx, nascent_idx)   # nascent–nascent
    f.addInteractionGroup(nascent_idx, ribo_idx)       # nascent–ribosome
    # {ribosome}×{ribosome} omitted -> never computed (ribosome is frozen anyway)
```

Same physics as the reference (no intra-ribosome nonbonded; nascent interacts with
everything), far cheaper. Build is called **without** `nb_exclusions`; the bonded
1-2/1-3 exclusions the builder adds stay uniform across the nonbonded forces.

## 7. Validation before hand-off

- Short smoke run (first ~5 residues, **constant 1000 steps/residue**) on GPU:
  confirm it builds, the frozen ribosome stays put, the chain grows at the PTC
  (newest at ~0.38 nm, earlier beads walled to x>0.38), and the interaction groups
  give sane energies (ribosome–ribosome not computed).
- Confirm the **trajectory contains only the nascent chain** (nascent-only PSF/DCD;
  no ribosome beads in the frames).
- Confirm `cosmo-elongate-movie -o <outdir>` stitches the per-length outputs.

## 8. Open questions — all resolved (decision log)

| # | Question | Resolution | Where |
|---|---|---|---|
| 1 | Per-iteration setup: reuse engine vs hand-rolled? | **Reuse `cosmo.engine`** | §5.1 |
| 2 | Output layout: `L_<L>/traj.*` vs `complex_l{N}.*`? | **`L_<L>/traj.*`** (movie tool works) | §5.2 |
| 3 | Scope: cosmo only, or topo too? | **cosmo only** (hps_kr) | §5.3 |
| 4 | Method: keep legacy, or re-architect onto `elongate.py`? | **Keep legacy** frozen-ribosome + origin-PTC | §5.4 |
| 5 | Inputs: coarse-grain / crop them? | **No** — `nascent.pdb` is CG (CA-only), `rib.pdb` already cropped; use as-is | §5.5 |
| 6 | `scale_factor` default for the test? | **Not used** — constant 1000 steps/residue; omitted from `md.ini`, returns in Phase 5 | §5.6 |
| 7 | Schedule for the first port? | **Constant 1000 steps/residue**; realistic codon kinetics deferred to Phase 5 | §1, §5.6 |
| 8 | Intra-ribosome exclusion: O(n²) `nb_exclusions` vs interaction groups? | **Interaction groups** (post-build, as `append_ribosome`) | §6 |
| 9 | What goes into the trajectory? | **Nascent chain only** (nascent-only PSF+DCD); sim + final PDB stay full-complex | §4.1 |

No open questions remain — the plan is ready to implement (Phases 0–4 of `TASKS.md`;
codon kinetics is Phase 5).
