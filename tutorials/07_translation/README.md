# Co-translational synthesis example (`cosmo.translation`, build step v1)

Grows a nascent protein chain **N→C, one residue at a time**, on the ribosome,
instead of folding the full-length chain in bulk. This is **build step v1**: the
simulated system is the **nascent chain only**; the truncated ribosome supplies
just two fixed points — the **P-anchor** (P-site tRNA residue-76 `R` bead) and the
**A-anchor** (A-site tRNA residue-76 `R` bead) — used to place each new residue and
as the C-terminus restraint target. No ribosome forces yet (that is build step v2;
see [`cosmo/translation/PLAN.md`](../../cosmo/translation/PLAN.md)).

## Run

```bash
cd tutorials/07_translation
cosmo-elongate -f elongate.ini
# equivalently:
python -m cosmo.translation -f elongate.ini
```

Handy overrides: `-o my_run` (output dir), `--device GPU`.

## How it works (per length `L = L0 .. L_max`)

1. Build the length-`L` cosmo model on the **first `L` residues** of `pdb_file`
   (cosmo is sequence-based — no STRIDE / native contacts; PLAN.md §2).
2. Seed coordinates: `L = L0` → an extended chain laid along the tunnel axis (+x)
   from the P-anchor; `L > L0` → residues `1..L-1` from the previous length's final
   structure, plus the new residue `L` placed at the A-anchor + `buffer`.
3. Restrain only residue `L` (the current C-terminus) to the P-anchor, so the new
   residue migrates **A→P**.
4. Minimize → run `n_steps` → save. Each length writes its own folder and its final
   structure seeds the next length.

## Outputs

One self-contained folder per length, `synth_out/L_<L>/`:

| file | what it is |
|------|------------|
| `traj.dcd` | trajectory for this length |
| `traj.log` | energy / temperature log |
| `traj.psf` | topology |
| `traj.chk` | OpenMM checkpoint |
| `traj_final.pdb` | final structure (seeds the next length) |
| `traj_runinfo.log` | run provenance |
| `native_1_<L>.pdb` | length-`L` native CA structure used for the build |
| `seed.pdb` | seeded starting coordinates for this length |

## Build step v2 — the rigid ribosome

Set `rigid_ribosome = yes` in `elongate.ini` to append the truncated ribosome as
**fixed (mass-0) scenery** and turn on the two ribosome↔nascent cross-interactions:

- **Excluded volume** — a dedicated `CustomNonbondedForce`, pure `ε·(σ_ij/r)¹²`,
  restricted to `{nascent}×{ribosome}`.
- **Electrostatics** — the existing Yukawa force extended over the ribosome charges
  (rRNA phosphate −1e, charged residues), `{nascent}×{nascent}` + `{nascent}×{ribosome}`.

The ribosome is held rigid (mass 0 → ~0 drift; no intra-ribosome forces). By
default v2 writes **nascent-only** trajectories/PSF/final (`nascent_only_output =
yes`); the checkpoint still holds the full system.

**C-terminus hold (`trna_tether`, default on)** — instead of a generic position
restraint, the C-terminus is tethered to the P-site tRNA the O'Brien way: a harmonic
**bond** `CA(L)–tRNA` plus an **orienting angle** `CA(L-1)–CA(L)–tRNA` (double-Gaussian
basins 91.7°/130°) that aims the chain down the tunnel. `trna_tether = no` falls back
to the position restraint.

**Tunnel wall (`tunnel_wall`, default on)** — a one-sided planar restraint
`U = k·min(x−x0,0)²` on every nascent bead keeps the chain at `x ≥ x0` so it can only
extrude **forward** (+x). `x0` defaults to the C-terminal-AA addition plane (PTC) =
P-anchor x + tether bond length (set `tunnel_wall_x0` to override).

**Excluded-volume strength & rRNA bead size (`ribo_eps`, `ribo_rna_radii`).**
O'Brien's default (rRNA radii **0.71 nm**, `ribo_eps` **0.00055 kJ/mol**) is a very
*soft* wall: the topo rRNA bead radius is large vs cosmo's CG tunnel bore, so a hard
wall would jam the chain — hence the soft wall, which lets the chain **thread** the
bore but also lets it visibly **penetrate** the ribosome. Shrinking the rRNA beads
opens the bore so a *real* wall fits: e.g. `ribo_rna_radii = 0.40`, `ribo_eps = 0.2`
cuts bead-into-ribosome overlap ~6× (0.22 → 0.04 nm) **without** jamming. The tutorial
uses these confining values; blank either key to fall back to O'Brien's soft defaults.

> **New residues are seeded at their tether rest position** (not ~1 nm away at the
> A-site), so the tether bond starts at equilibrium — this avoids a large placement
> strain and keeps the per-step temperature controlled.

`rna_model` chooses the ribosome RNA representation: **`topo`** (3/4-bead P/R/BR,
default) or **`cosmo`** (1 bead/nucleotide). Each has a ready truncated CG ribosome
in `cosmo/translation/structures/` (`..._cg_trunc.pdb` and `..._cg_cosmo_trunc.pdb`).
Regenerate them from the all-atom source yourself with:

```bash
python -m cosmo.translation.cg_ribosome -i <all_atom.pdb> -o cg.pdb --rna-model topo
python -m cosmo.translation.truncate_ribosome -i cg.pdb -o cg_trunc.pdb
```

## Visualize the growth in VMD

Stitch the per-length trajectories into one fixed-width "growing chain" movie:

```bash
cosmo-elongate-movie -o synth_out
vmd -e synth_out/movie.tcl
```

For a v2 run, pass `--ribosome <trunc.pdb>` to also load the static ribosome as
scenery the chain grows inside.

## Notes

- **No ribosome confinement in v1**, so the short nascent chain folds/collapses
  freely rather than extruding straight down the tunnel — turn on v2 for the tunnel.
- The `model` key selects the nascent-chain force field (`hps_urry` default).
  Choosing `mpipi` triggers the tabulated-potential dummy-id handling for the
  ribosome beads (PLAN.md §2.1) — handled automatically.
- Programmatic use: `from cosmo.translation import run_elongation, ElongationParams`
  (or `read_elongate_config`).
