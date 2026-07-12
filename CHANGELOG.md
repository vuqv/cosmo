# Changelog

All notable changes to COSMO are documented here. The format is loosely based on
[Keep a Changelog](https://keepachangelog.com/); releases correspond to git tags.
Releases use CalVer `YEAR.N` (the Nth release of the year). For releases before
`2026.1` (the old `vX.Y` tags), see the [git tags](https://github.com/vuqv/cosmo/tags).

## [2026.1] — 2026-07-12

Large release: the project was renamed to **COSMO**, RNA gained `mpipi` support, and a
new co-translational protein-synthesis subsystem (`cosmo.csp`) was added.
**87 commits, 324 files, +25,107 / −71,638** since v1.3.1.

### ⚠️ Breaking
- **Renamed the project and package `hpsOpenMM` / `hps` → `COSMO` / `cosmo`.** All
  import paths, class names, console commands, and docs moved from `hps*` to `cosmo*`.
  This release is **not import-compatible** with the previous release (`v1.3.1`).

### Added
- **`cosmo.csp` — co-translational protein synthesis subsystem** (mirrors the sibling
  `topo.csp`, on cosmo's sequence-based HPS/mpipi force field — no STRIDE, native
  contacts, or `domain.yaml`). Two confinement models:
  - **Cylinder model** (`cosmo-cylinder`) — an analytic cylindrical exit tunnel (no
    ribosome beads; fast, never jams).
  - **CG-ribosome model** (`cosmo-csp`) — an explicit coarse-grained ribosome as rigid
    scenery, with the 3-stage elongation cycle (peptidyl-transfer → translocation →
    tRNA-binding) and the A→P C-terminus switch. Ribosome structures (borrowed from
    topo's ribosome-preparation pipeline) are provided for **four organisms — *E. coli*,
    yeast, *N. crassa*, and human** — swapped by pointing the `ribosome` key at the
    organism's PDB.
- **Per-codon kinetics** (`cosmo.csp.kinetics`) with shipped **codon dwell-time
  tables** for *E. coli*, human, yeast, and *N. crassa*.
- **mRNA generator** (`cosmo-make-mrna`) — synthetic transcripts in **fastest**,
  **slowest**, and **median** synonymous-codon modes.
- **tRNA tether** (`trna_tether`) — per-stage A/P covalent attachment (bond + two
  orienting angles + improper); a `prev-N-R` backbone tether angle for all IDP models.
- **Resume** for synthesis runs — continue from the last completed residue via
  `progress.log` and a consolidated on-disk layout (HPC-requeue safe).
- **Rigid-bond (`AllBonds`) support** and model-independent CSP steric radii.
- **Post-synthesis phases** — optional `ejection` then `dissociation` free runs.
- **Movie tool** (`cosmo-csp-movie`) — stitch per-residue/-stage trajectories into a
  VMD movie, with `--tunnel` (analytic tunnel) or `--ribosome` (CG ribosome) scenery.
- **RNA support**: CA/P coarse-grained representation (retain the phosphate **P** bead
  alongside each Cα); **`mpipi` RNA support** (validated).
- New **console commands**: `cosmo-csp`, `cosmo-cylinder`, `cosmo-csp-movie`,
  `cosmo-make-mrna`.
- New **`utils`** module.
- New tutorials **7 (cylinder)** and **8 (CG ribosome)**; new usage pages for
  continuous synthesis, cylinder synthesis, ribosome preparation, codon dwell-times,
  resuming long runs, visualizing the synthesis process, and a consolidated synthesis
  control-options reference; new bead-chain logo + favicon.

### Changed
- `mpipi` uses a **3σ cutoff** for every interaction pair.
- Upgraded to **OpenMM 8.0 / ParmEd 4.0**; print the OpenMM version and the nonbonded
  cutoff distance at runtime; `openmm` imported as `mm`.
- Intermediate-scattering-function (ISF) analysis reimplemented in **Julia** (replacing
  the Python version), with curve-fit bounds.
- More intuitive `runinfo.log` format; simulation-control options grouped into ordered
  categories; API toctree wires in `cosmo.mdrun` and `cosmo.utils`.
- Tutorials list both the console command and the `python -m …` module form; CSP
  runners shown console-first; **cylinder ordered before the CG-ribosome docs** in the
  Protein-synthesis sidebar.

### Fixed
- Convert `eps_di` parameters from **kcal/mol → kJ/mol** and document units (Issue #11).
- Fix phosphorylation-residue handling.
- Check residues are consecutive before adding harmonic bonds.
- CSP: tether ½-factor, always-on PTC optimization, retired-key guard, stale
  codon-table docs.

### Removed
- Retired `cosmo.translation`; dropped the single-bead ribosome + prep tooling
  (`crop_ribosome`, `ctf_utils`, `filter_pdb`); stopped tracking `CLAUDE.md`.
