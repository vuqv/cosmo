# The ribosome structure (get one, or build your own)

The explicit-ribosome runner ({doc}`continuous_synthesis`) needs a **rigid, oriented,
truncated, coarse-grained large subunit** as its `ribosome` input — the
`4v9d_50S_PtR_5jte_AtR_model_cg_trunc.pdb` you see in `tutorials/08_csp_cg_ribosome/`.
(The cylinder runner needs none of this.)

Every such structure is oriented in the **tunnel frame**: PTC at the origin, exit tunnel
on **+x**, tRNA tails on **+y** — so the tunnel axis *is* the X-axis and the radial
distance to it is `sqrt(y²+z²)`.

## The canonical ribosome

The reference structure is the *E. coli* **4V9D** 50S subunit (with the 5JTE A-site tRNA
grafted in), coarse-grained to the O'Brien **P/R/BR** representation (4 576 beads) and
truncated around the exit tunnel: `4v9d_50S_PtR_5jte_AtR_model_cg_trunc.pdb`. It is
bundled in `tutorials/08_csp_cg_ribosome/`.

Point the `ribosome` key of a `csp.ini` at any such structure — it is a plain path:

```ini
pdb_file  = my_protein.pdb
ribosome  = /path/to/4v9d_50S_PtR_5jte_AtR_model_cg_trunc.pdb
```

## Other organisms / building your own

cosmo does **not** ship a structure-preparation pipeline (or bundle ribosome structures
in the package itself). The synthesis subsystem is a port of the sibling **`topo`**
package, and the two share the same tunnel frame, truncation rule, and CG mapping — so to
prepare a ribosome from any PDB deposition (or to use a ready-made yeast / *N. crassa* /
human structure), use topo's scripted pipeline and its guide:

- **topo's `assets/csp/prepare_ribosome/`** — the fetch → orient → coarse-grain →
  truncate pipeline plus four ready-made organisms (*E. coli*, *S. cerevisiae*,
  *N. crassa*, *H. sapiens*).
- **topo's `usage/ribosome_preparation`** doc page — the ready-made table, per-organism
  provenance/caveats, and how to re-truncate with your own criteria.

A structure prepared by that pipeline (P/R/BR representation) drops straight into
`cosmo-csp` — the same 4 576-bead file bundled in the tutorial.

```{note}
Only the **P/R/BR** rRNA representation is supported. An earlier single-bead-per-nucleotide
(`cosmo` 1-bead) rep existed for the retired `cosmo.translation` runner; it is not usable
by `cosmo-csp` (the PTC-geometry optimizer needs the P/R/BR tRNA beads) and has been removed.
```
