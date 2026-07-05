# The ribosome structure (get one, or build your own)

The explicit-ribosome runner ({doc}`continuous_synthesis`) needs a **rigid, oriented,
truncated, coarse-grained large subunit** as its `ribosome` input — the
`ribosome_trunc.pdb` you see in the tutorials. (The cylinder runner needs none of this.)

Every such structure is oriented in the **tunnel frame**: PTC at the origin, exit tunnel
on **+x**, tRNA tails on **+y** — so the tunnel axis *is* the X-axis and the radial
distance to it is `sqrt(y²+z²)`.

## What cosmo ships

cosmo bundles the *E. coli* **4V9D** 50S subunit (with the 5JTE A-site tRNA grafted in),
in `cosmo/csp/structures/`, in two coarse-grain representations:

| File | Beads | rRNA rep | Used by |
|------|------:|----------|---------|
| `4v9d_50S_PtR_5jte_AtR_model_cg_trunc.pdb` | 4 576 | topo P/R/BR | `cosmo-csp` (this is the one copied into the tutorials/sandbox as `ribosome_trunc.pdb`) |
| `4v9d_50S_PtR_5jte_AtR_model_cg_cosmo_trunc.pdb` | 1 706 | cosmo rep | the single-stage `cosmo-elongate` runner |

Point the `ribosome` key of a `csp.ini` at a truncated structure — it is a plain path:

```ini
pdb_file  = my_protein.pdb
ribosome  = /path/to/4v9d_50S_PtR_5jte_AtR_model_cg_trunc.pdb
```

## Other organisms / building your own

cosmo does **not** ship the structure-preparation pipeline itself. The synthesis
subsystem is a port of the sibling **`topo`** package, and the two share the same tunnel
frame, truncation rule, and CG mapping — so to prepare a ribosome from any PDB deposition
(or to use a ready-made yeast / *N. crassa* / human structure), use topo's scripted
pipeline and its guide:

- **topo's `assets/csp/prepare_ribosome/`** — the fetch → orient → coarse-grain →
  truncate pipeline plus four ready-made organisms (*E. coli*, *S. cerevisiae*,
  *N. crassa*, *H. sapiens*).
- **topo's `usage/ribosome_preparation`** doc page — the ready-made table, per-organism
  provenance/caveats, and how to re-truncate with your own criteria.

A structure prepared by that pipeline (topo P/R/BR representation) drops straight into
`cosmo-csp` — it is exactly the 4 576-bead file bundled above.
