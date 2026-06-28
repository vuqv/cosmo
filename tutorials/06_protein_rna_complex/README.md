# Tutorial 6 — Protein–RNA complex

**Goal:** build and simulate a **multi-component** system — a protein chain
**plus** an RNA chain — the minimal ingredients of a ribonucleoprotein (RNP)
condensate. You'll learn how COSMO handles mixed biomolecule types and why the
model choice is constrained when nucleic acids are involved.

**Time:** a few seconds on a CPU.

---

## Files in this folder

| File | Role |
|------|------|
| `protein.pdb` | The protein component (all-atom). |
| `rna.pdb` | The RNA component (all-atom). |
| `combine.py` | Merges the two into one PDB (RNA → chain A, protein → chain B). |
| `combined.pdb` | The pre-merged input used by `md.ini` (ready to run). |
| `md.ini` | Configuration — note `model = hps_kr`. |
| `run_simulation.py` | Thin runner wrapper. |

## Background concepts

- **One bead per residue, two chemistries.** COSMO coarse-grains a protein at the
  Cα and a nucleotide at the phosphate (P). A mixed system is just a topology
  that contains both kinds of bead, with the appropriate per-type parameters.
- **Model choice is constrained.** Only models that carry **nucleic-acid
  parameters** can score protein–RNA and RNA–RNA interactions. In COSMO that is
  **`hps_kr`** (Kapcha–Rossy scale, with RNA and phospho-protein parameters), so
  this tutorial uses `model = hps_kr`. Using a protein-only model (`hps_urry`,
  `hps_ss`) here would lack RNA parameters. (`mpipi` also has RNA parameters; see
  the model table in Tutorial 2.)
- **Everything is in one PDB.** COSMO reads a single `pdb_file`. Multi-component
  and multi-chain systems are assembled **upstream** into one PDB whose chains
  are distinct — exactly what `combine.py` does here.

## Step-by-step

### 1. (Optional) Rebuild the combined input
`combined.pdb` is already provided, but you can regenerate it to see how a
multi-component input is assembled:
```bash
python combine.py        # reads rna.pdb + protein.pdb -> writes combined.pdb
```
It puts the RNA on chain A and the protein on chain B and renumbers residues so
the merged file is clean.

### 2. Run the complex
```bash
python run_simulation.py -f md.ini
```
Outputs go to `traj/complex.*` (see Tutorial 1/4 for the full file list). The
build log reports **two chains** and uses the `hps_kr` force field, which scores
the protein–protein, protein–RNA, and RNA–RNA contacts together.

### 3. Watch them interact
```bash
vmd traj/complex.psf traj/complex.dcd
```
With attractive protein–RNA interactions the two chains associate rather than
drifting apart — the elementary event behind RNP condensate formation.

## Scaling up to a condensate

This two-chain system is the building block. To study an actual RNP condensate
you would:

1. replicate many protein and RNA chains into one PDB (the upstream assembly
   step, generalizing `combine.py`),
2. put them in a periodic box (Tutorial 3) and use the **slab protocol**
   (Tutorial 5) to measure coexisting dense/dilute phases,
3. scan the **protein:RNA stoichiometry** and temperature — RNA can both promote
   and (in excess) dissolve condensates, the well-known reentrant behavior.

## Try next

- Swap `model = mpipi` (also nucleic-acid capable) and compare the strength of
  association.
- Put the complex in a box (`pbc = yes`, `box_dimension = …`) and run NVT, as a
  bridge toward the slab setup in Tutorial 5.
