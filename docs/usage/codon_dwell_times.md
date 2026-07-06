# Codon dwell-time tables (per-codon timing)

Per-codon timing ({doc}`continuous_synthesis`) needs a **codon dwell-time table** — the
`codon_times` key of a `csp.ini` / `cylinder.ini`. It maps each codon to its **mean
in-vivo per-codon translation (dwell) time in seconds** (the codon's mean
first-passage time), which is what sets each residue's MD-segment length. There is **no
bundled default**: you pick a table for your organism (and, where relevant, elongation
rate). cosmo ships a small library of them.

## Where they live & the format

The tables are under **`assets/csp/codon_dwell_times/<organism>/`**. Each is a plain
text file, one row per codon:

```
# codon   dwell_time_seconds   amino_acid
UUU       0.068164             PHE
...
UAA       0.005513             STOP
```

- **RNA alphabet** (`U`), tab-separated, `#` comments ignored; **64 rows** (61 sense +
  3 stop). Larger time = slower codon; stop codons carry the amino acid `STOP`.
- The **third column (amino acid)** is what `mrna = fastest`/`slowest` uses to pick each
  residue's extreme synonymous codon — see [Fastest / slowest mRNA](#fastest-slowest-mrna).
- Every organism directory also has a **`*_dwell_time_methods.md`** with the full
  provenance (source dataset, references, how relative data were scaled to seconds).

## Shipped tables

| Organism | File(s) | Basis (see the methods note for full provenance) |
|----------|---------|--------------------------------------------------|
| *E. coli* | `ecoli_codon_dwell_times_310K.txt` | Fluitt, Pienaar & Viljoen (2007) mechanistic decoding times, O'Brien-rescaled so the codon-usage-weighted mean = 0.061 s/codon (16.5 aa/s) at **310 K**. Absolute; the most-validated table. |
| *S. cerevisiae* | `yeast_codon_dwell_times.txt` | Gardin et al. (2014) relative decoding rates (RRT), frequency-weighted → seconds (~0.108 s/codon; ~9.3 codons/s). |
| | `yeast_codon_dwell_times_OBrien.txt` | O'Brien-lab chemical-kinetic model, **absolute** (~0.233 s/codon; ~4.3 codons/s) — the most pipeline-native yeast table. |
| | `yeast_codon_dwell_times_OBrien_rescaled_9.3.txt` | The O'Brien table rescaled to 9.3 aa/s (matches the Gardin rate). |
| *H. sapiens* | `human_codon_dwell_times.txt` | Gobet et al. (2020) liver ribosome-profiling dwell times (Lintner dataset) → seconds. |
| | `human_codon_dwell_times_v2_freqweighted.txt` | Same, **frequency-weighted** to an elongation-rate anchor R = 5.6 codons/s. |
| *N. crassa* | `ncrassa_codon_dwell_times_R{5.35_low,6.7_mid,8.02_high}.txt` | Yang et al. (2019) relative codon decoding times (RCDT) → seconds at three **measured elongation-rate anchors** (5.35 / 6.7 / 8.02 codons/s). Least-certain table (relative values digitized from a figure). |
| | `..._v2_freqweighted.txt` (× 3) | Frequency-weighted variants of the three anchors. |

**Frequency-weighted (`v2`) variants.** These rescale the per-codon times so the
*codon-usage-weighted* mean over the organism's transcriptome hits the target
elongation rate, instead of an unweighted mean. Use them when you want the genome-wide
average speed to match the anchor; use the plain (`v1`) tables for the raw per-codon
values.

```{note}
References — *E. coli*: Fluitt, Pienaar & Viljoen, *Comput. Biol. Chem.* **31**:335–346
(2007). Yeast: Gardin et al., *eLife* **3**:e03735 (2014). Human: Gobet et al., *PNAS*
**117**(17):9630–9641 (2020). *N. crassa*: Yang et al., *Nucleic Acids Res.*
**47**(17):9243–9258 (2019).
```

## Use one in a run

Point the `codon_times` key at a table (a **path** selects per-codon timing; a plain
**number** would instead mean a uniform time for every codon):

```ini
pdb_file    = my_protein.pdb
mrna        = my_protein_mrna.txt
codon_times = assets/csp/codon_dwell_times/ecoli/ecoli_codon_dwell_times_310K.txt
```

Keep the thermostat consistent with the table's temperature (the *E. coli* table is
310 K); set `ref_t` to your table's temperature. See {doc}`synthesis_control` for the
full key reference.

(fastest-slowest-mrna)=
## Fastest / slowest mRNA

Because each table carries the codon→amino-acid mapping (column 3), you can set
`mrna = fastest` or `slowest` (instead of a sequence file) to auto-build a
synonymous-codon mRNA that encodes your protein with every residue's fastest/slowest
codon **for the chosen table** — the extremes of the codon-optimization axis. The runner
reads the protein sequence from `pdb_file`, picks each residue's extreme synonymous codon
(plus a terminating stop), and writes `mrna_fastest.txt` / `mrna_slowest.txt` next to the
PDB. This needs a `codon_times` **table** (it defines fast/slow); a uniform numeric
`codon_times` is rejected. To pre-generate one standalone:

```bash
cosmo-make-mrna --pdb my_protein.pdb \
    --codon-times assets/csp/codon_dwell_times/yeast/yeast_codon_dwell_times_OBrien.txt \
    --mode slowest
```
