# S. cerevisiae (yeast) — per-codon translation times (methods & provenance)

There are **two independent yeast tables**, from two different methods:

| file | method | data type | mean s/codon |
|---|---|---|---|
| `yeast_codon_dwell_times.txt` | Gardin RRT × R (frequency-weighted) | relative → scaled | ~0.108 |
| `yeast_codon_dwell_times_OBrien.txt` | Sharma/O'Brien chemical-kinetic | absolute (model) | ~0.233 |

They differ ~2× in absolute scale because they encode different rate
assumptions (9.3 vs ~4.3 codons/s). Both are internally consistent; the O'Brien
table is the more pipeline-native one (same chemical-kinetic framework).

---

## TABLE A — Gardin RRT → seconds  (`yeast_codon_dwell_times.txt`)

### A.1 Source
**Gardin et al. (2014)** "Measurement of average decoding rates of the 61 sense
codons in vivo." *eLife* 3:e03735. doi:10.7554/eLife.03735

- **Table 2A** gives, for all 61 sense codons: **RRT at footprint position 6**
  (the ribosomal A-site) *and* the genomic codon **usage per 1000**. RRT is
  normalized so the genome-average codon = 1.0; larger RRT = slower.
- This table is embedded in the eLife article HTML (retrieved directly), so no
  supplement download was needed.
- **Data fix:** as published, Table 2A lists codon `TCT` twice as "TCC"; the
  second entry (usage 23.5, RRT 0.98) is `TCT` and is corrected here.

### A.2 Elongation-rate anchor
**R = 9.3 codons/s** (mean; Chu et al. 2014, *Cell* 153:1589, validated vs
Arava et al. 2003 polysome data). `τ = 1/9.3 = 0.107527 s`.

### A.3 Conversion (frequency-weighted, Option B — this is what the file uses)
Usage frequencies `f_c` come from Table 2A's own usage column (Σ over 61 = 1).

```
r̄_w = Σ_c f_c · RRT_c            = 1.03858
t_c = RRT_c / ( R · r̄_w )
```

Verified `Σ_c f_c · t_c = τ = 0.107527 s`. T→U, 6 dp. Stops (UAA/UAG/UGA) =
mean = τ. Range: fastest **ACC 0.072473 s**, slowest **CUC 0.195676 s** (brackets
the E. coli mean 0.068 s).

### A.4 Unweighted alternative (Option A)
If you instead want the unweighted convention (as used for the human/N. crassa
v1 files): `t_c = RRT_c · τ / mean_over_61(RRT)`. Because Gardin RRT is already
≈mean-normalized, A and B differ by only ~1 % here.

### A.5 Reproduce from scratch
1. Read Table 2A (RRT@pos6 + usage/1000) for 61 codons; correct the TCT/TCC typo.
2. Choose R (9.3). τ = 1/R.
3. Normalize usage to f_c (Σ = 1); compute r̄_w = Σ f_c·RRT_c.
4. t_c = RRT_c/(R·r̄_w). Convert T→U. Stops = τ. Write 6 dp.

---

## TABLE B — O'Brien chemical-kinetic absolute times  (`yeast_codon_dwell_times_OBrien.txt`)

### B.1 Source
**Sharma, Sormanni, Ahmed, Ciryam, Friedrich, Kramer, O'Brien (2019)**
"A chemical kinetic basis for measuring translation initiation and elongation
rates from ribosome profiling data." *PLoS Comput Biol* 15(5):e1007070.
doi:10.1371/journal.pcbi.1007070

- **S1 File** (`pcbi.1007070.s016.xlsx`) — sheets **"Nissley data translation
  rates"** and **"Weissman data translation rates"**. Each has, per codon:
  count, **mean translation time (ms)**, median, SD, etc., for all 64 codons
  (stop codons included, measured).
- Absolute per-codon times from applying the authors' chemical-kinetic model
  (their Eq. 12) to yeast ribosome-profiling data — the same continuous-synthesis
  framework as this project. No external R is needed; the model sets absolute
  time (transcriptome-wide average anchored to ~200 ms → ~5 codons/s).

### B.2 Which numbers the file uses
- **MEAN** translation time (not median) — matches the mean first-passage
  convention of the Fluitt E. coli table and O'Brien simulations.
- **Weissman/Williams** dataset (364 transcripts), chosen over Nissley (117)
  because its per-codon means rest on ~3× more observations. The two agree
  closely (sense-codon average 233 vs 238 ms) and are well-correlated.
- Values converted ms→s. Range: fastest **GUU 0.145 s**, slowest **CCG 0.504 s**.
- The file header also carries all four series (both datasets × mean/median, in
  ms) so any can be substituted.

### B.3 Reproduce from scratch
1. Open `pcbi.1007070.s016.xlsx`; read sheet "Weissman data translation rates".
2. Take the "Mean translation time (ms)" column for all 64 codons; ÷1000 → s.
3. (Codons already in U alphabet.) Write 6 dp. Stops are measured — keep as given.
No R, no weighting step: the values are already absolute.

### B.4 Caveats
- Per-codon times vary up to ~16–26 fold with sequence context; the tabulated
  value is the mean over all instances.
- Stop-codon times differ markedly between the two datasets (termination is noisy).

---

## References
- Gardin et al. *eLife* 3:e03735 (2014) — RRT, Table 2A.
- Chu et al. *Cell* 153:1589 (2013/2014) — R = 9.3 aa/s.
- Sharma … O'Brien. *PLoS Comput Biol* 15:e1007070 (2019) — S1 File.
- Kazusa Codon Usage Database (if using an external usage table).

---

## Codon usage f_c used for the frequency-weighting

Source: **Gardin et al. 2014, Table 2A** codon usage (per 1000; 61 sense codons only — the RRT table has no stop codons). Used to weight BOTH the Gardin table and the O'Brien-rescaled-9.3 table. Normalized to Σ = 1 over the 61 sense codons to give f_c. Layout: `codon per-1000` (stops n/a).

```
UUU  26.1     UCU  23.5     UAU  18.8     UGU   8.1   
UUC  18.4     UCC  14.2     UAC  14.8     UGC   4.8   
UUA  26.2     UCA  18.7     UAA   n/a     UGA   n/a   
UUG  27.2     UCG   8.6     UAG   n/a     UGG  10.4   

CUU  12.3     CCU  13.5     CAU  13.6     CGU   6.4   
CUC   5.4     CCC   6.8     CAC   7.8     CGC   2.6   
CUA  13.4     CCA  18.3     CAA  27.3     CGA   3.0   
CUG  10.5     CCG   5.3     CAG  12.1     CGG   1.7   

AUU  30.1     ACU  20.3     AAU  35.7     AGU  14.2   
AUC  17.2     ACC  12.7     AAC  24.8     AGC   9.8   
AUA  17.8     ACA  17.8     AAA  41.9     AGA  21.3   
AUG  20.9     ACG   8.0     AAG  30.8     AGG   9.2   

GUU  22.1     GCU  21.2     GAU  37.6     GGU  23.9   
GUC  11.8     GCC  12.6     GAC  20.2     GGC   9.8   
GUA  11.8     GCA  16.2     GAA  45.6     GGA  10.9   
GUG  10.8     GCG   6.2     GAG  19.2     GGG   6.0
```

---

## Data files (in this folder)

The actual `CODON<TAB>seconds` tables are stored alongside this document.
(Identical copies also live in the project root.)

| file | version | mean s/codon | range (sense) | description |
|---|---|---|---|---|
| `yeast_codon_dwell_times.txt` | PRIMARY | 0.1164 | 0.0725–0.1957 | Gardin RRT x R=9.3, frequency-weighted. |
| `yeast_codon_dwell_times_OBrien.txt` | alt | 0.2334 | 0.1450–0.5040 | O'Brien chemical-kinetic ABSOLUTE times (Weissman mean); native ~4.3 aa/s scale. |
| `yeast_codon_dwell_times_OBrien_rescaled_9.3.txt` | alt | 0.1206 | 0.0749–0.2604 | O'Brien pattern uniformly rescaled to 9.3 aa/s (comparable to the Gardin table). |

All tables: RNA alphabet (U), 6 dp, 64 rows (61 sense + UAA/UAG/UGA).

---

## Cross-check: do the two independent yeast tables agree?

The Gardin (RRT) and O'Brien (chemical-kinetic) tables come from different data
and methods, so their *relative* codon pattern is an independent validation.
Over the 61 sense codons:

- **Pearson r = 0.68**, **Spearman ρ = 0.78** (Gardin time vs O'Brien time).

They agree well on which codons are fast/slow — the ~2× difference between them
was purely the absolute-rate anchor (9.3 vs ~4.3–4.8 codons/s), not the pattern.
`yeast_codon_dwell_times_OBrien_rescaled_9.3.txt` removes that scale difference
(uniform ×0.517 so the usage-weighted sense mean = 1/9.3 s), leaving only the
genuine per-codon methodological differences.
