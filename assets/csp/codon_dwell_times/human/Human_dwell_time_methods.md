# H. sapiens (human) — per-codon translation times (methods & provenance)

**Final files:**
- v1 (unweighted): `human_codon_dwell_times.txt`
- v2 (frequency-weighted): `human_codon_dwell_times_v2_freqweighted.txt`

**Type:** relative (ribosome-profiling dwell time) → scaled to seconds by R.

---

## 1. Source of the relative data
**Gobet, Weger, Marquis, Martin, Neelagandan, Gachon, Naef (2020)** "Robust
landscapes of ribosome dwell times and aminoacyl-tRNAs in response to nutrient
stress in liver." *PNAS* 117(17):9630–9641. doi:10.1073/pnas.1918145117
(PMC7196831).

- **Supplementary Dataset S1** (`pnas.1918145117.sd01.xlsx`), **sheet 2**,
  column **`Human_Lintner_DT`** — the single-codon dwell time for the human
  (Lintner et al.) ribosome-profiling dataset.
- The supplement is a binary xlsx behind a PMC CAPTCHA; it was retrieved by
  fetching the file in-browser (in-page `fetch` + JSZip unzip of the xlsx,
  parsing `sheet2.xml` against `sharedStrings.xml`).
- `DT` is reported by the authors' GLM as **log2 dwell time, mean-centered per
  ribosomal site**, as a function of codon position relative to the E site.
  Relative dwell factor: **`r_c = 2^DT_c`**. Higher DT = slower.

## 2. Choosing the A-site position (important modeling decision)
The dataset gives DT at many positions; the core E/P/A positions (0,1,2) are
blank for this column. We used the position of **maximal single-codon dwell
signal**, which is **position −3** (spread across codons ≈0.62 in log2 vs <0.4
elsewhere). Two checks supported this being the decoding (A-site) signal:
- It shows the expected **third-position / wobble signature** (NNU codons slow,
  NNG codons fast) — a hallmark of A-site decoding.
- It is the human analogue of the yeast A-site (Gardin position 6).

If a future reading of the authors' coordinate assigns the A site to a different
offset, re-extract `Human_Lintner_DT` at that position and re-run.

## 3. Elongation-rate anchor
**R = 5.6 codons/s** (mammalian; Ingolia, Lareau & Weissman 2011, *Cell* 147:789).
`τ = 1/5.6 = 0.178571 s`. (Human ≈ mammalian: codon dwell-time shapes are
conserved across mammals, so a mammalian R with the human relative table is
legitimate. Alternatives: 5.2 aa/s O'Brien-framework; 4.1 codons/s Tomuro 2024
single-molecule.)

## 4. Conversion

`r_c = 2^DT_c` for the 61 sense codons (undo the log2). Then:

### v1 — unweighted (Option A)
```
t_c = r_c · τ / mean_over_61(r)        mean_over_61(r) = 1.0955
```
Verified unweighted mean of t = τ. Range: fastest **GGG 0.053541 s**, slowest
**CGU 0.366171 s**.

### v2 — frequency-weighted (Option B)
Codon usage `f_c` = **Kazusa H. sapiens** (93,487 CDS; per-thousand), sense
codons renormalized to Σ = 1.
```
r̄_w = Σ_c f_c · r_c       = 1.0397
t_c = r_c · τ / r̄_w
```
Verified `Σ f_c t_c = τ`. Range: 0.056416 … 0.385833 s.

### v1 → v2 relationship
`v2 = v1 × (1.0955 / 1.0397) = ×1.0537` (**+5.4 %**), a single constant; the
relative codon pattern is identical. v2 is preferred: a usage-average human
transcript then translates at exactly R (consistent with the E. coli / yeast
tables). Under v1 such a transcript would translate ~5.4 % too fast.

**Stop codons** (UAA/UAG/UGA): assigned mean = τ (placeholder; the source has no
termination data).

## 5. Reproduce from scratch
1. Retrieve `pnas.1918145117.sd01.xlsx`; parse sheet 2; take `Human_Lintner_DT`
   at **position −3** for the 61 sense codons.
2. `r_c = 2^DT_c`.
3. R = 5.6 → τ = 0.178571.
4. **v1:** `t_c = r_c·τ / mean₆₁(r)`. **v2:** get Kazusa human usage → f_c →
   `r̄_w = Σ f_c r_c` → `t_c = r_c·τ / r̄_w`.
5. Stops = τ. Write U alphabet, 6 dp.

## 6. Caveats
- `Human_Lintner_DT` is a single human cell-line dataset; the mouse-liver columns
  in the same sheet (`Liver_Gobet_*`) are cleaner and, given mammalian conservation,
  a legitimate alternative "human/mammalian" table at the same R.
- The −3 wobble pattern may partly reflect ribosome-profiling ligation bias
  (discussed in the source paper).
- Occupancy/model-derived, not mechanistic; single-codon DTs only (codon-pair
  effects are in Dataset S2 and are not included).

## 7. References
- Gobet et al. *PNAS* 117:9630 (2020) — Dataset S1, `Human_Lintner_DT`.
- Ingolia, Lareau & Weissman. *Cell* 147:789 (2011) — R = 5.6.
- Kazusa Codon Usage Database, *H. sapiens* (93,487 CDS).

---

## Codon usage f_c used for the frequency-weighting (v2)

Source: **Kazusa** Codon Usage Database, *Homo sapiens* [gbpri], **93,487 CDS (40,662,582 codons)**, species 9606, https://www.kazusa.or.jp/codon/ . The 61 sense-codon per-thousand values are renormalized to Σ = 1 to give f_c. Layout: `codon per-thousand(absolute count)`.

```
UUU  17.6(714298)     UCU  15.2(618711)     UAU  12.2(495699)     UGU  10.6(430311)   
UUC  20.3(824692)     UCC  17.7(718892)     UAC  15.3(622407)     UGC  12.6(513028)   
UUA   7.7(311881)     UCA  12.2(496448)     UAA   1.0(40285)      UGA   1.6(63237)    
UUG  12.9(525688)     UCG   4.4(179419)     UAG   0.8(32109)      UGG  13.2(535595)   

CUU  13.2(536515)     CCU  17.5(713233)     CAU  10.9(441711)     CGU   4.5(184609)   
CUC  19.6(796638)     CCC  19.8(804620)     CAC  15.1(613713)     CGC  10.4(423516)   
CUA   7.2(290751)     CCA  16.9(688038)     CAA  12.3(501911)     CGA   6.2(250760)   
CUG  39.6(1611801)    CCG   6.9(281570)     CAG  34.2(1391973)    CGG  11.4(464485)   

AUU  16.0(650473)     ACU  13.1(533609)     AAU  17.0(689701)     AGU  12.1(493429)   
AUC  20.8(846466)     ACC  18.9(768147)     AAC  19.1(776603)     AGC  19.5(791383)   
AUA   7.5(304565)     ACA  15.1(614523)     AAA  24.4(993621)     AGA  12.2(494682)   
AUG  22.0(896005)     ACG   6.1(246105)     AAG  31.9(1295568)    AGG  12.0(486463)   

GUU  11.0(448607)     GCU  18.4(750096)     GAU  21.8(885429)     GGU  10.8(437126)   
GUC  14.5(588138)     GCC  27.7(1127679)    GAC  25.1(1020595)    GGC  22.2(903565)   
GUA   7.1(287712)     GCA  15.8(643471)     GAA  29.0(1177632)    GGA  16.5(669873)   
GUG  28.1(1143534)    GCG   7.4(299495)     GAG  39.6(1609975)    GGG  16.5(669768)
```

---

## Data files (in this folder)

The actual `CODON<TAB>seconds` tables are stored alongside this document.
(Identical copies also live in the project root.)

| file | version | mean s/codon | range (sense) | description |
|---|---|---|---|---|
| `human_codon_dwell_times_v2_freqweighted.txt` | RECOMMENDED (v2) | 0.1882 | 0.0564–0.3858 | Frequency-weighted, R=5.6. Consistent with E. coli/yeast convention. |
| `human_codon_dwell_times.txt` | v1 | 0.1786 | 0.0535–0.3662 | Unweighted, R=5.6. |

All tables: RNA alphabet (U), 6 dp, 64 rows (61 sense + UAA/UAG/UGA).
