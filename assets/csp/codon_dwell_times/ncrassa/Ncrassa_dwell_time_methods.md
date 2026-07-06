# N. crassa — per-codon translation times (methods & provenance)

**Final files** (3 rate anchors × 2 versions):
- v1 (unweighted): `ncrassa_codon_dwell_times_R5.35_low.txt`,
  `..._R6.7_mid.txt`, `..._R8.02_high.txt`
- v2 (frequency-weighted): the same three with `_v2_freqweighted` suffix.

**Type:** relative (ribosome-profiling decoding time) → scaled to seconds by R.
This is the least-certain table (relative values digitized from a figure; and
until recently R was unmeasured — now anchored to a measured range).

---

## 1. Source of the relative data
**Yang, Yu, Zhao, Dang, Wu, Xie, Sachs, Liu (2019)** "eRF1 mediates codon usage
effects on mRNA translation efficiency through premature termination at rare
codons." *Nucleic Acids Research* 47(17):9243–9258. doi:10.1093/nar/gkz710
(PMC6755126).

- The metric is **RCDT** (Relative Codon Decoding Time): the average relative
  A-site (RPF nucleotides 16–18) ribosome dwell time per codon, computed as
  CDT normalized to the **highest-occupied codon CCA (Pro)** → **CCA = 1.0** is
  the slowest reference; all other codons ≤ 1.0. Higher = slower. RCDT is a
  decoding *time*, directly proportional to dwell time (this supersedes the
  earlier RCDR of Yu et al. 2015 from the same lab).
- **The 61 RCDT values are NOT tabulated** anywhere in the paper or its
  supplements — they exist only as the **black dots in Figure 1** (RCDT on the
  right axis; CAI bars on the left; codons grouped by amino-acid family).

## 2. How the values were extracted (figure digitization)
Because the numbers are figure-only, they were pixel-digitized:
1. Rendered the main-text PDF page with Figure 1 at **300 dpi**.
2. Every panel shares the **same right-axis scale**: 0.5 at the plot-box bottom,
   0.75 mid, 1.0 at the top. Calibrated each row's box-bottom (RCDT 0.5) and
   box-top (RCDT 1.0) in pixels.
3. Read each black dot's height against that scale: `RCDT = 0.5 + 0.5 × (dot
   height on the 0–1 box)`.
4. **Consistency check:** CCA lands exactly on 1.0 (the definitional maximum),
   confirming calibration.
- **Estimated read precision: ± ~0.03 in RCDT.** Treat as figure-digitized
  approximations, not exact source values.
- Sanity check on the result: optimal/frequent codons come out fast, rare codons
  (CGA, AUA, CCA) come out slow — matching the paper's central finding.

The digitized RCDT values (10/row) are stored in each output file's header.

## 3. Elongation-rate anchor (now measured, previously an assumption)
Originally there was **no** published absolute elongation rate for N. crassa, so
an early version used an assumed ~5 codons/s. A measured value was then found:

- **Peptide chain elongation rate = 5.35–8.02 amino acids/s per ribosome**
  (growth-condition dependent; a genome-averaged rate).
- Source: **Karpinets et al. 2006**, *BMC Biol* 4:30 (PMID 16953894), Table 3,
  derived from **Alberghina et al. 1975**, *JBC* 250(12):4381 (PMID 124730);
  also **BNID 107872**.
- This is the right kind of anchor (a bulk, genome-averaged rate = the τ = 1/R
  a mean-composition transcript should reproduce).

Three files span the measured range: **R = 5.35 (low), 6.7 (mid, best single
estimate), 8.02 (high)**. τ = 0.186916 / 0.149254 / 0.124688 s respectively.

## 4. Conversion
`r_c = RCDT_c` (already proportional to dwell time; do not invert).

### v1 — unweighted (Option A)
```
t_c = RCDT_c · τ / mean_over_61(RCDT)      mean_over_61(RCDT) = 0.7213
```
Verified unweighted mean of t = τ.
Example (R = 6.7): fastest **ACC 0.1097 s**, slowest **CCA 0.2069 s**.

### v2 — frequency-weighted (Option B)
Codon usage `f_c` = **Kazusa N. crassa** (3,953 CDS; per-thousand), sense codons
renormalized to Σ = 1.
```
r̄_w = Σ_c f_c · RCDT_c      = 0.6818
t_c = RCDT_c · τ / r̄_w
```
Verified `Σ f_c t_c = τ`.

### v1 → v2 relationship
`v2 = v1 × (0.7213 / 0.6818) = ×1.0579` (**+5.8 %**), a single constant per file;
the relative pattern is identical across v1, v2, and all three R anchors. v2 is
preferred for consistency with the E. coli / yeast tables (a usage-average
N. crassa transcript then translates at exactly R).

**Stop codons** (UAA/UAG/UGA): assigned mean = τ (placeholder).

## 5. Reproduce from scratch
1. Render Yang 2019 Figure 1 at ≥300 dpi; digitize the 61 RCDT dots against the
   0.5–1.0 right axis (verify CCA = 1.0).
2. Pick R from the measured range (5.35 / 6.7 / 8.02); τ = 1/R.
3. **v1:** `t_c = RCDT_c·τ / mean₆₁(RCDT)`.
   **v2:** Kazusa N. crassa usage → f_c → `r̄_w = Σ f_c·RCDT_c` →
   `t_c = RCDT_c·τ / r̄_w`.
4. Stops = τ. Write U alphabet, 6 dp.

## 6. Caveats
- Relative values are figure-digitized (± ~0.03 RCDT) — the largest source of
  uncertainty besides R.
- R is a bulk, growth-averaged rate; the range itself (5.35–8.02) is real
  biological/condition variation.
- Profiling / in-vitro-translation temperature (~25–32 °C) is below the 37 °C
  E. coli reference; elongation is temperature-dependent.

## 7. References
- Yang et al. *Nucleic Acids Res* 47:9243 (2019) — RCDT, Figure 1.
- Karpinets et al. *BMC Biol* 4:30 (2006); Alberghina et al. *JBC* 250:4381 (1975); BNID 107872 — measured elongation rate.
- Kazusa Codon Usage Database, *N. crassa* (3,953 CDS).

---

## Codon usage f_c used for the frequency-weighting (v2)

Source: **Kazusa** Codon Usage Database, *Neurospora crassa* [gbpln], **3,953 CDS (2,048,035 codons)**, species 5141, https://www.kazusa.or.jp/codon/ . The 61 sense-codon per-thousand values are renormalized to Σ = 1 to give f_c. Layout: `codon per-thousand(absolute count)`.

```
UUU  11.8(24098)      UCU  11.9(24466)      UAU   8.5(17342)      UGU   3.4(6864)     
UUC  22.1(45213)      UCC  20.0(40934)      UAC  17.5(35768)      UGC   7.7(15795)    
UUA   2.7(5596)       UCA   9.2(18878)      UAA   0.6(1305)       UGA   0.8(1606)     
UUG  14.9(30615)      UCG  14.5(29720)      UAG   0.5(1105)       UGG  13.1(26845)    

CUU  14.2(29177)      CCU  15.1(30898)      CAU   9.5(19363)      CGU   8.9(18193)    
CUC  26.8(54859)      CCC  22.4(45919)      CAC  14.8(30271)      CGC  17.6(36131)    
CUA   6.0(12196)      CCA  12.4(25320)      CAA  17.0(34722)      CGA   7.1(14442)    
CUG  18.3(37407)      CCG  14.6(29827)      CAG  26.0(53342)      CGG   8.5(17492)    

AUU  14.0(28668)      ACU  11.2(22851)      AAU  10.3(21140)      AGU   8.7(17740)    
AUC  26.5(54226)      ACC  24.7(50609)      AAC  27.0(55293)      AGC  17.4(35699)    
AUA   4.1(8373)       ACA  10.7(22008)      AAA  11.7(23932)      AGA   7.9(16197)    
AUG  21.8(44655)      ACG  13.5(27722)      AAG  40.4(82727)      AGG  11.8(24252)    

GUU  13.8(28345)      GCU  21.1(43267)      GAU  24.0(49136)      GGU  18.3(37439)    
GUC  24.8(50857)      GCC  36.0(73666)      GAC  32.5(66663)      GGC  29.0(59424)    
GUA   5.4(11060)      GCA  12.6(25733)      GAA  22.4(45953)      GGA  13.6(27770)    
GUG  15.5(31763)      GCG  17.3(35344)      GAG  42.7(87418)      GGG  10.9(22396)
```

---

## Data files (in this folder)

The actual `CODON<TAB>seconds` tables are stored alongside this document.
(Identical copies also live in the project root.)

| file | version | mean s/codon | range (sense) | description |
|---|---|---|---|---|
| `ncrassa_codon_dwell_times_R6.7_mid_v2_freqweighted.txt` | RECOMMENDED (v2) | 0.1579 | 0.1160–0.2189 | Freq-weighted, R=6.7 (mid of measured 5.35-8.02 range). |
| `ncrassa_codon_dwell_times_R5.35_low_v2_freqweighted.txt` | v2 | 0.1977 | 0.1453–0.2741 | Freq-weighted, R=5.35 (slow end). |
| `ncrassa_codon_dwell_times_R8.02_high_v2_freqweighted.txt` | v2 | 0.1319 | 0.0969–0.1829 | Freq-weighted, R=8.02 (fast end). |
| `ncrassa_codon_dwell_times_R6.7_mid.txt` | v1 | 0.1493 | 0.1097–0.2069 | Unweighted, R=6.7. |
| `ncrassa_codon_dwell_times_R5.35_low.txt` | v1 | 0.1869 | 0.1373–0.2591 | Unweighted, R=5.35. |
| `ncrassa_codon_dwell_times_R8.02_high.txt` | v1 | 0.1247 | 0.0916–0.1729 | Unweighted, R=8.02. |

All tables: RNA alphabet (U), 6 dp, 64 rows (61 sense + UAA/UAG/UGA).
