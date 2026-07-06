# E. coli — per-codon translation times (methods & provenance)

**Final file:** `ecoli_codon_dwell_times_310K.txt`
**Type:** mechanistic/absolute, frequency-weighted rescaling.
**One-line summary:** Fluitt et al. 2007 per-codon aa-tRNA *insertion (decoding)*
times, uniformly rescaled by the O'Brien group so the codon-usage-weighted mean
equals the experimental average of 0.061 s/codon (16.5 aa/s at 310 K).

---

## 1. Original source of the per-codon numbers

**Fluitt, A., Pienaar, E., Viljoen, H. (2007)** "Ribosome kinetics and aa-tRNA
competition determine rate and fidelity of peptide synthesis." *Computational
Biology and Chemistry* 31:335–346.

- The per-codon values are in **Table 5**, "Average insertion time (ms)"
  (evaluated at 37 °C / 310 K), one row per codon, all 64 codons.
- This is a bottom-up kinetic model: aa-tRNA diffusion + competition among
  cognate / near-cognate / non-cognate tRNAs → mean time for correct aa-tRNA
  insertion at the A-site. It is the **codon-dependent decoding step**, not the
  full elongation cycle (peptidyl transfer and translocation are treated as
  codon-independent constants elsewhere in the model).
- Example values (ms): AGG 461, GCC 415, CGG 397 (slowest) … GUU 26, UGA 12,
  UAA 11 (fastest).

## 2. How the numbers reached the project file (the rescaling)

The file used by the pipeline is the O'Brien-group redistribution:

**Nature Chemistry (2022), doi:10.1038/s41557-022-01091-z**, Supplementary
Information — Methods Eq. (12) and **Supplementary Table 8** ("Original and
rescaled … data").

They took Fluitt's intrinsic codon times `τ_Fluitt(i)` and rescaled so the
codon-usage-weighted average equals the experimental average `τ_exp`:

```
τ_trans(i) = [ τ_exp / Σ_k w_k · τ_Fluitt(k) ] · τ_Fluitt(i)          (Eq. 12)
```

- `w_k` = normalized E. coli codon-usage frequency.
- `τ_exp = 0.061 s`. Justification (from the SI): the experimental E. coli codon
  translation rate "varies from 12 to 21 s⁻¹ … average of 16.5 s⁻¹" at 310 K, so
  `τ_exp = 1/16.5 ≈ 0.061 s`.
- The bracket is a single constant, so `τ_trans(i) = C · τ_Fluitt(i)` with, for
  this file, **C = 0.501207** (i.e. Fluitt's insertion times ÷ ~2, then ms→s).
- Supplementary Table 8 lists both Fluitt's original values (e.g. UUU 0.136 s =
  136 ms) and the rescaled column — the rescaled column **is** this file.

This is Option B (frequency-weighted) of the general framework (see
`../README.md` §3), applied to already-absolute data.

## 3. Verification performed

Regressing this file's values against Fluitt Table 5:
`your_time(ms) = 0.501207 × insertion_time − 0.00008`, **R² = 1.00000000**, no
codon deviating >0.01 %. The only scatter (ratio 0.50117–0.50121) is Table 5
being rounded to whole ms while the file carries full precision. Means:
all-64 = 0.067608 s, sense-only = 0.070588 s (≈ the stated ~0.068 s ≈ 15 aa/s).

## 4. Recipe to reproduce from scratch

1. Take Fluitt 2007 Table 5 "average insertion time" (ms) for all 64 codons; ÷1000 → s. Convert T→U.
2. Get E. coli codon-usage frequencies `w_k` (normalized, Σ = 1). (In the source, from ref. 31 of the Nat Chem paper; also listed in its Supplementary Table 8.)
3. Compute `Σ_k w_k · τ_Fluitt(k)` (the usage-weighted mean of the Fluitt times).
4. `C = 0.061 / Σ_k w_k · τ_Fluitt(k)` (here C ≈ 0.501207).
5. `τ_trans(i) = C · τ_Fluitt(i)`; write to 6 dp.

**Unweighted alternative (Option A, for consistency with an unweighted convention):**
replace step 3–4 with `C = 0.061 / mean_over_64(τ_Fluitt)`. Not what the shipped
file uses — the shipped file is frequency-weighted.

## 5. Notes / caveats

- The tabulated quantity is the **decoding/insertion** time. In the O'Brien
  model the codon-independent peptidyl-transfer ⟨τ_PT⟩ and translocation ⟨τ_TL⟩
  dwell times (also from Fluitt) are added separately. If your header describes
  the file as "total per-codon time," tighten it to "rescaled decoding time."
- 310 K / 37 °C. This is the warm, fast reference; eukaryote tables are slower
  and often at lower temperature.
- Stop codons carry Fluitt's measured insertion times (rescaled), not placeholders.

## 6. References
- Fluitt, Pienaar, Viljoen. *Comput Biol Chem* 31:335 (2007) — Table 5.
- O'Brien group. *Nat Chem* (2022), doi:10.1038/s41557-022-01091-z — SI Eq. 12, Suppl. Table 8.

---

## Codon usage f_c used for the frequency-weighting

Source: normalized E. coli codon-usage frequency `w_i` as used in the O'Brien rescaling (Nature Chem 2022 doi:10.1038/s41557-022-01091-z, Supplementary Table 8; originally ref. 31 therein). These are the exact weights behind Eq. 12. Layout: `codon w_i` (dimensionless frequency, Σ over 64 ≈ 1), in the standard codon-table order.

```
UUU 0.022499    UCU 0.008660    UAU 0.016380    UGU 0.005260  
UUC 0.016260    UCC 0.008840    UAC 0.012160    UGC 0.006400  
UUA 0.013880    UCA 0.007640    UAA 0.002060    UGA 0.001040  
UUG 0.013430    UCG 0.008810    UAG 0.000250    UGG 0.015170  

CUU 0.011470    CCU 0.007250    CAU 0.012850    CGU 0.020599  
CUC 0.010930    CCC 0.005580    CAC 0.009460    CGC 0.021459  
CUA 0.003930    CCA 0.008480    CAA 0.015100    CGA 0.003700  
CUG 0.051888    CCG 0.022659    CAG 0.029139    CGG 0.005700  

AUU 0.030229    ACU 0.009040    AAU 0.018269    AGU 0.009080  
AUC 0.024589    ACC 0.022869    AAC 0.021479    AGC 0.015900  
AUA 0.004910    ACA 0.007670    AAA 0.033879    AGA 0.002470  
AUG 0.027589    ACG 0.014470    AAG 0.010700    AGG 0.001510  

GUU 0.018409    GCU 0.015530    GAU 0.032319    GGU 0.024439  
GUC 0.015070    GCC 0.025409    GAC 0.019109    GGC 0.028619  
GUA 0.010980    GCA 0.020579    GAA 0.039419    GGA 0.008440  
GUG 0.025879    GCG 0.032749    GAG 0.018199    GGG 0.011270
```

---

## Data files (in this folder)

The actual `CODON<TAB>seconds` tables are stored alongside this document.
(Identical copies also live in the project root.)

| file | version | mean s/codon | range (sense) | description |
|---|---|---|---|---|
| `ecoli_codon_dwell_times_310K.txt` | REFERENCE | 0.0706 | 0.0130–0.2311 | Fluitt insertion times rescaled (O'Brien Eq.12) to a freq-weighted mean of 0.061 s (16.5 aa/s, 310 K). |

All tables: RNA alphabet (U), 6 dp, 64 rows (61 sense + UAA/UAG/UGA).
