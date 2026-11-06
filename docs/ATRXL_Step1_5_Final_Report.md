# A-TRXL — Phase 1 Step 1.5: Near-edge Uncertainty Test

**Report ID:** ATRXL-P1-S15-20261105  
**Date:** 2026-11-05  
**Status:** **NUMERICAL PASS — scoped Step 1.5**  
**Scope:** Br K pre-edge **neutral-atom scattering-factor and pair-weight sensitivity** to tabulation choice, interpolation rule, and explicitly assumed energy-axis displacement. No measured chemical shift, real photochemical trajectory, experimental uncertainty budget, or instrument SNR is claimed.

## 1. Executive findings

1. **Source robustness at the original benchmark:** at physical **q = 4 Å⁻¹** and **13.45–12.00 keV**, the Br–Br pair-weight difference is **−155.275847 e²** using native XrayDB/Chantler, **−159.315760 e²** with a Henke **Br-only swap**, and **−164.812430 e²** with a NIST archived **Br-only swap**. The maximum magnitude difference relative to the primary value is **6.14%**, while **all four Br pairs keep negative contrast over q = 0.2–6 Å⁻¹** for the three source channels.
2. **The very near-edge window is markedly more source sensitive:** at **13.460–12.000 keV**, the Br–Br q=4 values are **−174.889**, **−237.130**, and **−175.404 e²**, respectively. Henke differs by **35.59%** from primary XrayDB at this closer-to-edge energy; this is a tabulation-dependent prediction, **not** a measured edge shift.
3. **Interpolation is non-negligible and method dependent:** at 13.45–12.00 keV, replacing linear f′ interpolation of the *original sparse Br raw knots* with pre-edge PCHIP f′ changes q=4 Br–Br ΔW by **1.35% (Henke)** and **0.43% (NIST archived)**. The corresponding change for **13.45–13.24 keV** is **1.21% and 2.79%**, respectively. Unconstrained cubic splines can produce substantially larger excursions on sparse knots; these are stress tests, **not estimates of truth or statistical errors**.
4. **Hypothetical spectral displacement can matter:** applying **f_Br(E) → f_Br(E−δ)**, δ ∈ {−10, −5, −2, 0, +2, +5, +10} eV, while leaving non-Br factors unchanged, changes the q=4 primary Br–Br ΔW(13.45–12.00 keV) from **−174.638** (δ=−10 eV) to **−142.745 e²** (δ=+10 eV), compared with **−155.276 e²** at δ=0. The maximum excursion is **12.47%**; for the 13.45–13.24 keV contrast it reaches **21.01%**. **These shifts are deliberately prescribed numerical scenarios, not measured values for CBr₄, Br₂, atomic Br, Br⁻, or HgBr₂.**
5. **Step 1.4's coefficient-level tradeoff remains in the tabulation swap:** at q=4, the per-pair |ΔW_OO|/|ΔW_BrBr| is **0.02383%** for 13.45–12.00 keV and **0.005512%** for 13.45–13.24 keV using primary data. All three Br-only sources preserve this *relative coefficient ranking*. This **does not rank real experimental energy pairs** without transmission, bandwidth, photon counts, covariance and structural responses.

**Decision:** The Phase 1 source-based claim that sizable **Br-selective atomic pair contrast exists** is robust at 12.00/13.45 keV across the tested neutral-atom tabulations and explicitly declared stress tests. **Quantitative near-edge pair weights are not source invariant**, especially closer to the K-edge. The project can conclude its Phase 1 **atomic-factor numerical work**, with explicit unresolved experimental and chemical-state uncertainties. No Phase 2 computation was performed.

## 2. Question, source lineage, and boundaries

The project's *A-TRXL Python Plan.pdf*, Phase 1 Step 1.5, requests **near-edge uncertainty testing**: quantify the effect of Henke/Chantler selection and interpolation on anomalous contrast. The Phase 0 analytical derivation established

\[
 f_a(q,E)=f_{0,a}(q)+f'_a(E)+if''_a(E),\quad
 W_{ab}(q,E)=\operatorname{Re}[f_a(q,E)f_b(q,E)^*],\quad
 \Delta W_{ab}=W_{ab}(q,E_n)-W_{ab}(q,E_f).
\]

**Native primary input:** `data/reference/step12_native_xraydb_4_5_8.nc`, a byte-preserved copy of user-generated **XrayDB 4.5.8** Step 1.2 output from Spyder/Anaconda. It has 404 pre-edge energies (11,830–13,460 eV), q = 0 plus 0.2–6.0 Å⁻¹ in 0.02 Å⁻¹ steps, f₀/f′/f″ for Br, C, O, H, and three Br source channels. Its matching user metrics are preserved. **No new native XrayDB import was performed in Step 1.5.**

**Alternate Br source data:** the traceable Step 1.1 original *CXRO/Henke* Br raw-point excerpt (near-edge f1_total, f2) and *NIST FFAST.jl archived Chantler* Br raw-point excerpt, both copied verbatim from the Step 1.1 validated bundle. The Step 1.2 saved Henke and NIST source-swap channels are independently reproduced at nominal benchmark energies from these excerpts.

**Source-swap principle:** replace **only Br** anomalous f′/f″; retain **identical XrayDB C/O/H** factors and f₀ across all channels. Thus Henke vs NIST archive vs XrayDB is a Br-tabulation sensitivity comparison, **not three independent complete all-element theoretical models**. NIST archived Chantler and XrayDB's refined Chantler are also **not statistically independent physical experiments**.

All numerical values are **atomic neutral-element coefficients** in electron/electron² units, **per one unordered distinct pair**, excluding multiplicity 2, geometry sinc, photoproduct populations, interatomic correlations, sample concentration, instrument effects or photon noise. q = 0 is a forward-limit diagnostic, **not a measured angular point**. Energies are in eV internally; XrayDB f₀ expects s=q/(4π) in Å⁻¹, already applied upstream.

## 3. Executed methodology

### 3.1 Source selection

At E_n=13,450 and E_f=12,000 eV, and separately E_f=13,240 eV, compute ΔW on all 292 q points for pairs Br–Br, Br–C, Br–O, Br–H from the original native Step 1.2 W values. Test all three sources at identical physical q with all complex factors retained; compare source spread to the primary. The full raw input data are saved unchanged and linked by SHA-256.

### 3.2 Original-knot interpolation test

The exact Henke raw pre-edge points extend through 13,473.6 eV (first post-edge point 13,473.8). NIST archival pre-edge points extend through 13,460.23 eV (next post-edge 13,476.97). **Interpolation across these tabulated discontinuity brackets is forbidden.**

Reference processing replicates Step 1.1: **linear** interpolation of full real f1 versus E, **log-log** interpolation of positive imaginary f2. Alternative processing replaces only f1 with **shape-preserving PCHIP**, leaving the f2 rule unchanged. An additional **cubic spline** result is saved to expose unconstrained sparse-grid overshoots but is not recommended as an accurate near-edge model. Raw full f1 is converted as f′=f1−35; XrayDB's *native* f1_chantler is already the correction f′. The raw-knot analysis tests **algorithmic interpolation sensitivity**, not independent empirical correctness.

### 3.3 Synthetic energy-axis stress test

A seven-value predeclared δ grid (−10 to +10 eV) is applied only to the **Br** energy-dependent factor:

\[
 f_{\mathrm{Br}}^{(\delta)}(q,E)=f_{0,\mathrm{Br}}(q)+f'_{\mathrm{Br,primary}}(E-\delta)+if''_{\mathrm{Br,primary}}(E-\delta).
\]

Non-Br C/O/H factors are evaluated at nominal E. For off-grid E−δ, the script uses PCHIP **on the dense, already evaluated native XrayDB output curve**, which is **not the same thing as rerunning XrayDB on its fine raw knots**. The second-order sensitivity to interpolation of the stored curve is separately compared between PCHIP, linear and cubic. All shifted evaluations stay within the saved **11,830–13,460 eV** domain; no extrapolation above the Br K-edge is allowed. This intentionally simplified shift cannot model real state-specific f′/f″ line-shapes, XANES resonances or an on/off electronic response.

## 4. Numerical results

### 4.1 Tabulation choice at q=4 Å⁻¹

| Energy contrast | Pair | XrayDB ΔW (e²) | Henke Br-only ΔW (e²) | Archived NIST Br-only ΔW (e²) | Max source difference vs XrayDB |
|:--|:--|--:|--:|--:|--:|
| 13.45−12.00 | Br–Br | −155.275847 | −159.315760 | −164.812430 | 6.142% |
| 13.45−12.00 | Br–C | −9.514494 | −9.808143 | −10.096150 | 6.113% |
| 13.45−12.00 | Br–O | −15.583673 | −16.063769 | −16.534439 | 6.101% |
| 13.45−12.00 | Br–H | −0.887066 | −0.914543 | −0.941484 | 6.135% |
| 13.45−13.24 | Br–Br | −85.704209 | −87.540460 | −92.128586 | 7.496% |
| 13.45−13.24 | Br–C | −5.477104 | −5.630839 | −5.896362 | 7.655% |
| 13.45−13.24 | Br–O | −8.952139 | −9.203585 | −9.637373 | 7.654% |
| 13.45−13.24 | Br–H | −0.512333 | −0.526713 | −0.551568 | 7.658% |

The signs remain negative for all tested physical q points for both contrast pairs and source channels. **A per-pair ΔW sign or magnitude is not a double-difference intensity, inverse-problem result, or SNR estimate.**

### 4.2 Closer near-edge point: 13.460 keV

Holding E_far=12.00 keV, q=4 Å⁻¹:

| E_near | XrayDB Br–Br ΔW | Henke Br–Br ΔW | Archived NIST Br–Br ΔW |
|--:|--:|--:|--:|
| 13.450 keV | −155.276 | −159.316 | −164.812 |
| 13.460 keV | −174.889 | −237.130 | −175.404 |

Henke Br f′ is **−8.4469 e** at 13.460 keV under its reference interpolation; native XrayDB is **−6.4301 e**, archived NIST is **−6.2929 e**. These are source-specific modeled values, **not experimental chemical-state shifts**. This steep divergence reinforces the caution in the original project idea about measuring too close to a sharp edge.

### 4.3 Original-knot interpolation variation

| Br source | ΔW_BrBr(13.45−12.00) change PCHIP vs linear at q=4 | ΔW_BrBr(13.45−13.24) change PCHIP vs linear at q=4 |
|:--|--:|--:|
| Henke | +2.144813 e² (**1.346%**) | +1.061521 e² (**1.213%**) |
| NIST archived | +0.704914 e² (**0.428%**) | −2.571960 e² (**2.792%**) |

For the sparse original knots, cubic spline behavior can differ much more significantly (values exported). **Do not classify PCHIP or cubic deviations as calibrated 1σ uncertainty or as evidence that one method is true.** The archived NIST pre-edge contains very few points close to the jump. The native XrayDB fine source uses its own original grid and spline/log-log rules; the original raw-source interpolation experiment here is a distinct sensitivity test.

### 4.4 Energy-axis shift stress test

At q=4 Å⁻¹, 13.45−12.00 keV, native XrayDB Br-only shifted factors:

| Assumed δ (eV) | ΔW_BrBr (e²) | Change versus δ=0 |
|--:|--:|--:|
| −10 | −174.638 | +12.469% in magnitude |
| −5 | −163.767 | +5.468% |
| −2 | −158.466 | +2.055% |
| 0 | −155.276 | 0 |
| +2 | −152.346 | −1.887% |
| +5 | −148.400 | −4.428% |
| +10 | −142.745 | −8.070% |

At the alternative 13.45−13.24 keV pair, the largest tested change in Br–Br ΔW magnitude is **21.01%**. Over the tested ±10 eV scenarios, all four Br-related ΔW retain negative signs on all physical q points. **This is a parametric stress test, not an asserted physical range for any sample**.

### 4.5 Non-Br leakage context

From the same native source, q=4 Å⁻¹, O–O coefficient residual:

| Contrast | |ΔW_OO| / |ΔW_BrBr| primary | Henke Br-only | NIST Br-only |
|:--|--:|--:|--:|
| 13.45−12.00 | 0.023829% | 0.023224% | 0.022450% |
| 13.45−13.24 | 0.005512% | 0.005396% | 0.005128% |

The *coefficient-level* leakage ranking remains stable under source swaps; **it does not guarantee real solution heating suppression**, because solvent pair counts, geometry response, background and normalization are not modeled here.

## 5. Validation and audit

1. **132 pytest tests passed, 0 failed**, using the actual saved primary input and source raw-knot excerpts. They include all source/pair/contrast direct-difference checks, nominal zero-shift closure, chemical-swap invariance of non-Br factors, exact raw-source reproduction of historical benchmark values, no-interpolation-across-edge guards, q-range and finite data, exact complex-conjugation algebra and energy-contrast additivity. Two initial test failures were due to mistakenly choosing *within-range* shifts as "must reject" test cases; test input values were corrected and then all tests passed. **No numerical production calculation was altered to hide these test failures.**
2. **175 / 175 read-only independent numerical/provenance audit checks passed.** These reconstruct stored ΔW from original saved f₀/f′/f″ with independent complex arithmetic, verify source-swap invariance, original input hashes, zero-shift comparison, selected off-grid shifts and the archived interpolation CSV, with **0 e² maximum direct complex-formula residual** in the saved output precision.
3. Generated a new NetCDF (`step15_sensitivity.nc`) and reloaded it with xarray equality; all arrays finite. Source- and shift-dependent ΔW are supplied on the full q grid; original-knot interpolation summaries are exported to CSV at q=4 Å⁻¹.
4. **Provenance:** native XrayDB **4.5.8 was executed in Step 1.2 on the user's PC**, not in this Step. Step 1.5 calculations were **actually executed using its frozen factor arrays** in this Python environment. No XAS/XANES measurement, calibrated chemical shift or spectroscopy-specific line-shape fitting has been performed.

## 6. Scientific interpretation and decision gate

**Robust conclusion:** the existence and sign of Br anomalous pair-weight contrast at the proposed 12.00/13.45 keV pre-edge pair survive all source choices and tested sensitivity cases. This is a **numerical statement about neutral-atom coherent pair coefficients**; it does **not** demonstrate observable photoexcited Br-selective ΔΔI or arbitrary chemical-state invariance.

**Non-robust quantitative conclusions:** the precise near-edge coefficient magnitude (especially approaching 13.460 keV) and any declaration of an *optimal* energy pair. The source envelope is a set of model predictions, **not** a probability distribution or uncertainty confidence interval. No energy resolution/fluence/absorption/measurement covariance has been supplied for instrument-aware optimization.

**Chemical caveat:** published Br K-edge X-ray absorption spectroscopy demonstrates dependence on chemical bonding or compound identity, but these sources do not provide a measured f′/f″ curve or pump-on/off contrast for the actual proposed CBr₄/methanol or HgBr₂/methanol photoproducts. Replacing a real chemical-state spectrum with a rigid ±δ translation is merely a failure-mode **proxy**.

**Recommendation for subsequent phases:** use `XrayDB/Chantler` as primary and `CXRO/Henke` as independent Br comparison, `NIST/Chantler archived` as additional sensitivity reference. Preserve source- and interpolation-dependent envelopes and explicitly vary chemical-state assumptions in future forward/inverse calculations. Before a near-edge experimental proposal, obtain chemically matched Br XANES/fluorescence or independently calibrated anomalous factors and instrument parameters; do not claim instrument feasibility from the present numbers.

**Go/No-Go:** **Step 1.5 NUMERICAL PASS.** Phase 1's required factor/convention/contrast/residual/sensitivity analyses (Steps 1.1–1.5) are complete at the **neutral-atom computational level**, subject to the limitations above. **No Phase 2 forward scattering simulation was executed.**

## 7. Deliverables and reproducibility

- `src/step15_nearedge_uncertainty.py` — executed Step 1.5 numerical calculator, declared source/shift conventions, figures and tables.
- `src/audit_step15.py` — independently executed read-only artifact audit.
- `tests/test_step15_nearedge_uncertainty.py` — 132 passed pytest checks.
- `data/reference/step12_native_xraydb_4_5_8.nc` and `step12_calculation_metrics.json` — original frozen primary source and run metadata, including XrayDB 4.5.8.
- `data/reference/henke_br_preedge_excerpt.csv`, `nist_ffast_chantler_br_excerpt.csv` — traceable original Br tabulated excerpts.
- `results/tables/step15_sensitivity.nc` — source-dependent and synthetic-shift-dependent Br pair-weight sensitivity over q.
- `results/tables/source_sensitivity_q4.csv`, `raw_knot_interpolation_comparison_q4.csv`, `near_edge_source_divergence_q4.csv`, `br_shift_stress_test_q4.csv`, `source_dependent_leakage_ratio_q4.csv`, `shift_densegrid_interpolation_check_q4.csv`, `calculation_metrics.json`.
- `results/figures/01–05` and execution/pytest/independent audit logs, with SHA256SUMS.

Rerun from bundle root (Python 3.9+ with numpy, scipy, pandas, xarray, matplotlib, pytest):

```bash
python src/step15_nearedge_uncertainty.py
python -m pytest -q -rs tests/
python src/audit_step15.py
```

No XrayDB install is needed to **reproduce this Step 1.5 from preserved native outputs**; a fresh live-XrayDB extraction requires that package and explicit revalidation of data/package versions.

## 8. Authoritative sources

1. **XrayDB official Python API:** https://xraypy.github.io/XrayDB/python.html — `f0(s)`, `f1_chantler(E)`, `f2_chantler(E)` conventions; documentation identifies release 4.5.8 and SQLite version 6.
2. **XrayDB database scheme:** https://xraypy.github.io/XrayDB/dbschema.html — Chantler table stores correction f′ and f″.
3. **CXRO Henke Br raw table:** https://henke.lbl.gov/optical_constants/sf/br.nff ; Henke, Gullikson & Davis, *At. Data Nucl. Data Tables* **54**, 181–342 (1993), DOI https://doi.org/10.1006/adnd.1993.1013 .
4. **NIST archived Chantler source:** https://github.com/usnistgov/FFAST.jl/blob/master/data/mac%5B35%5D.csv . Original Chantler, *J. Phys. Chem. Ref. Data* **24**, 71–643 (1995), DOI https://doi.org/10.1063/1.555974 . NIST history https://www.nist.gov/pml/x-ray-form-factor-attenuation-and-scattering-tables/x-ray-form-factor-attenuation-and-scattering (v2.2, 2026 backend rewrite).
5. **Saito & Waseda (2000)**, *J. Synchrotron Radiation* **7**, 152–159, DOI https://doi.org/10.1107/S0909049500002971 — static AXS context for liquid partial functions, not time-resolved observability proof.
6. **Stephenson (1936)**, *Phys. Rev.* **50**, 790, DOI https://doi.org/10.1103/PhysRev.50.790 — Br K-edge shifts/widths and compound-dependent absorption; does **not** determine shifts for the proposed TRXL photoproducts.
7. **Feiters, Küpper & Meyer-Klaucke (2005)**, *J. Synchrotron Radiation* **12**, 85–93, DOI https://doi.org/10.1107/S0909049504027815 — organic Br/iodine-compound XAS context. Not a measured f′/f″ table for CBr₄ photoproducts.

**Integrity note:** All links, claims and method limitations above are recorded for reproducibility. Results are numerical predictions under specified isolated-atom assumptions, not experimental discovery.
