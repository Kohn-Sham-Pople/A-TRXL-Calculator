# A-TRXL — Phase 1 Step 1.2: Energy-dependent Pair Weights

**Report ID:** ATRXL-P1-S12-NATIVE-20261105  
**Date:** 2026-11-05  
**Scope:** Phase 1 **Step 1.2 only**.  
**Status:** **NUMERICAL PASS** — native XrayDB 4.5.8 output supplied by the user, independently audited against the prior frozen-upstream dataset. User-side full pytest (`24 passed`) output was **not** included in this submission and must not be represented as observed.  
**Primary provenance:** User-executed Spyder/Anaconda Python, `--engine xraydb`, recorded in `calculation_metrics.json` and the uploaded NetCDF. The present audit was executed in a separate Python environment by reading the actual uploaded files; it did **not** independently invoke the XrayDB API.

## 1. Research objective and strict scope

The *A-TRXL Python Plan.pdf*, Phase 1 Step 1.2, requires evaluating complex Br-containing pair weights on a physical-q grid using the actual Br absorption-edge energy contrast:

\[
 f_a(q,E) = f_{0,a}(q)+ f'_a(E)+if''_a(E),\qquad
 W_{ab}(q,E)=\operatorname{Re}\left[f_a(q,E) f_b(q,E)^*\right].
\]

For identically tabulated Br scatterers, \(W_{\mathrm{Br-Br}}=|f_{\rm Br}|^2\). For other pairs \(W_{\mathrm{Br-X}}=\operatorname{Re}(f_{\mathrm{Br}}f_X^*)\), with \(X\in\{\mathrm{C,O,H}\}\). Each quantity describes one **unordered distinct pair**, excluding both the Debye factor of two and the structural sinc kernel. These pair weights **must not** be equated with total measured intensity, time-resolved double-difference intensity, or experimental SNR.

Energy difference convention:

\[
\Delta W_{ab}(q;E_{\rm near},E_{\rm far}) = W_{ab}(q,E_{\rm near})-W_{ab}(q,E_{\rm far}).
\]

**Units:** E in eV for API calls (keV on plots), physical \(q\) in Å⁻¹, \(f\) in electrons, \(W\) and \(\Delta W\) in electrons². XrayDB `f0()` takes \(s=q/(4\pi)=\sin\theta/\lambda\), not the physical \(q\). Include both \(f'\) and \(f''\) in the complex product.

## 2. Data lineage and execution evidence

- **Primary:** actual locally executed `xraydb.f0`, `xraydb.f1_chantler`, `xraydb.f2_chantler` for Br, C, O, H; user-uploaded `contrast_weights(1).nc` and `calculation_metrics(1).json`. XrayDB **4.5.8** is explicitly recorded; Python **3.14.6**, NumPy **2.4.6**, SciPy **1.18.0**, xarray **2026.7.0**, pandas **3.0.3**, Matplotlib **3.11.0**.
- **Independent implementation comparison:** prior Step 1.2 `frozen-upstream` output, which uses Step 1.1 actual XrayDB Br factor output plus upstream XrayDB Waasmaier–Kirfel coefficients and Chantler C/O/H reference excerpts reconstructed with specified interpolation. This is an **implementation crosscheck**, not an independent theory of the underlying atomic data.
- **Henke_Br_swap:** replace *only Br* anomalous factors with CXRO/Henke Br values; retain native XrayDB C/O/H factors.
- **NIST_Br_swap:** replace *only Br* anomalous factors with archived NIST/Chantler Br values; retain native XrayDB C/O/H factors. It is **not** a third fully independent all-element tabulation.
- Five frozen reference-input files have SHA-256 hashes stored in the supplied metrics; the audit verified that **all five matched** the bundled reference files.
- The actual XrayDB installation was exercised on the user's PC; the results and metadata were then checked **in a separate runtime**. This report does not infer a user's pytest result absent a screenshot or execution log.

The source calculation covers 404 pre-edge energies **11.830–13.460 keV**; benchmark energies are **12.00, 13.24, 13.45 keV**. The physical trial grid is 291 q-values **0.2–6.0 Å⁻¹**, step **0.02 Å⁻¹**, with q=0 included separately for a diagnostic limit. No interpolation crosses the Br K-edge. The three listed energy contrasts are 13.24−12.00, 13.45−12.00, and 13.45−13.24 keV; only two are linearly independent.

## 3. Results — native XrayDB/Chantler

### 3.1 Benchmark at q = 4 Å⁻¹; 13.45 − 12.00 keV

| Pair | \(W_{\rm far}\) (e²) | \(W_{\rm near}\) (e²) | \(\Delta W\) (e²) | Relative \(\Delta W/W_{\rm far}\) |
|:---|---:|---:|---:|---:|
| Br–Br | 458.712626 | 303.436780 | **−155.275847** | **−33.850%** |
| Br–C | 50.778979 | 41.264486 | **−9.514494** | **−18.737%** |
| Br–O | 82.998096 | 67.414422 | **−15.583673** | **−18.776%** |
| Br–H | 4.749953 | 3.862886 | **−0.887066** | **−18.675%** |

The Br anomalous correction difference \(f'_{\rm Br}(13.45)-f'_{\rm Br}(12.00)\) is **−3.998119896 e**. For Br–Br at q=4 Å⁻¹, the \(f''\)-squared term contributes **−0.068447549 e²** to \(\Delta W\), retained despite being a small fraction of the full contrast.

### 3.2 q dependence

For 13.45−12.00 keV, all **Br–Br, Br–C, Br–O, Br–H** weights have **negative \(\Delta W\)** across the entire trial physical q-grid in all three Br-tabulation channels. For example, Br–Br relative contrast is approximately **−22.72% at q=0.2**, **−33.85% at q=4**, and **−41.48% at q=6 Å⁻¹**. This result describes only the **weight**, not the experimental scattering signal or signal-to-noise ratio.

At q=4 Å⁻¹, the Henke Br-only swap gives \(\Delta W_{\rm Br-Br}\approx-159.31576\) e², and NIST archived Br-only swap \(\approx-164.81243\) e², compared with primary **−155.27585 e²**. Thus the source choice changes the quantitative result but not the sign here. Chemical-state dependence and near-edge uncertainty remain open.

## 4. Independent numerical audit of native artifacts

**Actual input files:** `contrast_weights(1).nc` and `calculation_metrics(1).json`, supplied by the user; separate frozen reference: the Step 1.2 source-traceable NetCDF. Rebuilt verification script: `src/verify_native_artifacts.py` in the final bundle.

### 4.1 Native vs frozen-upstream over *entire* multidimensional grids

| Comparison | Maximum absolute difference |
|:---|---:|
| \(f_0(q)\) | \(7.106\times 10^{-15}\) electrons |
| \(f'(E)\) | \(1.714\times 10^{-10}\) electrons |
| \(f''(E)\) | \(1.735\times 10^{-18}\) electrons |
| \(W(q,E)\) | **\(4.921\times 10^{-9}\) e²** |
| \(\Delta W(q)\) | **\(4.437\times 10^{-9}\) e²** |

All 7 coordinate grids match exactly between both NetCDF files. The differences in C/O/H \(f'\), and hence \(W\), are compatible with the limited-knot reconstructed interpolation versus the native full XrayDB implementation. Their small magnitude is confirmed; **causal attribution solely to spline knots has not been proven by a separate upstream-grid ablation**. No difference has been silently rounded to zero in the audit.

### 4.2 Independently executed checks

**124 checks passed** in a fresh Python read-and-verify audit on the actual uploaded arrays:

- Coordinate agreement, sizes, 404 ascending energy points, 291 physical-q trial points and presence of benchmark energies;
- All stored arrays finite; exact input reference SHA-256 matches for all 5 source inputs;
- \(f_0(q\to0)\) within 0.01 electron of atomic Z for each Br/C/O/H; numerical values: **Br 34.995982, C 5.997198, O 7.999706, H 0.999978**;
- Native C/O/H \(f'\), \(f''\) identical across Br-only source-swap channels;
- Independently rebuilt \(\operatorname{Re}[f_a\operatorname{conj}(f_b)]\) matches the full saved \(W\) arrays for all source channels and 4 pair types, **maximum absolute residual 0 e²**;
- Each of the 3 stored \(\Delta W\) channels reproduces direct subtraction of the stored \(W(E)\) arrays, with maximum residual under \(10^{-11}\) e²;
- Three-energy closure, q=0 two-equivalent-Br closure: both **0 e²** within saved-number precision;
- q=4 Å⁻¹ numerical benchmarks reproduce the user-supplied metrics; full-grid primary source contrast has expected sign.

**Other execution evidence:** The patched source's `pytest` suite ran in this *audit environment* with **23 passed, 1 skipped**; the skip is because `xraydb` is not installed here. In earlier development, a test-only mock native path had passed 24 tests, but **this is not a user's actual native pytest run**. The user's new upload proves native calculation and the dataset-level consistency checks, **not** that their PC's pytest run returned `24 passed`. Record this transparently.

## 5. Validation decision

**Step 1.2: NUMERICAL PASS (native runtime confirmed through saved output; independent full-data audit successful).** The source-traceable q-dependent pair-weight result is robust to the tested XrayDB implementation method and preserves its contrast direction under Br-only source substitution.

**Outstanding, not grounds for claiming new results:**

1. User's **native pytest execution log** not supplied; 24/24 must not be reported as observed.
2. End-to-end independent all-element Henke anomalous calculation not carried out; the alternate sources are **Br-only swaps**.
3. Realistic **non-Br-only residual** calculations (C–C, C–O, O–O, O–H) belong to **Step 1.3**, which has **not** been executed.
4. Energy-pair/q optimization and full source/chemical-state sensitivity belong to Steps 1.4 and 1.5.
5. Geometry changes, actual \(\Delta I(q,t)\), structural identifiability and experimental SNR are **not verified** by these pair-weight calculations.

## 6. Reproducibility & deliverables

Final bundle contains:

- `src/step12_pair_weights.py` — native and frozen engines; `--engine xraydb` for a fresh native calculation;
- `src/verify_native_artifacts.py` — **executed independent array audit**, requiring only NumPy, xarray and a SciPy NetCDF reader;
- `tests/test_step12_pair_weights.py` — patched test suite;
- `data/reference/` — five frozen reference inputs with checksums in native metrics;
- `results/tables/contrast_weights.nc` — **user-executed native XrayDB results**;
- `results/reference/frozen_upstream_contrast_weights.nc` — independently calculated prior frozen comparison;
- `results/tables/calculation_metrics.json`, `numeric_audit.json` and benchmark CSV tables;
- `results/figures/01–06` — figures regenerated from the **native** dataset, with figure 06 showing native−frozen numerical discrepancies;
- `results/logs/independent_native_audit_run.txt`, source/README, and `SHA256SUMS.txt`.

For local rerun from the project root:

```bash
python src/verify_native_artifacts.py
python -m pytest -q -rs tests/
# Optional: creates / overwrites primary results after backup:
python src/step12_pair_weights.py --engine xraydb
```

Keep the uploaded native NetCDF unmodified as a reference; the bundle includes the file and its checksum. The old standalone report describes historical `frozen-upstream` work; **use this report as the final Step 1.2 handoff**. Step 1.3 is not automatically authorized.

## 7. Sources

- XrayDB official Python API: https://xraypy.github.io/XrayDB/python.html . `f0` uses \(\sin\theta/\lambda\); `f1_chantler`, `f2_chantler` use eV and electron units.
- D. Waasmaier and A. Kirfel (1995), *Acta Crystallographica A* **51**, 416–431, DOI https://doi.org/10.1107/S0108767394013292 — five-Gaussian form-factor representation.
- M. Saito and Y. Waseda (2000), *J. Synchrotron Radiation* **7**, 152–159, DOI https://doi.org/10.1107/S0909049500002971 — physical basis of anomalous contrast in liquids.
- XrayDB database construction source: https://github.com/xraypy/XrayDB/blob/master/data_sources/create_db.py .
- Reference lineages and hashes are explicitly recorded in `data/reference/provenance_step12.json` and the uploaded `calculation_metrics.json`.

**Scope boundary:** **Step 1.2 ends here.** No Step 1.3 residual, experiment, MD, or photochemistry simulation has been carried out.
