# A-TRXL — Phase 1 Step 1.3: Non-Br Residual Calculation

**Report ID:** ATRXL-P1-S13-20261105  
**Status:** **NUMERICAL PASS — Step 1.3 only**  
**Date:** 2026-11-05  
**Scope:** C, O, H non-Br-only *atomic coherent pair weights* over the Step 1.2 pre-edge energy/physical-q grid. This is **not** a solution-scattering or experimental cancellation demonstration.

## 1. Scientific question and source lineage

From the user-provided *A-TRXL Python Plan.pdf*, Phase 1 Step 1.3: calculate realistic non-Br residual weights for C–C, C–O, O–O, O–H and related pairs rather than assuming the non-Br-only contribution is mathematically zero. The underlying testable Phase 0 hypothesis is **P0-009** of `phase0_hypothesis_register.csv`.

- **Primary input:** `data/reference/step12_native_xraydb_4_5_8.nc`, a byte-preserved copy of the actual Step 1.2 NetCDF generated on the user's Spyder/Anaconda installation using **XrayDB 4.5.8** (`--engine xraydb`). Its provenance is independently documented in the Step 1.2 final report and the paired metrics JSON. Step 1.3 **reuses these native values**; it does **not** claim a fresh live `xraydb` import in the present calculation runtime.
- **Atomic factors:** the real user-generated `f0(q)`, `f1_chantler` (anomalous correction `f′(E)`) and `f2_chantler` (`f″(E)`) outputs for **C, O, H**. Br pair weights are used **only as a benchmark reference** from the same Step 1.2 dataset.
- **Independent references:** official XrayDB documentation https://xraypy.github.io/XrayDB/python.html ; XrayDB SQLite convention https://xraypy.github.io/XrayDB/dbschema.html ; Saito and Waseda, *J. Synchrotron Radiat.* **7**, 152–159 (2000), DOI https://doi.org/10.1107/S0909049500002971 ; Waasmaier and Kirfel, *Acta Cryst. A* **51**, 416–431 (1995), DOI https://doi.org/10.1107/S0108767394013292 .

### Numerical scope and units

- **Elements:** C, O, H; unordered distinct-atom pair types **C–C, C–O, C–H, O–O, O–H, H–H**. An X–X pair means **two separate atoms**, not a self term.
- **Energy domain:** **404 values, 11,830–13,460 eV** taken from the actual Step 1.2 source file; analyzed benchmark contrasts: **13.24−12.00, 13.45−12.00, 13.45−13.24 keV**. No interpolation crosses Br K-edge.
- **Physical q:** **0.2–6.0 Å⁻¹ in 0.02 Å⁻¹ increments (291 points)**, plus `q=0` (analytical diagnostic, not part of experimental window).
- **Units:** energy eV, q Å⁻¹, factors electrons, coherent pair weights and pair-weight residuals **electrons² per distinct unordered pair**. The Debye multiplicity factor 2, geometry sinc factor, atomic counts, normalization, flux, detector effects and noise are **not included**.
- **XrayDB convention:** `f0` takes `s = q/(4π)`, not physical q directly. The input Step 1.2 dataset had already applied this convention and was checked against neutral-atom q→0 limits.

## 2. Methodology

For each atom type `a` at a common physical q and photon energy E,

\[
f_a(q,E)=f_{0,a}(q)+f'_a(E)+i f''_a(E)=a_a(q,E)+i b_a(E).
\]

The *real coherent* weight (excluding the Debye double-counting factor) is

\[
W_{XY}(q,E)=\operatorname{Re}\left[f_X(q,E)f_Y^*(q,E)\right]
=a_Xa_Y+b_Xb_Y .
\]

The exact finite difference is

\[
\Delta W_{XY}(q;n,f)=W_{XY}(q,E_n)-W_{XY}(q,E_f).
\]

A **direct complex-amplitude implementation** and **independent real/imaginary finite-difference expansion** were compared:

\[
\Delta W_{XY}=\Delta a_X a_{Y,f}+a_{X,f}\Delta a_Y+\Delta a_X\Delta a_Y
+\Delta b_X b_{Y,f}+b_{X,f}\Delta b_Y+\Delta b_X\Delta b_Y .
\]

The **ideal cancellation control** freezes the *entire complex C/O/H factors* at their 12.00-keV values for all energies. Consequently all six ideal non-Br-only `ΔW` are **exactly zero**. The actual XrayDB factors remain energy dependent, so a finite residual is expected. This is a *coefficient-level* negative control, not a simulation of bulk solvent heating.

## 3. Executed numerical findings

### 3.1 Primary XrayDB benchmark: q = 4 Å⁻¹, 13.45−12.00 keV

| Pair | W(12.00 keV), e² | W(13.45 keV), e² | ΔW, e² | ΔW/W(12.00), % |
|:--|--:|--:|--:|--:|
| C–C | 5.625143852 | 5.616605169 | **−0.008538683** | **−0.151795** |
| C–O | 9.193803591 | 9.175512467 | **−0.018291124** | **−0.198951** |
| C–H | 0.526207437 | 0.525805490 | **−0.000401948** | **−0.076386** |
| O–O | 15.026522037 | 14.989522070 | **−0.036999967** | **−0.246231** |
| O–H | 0.860037213 | 0.858975436 | **−0.001061777** | **−0.123457** |
| H–H | 0.049224509 | 0.049224008 | **−5.0083243×10⁻⁷** | **−0.001017** |

**Conclusion A — strict invariance is false:** all six actual non-Br `ΔW` values are finite, although small on an absolute per-pair scale. Hence **perfect non-Br-only cancellation is an ideal approximation**, not an exact consequence of a Br-free solvent under physical tabulations.

**Conclusion B — residual size depends on the pair:** at q = 4, the largest tested non-Br absolute residual is O–O (**0.0370 e²**). Across the physical q range, the largest observed non-Br |ΔW| is O–O (**0.0763049 e²**). The relative O–O weight change spans approximately **−0.119% to −0.385%** over q=0.2–6 Å⁻¹. These figures are neither universal error bounds nor SNR estimates.

### 3.2 Comparison to Br-related coefficients

Step 1.2 native XrayDB reference at q=4 Å⁻¹ (13.45−12.00 keV):

| Reference pair | ΔW, e² |
|:--|--:|
| Br–Br | −155.275846659 |
| Br–C | −9.514493800 |
| Br–O | −15.583673439 |
| Br–H | −0.887066479 |

For O–O vs Br–Br at q=4,

\[
\frac{|\Delta W_{O-O}|}{|\Delta W_{Br-Br}|} = \mathbf{0.0002382854} = \mathbf{0.02383\%}.
\]

Over the **entire 0.2–6 Å⁻¹ q grid**, the greatest per-pair `|ΔW_nonBr|/|ΔW_BrBr|` among the six non-Br types is about **0.03071%** (O–O). **However**, O–O relative to the much smaller Br–H coefficient reaches **4.17% at q=4** and up to **7.40%** across the trial grid.

**Interpretation boundary:** These comparisons are *ratios of per-pair weight coefficients*, not predicted bulk-solvent suppression. A sample contains vastly different numbers of solvent–solvent, solvent–solute and intramolecular pairs; their energy-dependent geometry kernels `ΔK_XY(q,t)` can also be different, have either sign and undergo cancellation. An experimental `D(q,t)` must be computed as `2Σ_pairs ΔW_ij(q)ΔK_ij(q,t)` with the appropriate structural ensemble and calibrated measurements. **No measured or physically simulated solvent heating residual, signal-to-noise ratio or suppression percentage is established here.**

### 3.3 Complex imaginary-factor contribution

The `f″_X f″_Y` contribution was explicitly retained. For the O–O benchmark at q=4 Å⁻¹ its finite energy difference is **−7.33048×10⁻⁵ e²**, approximately **0.198%** of the total O–O `ΔW`. Omitting `f″` would not change the qualitative result here but would violate the specified complex scattering formalism.

### 3.4 Source independence limitation

The Step 1.2 `Henke_Br_swap` and `NIST_Br_swap` **replace only Br** while holding the native XrayDB C/O/H factors identical. Therefore they cannot provide an independent C/O/H tabulation comparison for Step 1.3; agreement of non-Br residual across those channels is **by construction**. Independent all-element source sensitivity, if required, belongs to **Step 1.5**. No scientific conclusion here assumes that three independent C/O/H models were tested.

## 4. Validation and reproducibility

| Check | Observed result | Status |
|:--|:--|:--|
| Native input provenance | Step 1.2 native XrayDB 4.5.8 NetCDF, preserved SHA-256 | PASS |
| Source energy/q coverage and units | 404 E, 292 q including zero, 6 non-Br pairs, 3 contrasts | PASS |
| q→0 neutral f₀ | C 5.997198, O 7.999706, H 0.999978 electrons | PASS within 0.01 e |
| Real/imag decomposition | max absolute discrepancy **1.253×10⁻¹⁴ e²** | PASS |
| Independent finite-difference expansion | max absolute discrepancy **8.105×10⁻¹⁵ e²** | PASS |
| Three-contrast additivity | saved numerical residual **0 e²** | PASS |
| Ideal non-Br `ΔW=0` | exact 0 for every pair/contrast/q | PASS |
| C/O/H invariant under Br-only source swaps | exact 0 source-factor difference | PASS |
| NetCDF save/reopen | `xr.testing.assert_equal` | PASS |
| Actual pytest run | **40 passed**, no skips | PASS |
| Independent saved-artifact script | **54 checks passed**, no failures | PASS |
| Actual native XrayDB *new* API invocation in Step 1.3 | Not performed; **native factor arrays** are an audited Step 1.2 input | Provenance boundary |

**Execution runtime of Step 1.3:** Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0, xarray 2026.4.0, pandas 2.2.3, Matplotlib 3.10.8, pytest 9.0.2. **Original source runtime** recorded in Step 1.2 metrics: XrayDB 4.5.8, Python 3.14.6. Different runtimes are explicitly distinguished; no new xraydb package import is invented.

Input SHA-256, full q/E arrays, numerical metrics and independent audit checklist are preserved in the bundle. Reproduction from the bundle root:

```bash
python src/step13_nonbr_residual.py
python -m pytest -q -rs tests/
python src/audit_step13.py
```

The script writes NetCDF via an in-memory serialization buffer then normal file I/O to avoid the Windows SciPy-backend direct-file write error previously observed with Step 1.2.

## 5. Decision gate and unresolved issues

**Step 1.3 — NUMERICAL PASS.** Realistic neutral-atom XrayDB C/O/H factors give **small but finite** non-Br-only pair-weight residuals. Ideal energy-invariant non-Br factors give exact coefficient-level cancellation. This quantitatively tests the pair-weight part of Phase 0 **P0-009** and confirms that the approximation **ΔW_XY≈0** is *sometimes useful*, but not the stronger identity **ΔW_XY=0** for real tabulations.

**Unresolved / not claimed:** (i) sample pair counts and geometry kernels, (ii) actual photochemical or solvent dynamics, (iii) detector/normalization errors, (iv) experimental background suppression and photon SNR, (v) independent all-element C/O/H tabulation uncertainty. A coefficient ratio cannot resolve these questions. No `g(r)`, MD trajectory, TRXL `ΔI`, `ΔΔI`, inference, or instrument model was generated.

**Scope closure:** No Phase 1 Step 1.4 energy sensitivity optimization has been carried out in this work.

## 6. Files

- `src/step13_nonbr_residual.py` — executed Step 1.3 calculator; all source handling explicit.
- `src/audit_step13.py` — independently executed saved-artifact audit.
- `tests/test_step13_nonbr_residual.py` — 40 executed unit/regression tests.
- `data/reference/step12_native_xraydb_4_5_8.nc` and `step12_calculation_metrics.json` — frozen user-generated upstream source and provenance.
- `results/tables/nonbr_residual.nc` — 404 × 292 full E/q dataset, 6 pairs, 3 contrasts, ideal control and Br coefficient benchmark.
- `results/tables/nonbr_residual_all_q.csv`, `nonbr_benchmarks.csv`, `per_pair_summary_13p45_vs_12.csv`, `calculation_metrics.json`.
- `results/figures/01–05` — PNG and vector PDF.
- `results/logs/` — actual execution, pytest and independent numerical audit outputs.
