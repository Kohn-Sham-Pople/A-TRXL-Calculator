# A-TRXL — Phase 1, Step 1.1: Br Scattering Factor Verification

**Report ID:** ATRXL-P1-S11-20261103  
**Date:** 2026-11-03  
**Scope:** Phase 1, **Step 1.1 only**; Br K-edge original tabulations; independent numerical comparison.  
**Execution status:** **UPDATED 2026-11-03** — original CXRO/NIST calculations executed locally; user provided XrayDB runtime-generated numerical CSV, and arithmetic/grid consistency reverified independently. See §10 addendum. Original §1–9 statements document the *initial* offline run.  
**Assessment (updated):** Step 1.1 **PASS**, with quantified source discrepancies, no independently audited installed XrayDB binary in the reporting environment and unrecorded XrayDB version. Initial assessment appears in §7 for historical traceability.

## 1. Objectives and provenance

The project plan *A-TRXL Python Plan.pdf* specifies Step 1.1: extract CXRO `f1,f2`, compare to Chantler `f′,f″`, record conventions and interpolation rules. The *Phase 0 analytical model*, *pair-weight derivation*, and *hypothesis register* motivate checking that a real Br anomalous contrast exists while avoiding mixed `f1`/`f′` conventions, incorrect eV/keV units, and blind near-edge interpolation. The primary relevant Phase 0 hypothesis is **P0-010** (source convention), and q-input bookkeeping is **P0-011**. Numeric Br/non-Br pair-weight leakage hypotheses P0-008/009 and complete W(q,E) are **not evaluated here** because they belong to later steps.

**Original, independent source tables:**

1. **CXRO/Henke:** LBNL CXRO Br `br.nff` (`E` in eV, **full forward** `f1`, `f2`), https://henke.lbl.gov/optical_constants/sf/br.nff . Underlying paper: B. L. Henke, E. M. Gullikson, J. C. Davis, *Atomic Data and Nuclear Data Tables* **54**, 181–342 (1993), DOI **10.1006/adnd.1993.1013**.
2. **NIST/Chantler:** `usnistgov/FFAST.jl`, original archival CSV `data/mac[35].csv`, https://github.com/usnistgov/FFAST.jl/blob/master/data/mac%5B35%5D.csv , Git blob SHA **5c8932cc8a019b9178055246ed8926dd79733d0d**, `E` in keV, **full forward** `f1`, `f2`. Underlying paper: C. T. Chantler, *J. Phys. Chem. Ref. Data* **24**, 71–643 (1995), DOI **10.1063/1.555974**. Archive is not necessarily identical to the *live* May-2026 NIST FFAST v2.2 backend.
3. **XrayDB interface checked:** https://xraypy.github.io/XrayDB/python.html ; GitHub source implementation https://github.com/xraypy/XrayDB/blob/master/python/xraydb/xraydb.py . Docs explicitly define `f1_chantler` as anomalous **correction** `f′`, and `f2_chantler` as `f″`, with photon-energy input in **eV**; `f0` expects `sin(theta)/lambda`.

Raw source tabulations were accessed via the official web text and NIST organization GitHub source through a repository tool. The locally bundled excerpts faithfully reproduce the selected rows, but **a complete remote full-table file was not saved**, and the optional online full-table checker could not be executed in this environment. See `data/reference/provenance.json`.

## 2. Numerical conventions

At fixed physical momentum transfer q (Å⁻¹) and photon energy E (eV), the Phase 0 convention is

```text
f(q,E) = f0(q) + fprime(E) + i fdoubleprime(E)
q_phys = 4*pi*sin(theta)/lambda [inverse Angstrom]
s_for_xraydb_f0 = q_phys / (4*pi) [inverse Angstrom]
```

CXRO Henke and NIST FFAST original archival tabulations use **full** real forward factor `f1_total≈Z+f′` and imaginary `f2=f″`. Thus `fprime=f1_total−35` for Br. XrayDB `f1_chantler()` already returns `f′`; subtracting Z from its output would be wrong. For clarity, the dataset contains both `f1_total` and `fprime`, never interchangeably.

Interpolation for the **source-level comparison**: piecewise linear `f1(E)` and piecewise log-log `f2(E)`, using original independent tabulated values. No extrapolation. No interpolation across K-edge discontinuities. The actual XrayDB source uses a `UnivariateSpline(s=0)` for `f1` and log-log interpolation for `f2`. Consequently the NIST archival linear-interpolation calculations in this report are **not** presented as outputs of `xraydb.f1_chantler()`.

## 3. Numerical results — actual Python calculation

Benchmarks were calculated at **nominal** 12,000 / 13,240 / 13,450 eV, not confused with the slightly offset raw CXRO tabulated entries 12,020.8 / 13,235.6 / 13,449.7 eV. Units for all table scattering factors are electrons.

| E (keV) | Henke f′ | Chantler f′ | Henke f″ | Chantler f″ | Δf′ (Chantler−Henke) |
|---:|---:|---:|---:|---:|---:|
| 12.00 | -1.898335 | -1.740330 | 0.653927 | 0.601470 | 0.158005 |
| 13.24 | -3.646285 | -3.497684 | 0.547469 | 0.499599 | 0.148601 |
| 13.45 | -6.020298 | -5.983724 | 0.532027 | 0.485080 | 0.036574 |

**Main near–far comparison** (`E_near=13450 eV`, `E_far=12000 eV`):

- Henke `Δf′ = -4.121962885 e`; Chantler original `Δf′ = -4.243393667 e`.
- Henke `Δf″ = -0.121899932 e`; Chantler original `Δf″ = -0.116390917 e`.
- Difference between predicted `Δf′` magnitudes ≈ `2.95%` relative to Henke. **Both predict the same sign and order of magnitude**. This establishes *neutral-atom factor contrast*, not a molecular double-difference intensity or SNR.
- The largest `|f′_NIST−f′_Henke|` on the common dense grid 11.83–13.46 keV occurs at **13460.0 eV**, reaching **2.154020 e**. The maximum `|f″|` inter-source difference on that grid is **0.052976 e**. Agreement at 13.45 keV must **not** be extrapolated to energies arbitrarily closer to the edge.

**Critical near-edge observation:** Henke’s tabulated f′ changes steeply between 13,449.7 and 13,473.6 eV. Thus Henke at the **exact** nominal 13,450.0 eV is `f′≈−6.0203 e`, while the 13,449.7-eV **raw** point is `f′=−5.9475 e`. The original proposal quoted the raw point as an approximation to 13.45 keV; this report distinguishes them. Energies within a few eV of the edge are highly interpolation/edge-position sensitive.

## 4. K-edge behavior, discontinuity guard

The selected source tables locate a sudden f″ increase at different grid intervals:

| Source | Last tabulated pre-edge (eV) | First tabulated post-edge (eV) | f″ pre/post (electrons) |
|---|---:|---:|---|
| CXRO/Henke | 13,473.6 | 13,473.8 | 0.530333 / 3.862040 |
| NIST/Chantler archive | 13,460.23 | 13,476.97 | 0.484390 / 3.802300 |

These are **source-grid bracket locations**, not evidence of a measured chemical shift. NIST's 13,460.23→13,476.97 eV gap spans its edge and is not safe for naive single-interval interpolation. Our common pre-edge grid ends at 13,460.0 eV and all benchmarks lie below both source jumps. Caution: do not equate the Henke and NIST edge locations to an experimentally established chemical-state-dependent edge shift.

## 5. Validation criteria and executed tests

| Criterion | Result | Boundary |
|---|---|---|
| CXRO raw values agree with fixed checkpoint entries | **PASS** | 5 exact source checkpoints tested (including pre/post edge) |
| NIST Chantler archive checkpoints agree | **PASS** | 3 exact NIST source checkpoints |
| `f1_total−Z=f′`, `f2=f″` and eV/keV interpretation | **PASS** | Both tabulations explicitly normalized to common forward-factor convention |
| Source interpolation returns exact knots | **PASS** | Both f1 and f2 at original nodes |
| No pre-edge extrapolation across edge | **PASS** | Explicit exception for out-of-domain energies |
| Physical q → XrayDB `f0` argument conversion | **PASS (convention algebra)** | q/(4π) check; no actual XrayDB f0 call |
| q→0 complex 2-Br scattering identity / real intensity | **PASS (forward algebra)** | `|2f|² = 2|f|²+2Re[f f*]`; absolute closure residual **9.095e-13 e²** |
| Source-independent sign of `Δf′` (12.00→13.45 keV) | **PASS** | Both independent original tabulations predict large negative contrast |
| Actual installed XrayDB numerical values, Waasmaier–Kirfel f0(0) | **NOT RUN** | Package not installed; pytest appropriately skipped |
| Molecule-specific q-dependent pair-weight intensities | **OUT OF SCOPE** | Step 1.2 onward |

**Test runner:** pytest output is **11 passed, 1 skipped**. The skip is the XrayDB-dependent test; it is **not** counted as independently validated. The generated NetCDF and CSV files were written and subsequently inspected as files by Python. No simulation or experimental data are claimed.

## 6. Scientific findings and conclusions

1. **Br K pre-edge contrast exists at the factor level.** Both source data sets predict ~4.1–4.2 electron real anomalous-factor change between 12.0 and 13.45 keV, supporting the fundamental numerical premise of the original idea within an isolated neutral-atom tabulation model.
2. **Absolute anomalous factor values are source-dependent.** At 12.0 keV the `f′` source discrepancy is 0.158 e, and at 13.45 keV it is 0.037 e; the maximum near 13.46 keV is much larger (~2.154 e). Tabulation choice can be important extremely close to the edge.
3. **The Henke f1 versus XrayDB f′ convention is essential.** A mistaken extra Z would overwhelm the real anomalous contrast. q and energy units have also been fixed for future Python code.
4. **`f″` must be kept.** A real-amplitude-only treatment is not mathematically adequate; full complex amplitudes are the proper input for subsequent pair calculations.
5. **Contrast sign robustness is not experimental feasibility.** Nothing here demonstrates solvent cancellation, near-edge chemical-state invariance, structure identifiability, or detectability.

## 7. Unresolved issues and gate

- **XrayDB runtime comparison remains pending.** Environment network DNS was unavailable and `xraydb` import/pip install failed. Code includes `--xraydb` to invoke its real functions and save differences after installation. The existence of XrayDB documentation is **not** substituted for an actual XrayDB runtime test.
- **NIST archive versus live FFAST v2.2 (May 2026)** may differ in backend/format; the exact archived Git blob was recorded to prevent silent version drift. Current live NIST numeric download was not performed.
- **Input file snapshots are excerpts** of the full tabulations. Their rows are independently documented, and an optional online full-table checker is provided, but it was not executed.
- **`f0(q)` full Waasmaier–Kirfel database comparison not run**; q convention and q→0 complex intensity algebra checked, but numerical f0(0) API needs runtime validation.
- **Near-edge chemical-state dependence** and measured energy resolution require later research; a neutral-atom table cannot establish their absence.

**Decision:** Step 1.1 has a **completed source-based numerical demonstration** and passes all **locally executable** checks. Mark **CONDITIONAL (not unconditional full validation)** until XrayDB runtime and optional live-source verification are executed. This does not authorize proceeding to Step 1.2 in this request.

## 8. Deliverables and scope rule

Created: `data/reference/*.csv`, `data/reference/provenance.json`, `results/tables/atomic_factors.nc`, `factor_comparison.csv`, `benchmark_comparison.csv`, `calculation_metrics.json`, execution/test logs, `results/figures/01–04` in PNG/PDF, Python source, tests, and this report.

**Not created by design:** `contrast_weights.nc`, `Wab(q,E)` and `ΔWab(q)` graphs because the research plan assigns those to Step 1.2 onward. The phase-level required outputs are therefore *not* all completed in Step 1.1; that is a deliberate scope boundary, not a missing Step 1.1 calculation.

## 9. References

- B. L. Henke, E. M. Gullikson & J. C. Davis, *Atomic Data Nucl. Data Tables* **54**, 181–342 (1993), https://doi.org/10.1006/adnd.1993.1013 .
- C. T. Chantler, *J. Phys. Chem. Ref. Data* **24**, 71–643 (1995), https://doi.org/10.1063/1.555974 .
- C. T. Chantler, *J. Phys. Chem. Ref. Data* **29**, 597 (2000), https://doi.org/10.1063/1.1321055 (revised soft X-ray edge treatment; important provenance, not direct 13-keV raw data source).
- LBNL CXRO, Br original data: https://henke.lbl.gov/optical_constants/sf/br.nff .
- NIST archival FFAST.jl, Br raw data: https://github.com/usnistgov/FFAST.jl/blob/master/data/mac%5B35%5D.csv .
- XrayDB official documentation: https://xraypy.github.io/XrayDB/python.html ; implementation https://github.com/xraypy/XrayDB/blob/master/python/xraydb/xraydb.py .
- NIST FFAST Version History: https://www.nist.gov/pml/x-ray-form-factor-attenuation-and-scattering-tables/x-ray-form-factor-attenuation-and-scattering .


## 10. User-executed XrayDB runtime validation — addendum, 2026-11-03

**Evidence:** user supplied `xraydb_runtime_comparison.csv` (404 rows, 7 columns), updated `calculation_metrics.json`, and a Spyder screenshot showing `12 passed in 3.43s` / `Return code: 0`. The original `src/step11_br_factors.py --xraydb` invokes real `xraydb.f1_chantler('Br', energies)` and `xraydb.f2_chantler('Br', energies)` if installed. This report does **not** claim a separate XrayDB installation or API execution on the report-generation computer.

**Independent artifact checks carried out on the submitted files:**

- Exactly **404 unique, monotonically increasing energies**, 11,830–13,460 eV; all 7 numeric columns finite.
- Energy grid exactly matches the original `factor_comparison.csv`.
- Against the prior NIST reference dataset: max |difference| for `f′_NIST` <= 5e-12 e, for `f″_NIST` <= 4.99e-13 e (round-off attributable to original CSV precision).
- CSV-reported `XrayDB−NIST` columns recomputed independently: maximum subtraction errors 9.23e-16 e for `f′`, 1.54e-16 e for `f″`.
- User-generated metrics contain `xraydb_runtime_status=RUN`; all benchmark energies occur on the common grid.

### 10.1 Full three-source benchmarks (electron units)

| E (keV) | Henke f′ | NIST archive f′ | XrayDB f′ | Henke f″ | NIST archive f″ | XrayDB f″ |
|---:|---:|---:|---:|---:|---:|---:|
| 12.00 | -1.898335 | -1.740330 | -1.857902 | 0.653927 | 0.601470 | 0.602261 |
| 13.24 | -3.646285 | -3.497684 | -3.546821 | 0.547469 | 0.499599 | 0.505109 |
| 13.45 | -6.020298 | -5.983724 | -5.856022 | 0.532027 | 0.485080 | 0.542467 |

The near-minus-far (`13.45−12.00 keV`) real anomalous factor is **-4.121963 e** (Henke), **-4.243394 e** (NIST archive), **-3.998120 e** (actual XrayDB). XrayDB magnitude is **5.78% lower** than NIST archive and **3.00% lower** than Henke; **all predict negative ~4 e contrast**.

The imaginary contrast differs more: Henke **-0.121900 e**, NIST **-0.116391 e**, XrayDB **-0.059794 e**. At 13.45 keV, XrayDB `f″=0.542467` versus archived NIST `f″=0.485080` (**11.83%** difference relative to archived NIST). The XrayDB `f″` near edge is not simply an exact reproduction of archived FFAST; the difference cannot be attributed *solely* to f′ spline interpolation, since XrayDB uses log-log interpolation for f″. The precise origin (different underlying values/grids/version) is **not yet demonstrated**.

At 13.460 keV, the largest `|f″_XrayDB−f″_NIST|` in the current common domain is **0.101010 e**, about **20.85%** of archived NIST; largest `|f′_XrayDB−f′_NIST|` on the full grid is **0.181862 e** at 13439.625 eV. These source-sensitive differences must be retained as a warning, not silently treated as perfect numerical agreement.

### 10.2 Interpretation and final gate

**PASS for Step 1.1:** actual XrayDB values now exist in the saved dataset, all 12 reported tests passed, source conventions and interpolation boundaries documented, differences explicitly quantified, and the central conclusion (Br real anomalous contrast near −4 e for 12.00/13.45 keV) is robust to the three implementations tested. This is **not** a claim that the three tables match to arbitrary tolerance or that A-TRXL structural sensitivity/SNR is proven.

**Open provenance item:** The exact installed **XrayDB package/SQLite data version** was not recorded in the user's supplied metrics. Thus the runtime input files and output values are auditable, but a fully version-pinned independent rerun has not been demonstrated. The remaining near-edge source/chemical-state sensitivity is a known future Step 1.5 question; it does not invalidate the Step 1.1 conclusion.

**Added artifacts:** `results/tables/xraydb_runtime_comparison.csv`, `runtime_crosscheck_summary.json`, updated `atomic_factors.nc` with XrayDB factors and original reference outputs, `results/figures/05_xraydb_runtime_validation.png`. Initial offline metrics are retained as `calculation_metrics_initial_offline.json`.

**Hashes of unchanged user-supplied inputs:**

- `xraydb_runtime_comparison.csv`: SHA256 `f5c28a8f77cd99a6dc49d5c05c60238cf353711b85fc074dbc77e5eba6a6c106`
- `calculation_metrics.json`: SHA256 `6f45260c4911ca3898c12f91e619b3d4ec87737087317446ab566256ad8aabd8`

**Phase 1 scope:** Step 1.2 energy-dependent pair weights, non-Br leakage, and experimental feasibility are **not** calculated here.
