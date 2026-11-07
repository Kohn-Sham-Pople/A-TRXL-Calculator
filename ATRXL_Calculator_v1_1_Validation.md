# A-TRXL Calculator v1.1 — Scientific and GUI validation

Date: 2026-11-07 (Asia/Seoul). Base: `atrxl-calculator-v1`, commit
`a0d65fe6e94d122976fdbf1030cb91327e90834a`. Work branch: `atrxl-calculator-v1.1`.

## Outcome and scope

Implemented only the four requested improvements: dynamic source/pair energy domains,
distinct W/ΔW calculation modes, independent single-point/energy/q GUI sections, and
full-element source adapters. Calculation formulas, physical q conversion, complex
conjugation, imaginary term, neutral-atom model and near-minus-far sign are unchanged.
No Phase 2 engine, geometry, MD, dynamics, fitting, experimental intensity or SNR model.

Final executed suite: **611 passed, 0 failed, 0 skipped**, 29.97 seconds.
Seven warnings are pandas `float_precision` deprecations in round-trip CSV tests.
The tests retain the explicit round-trip parser to test float64 preservation; production
CSV export is unaffected. See `results/v1.1/logs/pytest.txt` for the actual runner output.

## Preserved Phase 1 / v1.0 compatibility

All 8 original reference/report SHA-256 checks pass. `data/reference/`, the five Phase 1
reports and `SHA256SUMS.txt` are unchanged. Historical validation reports, execution logs
and sample artifacts remain preserved in their v1.0 paths.

The original **83 scientific test cases** retain their numerical expected values,
mathematical identities, NetCDF coordinates, assertions and tolerances. Their source IDs
were explicitly migrated from `Henke`/`NIST` to `Henke_Br_swap`/`NIST_Br_swap`, because
full Henke/NIST now use source-specific factors for **both** atoms. The invalid-element
example `Fe-Br` was replaced by `Xx-Br`, since Fe is now intentionally supported.
The original 4 GUI test cases were adapted to keyed A/B controls, new mode-specific plots
and OFF-by-default independent scans, rather than preserving the old GUI defect.
All **87 compatibility cases passed**, not skipped. The separate compatibility run is
recorded in `results/v1.1/logs/legacy_87.txt`.

The legacy channels are explicit opt-in IDs, restricted to Phase 1 Br/C/O/H. They retain
native C/O/H anomalous factors by definition and are labeled Br-only in metadata and
Source Comparison. Full `Henke`/`NIST` do **not** perform such substitutions.

| Br–Br, physical q=4 Å⁻¹, Enear−Efar=13.45−12 keV | Actual ΔW (e²) |
|---|---:|
| XrayDB 4.5.8 / Chantler | −155.27584665917908 |
| Full Henke / legacy Henke Br-only | −159.31575985496892 |
| Full NIST archived / legacy NIST Br-only | −164.81242968413034 |

Primary reference remains −155.27584665917914 e². The independent metrics check retains
10⁻¹¹ e² absolute tolerance. Historical archived six-decimal benchmarks retain 5×10⁻⁷ e².
The full reference audit again compared **1,415,616 W values** with no failed comparisons;
`results/v1.1/logs/reference_audit.json` retains all per-source/pair maximum residuals.
No archived table or reference value was altered to make tests pass.

The intentional full-source change is visible in Br–C at the same benchmark:

| Source | ΔW_BrC (e²) |
|---|---:|
| XrayDB | −9.5144938004816808 |
| Full Henke | −9.8199258066186914 |
| Legacy Henke Br-only | −9.8081426180117859 |
| Full NIST archived | −10.11109736762576 |
| Legacy NIST Br-only | −10.096150302654415 |

These are different source definitions, not an unexplained regression or silent fallback.

## Actual source acquisition, support and provenance

| Source | Discovered neutral species | Original inputs |
|---|---:|---|
| XrayDB/Chantler | 92, H–U | Installed 4.5.8 SQLite DB 9.2; neutral Waasmaier/Chantler row intersection |
| CXRO/Henke | 92, H–U | Official https://henke.lbl.gov/optical_constants/sf.tar.gz |
| NIST/Chantler archived | 92, H–U | Official https://github.com/usnistgov/FFAST.jl pinned snapshot |

Henke archive SHA-256: `e8fbbea874f23a9cef8e8feb8eea0ea7ace68874d8f8738377b1028e312189c6`.
NIST commit: `a3ccd0b60721a5557d6c77c167e64d89076b77ca`.
Its Br file still has the original Phase 1 Git blob
`5c8932cc8a019b9178055246ed8926dd79733d0d`.
NIST here means this archive, not an assertion of identity with the current live FFAST backend.

The 184 original tables, NIST shell data and license are stored locally under
`data/sources/full/`; identical package copies are in `src/atrxl/data_sources/full/`.
`manifest.json` records actual retrieval UTC, per-element URL/version/Z/SHA-256,
raw and valid knot counts, min/max energy, forbidden interpolation intervals and reasons.
`results/v1.1/source_inventory.csv` records all 276 source/element capabilities.
No network request occurs during factor evaluation or GUI calculation.

All sources use native XrayDB/Waasmaier–Kirfel **neutral f0(q/(4π))** with provenance.
Henke/NIST use **both atoms' own source tables** for anomalous factors. Independent tests
replace native f1/f2 functions with failures while evaluating full-source Br–C, proving
no native anomalous fallback is needed. Native f0 is intentionally common.

## Energy domains and edge-aware interpolation

Domain values come from installed native `chantler_energies()` or verified raw table
coverage, never from a manually chosen universal photon-energy limit. The pair domain
is the intersection of A and B. Envelope extrema differ by element/source. For example,
Fe/Au in these snapshots has native 0.00101–966.26674 keV, Henke 0.0293–30 keV and
NIST 0.01069–432.9451 keV. Domain displays are computed from those inputs on each
selection change; these examples are observations, not application constants.

An envelope is not a continuous interpolation permission. The GUI also supplies usable
bands, per-element coverage and edges/gaps. Native scan bands split at installed DB shell
edges. A scan from below to exactly a native edge is blocked (lower band's upper endpoint
is open for scans); a standalone native API point there remains a source model evaluation.
Raw gap endpoint knots are readable, but interpolation or contrasts traversing the gap
are rejected. No curve interpolation across shell discontinuities is introduced.

Full-source adapters apply linear **full f1(E)** interpolation and subtract the selected
atom's **Z**, and log-log interpolation of positive f2. Zero–zero intervals return zero.
Zero/positive transitions are excluded from interpolation. Native XrayDB returns f′
directly and is not subjected to an extra Z subtraction.

Actual source quality issues and conservative handling:

- Henke undefined f1 sentinel −9999 is excluded; the resulting undefined regions are not extrapolated.
- Henke Mg, Pt and Er contain duplicate energy rows with conflicting factors; all rows at
  those ambiguous energies are excluded and interpolation across their gaps is forbidden.
  Er also has an out-of-order row; a stable energy sort is applied to the evaluation view.
  Original bytes remain unchanged.
- Narrow adjacent Henke knot brackets are detected as discontinuities. Installed native
  shell locations provide a conservative supplementary locator: if no matching raw narrow
  bracket is present, its enclosing raw interval is blocked and labeled with that provenance.
  This policy can exclude extra intervals; it is not a measured edge-energy determination.
- NIST shell edge locations come from the **same archive's shelldata.csv**. Their enclosing
  intervals are blocked; exact valid raw knots can still be queried. The shell metadata,
  shell-data hash and both original/package table hashes are checked in tests.
- A native W curve or energy contrast spanning any selected atom's shell edge is rejected.
  Full-source curves/contrasts also reject table-gap crossings. No exploratory bypass exists.

Historical native Br/C/O/H f1 evaluation within 11.83–13.46 keV retains Phase 1 endpoint
anchors `[11830, E_eV, 13460]`. This is a direct native API protocol, not NetCDF interpolation.
It is necessary because native spline knot selection depends on request extent. New
species and energies outside that interval use scalar native API evaluation. Native
near-edge spline predictions do not establish chemical-state accuracy or experimental validity.

Default numerical q/E input limits remain as v1.0. Phase 1 physical/research validation
covers **Br/C/O/H only**. New species received source-level availability, convention,
interpolation and numerical identity tests. They are not claimed Phase 1 validated.

## Mode, scans and exports actually verified

- **W**: one E; no far/near widgets or ΔW computation. Point W is the main metric; optional
  W(q) and W(E). W CSV/JSON contains W and one-energy factors, not delta/far/near fields.
- **ΔW**: Efar/Enear; metrics Wfar/Wnear/ΔW. Optional ΔW(q), and energy scan
  ΔW(Enear; fixed Efar). Near-minus-far sign is tested against direct formulas.
- Three labeled GUI sections in required order, divided with spacing and rules.
- Both scan checkboxes default OFF and react without Calculate. Disabled inputs are ignored.
  Tests replace scan functions with failures to prove disabled scans do not execute, and
  replace ΔW evaluation with a failure to prove W mode does not calculate ΔW.
- All 8 W/ΔW × Energy Scan ON/OFF × q Scan ON/OFF combinations are exercised by AppTest.
- Data-source and A/B selection changes immediately update the Extended domain before
  any numerical result exists. Supported species selectors and missing-factor errors are tested.
- Invalid enabled scan inputs clear result/export state. Changing controls hides stale
  results until Calculate, rather than exporting old results under new labels.
- Actual mode-specific plotted trace names/values, CSV float64 round-trip, JSON precision,
  UUID linkage, full source provenance and PNG output are checked. No NetCDF writes.

## Executed checks / logs

| Operation | Outcome |
|---|---|
| `python -m pytest -q -ra` | 611 passed, 0 failed, 0 skipped; 29.97 s |
| Legacy scientific + adapted original GUI cases | 87 passed |
| Source discovery | 92 verified species per source |
| Per-element factor and arbitrary H-pair identities | Every one of 92 species × 3 sources executed |
| Raw-knot / off-knot interpolation / all registered forbidden gaps | Every archived element table executed |
| `python scripts/audit_reference.py` | 1,415,616 original W values compared |
| `sha256sum -c SHA256SUMS.txt` | 8/8 original files unchanged |
| `pip install --no-deps -e .`, `pip check` | Successful, no broken requirements |
| `pip wheel . --no-deps`, separate target install | All 184 tables packaged; 92 species/source discovered outside checkout |
| Separate wheel Fe–Au W and primary Br–Br benchmark | Passed; mode-specific W JSON verified |
| Streamlit on 8502 + `python scripts/browser_smoke.py` | Actual Chromium calculation, mode change, scan toggle and JSON download |

Logs and browser screenshots are under `results/v1.1/`. The full source inventory, new
Fe–Au W sample/figure, Br–C full/legacy comparison, and ΔW energy sample are also included.

Early debugging is not counted as PASS: an initial full-source scan identity test had
zero absolute tolerance at cancellation and differed by 5.68×10⁻¹⁴ e²; it now retains the
identity assertion with 2×10⁻¹² e² absolute floating-point tolerance. Production formulas
were not changed for that test. Initial browser probes needed visible-label/radio-group
selectors instead of hidden radio inputs or ambiguous metadata text. A long-running
Streamlit process was restarted after package code changes before final browser testing.

## Versions and limits

Runtime: Python 3.12.14 / Linux; calculator 1.1.0; XrayDB 4.5.8 / DB 9.2;
NumPy 2.5.4; SciPy 1.18.1; pandas 3.1.1; Streamlit 1.67.0; Plotly 7.1.0;
Matplotlib 3.11.2; pytest 9.2.0; xarray 2026.9.0.
Optional browser tooling: Playwright 1.64.1 and installed Chromium 156.0.8078.79.
Full actual version evidence: `results/v1.1/logs/versions.json`.

Windows and Python 3.14.6 runtime behavior have not been tested; Windows instructions
are supplied and the existing pinned dependency set is unchanged. There are no skipped
scientific tests hiding unavailable source tables. The plan PDF and original historical
sparse excerpts remain unavailable and have not been read. Source-level numerical
validation is not experimental validation, chemical-state modeling, uncertainty calibration,
or proof of physical accuracy at every mathematically evaluable q/E.
