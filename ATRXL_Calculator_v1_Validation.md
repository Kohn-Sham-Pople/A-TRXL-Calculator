# A-TRXL Calculator v1.0 — Scientific validation

검증 작성일: 2026-11-07 (Asia/Seoul). 로그의 UTC 실행 시각은 2026-11-06입니다.

## Result

실행 가능한 Streamlit local application, GUI-independent engine, 세 source adapter,
full-precision exports, plotting, tests 및 sample artifacts를 구현했습니다.
`python -m pytest -q -ra`: **87 passed, 0 failed, 0 skipped**, 14.01 s.
Pandas `float_precision` test argument에 대한 deprecation warning 1개가 있습니다.
scientific assertion을 삭제하거나 benchmark를 production 계산에 하드코딩하지 않았습니다.

## Inputs actually inspected

- 원본 `docs/ATRXL_Step1_1_Report.md`부터 `ATRXL_Step1_5_Final_Report.md`까지의 5개 보고서.
- `data/reference/contrast_weights.nc`: 404 energies × 292 q, 3 sources, 4 Br pairs, native factors.
- `data/reference/calculation_metrics.json`, `xraydb_runtime_comparison.csv`.
- `SHA256SUMS.txt`: 8개 원본 reference/report의 SHA-256 모두 유지, 재검증 PASS.
- 미제공: `A-TRXL Python Plan.pdf`, Phase 0 문서, 과거 Henke/NIST sparse source excerpts.
  이 파일들을 읽거나 byte 비교했다고 주장하지 않습니다.

새 alternate raw tables는 보고서에 기록된 official Henke URL 및 NIST Git blob에서
TLS 검증을 유지하며 확보했습니다. NIST Git blob SHA-1
`5c8932cc8a019b9178055246ed8926dd79733d0d`를 Git object encoding으로 독립 검증했습니다.
`data/sources/provenance.json`에는 URL, version, 실제 acquisition UTC, SHA-256,
허용 pre-edge domain 및 interpolation convention이 기록되어 있습니다.
Henke의 과거 파일과 byte-identical임은 확인할 수 없지만 nominal benchmarks 및 전체
저장 NetCDF의 source swap 결과와 수치적으로 일치합니다.

## Native evaluation protocol: important finding

XrayDB 4.5.8의 `_from_chantler`는 요청한 energy minimum/maximum으로 주변 spline knots를
선택합니다. `f1_chantler('Br', 13450)`의 scalar 값과
`f1_chantler('Br', [11830, 13450, 13460])[1]`은 다릅니다.

- scalar Br f′ at 13450 eV: −5.867239261730717 e.
- Phase 1 extent Br f′: −5.856022224264764 e (original NetCDF와 동일).
- unanchored scalar protocol ΔW: −155.6663118992356 e².

재현 가능한 계산 규칙: Phase 1 energy interval 내부에서는 각 E에 대해 native API에
`[11830, E_eV, 13460]`을 전달하여 중앙값을 선택합니다. 따라서 요청 배열 크기/순서에
의존하지 않고 원본 Phase 1의 spline extent를 유지합니다. NetCDF curve를 외삽하거나
저장 값을 production에 대입하지 않습니다. f0와 f2도 실제 native API를 사용합니다.
Phase 1 interval 밖 Extended E는 그 E의 scalar native API를 사용하며 source metadata에
이 차이를 명시합니다. 경계에서 spline-context에 따른 작은 차이가 가능하고, extended
결과의 physical accuracy 또는 spectroscopy accuracy는 검증하지 않았습니다.

## Benchmarks and independent agreement

| Br–Br, q=4 Å⁻¹, Enear−Efar=13.45−12.00 keV | Computed ΔW (e²) | Reference |
|---|---:|---:|
| XrayDB 4.5.8 / Chantler | −155.27584665917908 | −155.27584665917914 |
| Henke Br-only swap | −159.31575985496892 | ≈−159.315760 |
| NIST archived Br-only swap | −164.81242968413034 | ≈−164.812430 |

Primary benchmark agreement tolerance: 5×10⁻⁷ e² in the parameterized published benchmark
check (shared with the six-decimal archived benchmarks); primary exact metrics check separately
requires 10⁻¹¹ e². No newer XrayDB engine was substituted.

`python scripts/audit_reference.py` compares **1,415,616 W values**, all 404 energies,
all 292 q points, all 3 sources and 4 Br-containing pairs.
Maximum primary residual: 9.095×10⁻¹³ e².
Maximum Henke residual: 3.290×10⁻¹⁰ e².
Maximum NIST residual: 3.317×10⁻¹⁰ e².
Tolerance: absolute 10⁻⁹ e² + relative 2×10⁻¹³.

Initial stricter 3×10⁻¹¹ e² archived tests failed. Investigation found saved alternate
f′ factors differ from direct raw interpolation by at most 4.997×10⁻¹² e and f″ by
4.998×10⁻¹³ e, consistent with upstream serialization precision. A conservative
10⁻⁹ e² absolute bound covers the independently measured factor-rounding propagation.
Primary residuals remain near floating-point precision. No reference files were changed.
Full residuals: `results/logs/reference_audit.json`.

## Executed tests and operations

| Command / capability | Actual outcome |
|---|---|
| `python -m pip install -e .` | Editable package built and installed |
| `python -m pip check` | No broken requirements |
| `python -m pytest -q -ra` | 87 passed, 0 failed, 0 skipped |
| `python scripts/audit_reference.py` | 1,415,616 coefficient comparisons PASS |
| `sha256sum -c SHA256SUMS.txt` | 8/8 originals unchanged |
| `python scripts/sample_calculation.py` | CSV/JSON/PNG samples generated |
| `streamlit run app.py --server.headless true --server.address 127.0.0.1 --server.port 8501 --browser.gatherUsageStats false` | Server actually started |
| Local GET `/`, `/_stcore/health` | HTTP 200, Streamlit document / `ok` |
| `pip wheel . --no-deps` then isolated target install | Wheel import and all 3 source benchmarks PASS outside checkout |
| CPython 3.14 / win_amd64 `pip install --dry-run --ignore-installed --only-binary=:all: -r requirements.txt` | All pinned wheel dependencies resolve |

Tests include physical q/(4π) conversion, neutral forward f0 limit (0.01 e tolerance),
complex-conjugation and independent real/imag algebra, imaginary contribution,
self-pair |f|², pair symmetry, exact zero-energy contrast, near-minus-far identity,
non-Br swap invariance across all six non-Br pairs, arbitrary off-grid native q/E,
batch invariance, Default/Extended bounds, Br K-edge crossing, pre-edge source range,
exact raw knots and independently reconstructed interpolation, unavailable/hash-invalid
source rejection, invalid dtype/shape/q/energy, negative/zero scan step and oversized
broadcast/scan guards, full-precision CSV/JSON round-trip, metadata linkage, PNG decoding,
packaged raw-table identity, all original source hashes, and full NetCDF W agreement.

Four `streamlit.testing.v1.AppTest` tests actually execute GUI scientific calculations,
validate the displayed primary benchmark, W controls, single-point mode, download controls,
source comparison, navigation, invalid input, Extended warning and NIST range rejection.
An initial UI test incorrectly reran the form between field edit and submission; fixing
its interaction sequence retained the actual invalid-input assertion, which passes.

`results/logs/pytest.txt`, `http_smoke.json`, `install.txt`, `wheel_build.txt`,
`wheel_smoke.txt`, `pip_check.txt`, `original_reference_hash_check.txt`,
`windows_python314_dependency_check.txt`, `versions.json` preserve actual evidence.

## Runtime versions

| Component | Actual version |
|---|---|
| Python / platform | 3.12.14 / Linux |
| XrayDB / SQLite data | 4.5.8 / 9.2 |
| NumPy / SciPy | 2.5.4 / 1.18.1 |
| pandas / xarray | 3.1.1 / 2026.9.0 |
| Streamlit / Plotly | 1.67.0 / 7.1.0 |
| Matplotlib / pytest | 3.11.2 / 9.2.0 |

Full pinned environment: `requirements.txt`. Database file SHA-256 is emitted in each
primary calculation's source metadata. Native package/database versions remain present
for Br-only swaps as well.

## Delivered behavior and limitations

- GUI provides 10 pairs, W/ΔW output, single-point values, q scans, Wnear/Wfar and ΔW plots,
  source comparison and requested-domain energy curves; CSV/JSON/metadata and PNG exports.
- q is physical Å⁻¹, energy external keV/internal eV. float64/complex128 calculations are
  independent of UI rounding; CSV uses 17 significant digits, JSON retains Python float precision.
- All source swaps are Br-only; no all-element Henke or NIST independence is asserted.
- Br K-edge native DB value is 13.474 keV. Cross-edge intervals are blocked; there is no
  exploratory override. Archived adapters intentionally support only registered pre-edge
  ranges, even though retained full raw tables extend further.
- Extended mode is numerical exploration. Neutral-atom tables do not establish chemical-state
  line shapes, experimental sensitivity, solution heating cancellation, or SNR.
- No Phase 2 engine, geometry, MD, dynamics, fitting, experimental intensity or SNR was implemented.
- Windows/Python 3.14.6 execution and automatic browser launch were **not** tested on Windows.
  Dependency resolution alone does not prove runtime behavior. Python 3.12 is the tested workflow.
- There are no skipped tests hiding unavailable reference data. Missing PDF and historical sparse
  excerpts are explicitly listed above; raw-source lineage and numerical agreement are documented
  without claiming historical byte identity.
- NetCDF is read-only reference input. Application writes no NetCDF; exports are generated in
  memory and paths use pathlib to avoid the earlier Windows direct NetCDF-file writing failure.
- No environment publication or remote Git push was performed.

## Archive / Unicode path check

전체 ZIP의 CRC 검사를 통과한 후 `/tmp/atrxl-archive-check/연구 OneDrive/` 아래에
압축을 풀고 새 checkout 경로에서 동일한 `python -m pytest -q -ra`를 실행했습니다.
**87 passed, 0 failed, 0 skipped**, 13.25 s. 실제 Linux Unicode/space 경로 검증이며
Windows/OneDrive sync service 실행을 검증한 것은 아닙니다.
로그: `results/logs/archive_unicode_smoke.txt`.
