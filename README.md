# A-TRXL Scattering Factor & Contrast Calculator v1.1

A local Streamlit scientific calculator with a standalone NumPy/XrayDB engine.
Scope: atomic pair coefficients only. No Debye engine, molecular geometry, experimental
intensity, pump–probe double difference, or SNR model is included.

## Windows / Anaconda Prompt: 설치 및 실행

1. GitHub의 **atrxl-calculator-v1.1** 브랜치를 선택하고 Code → Download ZIP으로 내려받아
   압축을 풀거나 해당 브랜치를 clone합니다. 기존 v1.0은 `atrxl-calculator-v1`에 보존되어 있습니다.
2. Anaconda Prompt에서 `app.py`가 보이는 폴더로 이동합니다. Unicode/OneDrive 경로는 따옴표를 사용합니다.
   ```bat
   cd /d "C:\Users\사용자\OneDrive\연구\A-TRXL-Calculator"
   conda create -n atrxl-v11 python=3.12 pip -y
   conda activate atrxl-v11
   python -m pip install -r requirements.txt
   python -m pip install --no-deps -e .
   python -m pip check
   streamlit run app.py
   ```
3. 브라우저가 자동으로 열리지 않으면 Streamlit이 터미널에 표시한 Local URL을 열어주세요.
   종료는 Ctrl+C입니다. 포트 충돌 시 `streamlit run app.py --server.port 8503`을 사용합니다.
4. Tests:
   ```bat
   python -m pytest -q -ra
   python scripts/audit_reference.py
   ```

실제 실행 검증은 Linux / Python 3.12.14에서 수행했습니다. Windows/Python 3.14.6의 실제
실행 검증을 주장하지 않습니다. 기존 고정 requirements의 Windows CPython 3.14 wheel 해석
결과는 v1.0 로그에 보존되어 있지만 v1.1 Windows 실행과는 구분해야 합니다.
XrayDB 4.5.8을 유지했고 runtime dependency는 v1.0과 같습니다.
Source 원본 데이터가 프로젝트와 설치 wheel에 포함되므로 계산 중 인터넷 다운로드는 없습니다.

## Four v1.1 improvements

### 1. Dynamic Extended energy domains

Extended를 선택하면 현재 source와 Element A/B에 대한 실제 energy-domain intersection을
keV로 표시합니다. Source 또는 원소를 바꾸면 즉시 갱신됩니다. 계산은 Calculate를 눌러야 실행됩니다.

최소/최대 값은 **전체 envelope**입니다. 그 사이에 absorption edges 또는 interpolation gaps가
있으므로 펼침 영역의 `usable_bands_keV`, 원소별 domain, edge 및 금지 구간도 확인하세요.
Native source의 scan bands는 설치 DB의 shell edges에서 분리합니다.
확인한 공통 domain 밖의 에너지, 금지된 구간 보간 또는 edge를 건너는 scan/contrast는
오류로 차단하며 extrapolation과 source fallback을 하지 않습니다.

Native DB는 원소별 `chantler_energies()`의 실제 양끝을 사용하고,
Henke/NIST는 locally validated raw tables의 양끝과 manifest의 금지 구간을 사용합니다.
예를 들어 Fe/Au는 이번 데이터에서 native 0.00101–966.26674 keV,
Henke 0.0293–30 keV, NIST archived 0.01069–432.9451 keV입니다.
이 값은 구현에서 hard-code한 한계가 아니며 pair/source마다 다시 계산됩니다.

### 2. W / ΔW are distinct calculations

| Mode | Inputs | Main result | q scan | Energy scan |
|---|---|---|---|---|
| W | q, **E** | W(q,E) | W(q) at fixed E | W(E) at fixed q |
| ΔW | q, **Efar/Enear** | Wfar, Wnear, ΔW | ΔW(q) | ΔW(Enear; fixed Efar) |

W에서는 Efar/Enear 입력과 ΔW 계산·결과가 없습니다. W CSV/JSON에는 `W_e2`, `energy_keV`와
단일 에너지의 factors만 들어 있습니다. ΔW에는 `W_far_e2`, `W_near_e2`, `delta_W_e2`,
relative change 및 두 에너지의 factors가 포함됩니다. Relative denominator |Wfar|≤1e−12이면
CSV 빈값 / JSON null입니다. PNG는 해당 mode의 동일 numerical arrays로 생성합니다.

### 3. Three independent sections

1. **Single-point Calculation**: Element A/B, source, domain mode, q, W/ΔW와 해당 에너지 입력.
2. **Energy Scan**: checkbox와 energy min/max. 꺼져 있으면 입력이 비활성화되고 scan을 계산하지 않습니다.
3. **q Scan**: checkbox와 qmin/qmax/step. Energy Scan과 독립적으로 동작합니다.

두 scan은 기본 OFF입니다. Checkbox와 mode는 즉시 UI에 반영됩니다. 입력 변경 후 기존 결과는
숨겨지고 Calculate 재실행 안내가 나옵니다. Calculate는 모든 활성화 입력을 검증한 뒤 결과와
다운로드를 함께 표시하므로 오류가 나면 일부 성공 결과를 새 결과로 export하지 않습니다.
q grid는 qmin+n·step≤qmax로 구성하며 qmax는 step에 맞는 경우에만 포함됩니다. 최대 10000 points.

### 4. Full-element, source-specific anomalous factors

이번 설치 데이터에서 세 source 모두 **H–U의 92 neutral species**를 지원합니다.
지원 목록은 실제 neutral f0 및 factor rows/raw files를 확인하여 생성합니다.
A/B는 독립 선택하며 same-element/different-element pair가 모두 가능합니다.
손상·누락된 factor는 지원 목록에서 제외하고 사유를 표시합니다. 다른 source로 대체하지 않습니다.

- **XrayDB/Chantler**: neutral Waasmaier f0 및 Chantler f1/f2 row의 실제 intersection.
- **CXRO/Henke**: 공식 `sf.tar.gz`의 92개 original `.nff` tables.
- **NIST/Chantler archived**: 공식 usnistgov/FFAST.jl의 pinned snapshot, 92개 `mac[Z].csv`와 `shelldata.csv`.
  Live NIST FFAST backend의 최신 table이라고 주장하지 않습니다.
- 모든 source의 f0는 공통 XrayDB/Waasmaier–Kirfel neutral model입니다.
- Henke/NIST에서는 **선택한 두 원소 모두** 해당 source의 anomalous factors를 사용합니다.
  full f1에서 각 원소의 Z를 빼서 f′로 변환하며, f2는 f″입니다.

Provenance: `data/sources/full/manifest.json`에 table별 URL, version/commit, retrieval UTC,
SHA-256, 실제 domain, excluded intervals, interpolation 및 common-f0 provenance가 있습니다.
원본 표는 수정하지 않았고 동일 package copies는 `src/atrxl/data_sources/full/`에 포함됩니다.

**Edge / interpolation policy**:

- f1: piecewise linear in E, then subtract Z. Positive f2: log-log; zero–zero intervals: exact zero.
- Henke −9999 undefined f1은 사용할 수 없습니다. Mg/Pt/Er 등의 상충 duplicate-energy rows는
  모두 제외하고 그 빈 구간 보간을 막습니다. 원본 raw file은 보존합니다. Er의 역순 행은 안정 정렬합니다.
- Henke narrow raw-knot edge brackets를 검출하고, native shell locator 주변에 bracket이 없으면
  해당 enclosing interval을 보수적으로 금지합니다. 이 보조 locator의 출처도 metadata에 명시합니다.
- NIST는 **동일 archive의 shelldata**에서 각 shell edge를 확인하여 enclosing interval 보간을 금지합니다.
- Zero/positive f2 경계는 log-log가 정의되지 않으므로 금지합니다.
- 보수적으로 금지한 넓은 bracket 안의 **정확한 유효 원본 knot**는 직접 읽을 수 있습니다.
  그 사이를 건너는 scan/ΔW는 허용하지 않습니다. 상충 duplicate knot 자체는 사용할 수 없습니다.
- Native source는 직접 XrayDB API를 호출합니다. Near-edge native spline 및 raw tabulation의
  numerical prediction은 chemical-state 또는 experimental validity 검증을 의미하지 않습니다.

## Scientific conventions / Phase 1 compatibility

`f=f0(q)+f′(E)+i f″(E)`, `W=Re[f_a conj(f_b)]`, `ΔW=Wnear−Wfar`.
외부 q는 physical Å⁻¹, energy는 keV, 내부 XrayDB는 q/(4π)와 eV입니다.
Factors는 electrons, W/ΔW는 electrons² per one unordered pair입니다.
Complex conjugation, imaginary term, neutral model을 유지하며 multiplicity 2와 geometry는 포함하지 않습니다.

Default numerical input range: q=0.2–6 Å⁻¹, E=11.83–13.46 keV, q=0은 forward diagnostic.
**Phase 1 research validation은 Br/C/O/H에 한정**됩니다. 새 원소는 source-level numerical tests이며
Phase 1에서 모두 검증되었다고 주장하지 않습니다. Extended f0 evaluation은 physical accuracy 보장이 아닙니다.

XrayDB 4.5.8 f1 spline은 요청 energy extent에 의존합니다. 기존 Br/C/O/H의 11.83–13.46 keV
계산은 `[11830, E_eV, 13460]` native extent anchors를 유지하여 Phase 1을 재현합니다.
새 원소와 그 밖의 E는 scalar native API를 직접 호출합니다. 저장 NetCDF의 외삽은 없습니다.

**Source-ID migration**: v1.1에서 `Henke`/`NIST`는 두 원소 모두 해당 source를 뜻합니다.
역사적 Br-only 값은 `Henke_Br_swap`/`NIST_Br_swap`으로 명시적으로 보존했습니다.
Source Comparison에서 full 및 labeled Phase 1 legacy channels를 함께 선택할 수 있습니다.
기존 Br–Br benchmark는 full/legacy 모두 동일하며, Br–C 등의 full-source 결과는 non-Br도 바뀌므로
legacy Br-only 값과 다를 수 있습니다. 이는 fallback이 아니라 정의를 명확히 분리한 것입니다.

## Python API / exports

```python
import numpy as np
from atrxl import atomic_factor, pair_weight, delta_pair_weight
from atrxl.sources import supported_elements, pair_domain
from atrxl.export import calculate_w, calculate, energy_scan

print(supported_elements('Henke'))
print(pair_domain(('Fe', 'Au'), 'Henke'))
f = atomic_factor('Fe', 4., 12., 'Henke')
w = pair_weight(('Fe', 'Au'), 4., 12., 'Henke')
dw = delta_pair_weight('Br-Br', 4., 13.45, 12., 'XrayDB')
point, metadata = calculate_w('Fe-Au', 4., 12., 'NIST')
curve, metadata = energy_scan('Br-Br', 4., np.linspace(12, 13.45, 50),
                              calculation_mode='ΔW', e_far=12.)
legacy = delta_pair_weight('Br-C', 4., 13.45, 12., 'Henke_Br_swap')
```

Scalar/NumPy arrays use float64/complex128 and explicit shape/finite checks.
CSV uses 17 significant digits; JSON retains Python float precision. UUID links numerical results to
metadata including species, q/E, mode, source tables/hash/version, native engine version, warnings and UTC.
No NetCDF is written: original NetCDF is read-only regression reference. Downloads are generated in memory
and pathlib handles Unicode paths; browser handles final file saving under Windows/OneDrive.

## Validation and troubleshooting

See `ATRXL_Calculator_v1_1_Validation.md` and `results/v1.1/logs/` for executed checks and limitations.
Historical v1.0 validation/report/logs/sample outputs remain in their original paths.
`data/reference/`, Phase 1 reports and `SHA256SUMS.txt` remain unchanged.
The plan PDF and historical sparse excerpts were not supplied and have not been read.

- Module missing: activate the intended conda environment and repeat installation commands.
- Source unavailable: inspect the reported SHA-256/missing-factor error and restore the recorded table;
  changing sources is an explicit user action.
- Domain/edge error: use Extended display's actual bands and split scans into one band.
- Optional browser validation: `pip install -e ".[browser]"`; install Chromium and set
  `ATRXL_CHROMIUM` to its executable path, then `python scripts/browser_smoke.py` against a local server on 8502.
- Developer-only dataset refresh: `scripts/prepare_full_sources.py HENKE_TAR FFAST_CHECKOUT` prepares
  already-downloaded authoritative inputs. Pin/verify the recorded source snapshots and rerun all tests;
  this is never run by the GUI. See manifest and validation report for exact acquisition URLs and versions.
- 소스나 패키지 업데이트 후에는 실행 중인 Streamlit을 Ctrl+C로 종료하고 다시 시작하세요.
  장시간 실행 중인 Python 프로세스가 이전 import module을 유지할 수 있습니다.
