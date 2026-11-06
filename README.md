# A-TRXL Scattering Factor & Contrast Calculator v1.0

실제로 설치한 XrayDB 4.5.8을 호출하는 neutral-atom scientific calculator입니다.
GUI 없이도 `atrxl` Python API로 사용할 수 있습니다. **계산값은 원자 pair coefficient이며 실험 intensity 또는 SNR가 아닙니다.**

## Windows / Anaconda Prompt 설치

1. 프로젝트 ZIP을 압축 해제하거나 이 저장소를 clone합니다. `app.py`가 보이는 폴더로 이동합니다.
   Unicode 및 OneDrive 경로는 따옴표로 묶으세요.
   ```bat
   cd /d "C:\Users\사용자\OneDrive\연구\A-TRXL-Calculator"
   ```
2. 격리된 환경을 생성하고 활성화합니다. 실제 실행 검증 환경은 Python 3.12.14입니다.
   ```bat
   conda create -n atrxl python=3.12 pip -y
   conda activate atrxl
   python -m pip install -r requirements.txt
   python -m pip install --no-deps -e .
   python -m pip check
   ```
3. 프로그램을 실행합니다. 일반 local desktop 환경에서는 Streamlit이 브라우저를 엽니다.
   ```bat
   streamlit run app.py
   ```
   브라우저가 자동으로 열리지 않으면 터미널에 표시된 Local URL을 브라우저에 입력합니다.
   종료는 Prompt에서 Ctrl+C입니다. 클라우드 검증은 headless로 수행했으며 브라우저 자동 실행은 Windows에서 실행 검증하지 않았습니다.
4. 검증 및 sample 생성:
   ```bat
   python -m pytest -q -ra
   python scripts/sample_calculation.py
   ```

Python 3.14.6 / Windows는 `results/logs/windows_python314_dependency_check.txt`의
공식 win_amd64 CPython 3.14 wheel dependency dry-run으로 호환성만 확인했습니다.
Windows 또는 Python 3.14.6에서 실제 scientific/UI tests를 실행한 것은 아닙니다.
Python 3.12 환경을 우선 권장하며 3.14 사용 시 별도 conda 환경에서 동일한 tests를 실행하세요.
`requirements.txt`는 실제 검증 환경의 전체 고정 버전이며 XrayDB는 4.5.8입니다.

## 사용 방법

- **Home**: 공식과 기본 benchmark 안내.
- **Calculator**: 10개 pair, source, W/ΔW, q, Efar/Enear 및 Default/Extended 선택 후 **Calculate**.
- 단일점 출력에는 각 원소의 f0, f′, f″, complex factor의 실수/허수부, Wfar/Wnear/ΔW와 상대변화가 포함됩니다.
- q scan은 qmin + n·step ≤ qmax인 점을 사용합니다. qmax는 step에 맞는 경우에만 포함됩니다.
- q scan 2개 그래프, energy-dependent W와 PNG/CSV/JSON 다운로드를 제공합니다.
- **Source Comparison**: 원하는 source들을 선택해 W(Enear) 또는 near-minus-far ΔW를 비교합니다.
- **About / Scientific Notes**: 해석 한계, source metadata 및 SHA-256.
- 입력 변경 후 Calculate를 다시 눌러야 결과가 갱신됩니다. 다운로드에는 마지막 성공한 계산의 입력 및 UUID가 연결됩니다.
- PNG는 같은 full-precision arrays를 Matplotlib으로 렌더링하며 Plotly 화면은 interactive입니다.

## Scientific conventions

`f=f0(q)+f′(E)+i f″(E)`, `W=Re[f_a conj(f_b)]`, `ΔW=Wnear−Wfar`.
Physical q는 Å⁻¹, 외부 API/GUI energy는 keV, XrayDB energy는 eV,
f0 argument는 `q/(4π)`입니다. Factors는 electrons, weights는 electrons²입니다.
단일 unordered pair coefficient이며 Debye multiplicity 2, geometry sinc,
분자 구조, sample pair counts 및 실험 intensity를 포함하지 않습니다.

Default: q=0.2–6.0 Å⁻¹, E=11.83–13.46 keV. q=0은 별도 forward diagnostic입니다.
Extended: q≥0 및 실제 native table domain 안에서만 계산하고 경고합니다.
f0 parameterization이 수치적으로 평가된다는 사실은 그 q의 physical accuracy 검증이 아닙니다.
Br K-edge는 native DB에서 13.474 keV로 확인하며 그 경계를 건너는 energy interval은 차단합니다.
입력은 real finite scalar 또는 1/2-D NumPy-broadcastable arrays; 최종 최대 10000 points입니다.

### 중요한 native API 재현 규칙

XrayDB 4.5.8은 `f1_chantler`에 전달한 energy extent에 따라 spline knots를 선택합니다.
동일 버전의 unanchored scalar 호출만으로는 Phase 1 benchmark를 재현하지 못합니다
(ΔW≈−155.6663118992356 e²).
따라서 11.83–13.46 keV 내의 각 요청에는 native API에 `[11830, E_eV, 13460]`을 전달하고
가운데 결과를 선택합니다. 이것은 저장 NetCDF 보간이 아니며 native API 직접 평가입니다.
이를 통해 arbitrary q/E에서도 batch/shape에 무관한 Phase 1-compatible 값을 얻습니다.
범위 밖 Extended E는 scalar native API를 직접 사용합니다. 두 규칙 경계에서의 작은
spline context 차이는 의도적으로 metadata에 명시하며 physical continuity를 보증하지 않습니다.

### Sources and references

Primary: neutral XrayDB/Chantler 4.5.8, database 9.2.
Henke 및 NIST는 **Br-only source swap**이며 모든 f0와 C/O/H는 동일 XrayDB입니다.
Non-Br pair는 swap에 의해 변하지 않습니다. Br raw full f1에서 Z=35를 빼고,
f1는 linear E, 양수 f2는 log-log interpolation을 적용합니다.
검증·등록된 pre-edge domain만 허용합니다: Henke 11.830–13.4736 keV,
NIST 11.830–13.46023 keV. 전체 원시 표가 더 넓더라도 이 adapter는 그 밖의 값을 제공하지 않습니다.
Extrapolation 및 absorption discontinuity를 건너는 보간은 없습니다.

`data/sources/provenance.json`에 URL/version/retrieval UTC/SHA-256/보간법이 있습니다.
NIST는 Phase 1에서 지목한 Git blob의 Git SHA-1까지 검증 후 저장했습니다.
Henke는 공식 URL의 현재 raw file이며 과거 retrieval과 byte-identical이라고 주장하지 않습니다.
두 표의 nominal benchmark와 전체 저장 NetCDF와의 일치를 tests에서 확인합니다.
wheel 설치를 위한 동일 raw files가 `src/atrxl/data_sources/`에도 들어 있습니다.

제공된 Phase 1 files는 `data/reference/`와 `docs/ATRXL_Step1_*`이며 원본 SHA-256을 tests로 검사합니다.
보고서에 언급된 `A-TRXL Python Plan.pdf`, 원래 sparse source excerpts, Phase 0 문서는 이 checkout에 없습니다.
이 파일들을 읽었다고 주장하지 않으며 현재 확보한 raw tables의 provenance를 구분합니다.

## GUI-independent Python API

```python
import numpy as np
from atrxl import atomic_factor, pair_weight, delta_pair_weight
from atrxl import pair_weight_curve, delta_weight_curve
f = atomic_factor('Br', 4., 13.45)  # complex electrons
w = pair_weight(('Br','C'), 4., 12.)
dw = delta_pair_weight('Br-Br', 4., 13.45, 12.)
curve = delta_weight_curve('Br-Br', np.linspace(.2, 6, 291), 13.45, 12.)
```

## Troubleshooting

- `ModuleNotFoundError`: 해당 conda 환경 활성화 후 `python -m pip install -r requirements.txt`와 editable install을 다시 실행하세요.
- Domain error: 메시지의 q/E 범위를 확인합니다. Extended도 source domain 밖 extrapolation을 허용하지 않습니다.
- Source unavailable: raw file과 provenance SHA-256을 검사하세요. 다른 source로 자동 대체하지 않습니다.
- OneDrive: 다운로드는 메모리에서 생성합니다. 파일 저장은 브라우저가 처리하므로 NetCDF backend가 Unicode 경로를 직접 열지 않습니다.
  잠긴 다운로드 파일은 닫은 뒤 다른 파일명으로 저장하세요.
- Relative change denominator |Wfar|≤1e−12이면 CSV 빈값 / JSON null입니다.
- 그래프 글꼴/cache 오류: 쓰기 가능한 local 폴더로 `MPLCONFIGDIR`을 지정할 수 있습니다.
- Port 충돌: `streamlit run app.py --server.port 8502`.

`ATRXL_Calculator_v1_Validation.md`, `results/logs/`, `results/sample_*`, `results/figures/`에
실제 검증 명령, 버전, 결과와 scientific limitations가 기록되어 있습니다.
