# A-TRXL Calculator

Use the existing checkout; do not create a worktree unless explicitly requested.
Install with `python -m pip install -r requirements.txt` then `python -m pip install --no-deps -e .`.
Run `python -m pytest -q -ra` and `streamlit run app.py` from this directory.
Do not modify `data/reference/`, existing Phase 1 reports, or SHA256SUMS.txt.
Scientific API uses physical q in Å⁻¹ and photon energy in keV; XrayDB receives q/(4π) and eV.
Retain complex factors, unordered coefficient normalization, and near-minus-far sign.
Native f1 spline extent matters: the Phase 1 extent anchors [11830,13460] eV are
required inside the validated energy interval. Do not replace with unanchored
scalar calls without reporting and testing numerical differences.
No silent extrapolation, source fallback, manufactured benchmark, or claims of experimental intensity/SNR.
Legacy Br-only tables retain pre-edge domains. Full adapters use verified source domains, edge/gap guards and SHA-256 provenance.
Use pathlib and in-memory exports. Do not add NetCDF writes to the GUI.
Regenerate samples with `python scripts/sample_calculation.py` only intentionally;
metadata timestamps and IDs will change. Tests should not rewrite samples.

## v1.1 source and UI contracts

`Henke` and `NIST` mean full source factors for BOTH atoms. Preserve Phase 1 comparisons
under explicit `Henke_Br_swap`/`NIST_Br_swap` IDs; do not reintroduce implicit fallback.
The historical extent anchors apply to historical Br/C/O/H only; new native species
use the native scalar API. New species have source-level tests, not Phase 1 research validation.
Do not hard-code displayed domains: use source table/DB intersection and excluded bands.
W uses one E and has no ΔW inputs/results; ΔW uses near-minus-far with fixed Efar in energy scans.
Energy/q scans are independent, OFF means no calculation, and controls rerun without
numerical work until Calculate. Keep metadata and CSV/JSON/PNG mode-consistent.
Full original tables and pinned source provenance are under data/sources/full with
packaged copies under src/atrxl/data_sources/full. No runtime downloads or source substitutions.
Run the entire suite after source changes; retain all original benchmark numbers/tolerances.
