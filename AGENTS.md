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
Raw alternate source tables have verified pre-edge domains and SHA-256 provenance.
Use pathlib and in-memory exports. Do not add NetCDF writes to the GUI.
Regenerate samples with `python scripts/sample_calculation.py` only intentionally;
metadata timestamps and IDs will change. Tests should not rewrite samples.
