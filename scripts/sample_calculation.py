"""Generate reproducible samples without writing NetCDF or altering reference inputs."""
from pathlib import Path
import sys,json,platform,importlib.metadata
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from atrxl.export import calculate,csv_bytes,json_bytes,results_json
from atrxl.plotting import png
from atrxl.weights import delta_pair_weight
root=Path(__file__).resolve().parents[1]/'results'
root.mkdir(exist_ok=True);(root/'figures').mkdir(exist_ok=True)
point,meta=calculate('Br-Br',4.,13.45,12.)
(root/'sample_point.json').write_bytes(results_json(point,meta))
(root/'sample_point.csv').write_bytes(csv_bytes(point))
curve,cm=calculate('Br-Br',np.linspace(.2,6,291),13.45,12.)
(root/'sample_curve.csv').write_bytes(csv_bytes(curve));(root/'sample_metadata.json').write_bytes(json_bytes(cm))
(root/'figures/benchmark.png').write_bytes(png(curve.q_Ainv,{'W far':curve.W_far_e2,'W near':curve.W_near_e2},'Physical q (Å⁻¹)','Br–Br · XrayDB/Chantler 4.5.8'))
series={s:delta_pair_weight('Br-Br',curve.q_Ainv.to_numpy(),13.45,12,s) for s in ('XrayDB','Henke','NIST')}
(root/'figures/source_contrast.png').write_bytes(png(curve.q_Ainv,series,'Physical q (Å⁻¹)','Br–Br ΔW · secondary channels: Br-only swap'))
versions={'python':platform.python_version(),**{p:importlib.metadata.version(p) for p in ('xraydb','numpy','scipy','pandas','streamlit','plotly','matplotlib','pytest','xarray')}}
(root/'logs/versions.json').write_bytes(json_bytes(versions))
print(point.to_string(index=False));print(versions)
