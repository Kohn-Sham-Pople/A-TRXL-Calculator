"""Read-only direct API vs original NetCDF audit; emits maximum residuals."""
from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import xarray as xr
import numpy as np
from atrxl import pair_weight
root=Path(__file__).resolve().parents[1]
report={};count=0
with xr.open_dataset(root/'data/reference/contrast_weights.nc') as ds:
    for source,channel in [('XrayDB','XrayDB_primary'),('Henke_Br_swap','Henke_Br_swap'),('NIST_Br_swap','NIST_Br_swap')]:
        for pair in ds.pair.values:
            max_error=0.
            for start in range(0,ds.sizes['energy_eV'],20):
                es=ds.energy_eV.values[start:start+20]/1000
                actual=pair_weight(str(pair),ds.q_Ainv.values[None,:],es[:,None],source)
                original=ds.W.sel(source=channel,pair=pair).isel(energy_eV=slice(start,start+20)).values
                residual=float(np.max(np.abs(actual-original)))
                max_error=max(max_error,residual);count+=actual.size
                assert np.allclose(actual,original,rtol=2e-13,atol=1e-9)
            report[f'{source}/{pair}']={'max_absolute_residual_e2':max_error}
report['compared_values']=count
report['tolerance']={'absolute_e2':1e-9,'relative':2e-13}
print(json.dumps(report,indent=2))
