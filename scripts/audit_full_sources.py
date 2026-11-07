"""Read-only source inventory and reproducible v1.1 mode-specific examples."""
from pathlib import Path
import sys,json,platform,importlib.metadata
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
import pandas as pd
from atrxl.sources import SOURCES,supported_elements,element_domain,full_table,source_metadata,source_availability
from atrxl.export import calculate_w,calculate,csv_bytes,results_json,energy_scan
from atrxl.plotting import png
root=Path(__file__).resolve().parents[1]/'results/v1.1';(root/'figures').mkdir(parents=True,exist_ok=True)
rows=[]
for source in SOURCES:
    for element in supported_elements(source):
        d=element_domain(source,element)
        m=source_metadata(source,(element,))
        row={'source':source,'element':element,'min_keV':d['min_keV'],'max_keV':d['max_keV'],'excluded_intervals':len(d['forbidden_intervals_keV'])}
        if source!='XrayDB':
            table=m['tables'][element]
            row.update({'Z':table['Z'],'valid_knots':table['valid_rows'],'raw_rows':table['raw_rows'],'sha256':table['sha256'],'url':table['url'],'version':table['version']})
        else: row.update({'sha256':m['sha256'],'version':m['database_version'],'url':m['url']})
        rows.append(row)
(root/'source_inventory.csv').write_bytes(csv_bytes(pd.DataFrame(rows)))
point,meta=calculate_w('Fe-Au',4,12,'Henke')
(root/'Fe_Au_Henke_W.json').write_bytes(results_json(point,meta))
q=np.linspace(.2,6,291);curve,cm=calculate_w('Fe-Au',q,12,'Henke')
(root/'Fe_Au_Henke_W_curve.csv').write_bytes(csv_bytes(curve))
(root/'Fe_Au_Henke_W_curve_metadata.json').write_text(json.dumps(cm,indent=2))
(root/'figures/Fe_Au_Henke_W.png').write_bytes(png(q,{'W(12 keV)':curve.W_e2},'Physical q (Å⁻¹)','Fe–Au · both Henke atoms · common native f0'))
series={};sample=[]
for source in (*SOURCES,'Henke_Br_swap','NIST_Br_swap'):
    frame,meta=calculate('Br-C',q,13.45,12,source)
    series[source]=frame.delta_W_e2
    sample.append({'source':source,'delta_W_BrC_q4_e2':float(calculate('Br-C',4,13.45,12,source)[0].delta_W_e2.iloc[0])})
(root/'figures/Br_C_full_and_legacy.png').write_bytes(png(q,series,'Physical q (Å⁻¹)','Br–C · full-source vs explicit legacy Br-only ΔW'))
(root/'Br_C_full_and_legacy.csv').write_bytes(csv_bytes(pd.DataFrame(sample)))
es,em=energy_scan('Br-Br',4,np.linspace(12,13.45,160),'XrayDB',calculation_mode='ΔW',e_far=12)
(root/'Br_Br_delta_energy.json').write_bytes(results_json(es,em))
versions={'python':platform.python_version(),**{p:importlib.metadata.version(p) for p in ('atrxl-calculator','xraydb','numpy','scipy','pandas','streamlit','plotly','matplotlib','pytest','xarray','playwright')}}
(root/'logs/versions.json').write_text(json.dumps(versions,indent=2))
print(json.dumps({'counts':{s:len(supported_elements(s)) for s in SOURCES},'native_unavailable':source_availability('XrayDB')['unavailable'],'samples':sample,'versions':versions},indent=2))
