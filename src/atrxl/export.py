"""Full precision values linked to calculation metadata by UUID."""
from datetime import datetime,timezone
from uuid import uuid4
import json
import numpy as np
import pandas as pd
import xraydb
from .weights import pair_weight,delta_pair_weight
from .factors import factor_components
from .validation import validate,parse_pair
from .sources import source_metadata

def json_bytes(value):
    return json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False).encode("utf-8")

def calculate(pair,q,e_near,e_far,source="XrayDB",mode="Default"):
    q=np.asarray(q,dtype=float)
    _,_,wn=validate(q,e_near,mode);_,_,wf=validate(q,e_far,mode)
    near=pair_weight(pair,q,e_near,source,mode=mode)
    far=pair_weight(pair,q,e_far,source,mode=mode)
    delta=delta_pair_weight(pair,q,e_near,e_far,source,mode=mode)
    relative=np.full(np.shape(delta),np.nan)
    np.divide(100*delta,far,out=relative,where=np.abs(far)>1e-12)
    uid=str(uuid4())
    frame=pd.DataFrame({"calculation_id":uid,"q_Ainv":q.ravel(),"W_far_e2":np.ravel(far),"W_near_e2":np.ravel(near),"delta_W_e2":np.ravel(delta),"relative_change_percent":np.ravel(relative)})
    for atom in dict.fromkeys(parse_pair(pair)):
        for label,energy in (("far",e_far),("near",e_near)):
            f0,fp,fpp=factor_components(atom,q,energy,source,mode=mode)
            for name,values in (("f0",f0),("fprime",fp),("fdoubleprime",fpp),("f_real",f0+fp),("f_imag",fpp)):
                frame[f"{atom}_{label}_{name}_e"]=np.ravel(values)
    warnings=list(dict.fromkeys(wn+wf))
    if source!="XrayDB" and "Br" not in parse_pair(pair): warnings.append("Br-only source swap does not change a non-Br pair.")
    warnings.append("Neutral-atom near-edge predictions do not model chemical-state spectra or experimental uncertainty.")
    metadata={"calculation_id":uid,"calculator_version":"1.0.0","timestamp_utc":datetime.now(timezone.utc).isoformat(),"pair":list(parse_pair(pair)),"q_Ainv":q.tolist(),"energy_near_keV":e_near,"energy_far_keV":e_far,"mode":mode,"source":source_metadata(source),"xraydb_version":xraydb.__version__,"database_version":xraydb.get_xraydb().get_version(),"Br_K_edge_eV":xraydb.xray_edge("Br","K").energy,"units":{"f":"electrons","W":"electrons²","q":"Å⁻¹","energy":"keV"},"warnings":warnings,"convention":"near-minus-far; one unordered pair; no multiplicity or geometry"}
    return frame,metadata

def csv_bytes(frame):
    return frame.to_csv(index=False,float_format="%.17g").encode("utf-8-sig")

def results_json(frame,metadata):
    # pandas converts undefined relative values to JSON null.
    return json_bytes({"metadata":metadata,"results":[{k: (None if isinstance(v,float) and not np.isfinite(v) else v) for k,v in row.items()} for row in frame.to_dict(orient="records")]})
