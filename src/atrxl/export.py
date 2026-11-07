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
from .sources import source_metadata, pair_domain

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
    if source.endswith("_Br_swap") and "Br" not in parse_pair(pair): warnings.append("Br-only source swap does not change a non-Br pair.")
    if any(atom not in ("Br","C","O","H") for atom in parse_pair(pair)): warnings.append("Expanded species: source-level numerical validation only; not Phase 1 research validation.")
    warnings.append("Neutral-atom near-edge predictions do not model chemical-state spectra or experimental uncertainty.")
    metadata={"calculation_id":uid,"calculator_version":"1.1.0","timestamp_utc":datetime.now(timezone.utc).isoformat(),"pair":list(parse_pair(pair)),"q_Ainv":q.tolist(),"energy_near_keV":np.asarray(e_near).tolist(),"energy_far_keV":np.asarray(e_far).tolist(),"mode":mode,"source":source_metadata(source,parse_pair(pair)),"energy_domain":pair_domain(pair,source),"calculation_mode":"ΔW","xraydb_version":xraydb.__version__,"database_version":xraydb.get_xraydb().get_version(),"Br_K_edge_eV":xraydb.xray_edge("Br","K").energy,"units":{"f":"electrons","W":"electrons²","q":"Å⁻¹","energy":"keV"},"warnings":warnings,"convention":"near-minus-far; one unordered pair; no multiplicity or geometry"}
    return frame,metadata

def csv_bytes(frame):
    return frame.to_csv(index=False,float_format="%.17g").encode("utf-8-sig")

def results_json(frame,metadata):
    # pandas converts undefined relative values to JSON null.
    return json_bytes({"metadata":metadata,"results":[{k: (None if isinstance(v,float) and not np.isfinite(v) else v) for k,v in row.items()} for row in frame.to_dict(orient="records")]})

def calculate_w(pair,q,energy,source='XrayDB',mode='Default'):
    """Single-energy W results: no near/far inputs or ΔW fields."""
    q,e,warnings=validate(q,energy,mode)
    q,e=np.broadcast_arrays(q,e)
    w=pair_weight(pair,q,e,source,mode=mode)
    uid=str(uuid4())
    frame=pd.DataFrame({'calculation_id':uid,'q_Ainv':q.ravel(),'energy_keV':e.ravel(),'W_e2':np.ravel(w)})
    for atom in dict.fromkeys(parse_pair(pair)):
        f0,fp,fpp=factor_components(atom,q,e,source,mode=mode)
        for name,value in [('f0',f0),('fprime',fp),('fdoubleprime',fpp),('f_real',f0+fp),('f_imag',fpp)]:
            frame[f'{atom}_{name}_e']=np.ravel(value)
    if any(atom not in ('Br','C','O','H') for atom in parse_pair(pair)):
        warnings.append('Expanded species: source-level numerical validation only; not Phase 1 research validation.')
    warnings.append('Neutral-atom coefficients; near-edge chemical-state and experimental validity are not established.')
    meta={'calculation_id':uid,'calculator_version':'1.1.0','timestamp_utc':datetime.now(timezone.utc).isoformat(),
          'pair':list(parse_pair(pair)),'q_Ainv':q.tolist(),'energy_keV':e.tolist(),'calculation_mode':'W','mode':mode,
          'source':source_metadata(source,parse_pair(pair)),'energy_domain':pair_domain(pair,source),
          'xraydb_version':xraydb.__version__,'database_version':xraydb.get_xraydb().get_version(),
          'units':{'f':'electrons','W':'electrons²','q':'Å⁻¹','energy':'keV'},'warnings':list(dict.fromkeys(warnings)),
          'convention':'one unordered pair; no multiplicity or geometry'}
    return frame,meta

def energy_scan(pair,q,energies,source='XrayDB',mode='Default',*,calculation_mode='W',e_far=None):
    """Full-precision energy curve for either W(E) or ΔW(Enear; fixed Efar)."""
    from .sources import guard_source_interval
    from .validation import numeric
    energies=numeric(energies,'energy scan')
    if energies.ndim!=1: raise ValueError('energy scan must be one-dimensional')
    if calculation_mode not in ('W','ΔW'): raise ValueError('calculation mode must be W or ΔW')
    guard_source_interval(pair,energies,source)
    if calculation_mode=='W':
        frame,meta=calculate_w(pair,q,energies,source,mode)
    else:
        if e_far is None: raise ValueError('ΔW energy scan requires fixed Efar')
        frame,meta=calculate(pair,np.full_like(energies,float(q)),energies,e_far,source,mode)
        frame['energy_near_keV']=energies
    meta['curve_kind']='energy_scan'
    return frame,meta
