import numpy as np
import xraydb
from .validation import validate, ELEMENTS, CalculationError
from .sources import anomalous, native_elements

def factor_components(element, q, energy, source="XrayDB", *, mode="Default"):
    """Return float64 f0, f′, f″ in electrons, broadcast over physical q and keV."""
    if element not in native_elements(): raise CalculationError(f"{element}: neutral f0 or source anomalous factors unavailable")
    q,e,_=validate(q,energy,mode)
    q,e=np.broadcast_arrays(q,e)
    with np.errstate(over="ignore",invalid="ignore"):
        f0=np.asarray(xraydb.f0(element,q.ravel()/(4*np.pi))).reshape(q.shape)
        fp,fpp=anomalous(element,e,source)
    if not all(np.all(np.isfinite(v)) for v in (f0,fp,fpp)): raise CalculationError("Nonfinite native factor evaluation")
    return f0,fp,fpp

def atomic_factor(element,q,energy,source="XrayDB",*,mode="Default"):
    f0,fp,fpp=factor_components(element,q,energy,source,mode=mode)
    return f0+fp+1j*fpp
