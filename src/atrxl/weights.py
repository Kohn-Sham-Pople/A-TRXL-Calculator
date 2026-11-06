import numpy as np
from .factors import atomic_factor
from .validation import parse_pair,guard_edge,validate, CalculationError, MAX_POINTS

def pair_weight(pair,q,energy,source="XrayDB",*,mode="Default"):
    """Re[f_a f_b*], electrons² per unordered pair; q Å⁻¹, energy keV."""
    a,b=parse_pair(pair)
    guard_edge(energy)
    return np.real(atomic_factor(a,q,energy,source,mode=mode)*np.conj(atomic_factor(b,q,energy,source,mode=mode)))

def delta_pair_weight(pair,q,e_near,e_far,source="XrayDB",*,mode="Default"):
    validate(q,e_near,mode);validate(q,e_far,mode)
    try:
        shape=np.broadcast_shapes(np.shape(q),np.shape(e_near),np.shape(e_far))
    except ValueError as exc: raise CalculationError("q, near and far shapes cannot broadcast") from exc
    if np.prod(shape)>MAX_POINTS: raise CalculationError("contrast allocation exceeds 10000 points")
    guard_edge(np.concatenate([np.asarray(e_near).ravel(),np.asarray(e_far).ravel()]))
    return pair_weight(pair,q,e_near,source,mode=mode)-pair_weight(pair,q,e_far,source,mode=mode)

pair_weight_curve = pair_weight
delta_weight_curve = delta_pair_weight
