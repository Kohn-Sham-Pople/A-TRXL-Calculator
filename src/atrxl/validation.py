"""Input policy shared by GUI and standalone numerical API."""
import numpy as np
import xraydb
ELEMENTS = ("Br", "C", "O", "H")
PAIRS = tuple(f"{a}-{b}" for i,a in enumerate(ELEMENTS) for b in ELEMENTS[i:])
MAX_POINTS = 10000
class CalculationError(ValueError):
    pass

def numeric(value, name):
    raw = np.asarray(value)
    if raw.dtype.kind not in "iuf" or raw.size == 0 or raw.size > MAX_POINTS or raw.ndim > 2:
        raise CalculationError(f"{name}: real numeric scalar/array required, at most {MAX_POINTS} points and 2 dimensions")
    a = raw.astype(np.float64)
    if not np.all(np.isfinite(a)):
        raise CalculationError(f"{name}: values must be finite")
    return a

def parse_pair(pair):
    p = pair.replace("–", "-").split("-") if isinstance(pair,str) else pair
    if not isinstance(p,(list,tuple)) or len(p)!=2 or any(x not in ELEMENTS for x in p):
        raise CalculationError("pair must contain two of Br, C, O, H")
    return tuple(p)

def validate(q, energy, mode="Default"):
    if mode not in ("Default", "Extended"):
        raise CalculationError("mode must be Default or Extended")
    q,e = numeric(q,"physical q"),numeric(energy,"energy (keV)")
    if np.any(q < 0): raise CalculationError("physical q must be nonnegative")
    if np.any(e <= 0): raise CalculationError("energy must be positive")
    warnings=[]
    outside_q = np.any((q!=0)&((q<0.2)|(q>6)))
    outside_e = np.any((e<11.83)|(e>13.46))
    if mode == "Default" and (outside_q or outside_e):
        raise CalculationError("Default domain: q = 0 (forward diagnostic) or 0.2–6.0 Å⁻¹; energy = 11.83–13.46 keV")
    if outside_q: warnings.append("Outside Phase 1 validated q-range; numerical f0 evaluation does not establish physical accuracy.")
    if outside_e: warnings.append("Outside Phase 1 validated energy-range; neutral-atom near-edge predictions have source and chemical-state uncertainty.")
    if np.any(q==0): warnings.append("q=0 is an analytical forward-limit diagnostic.")
    try:
        shape=np.broadcast_shapes(q.shape,e.shape)
    except ValueError as exc: raise CalculationError("q and energy shapes cannot broadcast") from exc
    if np.prod(shape)>MAX_POINTS: raise CalculationError("broadcast allocation exceeds 10000 points")
    return q,e,warnings

def guard_edge(energies):
    e=numeric(energies,"energies")*1000
    edge=xraydb.xray_edge("Br","K").energy
    if e.min() < edge <= e.max():
        raise CalculationError(f"Energy interval crosses Br K-edge ({edge/1000:g} keV); edge-crossing calculations are blocked.")

def q_grid(lo, hi, step):
    lo,hi,step=[float(numeric(x,"q scan input")) for x in (lo,hi,step)]
    if lo<0 or hi<lo or step<=0: raise CalculationError("q scan requires 0 ≤ qmin ≤ qmax and step > 0")
    count=(hi-lo)/step
    if not np.isfinite(count) or count>MAX_POINTS-1: raise CalculationError("q scan exceeds 10000 points")
    return lo+np.arange(int(np.floor(count+1e-10))+1)*step
