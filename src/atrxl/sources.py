"""Verified raw Br table adapters. No NetCDF interpolation and no fallback."""
from pathlib import Path
from functools import lru_cache
import hashlib,json
import numpy as np
import pandas as pd
import xraydb
from .validation import CalculationError, ELEMENTS
ROOT = Path(__file__).resolve().parent / "data_sources"
# Source checkout layout; packaged copies are included for wheel installs.
if not ROOT.exists(): ROOT = Path(__file__).resolve().parents[2] / "data" / "sources"
SOURCES = ("XrayDB", "Henke", "NIST")

@lru_cache
def load_table(source):
    if source not in ("Henke","NIST"): raise CalculationError(f"Unknown archived source: {source}")
    try:
        m=json.loads((ROOT/"provenance.json").read_text())[source]
        p=ROOT/m["file"]
        if hashlib.sha256(p.read_bytes()).hexdigest()!=m["sha256"]: raise ValueError("SHA-256 mismatch")
        if source=="Henke":
            a=np.loadtxt(p,skiprows=1)[:,:3]
        else:
            a=pd.read_csv(p)[["E","f1","f2"]].to_numpy(); a[:,0]*=1000
        lo,hi=m["domain_eV"]
        # Include bounding original knots below lower limit; stop at final verified pre-edge knot.
        a=a[(a[:,0]>=11000)&(a[:,0]<=hi+1e-6)]
        if len(a)<2 or not np.all(np.isfinite(a)) or np.any(np.diff(a[:,0])<=0) or np.any(a[:,2]<=0) or a[0,0]>lo or a[-1,0]<hi-1e-6:
            raise ValueError("invalid raw table coverage")
        return a,m
    except (OSError,ValueError,KeyError) as exc:
        raise CalculationError(f"{source} unavailable: {exc}") from exc

def source_metadata(source):
    if source=="XrayDB":
        db=Path(xraydb.get_xraydb().dbname)
        return {"source":"XrayDB/Chantler", "version":xraydb.__version__,"database_version":xraydb.get_xraydb().get_version(),"sha256":hashlib.sha256(db.read_bytes()).hexdigest(),"url":"https://xraypy.github.io/XrayDB/python.html","interpolation":"native API: f1 spline(s=0), Phase 1 extent anchors [11830,13460] eV within that interval; outside it scalar native knots; f2 log-log; f0 Waasmaier-Kirfel", "scope":"neutral atoms", "retrieved_utc":json.loads((ROOT/"provenance.json").read_text())["Henke"]["retrieved_utc"],"retrieval_note":"XrayDB 4.5.8 wheel installed during this session; timestamp records acquisition session"}
    return {"source":source,**load_table(source)[1]}

def anomalous(element, energy_keV, source):
    if source not in SOURCES: raise CalculationError(f"Unavailable data source: {source}")
    if source!="XrayDB": load_table(source)  # fail even for a non-Br pair if source unavailable
    e=np.asarray(energy_keV,dtype=float)*1000
    grid=xraydb.chantler_energies(element)
    if np.any(e<grid[0]) or np.any(e>grid[-1]):
        raise CalculationError(f"{element}: native tabulated energy domain is {grid[0]/1000:g}–{grid[-1]/1000:g} keV")
    if element=="Br" and source!="XrayDB":
        a,m=load_table(source);lo,hi=m["domain_eV"]
        if np.any(e<lo) or np.any(e>hi): raise CalculationError(f"{source} Br-only verified pre-edge domain: {lo/1000:g}–{hi/1000:g} keV; no extrapolation")
        return np.interp(e,a[:,0],a[:,1])-35, np.exp(np.interp(np.log(e),np.log(a[:,0]),np.log(a[:,2])))
    # Use a common native evaluation context covering each element's table.
    # XrayDB chooses spline knots from the energy extent, so deterministic batching
    # uses one invocation per unique energy, independent of array shape/order.
    flat=e.ravel();unique,inv=np.unique(flat,return_inverse=True)
    fp=np.array([native_fprime(element,float(v)) for v in unique])[inv].reshape(e.shape)
    fpp=np.array([xraydb.f2_chantler(element,float(v)) for v in unique])[inv].reshape(e.shape)
    return fp,fpp

@lru_cache(maxsize=40000)
def native_fprime(element, energy_eV):
    # Phase 1 batch extent is part of the reproducible native evaluation protocol.
    if 11830 <= energy_eV <= 13460:
        return float(xraydb.f1_chantler(element,np.array([11830.,energy_eV,13460.]))[1])
    return float(xraydb.f1_chantler(element,energy_eV))
