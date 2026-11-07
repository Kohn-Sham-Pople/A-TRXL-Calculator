"""Local full-element source adapters and explicit Phase 1 legacy channels.

Raw full f1 is converted to anomalous f′ by subtracting Z. No source fallback
or extrapolation is performed; source/gap provenance is retained in metadata.
"""
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
    source = source.removesuffix("_Br_swap")
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

def _legacy_metadata(source):
    if source=="XrayDB":
        db=Path(xraydb.get_xraydb().dbname)
        return {"source":"XrayDB/Chantler", "version":xraydb.__version__,"database_version":xraydb.get_xraydb().get_version(),"sha256":hashlib.sha256(db.read_bytes()).hexdigest(),"url":"https://xraypy.github.io/XrayDB/python.html","interpolation":"native API: f1 spline(s=0), Phase 1 extent anchors [11830,13460] eV within that interval; outside it scalar native knots; f2 log-log; f0 Waasmaier-Kirfel", "scope":"neutral atoms", "retrieved_utc":json.loads((ROOT/"provenance.json").read_text())["Henke"]["retrieved_utc"],"retrieval_note":"XrayDB 4.5.8 wheel installed during this session; timestamp records acquisition session"}
    return {"source":source,**load_table(source)[1]}

def _legacy_anomalous(element, energy_keV, source):
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
    if element in ELEMENTS and 11830 <= energy_eV <= 13460:
        return float(xraydb.f1_chantler(element,np.array([11830.,energy_eV,13460.]))[1])
    return float(xraydb.f1_chantler(element,energy_eV))

# Explicit legacy channels retain the Phase 1 coefficients; the primary source
# names now always mean both-elements anomalous factors from that source.
LEGACY_SOURCES = ("XrayDB", "Henke_Br_swap", "NIST_Br_swap")

@lru_cache
def native_elements():
    """Intersection of actual neutral Waasmaier rows and Chantler f1/f2 rows."""
    db=xraydb.get_xraydb()
    with db.engine.connect() as connection:
        c=db.tables['Chantler'];w=db.tables['Waasmaier']
        chantler={r['element'] for r in connection.execute(c.select()).mappings()
                  if len(json.loads(r['energy']))>=4 and len(json.loads(r['f1']))==len(json.loads(r['energy']))==len(json.loads(r['f2']))}
        neutral={r['element'] for r in connection.execute(w.select()).mappings() if r['ion']==r['element']}
    return tuple(sorted(chantler&neutral,key=xraydb.atomic_number))

@lru_cache
def full_manifest():
    try:
        return json.loads((ROOT/'full/manifest.json').read_text(encoding='utf-8'))
    except (OSError,ValueError) as exc:
        raise CalculationError(f'Full-element source manifest unavailable: {exc}') from exc

@lru_cache
def full_table(source,element):
    try:
        m=full_manifest()['sources'][source][element]
        path=ROOT/'full'/m['file']
        if hashlib.sha256(path.read_bytes()).hexdigest()!=m['sha256']:
            raise ValueError('raw table SHA-256 mismatch')
        if source=='Henke':
            raw=np.loadtxt(path,skiprows=1)[:,:3]
            raw=raw[np.argsort(raw[:,0],kind='stable')]
            _,inverse,counts=np.unique(raw[:,0],return_inverse=True,return_counts=True)
            valid=np.isfinite(raw).all(axis=1)&(raw[:,1]>-9999)&(raw[:,2]>=0)&(counts[inverse]==1)
        else:
            raw=pd.read_csv(path)[['E','f1','f2']].to_numpy();raw[:,0]*=1000
            valid=np.isfinite(raw).all(axis=1)&(raw[:,2]>=0)
        a=raw[valid]
        if len(a)!=m['valid_rows'] or np.any(np.diff(a[:,0])<=0) or a[0,0]!=m['domain_eV'][0] or a[-1,0]!=m['domain_eV'][1]:
            raise ValueError('invalid factor availability or raw domain')
        return a,m
    except (KeyError,OSError,ValueError) as exc:
        raise CalculationError(f'{source}/{element} unavailable: {exc}') from exc

def source_availability(source):
    """List verified species and explicit missing-factor/integrity reasons."""
    if source=='XrayDB':
        symbols={xraydb.atomic_symbol(z) for z in range(1,99)}
        neutral=[el for el in xraydb.f0_ions() if el in symbols]
        return {'supported':list(native_elements()),'unavailable':{el:'No native Chantler f1/f2 row' for el in neutral if el not in native_elements()}}
    if source.endswith('_Br_swap') and source in LEGACY_SOURCES:
        load_table(source)
        return {'supported':list(ELEMENTS),'unavailable':{}}
    if source not in SOURCES: raise CalculationError(f'Unknown source {source}')
    result=[];missing={}
    for element in full_manifest()['sources'][source]:
        if element not in native_elements():
            missing[element]='Common neutral f0 or native metadata unavailable'
            continue
        try: full_table(source,element)
        except CalculationError as exc: missing[element]=str(exc)
        else: result.append(element)
    return {'supported':sorted(result,key=xraydb.atomic_number),'unavailable':missing}

def supported_elements(source):
    return tuple(source_availability(source)['supported'])

@lru_cache
def native_domain(element):
    if element not in native_elements(): raise CalculationError(f'XrayDB/{element}: neutral f0 or anomalous factors unavailable')
    grid=xraydb.chantler_energies(element)
    return float(grid[0]),float(grid[-1])

def element_domain(source,element):
    if source=='XrayDB' or (source in LEGACY_SOURCES and element!='Br'):
        lo,hi=native_domain(element)
        gaps=[]
    elif source in LEGACY_SOURCES:
        _,m=load_table(source);lo,hi=m['domain_eV'];gaps=[]
    else:
        _,m=full_table(source,element);lo,hi=m['domain_eV'];gaps=m['forbidden_intervals_eV']
    if element not in native_elements(): raise CalculationError(f'{element}: common neutral f0 unavailable')
    return {'element':element,'min_keV':lo/1000,'max_keV':hi/1000,
            'forbidden_intervals_keV':[[l/1000,h/1000,reason] for l,h,reason in gaps]}

def pair_domain(pair,source):
    from .validation import parse_pair
    atoms=parse_pair(pair)
    domains=[element_domain(source,atom) for atom in atoms]
    lo=max(d['min_keV'] for d in domains);hi=min(d['max_keV'] for d in domains)
    if lo>hi: raise CalculationError('No common energy domain for selected elements/source')
    gaps=sorted(set((max(lo,l),min(hi,h),d['element']+': '+reason)
                   for d in domains for l,h,reason in d['forbidden_intervals_keV'] if h>lo and l<hi))
    # Represent the usable bands as well as the overall envelope. Knot endpoints
    # remain usable; an ambiguous duplicate knot inside a gap does not.
    edges=pair_edges(atoms,source)
    scan_breaks=[(edge['energy_keV'],edge['energy_keV'],'native shell edge') for edge in edges if 'energy_keV' in edge and lo<edge['energy_keV']<hi]
    bands=[];cursor=lo
    for l,h,_ in sorted(gaps+scan_breaks):
        if l>cursor: bands.append([cursor,l])
        cursor=max(cursor,h)
    if cursor<hi: bands.append([cursor,hi])
    return {'source':source,'elements':list(atoms),'min_keV':lo,'max_keV':hi,
            'element_domains':domains,'forbidden_intervals_keV':[list(g) for g in gaps],
            'usable_bands_keV':bands,'edges':edges,
            'warning':'Domain envelope contains excluded interpolation gaps. Use a single usable band; raw knots at a gap boundary remain usable. No extrapolation.'}

def pair_edges(pair,source):
    edges=[]
    for element in dict.fromkeys(pair):
        if source=='NIST':
            for label,e in full_table(source,element)[1]['edges_eV'].items():
                edges.append({'element':element,'shell':label,'energy_keV':e/1000,'provenance':'NIST archived shelldata.csv'})
        elif source=='Henke':
            for l,h,reason in full_table(source,element)[1]['forbidden_intervals_eV']:
                edges.append({'element':element,'shell':reason,'bracket_keV':[l/1000,h/1000],'provenance':'Henke raw-knot bracket or conservative shell locator; see manifest'})
        else:
            for label,edge in xraydb.xray_edges(element).items():
                edges.append({'element':element,'shell':label,'energy_keV':edge.energy/1000,'provenance':'installed XrayDB xray_edges'})
    return edges

def ensure_energy(element,energy_keV,source):
    from .validation import numeric
    e=numeric(energy_keV,'energy (keV)')
    d=element_domain(source,element)
    if np.any(e<d['min_keV']) or np.any(e>d['max_keV']):
        raise CalculationError(f"{source}/{element}: energy must lie in {d['min_keV']:.10g}–{d['max_keV']:.10g} keV; no extrapolation")
    # Check in original eV coordinates to avoid different floating-point rounding
    # of a knot converted keV -> eV -> keV. Exact raw knots are always retained.
    if source in ('Henke','NIST'):
        a,m=full_table(source,element);ev=e*1000
        knot=np.any(np.isclose(ev[...,None],a[:,0],rtol=0,atol=1e-9),axis=-1)
        for l,h,reason in m['forbidden_intervals_eV']:
            if np.any((ev>l+1e-9)&(ev<h-1e-9)&~knot):
                raise CalculationError(f'{source}/{element}: forbidden interpolation gap {l/1000:g}–{h/1000:g} keV ({reason})')

def guard_source_interval(pair,energies,source):
    """Reject curves/contrasts spanning a shell edge or tabulation gap."""
    from .validation import parse_pair,numeric
    atoms=parse_pair(pair);e=numeric(energies,'energy interval');lo=float(e.min());hi=float(e.max())
    if source in LEGACY_SOURCES and 'Br' in atoms:
        from .validation import guard_edge
        guard_edge(e)
    for atom in atoms: ensure_energy(atom,e,source)
    if lo==hi: return
    for edge in pair_edges(atoms,source):
        if 'energy_keV' in edge:
            crossing=lo<edge['energy_keV']<=hi
        else:
            a,b=edge['bracket_keV'];crossing=lo<b and hi>a
        if crossing:
            raise CalculationError(f"{source}/{edge['element']} {edge['shell']}: energy interval crosses an absorption edge/interpolation boundary; split into one usable band")
    for atom in atoms:
        for l,h,reason in element_domain(source,atom)['forbidden_intervals_keV']:
            if lo<h and hi>l:
                raise CalculationError(f'{source}/{atom}: energy interval crosses forbidden interpolation gap ({reason})')

def source_metadata(source,elements=None):
    if source in LEGACY_SOURCES:
        base=source.removesuffix('_Br_swap')
        result=_legacy_metadata(base)
        result['source']=source
        if source!='XrayDB': result['scope']='Br-only source swap (EXPLICIT LEGACY); non-Br anomalous factors are native XrayDB'
        result['f0_provenance']='XrayDB/Waasmaier–Kirfel neutral atom, physical q/(4π)'
        return result
    if source not in SOURCES: raise CalculationError(f'Unavailable source {source}')
    elements=elements or supported_elements(source)
    tables={el:full_table(source,el)[1] for el in dict.fromkeys(elements)}
    return {'source':source,'scope':'Both selected elements use source-specific anomalous factors; no Br-only swap',
            'version':full_manifest()['NIST_commit'] if source=='NIST' else 'official Henke current archive',
            'retrieved_utc':full_manifest()['retrieved_utc'],'tables':tables,
            'f0_provenance':_legacy_metadata('XrayDB'),
            'interpolation':'linear full f1 minus atomic Z; positive f2 log-log; zero-zero intervals constant; edge/gap interpolation forbidden'}

def anomalous(element,energy_keV,source):
    if source in LEGACY_SOURCES:
        if source!='XrayDB' and element not in ELEMENTS:
            raise CalculationError('Legacy Br-only source channels support only Phase 1 Br/C/O/H; select a full source for other elements')
        ensure_energy(element,energy_keV,source)
        return _legacy_anomalous(element,energy_keV,source.removesuffix('_Br_swap'))
    if source not in SOURCES: raise CalculationError(f'Unavailable data source: {source}')
    a,m=full_table(source,element);ensure_energy(element,energy_keV,source)
    e=np.asarray(energy_keV,dtype=float)*1000
    fp=np.interp(e,a[:,0],a[:,1])-m['Z']
    left=np.clip(np.searchsorted(a[:,0],e,side='right')-1,0,len(a)-2);right=left+1
    f2l,f2r=a[left,2],a[right,2]
    # Return exact knots including physically tabulated zero f2 values.
    fpp=np.zeros_like(e)
    positive=(f2l>0)&(f2r>0)
    t=(np.log(e)-np.log(a[left,0]))/(np.log(a[right,0])-np.log(a[left,0]))
    fpp[positive]=np.exp(np.log(f2l[positive])+t[positive]*(np.log(f2r[positive])-np.log(f2l[positive])))
    for i,v in enumerate(e.ravel()):
        j=np.searchsorted(a[:,0],v)
        if j<len(a) and abs(a[j,0]-v)<1e-9: fpp.ravel()[i]=a[j,2]
        elif j>0 and abs(a[j-1,0]-v)<1e-9: fpp.ravel()[i]=a[j-1,2]
    if not np.all(np.isfinite(fp)) or not np.all(np.isfinite(fpp)):
        raise CalculationError(f'{source}/{element}: nonfinite factor')
    return fp,fpp
