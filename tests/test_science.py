import hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
import xarray as xr
import xraydb
from atrxl import atomic_factor, pair_weight, delta_pair_weight
from atrxl.factors import factor_components
from atrxl.validation import CalculationError, PAIRS, ELEMENTS, q_grid
from atrxl.sources import LEGACY_SOURCES as SOURCES,load_table,anomalous
from atrxl.export import calculate,csv_bytes,results_json
ROOT=Path(__file__).resolve().parents[1]

@pytest.mark.parametrize('source,expected',[('XrayDB',-155.27584665917914),('Henke_Br_swap',-159.315760),('NIST_Br_swap',-164.812430)])
def test_published_benchmarks(source,expected):
    assert delta_pair_weight('Br-Br',4,13.45,12,source)==pytest.approx(expected,abs=5e-7)

@pytest.mark.parametrize('element',ELEMENTS)
def test_q_conversion_and_forward_limit(element):
    q=np.array([0.,0.2,1.234567,4.,6.])
    f0,_,_=factor_components(element,q,12)
    np.testing.assert_allclose(f0,xraydb.f0(element,q/(4*np.pi)),rtol=1e-14)
    assert f0[0]==pytest.approx(xraydb.atomic_number(element),abs=0.01)

@pytest.mark.parametrize('source',SOURCES)
@pytest.mark.parametrize('pair',PAIRS)
def test_independent_identities(pair,source):
    q=np.array([0.,.2,2.71828,4.,6.]);a,b=pair.split('-')
    fa=atomic_factor(a,q,13.45,source);fb=atomic_factor(b,q,13.45,source)
    w=pair_weight(pair,q,13.45,source)
    np.testing.assert_allclose(w,fa.real*fb.real+fa.imag*fb.imag,atol=1e-12)
    np.testing.assert_array_equal(w,pair_weight((b,a),q,13.45,source))
    if a==b: np.testing.assert_allclose(w,abs(fa)**2,atol=1e-12)
    np.testing.assert_array_equal(delta_pair_weight(pair,q,12,12,source),np.zeros(5))
    np.testing.assert_allclose(delta_pair_weight(pair,q,13.45,12,source),w-pair_weight(pair,q,12,source),atol=1e-12)
    assert np.any(w-fa.real*fb.real>0)
    if 'Br' not in (a,b): np.testing.assert_array_equal(w,pair_weight(pair,q,13.45,'XrayDB'))

def test_native_extent_and_batch_independence():
    q=3.14159265;e=13.41789
    f=atomic_factor('Br',q,e)
    expected=xraydb.f0('Br',q/(4*np.pi))[0]+xraydb.f1_chantler('Br',np.array([11830.,e*1000,13460.]))[1]+1j*xraydb.f2_chantler('Br',e*1000)
    assert f==pytest.approx(expected,abs=1e-13)
    np.testing.assert_array_equal(atomic_factor('Br',[q,q],[e,e]),[f,f])
    # This explicit difference is why the documented Phase 1 extent is necessary.
    assert abs(xraydb.f1_chantler('Br',13450.)-xraydb.f1_chantler('Br',[11830.,13450.,13460.])[1])>0.01

@pytest.mark.parametrize('source',SOURCES)
@pytest.mark.parametrize('pair',['Br-Br','Br-C','Br-O','Br-H'])
def test_full_netcdf_weight_agreement(source,pair):
    channel={'XrayDB':'XrayDB_primary','Henke_Br_swap':'Henke_Br_swap','NIST_Br_swap':'NIST_Br_swap'}[source]
    with xr.open_dataset(ROOT/'data/reference/contrast_weights.nc') as ds:
        q=ds.q_Ainv.values
        for start in range(0,404,20):
            e=ds.energy_eV.values[start:start+20]/1000
            actual=pair_weight(pair,q[None,:],e[:,None],source)
            expected=ds.W.sel(source=channel,pair=pair).isel(energy_eV=slice(start,start+20)).values
            # Archived swap factors were serialized to ~12 significant decimal digits
            # (max observed fprime residual 5e-12 e); 1e-9 e² bounds propagated rounding.
            np.testing.assert_allclose(actual,expected,rtol=2e-13,atol=1e-9)

def test_metrics_and_runtime_csv():
    m=json.loads((ROOT/'data/reference/calculation_metrics.json').read_text())
    for pair,values in m['q4_contrasts'].items():
        assert delta_pair_weight(pair,4,13.45,12)==pytest.approx(values['delta_q4_e2'],abs=1e-11)
    raw=pd.read_csv(ROOT/'data/reference/xraydb_runtime_comparison.csv')
    _,fp,fpp=factor_components('Br',4,raw.energy_eV.to_numpy()/1000)
    np.testing.assert_allclose(fp,raw.fprime_xraydb_e,atol=1e-12)
    np.testing.assert_allclose(fpp,raw.fdoubleprime_xraydb_e,atol=1e-12)

@pytest.mark.parametrize('source',['Henke_Br_swap','NIST_Br_swap'])
def test_original_knots_and_interpolation(source):
    a,m=load_table(source)
    knots=a[(a[:,0]>=m['domain_eV'][0])&(a[:,0]<=m['domain_eV'][1])]
    fp,fpp=anomalous('Br',knots[:,0]/1000,source)
    np.testing.assert_allclose(fp,knots[:,1]-35,atol=1e-12)
    np.testing.assert_allclose(fpp,knots[:,2],atol=1e-13)
    e=(a[-2,0]+a[-1,0])/2
    fp,fpp=anomalous('Br',e/1000,source)
    assert fp==pytest.approx((a[-2,1]+a[-1,1])/2-35)
    assert fpp==pytest.approx(np.exp(np.interp(np.log(e),np.log(a[-2:,0]),np.log(a[-2:,2]))))
    for e in [11.0,m['domain_eV'][1]/1000+.00001,14.]:
        with pytest.raises(CalculationError): pair_weight('Br-Br',4,e,source,mode='Extended')

@pytest.mark.parametrize('q,e',[(0,12),(.2,11.83),(6,13.46)])
def test_default_boundaries(q,e): assert np.isfinite(pair_weight('Br-Br',q,e))

@pytest.mark.parametrize('q,e',[(.1,12),(6.1,12),(4,11.82),(4,13.461)])
def test_extended_policy(q,e):
    with pytest.raises(CalculationError): pair_weight('Br-Br',q,e)
    assert np.isfinite(pair_weight('Br-Br',q,e,mode='Extended'))

@pytest.mark.parametrize('q',[-1,np.nan,np.inf,1+1j,'4',True,[],np.ones(10001)])
def test_invalid_q(q):
    with pytest.raises(CalculationError): pair_weight('Br-Br',q,12,mode='Extended')

@pytest.mark.parametrize('e',[-1,0,np.inf,np.nan,1e-6,1e6])
def test_invalid_energy(e):
    with pytest.raises(CalculationError): pair_weight('Br-Br',4,e,mode='Extended')

@pytest.mark.parametrize('source',SOURCES)
def test_edge_crossing(source):
    with pytest.raises(CalculationError,match='K-edge'): delta_pair_weight('Br-Br',4,13.5,12,source,mode='Extended')
    with pytest.raises(CalculationError,match='K-edge'): pair_weight('Br-Br',4,[13.4,13.5],source,mode='Extended')

def test_invalid_shapes_sources_and_scans():
    for pair in ['Xx-Br','Br','Br-Br-C',None]:
        with pytest.raises(CalculationError): pair_weight(pair,4,12)
    with pytest.raises(CalculationError): pair_weight('Br-Br',[1,2],[12,12,12])
    with pytest.raises(CalculationError): pair_weight('Br-Br',np.ones((101,1)),np.ones((1,101))*12)
    with pytest.raises(CalculationError): pair_weight('C-C',4,12,'unavailable')
    for args in [(0,6,0),(0,6,-1),(6,0,.1),(0,100,1e-8)]:
        with pytest.raises(CalculationError): q_grid(*args)
    np.testing.assert_allclose(q_grid(.2,6,.02)[[0,-1]],[.2,6])

def test_export_full_precision():
    frame,meta=calculate('Br-Br',[0.2,4.,6.],13.45,12)
    from io import BytesIO
    read=pd.read_csv(BytesIO(csv_bytes(frame)),float_precision='round_trip')
    np.testing.assert_array_equal(read.delta_W_e2,frame.delta_W_e2)
    payload=json.loads(results_json(frame,meta))
    assert payload['results'][1]['delta_W_e2']==frame.delta_W_e2.iloc[1]
    assert set(read.calculation_id)=={meta['calculation_id']}
    assert meta['source']['sha256'] and meta['timestamp_utc']

def test_reference_original_hashes():
    for line in (ROOT/'SHA256SUMS.txt').read_text().splitlines():
        sha,path=line.split(maxsplit=1)
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==sha

def test_contrast_combined_allocation_guard():
    with pytest.raises(CalculationError,match='allocation'):
        delta_pair_weight('Br-Br',4,np.ones((101,1))*13,np.ones((1,101))*12)

def test_source_integrity_unavailable_and_packaged_copies(tmp_path,monkeypatch):
    import atrxl.sources as sources
    import shutil
    for s in ('Henke_Br_swap','NIST_Br_swap'):
        _,m=load_table(s)
        assert (ROOT/'data/sources'/m['file']).read_bytes()==(ROOT/'src/atrxl/data_sources'/m['file']).read_bytes()
    shutil.copytree(sources.ROOT,tmp_path/'tables')
    monkeypatch.setattr(sources,'ROOT',tmp_path/'tables')
    sources.load_table.cache_clear()
    (tmp_path/'tables/henke_br.nff').write_text('corrupt')
    try:
        with pytest.raises(CalculationError,match='SHA-256'):
            pair_weight('C-C',4,12,'Henke_Br_swap')
    finally: sources.load_table.cache_clear()

def test_png_and_metadata_exports():
    from atrxl.plotting import png
    from PIL import Image
    from io import BytesIO
    frame,meta=calculate('Br-C',[.2,4,6],13.45,12,'Henke_Br_swap')
    data=png(frame.q_Ainv,{'ΔW':frame.delta_W_e2},'Physical q (Å⁻¹)','Henke Br-only swap')
    image=Image.open(BytesIO(data));assert image.format=='PNG' and image.width>=1000
    assert meta['source']['scope'].startswith('Br-only')
