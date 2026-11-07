"""v1.1 source-level tests, distinct from Phase 1 research validation."""
from io import BytesIO
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
import pytest
import xraydb
from atrxl import atomic_factor,pair_weight,delta_pair_weight
from atrxl.factors import factor_components
from atrxl.sources import (SOURCES,supported_elements,native_elements,full_manifest,full_table,
                          element_domain,pair_domain,guard_source_interval,anomalous,source_metadata)
from atrxl.validation import CalculationError
from atrxl.export import calculate,calculate_w,energy_scan,csv_bytes,results_json

ATOMS=tuple(xraydb.atomic_symbol(z) for z in range(1,93))

@pytest.mark.parametrize('source',SOURCES)
def test_discovery_matches_actual_data(source):
    expected=native_elements() if source=='XrayDB' else tuple(full_manifest()['sources'][source])
    assert set(supported_elements(source))==set(expected)==set(ATOMS)
    for element in expected:
        assert np.isfinite(xraydb.f0(element,0)).all()
        if source!='XrayDB':
            a,m=full_table(source,element)
            assert len(a)==m['valid_rows'] and np.isfinite(a).all()
            assert np.all(np.diff(a[:,0])>0)
        else:
            db=xraydb.get_xraydb();row=db.get_cache('Chantler',column='element',value=element)[0]
            assert np.isfinite(json.loads(row.f1)).all() and np.isfinite(json.loads(row.f2)).all()

@pytest.mark.parametrize('source',SOURCES)
@pytest.mark.parametrize('element',ATOMS)
def test_all_elements_factors_and_arbitrary_pairs(source,element):
    domain=pair_domain((element,'H'),source)
    lo,hi=max(domain['usable_bands_keV'],key=lambda interval:interval[1]-interval[0])
    energy=lo+.61*(hi-lo)
    q=np.array([0.,1.12345,4.])
    f0,fp,fpp=factor_components(element,q,energy,source,mode='Extended')
    f=atomic_factor(element,q,energy,source,mode='Extended')
    np.testing.assert_allclose(f0,xraydb.f0(element,q/(4*np.pi)),rtol=1e-14)
    np.testing.assert_allclose(f,f0+fp+1j*fpp,rtol=1e-14)
    assert np.isfinite(f).all() and np.all(fpp>=0)
    w=pair_weight((element,element),q,energy,source,mode='Extended')
    np.testing.assert_allclose(w,f.real**2+f.imag**2,atol=1e-10)
    mixed=pair_weight((element,'H'),q,energy,source,mode='Extended')
    np.testing.assert_array_equal(mixed,pair_weight(('H',element),q,energy,source,mode='Extended'))
    fh=atomic_factor('H',q,energy,source,mode='Extended')
    np.testing.assert_allclose(mixed,f.real*fh.real+f.imag*fh.imag,rtol=1e-14)
    if source in ('Henke','NIST'):
        a,m=full_table(source,element)
        # Verify full f1 -> correction f′ at original valid knots, including endpoints.
        indices=np.unique([0,len(a)//2,len(a)-1])
        node_fp,node_fpp=anomalous(element,a[indices,0]/1000,source)
        np.testing.assert_allclose(node_fp,a[indices,1]-xraydb.atomic_number(element),atol=2e-10,rtol=1e-12)
        np.testing.assert_allclose(node_fpp,a[indices,2],atol=2e-10,rtol=1e-12)

@pytest.mark.parametrize('source',['Henke','NIST'])
@pytest.mark.parametrize('element',ATOMS)
def test_each_table_interpolation_and_forbidden_brackets(source,element):
    a,m=full_table(source,element)
    # Independent reconstruction of one interpolation interval away from all guards.
    for i in range(len(a)-2,0,-1):
        mid=(a[i,0]+a[i+1,0])/2
        if any(lo<mid<hi for lo,hi,_ in m['forbidden_intervals_eV']): continue
        fp,fpp=anomalous(element,mid/1000,source)
        assert fp==pytest.approx((a[i,1]+a[i+1,1])/2-m['Z'],abs=1e-10)
        if a[i,2]>0 and a[i+1,2]>0:
            t=np.log(mid/a[i,0])/np.log(a[i+1,0]/a[i,0])
            expected=np.exp(np.log(a[i,2])+t*np.log(a[i+1,2]/a[i,2]))
            assert fpp==pytest.approx(expected,rel=1e-12,abs=1e-13)
        else: assert fpp==0
        break
    else: pytest.fail('No usable interpolation interval')
    for lo,hi,_ in m['forbidden_intervals_eV']:
        # Raw exact knots inside a conservative wide bracket remain readable;
        # test an off-knot point to ensure genuine rejection.
        energy=lo+.4321987*(hi-lo)
        if np.min(abs(a[:,0]-energy))<1e-8: continue
        with pytest.raises(CalculationError,match='gap'):
            anomalous(element,energy/1000,source)

@pytest.mark.parametrize('source',SOURCES)
@pytest.mark.parametrize('pair',[('H','He'),('Fe','Au'),('U','O'),('Br','C')])
def test_domain_intersection_and_no_extrapolation(source,pair):
    d=pair_domain(pair,source);a=element_domain(source,pair[0]);b=element_domain(source,pair[1])
    assert d['min_keV']==max(a['min_keV'],b['min_keV'])
    assert d['max_keV']==min(a['max_keV'],b['max_keV'])
    for e in [d['min_keV']*.99,d['max_keV']*1.01]:
        with pytest.raises(CalculationError,match='energy must lie'):
            pair_weight(pair,4,e,source,mode='Extended')

@pytest.mark.parametrize('source',SOURCES)
def test_fe_k_edge_and_br_k_edge(source):
    for pair,near,far in [('Fe-H',7.2,7.),('Br-Br',13.5,13.4)]:
        with pytest.raises(CalculationError,match='edge|boundary|gap'):
            delta_pair_weight(pair,4,near,far,source,mode='Extended')
        with pytest.raises(CalculationError,match='edge|boundary|gap'):
            guard_source_interval(pair,[far,near],source)

@pytest.mark.parametrize('source',['Henke','NIST'])
def test_both_atoms_use_selected_source_without_native_substitution(source,monkeypatch):
    import atrxl.sources as sources
    def fail(*args,**kwargs): raise AssertionError('Unexpected native anomalous fallback')
    monkeypatch.setattr(xraydb,'f1_chantler',fail);monkeypatch.setattr(xraydb,'f2_chantler',fail)
    factors=[]
    for element in ('Br','C'):
        a,m=full_table(source,element);e=12000.
        expected_fp=np.interp(e,a[:,0],a[:,1])-m['Z']
        expected_fpp=np.exp(np.interp(np.log(e),np.log(a[:,0]),np.log(np.maximum(a[:,2],1e-300))))
        f=atomic_factor(element,4,12,source)
        assert f.real==pytest.approx(xraydb.f0(element,4/(4*np.pi))[0]+expected_fp,abs=1e-12)
        assert f.imag==pytest.approx(expected_fpp,abs=1e-12)
        factors.append(f)
    assert pair_weight('Br-C',4,12,source)==pytest.approx((factors[0]*np.conj(factors[1])).real)

@pytest.mark.parametrize('element',['Xx','Br1-','Cval','Cf','Fe2+'])
@pytest.mark.parametrize('source',SOURCES)
def test_unsupported_species(element,source):
    with pytest.raises(CalculationError): atomic_factor(element,4,12,source)

@pytest.mark.parametrize('source',SOURCES)
def test_mode_specific_results_and_csv_json(source):
    w,wm=calculate_w('Fe-Au',[.2,4,6],12,source)
    assert {'W_e2','energy_keV'}<=set(w)
    assert not any('far' in key or 'near' in key or 'delta' in key for key in w)
    assert 'energy_near_keV' not in wm and wm['calculation_mode']=='W'
    d,dm=calculate('Br-Br',[.2,4,6],13.45,12,source)
    assert {'W_far_e2','W_near_e2','delta_W_e2'}<=set(d)
    assert dm['calculation_mode']=='ΔW'
    for frame,meta in [(w,wm),(d,dm)]:
        exported=pd.read_csv(BytesIO(csv_bytes(frame)),float_precision='round_trip')
        for column in frame.select_dtypes('number'): np.testing.assert_array_equal(exported[column],frame[column])
        payload=json.loads(results_json(frame,meta))
        assert payload['metadata']['calculation_id']==frame.calculation_id.iloc[0]
        assert payload['results'][1][('W_e2' if frame is w else 'delta_W_e2')]==frame.iloc[1][('W_e2' if frame is w else 'delta_W_e2')]

@pytest.mark.parametrize('source',SOURCES)
def test_distinct_energy_scan_semantics(source):
    e=np.array([12.,12.5,13.,13.45])
    w,wm=energy_scan('Br-Br',4,e,source,calculation_mode='W')
    d,dm=energy_scan('Br-Br',4,e,source,calculation_mode='ΔW',e_far=12.)
    np.testing.assert_allclose(w.W_e2,pair_weight('Br-Br',4,e,source))
    np.testing.assert_allclose(d.delta_W_e2,w.W_e2-pair_weight('Br-Br',4,12,source),atol=2e-12)
    assert d.delta_W_e2.iloc[0]==0
    assert dm['energy_far_keV']==12 and dm['energy_near_keV']==e.tolist()
    assert 'energy_far_keV' not in wm
    for frame,meta in [(w,wm),(d,dm)]: assert json.loads(results_json(frame,meta))

def test_conflicting_duplicate_knots_are_unavailable():
    for element,energy in [('Mg',.0103),('Pt',.006199),('Er',1.399)]:
        with pytest.raises(CalculationError): atomic_factor(element,4,energy,'Henke',mode='Extended')

def test_full_metadata_and_integrity(tmp_path,monkeypatch):
    import atrxl.sources as s
    m=source_metadata('NIST',('Fe','Au'))
    assert set(m['tables'])=={'Fe','Au'}
    for element in m['tables']:
        row=m['tables'][element]
        assert row['Z']==xraydb.atomic_number(element) and row['sha256'] and row['url']
    path=tmp_path/'full';path.mkdir()
    (path/'manifest.json').write_text(json.dumps(full_manifest()))
    (path/'Henke').mkdir();(path/'Henke/Fe.nff').write_text('corrupt')
    monkeypatch.setattr(s,'ROOT',tmp_path);s.full_manifest.cache_clear();s.full_table.cache_clear()
    try:
        with pytest.raises(CalculationError,match='SHA-256'):
            pair_weight('Fe-Au',4,12,'Henke')
    finally: s.full_manifest.cache_clear();s.full_table.cache_clear()

@pytest.mark.parametrize('source,expected',[('XrayDB',-155.27584665917914),('Henke',-159.315760),('NIST',-164.812430)])
def test_full_source_primary_and_brbr_benchmarks(source,expected):
    assert delta_pair_weight('Br-Br',4,13.45,12,source)==pytest.approx(expected,abs=5e-7)

def test_nist_edge_metadata_and_all_raw_checksums():
    import atrxl.sources as s
    manifest=full_manifest();shellpath=s.ROOT/'full/NIST/shelldata.csv'
    assert hashlib.sha256(shellpath.read_bytes()).hexdigest()==manifest['NIST_shelldata_sha256']
    shell=pd.read_csv(shellpath)
    for source in ('Henke','NIST'):
        for element,m in manifest['sources'][source].items():
            for base in [s.ROOT/'full',Path(__file__).resolve().parents[1]/'data/sources/full']:
                assert hashlib.sha256((base/m['file']).read_bytes()).hexdigest()==m['sha256']
            if source=='NIST':
                for label,e in m['edges_eV'].items(): assert shell.iloc[m['Z']-1][label]==e

def test_availability_reports_missing_factor_without_fallback(tmp_path,monkeypatch):
    import atrxl.sources as s
    manifest=json.loads(json.dumps(full_manifest()))
    manifest['sources']['Henke']={'Fe':manifest['sources']['Henke']['Fe']}
    (tmp_path/'full/Henke').mkdir(parents=True)
    (tmp_path/'full/manifest.json').write_text(json.dumps(manifest))
    monkeypatch.setattr(s,'ROOT',tmp_path);s.full_manifest.cache_clear();s.full_table.cache_clear()
    try:
        available=s.source_availability('Henke')
        assert available['supported']==[] and 'Fe' in available['unavailable']
        with pytest.raises(CalculationError,match='unavailable'): atomic_factor('Fe',4,12,'Henke')
    finally: s.full_manifest.cache_clear();s.full_table.cache_clear()
