import json
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pytest
from streamlit.testing.v1 import AppTest
from atrxl.sources import pair_domain
APP=Path(__file__).resolve().parents[1]/'app.py'
def app(): return AppTest.from_file(str(APP),default_timeout=90).run()

@pytest.mark.parametrize('source',['XrayDB','Henke','NIST'])
def test_dynamic_source_element_domains_before_calculate(source):
    at=app()
    at.selectbox(key='data_source').set_value(source).run()
    at.selectbox(key='domain_mode').set_value('Extended').run()
    for element in ['H','Fe','U']:
        at.selectbox(key='element_a').set_value(element).run()
        at.selectbox(key='element_b').set_value('H').run()
        expected=pair_domain((element,'H'),source)
        assert any(f"{expected['min_keV']:.10g}–{expected['max_keV']:.10g} keV" in info.value for info in at.info)
        assert 'result' not in at.session_state
        assert len(at.selectbox(key='element_a').options)==92
        assert not at.exception

@pytest.mark.parametrize('energy_on,q_on',[(False,False),(True,False),(False,True),(True,True)])
@pytest.mark.parametrize('mode',['W','ΔW'])
def test_independent_toggles_and_mode_specific_results(mode,energy_on,q_on):
    at=app();at.radio(key='calculation_mode').set_value(mode).run()
    assert ('Efar (keV)' in [x.label for x in at.number_input])==(mode=='ΔW')
    assert ('Enear (keV)' in [x.label for x in at.number_input])==(mode=='ΔW')
    assert ('E (keV)' in [x.label for x in at.number_input])==(mode=='W')
    at.checkbox(key='enable_energy').set_value(energy_on).run()
    at.checkbox(key='enable_q').set_value(q_on).run()
    assert at.number_input(key='emin').disabled==(not energy_on)
    assert at.number_input(key='qlo').disabled==(not q_on)
    at.button(key='calculate').click().run()
    assert not at.exception and not at.error
    result=at.session_state['result']
    assert ('energy' in result)==energy_on and ('q' in result)==q_on
    assert len(at.get('plotly_chart'))==int(energy_on)+int(q_on)
    point,meta=result['point']
    assert meta['calculation_mode']==mode
    if mode=='W':
        assert 'W_e2' in point and 'delta_W_e2' not in point
        assert len(at.metric)==1 and 'energy_far_keV' not in meta
    else:
        assert 'delta_W_e2' in point and len(at.metric)==3
    for key in ['q','energy']:
        if key in result:
            frame,m=result[key]
            assert m['calculation_mode']==mode
            assert ('W_e2' if mode=='W' else 'delta_W_e2') in frame
    if energy_on:
        frame,m=result['energy']
        if mode=='ΔW':
            from atrxl import delta_pair_weight
            np.testing.assert_allclose(frame.delta_W_e2,delta_pair_weight('Br-Br',4,frame.energy_near_keV.to_numpy(),12),atol=1e-12)
    # Assert actual plotted ordinates correspond to the requested mode.
    for chart in at.get('plotly_chart'):
        spec=json.loads(chart.proto.spec)
        assert all(('ΔW' in trace['name'])==(mode=='ΔW') for trace in spec['data'])

def test_disabled_scans_never_execute_and_invalid_disabled_inputs_ignored():
    at=app()
    at.checkbox(key='enable_energy').set_value(True).run()
    at.number_input(key='emax').set_value(1.)
    at.checkbox(key='enable_energy').set_value(False).run()
    at.checkbox(key='enable_q').set_value(True).run()
    at.number_input(key='qstep').set_value(-1.)
    at.checkbox(key='enable_q').set_value(False).run()
    with patch('atrxl.export.energy_scan',side_effect=AssertionError('disabled energy scan ran')),patch('atrxl.validation.q_grid',side_effect=AssertionError('disabled q scan ran')):
        at.button(key='calculate').click().run()
    assert not at.exception and not at.error
    assert not any(key in at.session_state['result'] for key in ['energy','q'])

def test_w_mode_never_calls_delta_calculation():
    at=app();at.radio(key='calculation_mode').set_value('W').run()
    at.checkbox(key='enable_energy').set_value(True).run()
    at.checkbox(key='enable_q').set_value(True).run()
    with patch('atrxl.export.delta_pair_weight',side_effect=AssertionError('W called ΔW')):
        at.button(key='calculate').click().run()
    assert not at.exception and not at.error

def test_inputs_react_immediately_but_no_auto_calculation():
    at=app();at.button(key='calculate').click().run()
    uid=at.session_state['result']['point'][1]['calculation_id']
    at.checkbox(key='enable_energy').set_value(True).run()
    assert not at.number_input(key='emin').disabled
    assert at.session_state['result']['point'][1]['calculation_id']==uid
    assert len(at.metric)==0 and any('Click Calculate' in msg.value for msg in at.info)
    at.button(key='calculate').click().run()
    assert at.session_state['result']['point'][1]['calculation_id']!=uid
    assert 'energy' in at.session_state['result']

def test_invalid_enabled_scan_rejects_and_leaves_no_partial_export():
    at=app();at.checkbox(key='enable_energy').set_value(True).run()
    at.number_input(key='emax').set_value(10.)
    at.button(key='calculate').click().run()
    assert at.error and not at.exception and len(at.get('download_button'))==0
    assert 'result' not in at.session_state

def test_full_source_and_explicit_phase1_legacy_comparison():
    at=app();at.sidebar.radio(key='page').set_value('Source Comparison').run()
    at.multiselect(key='comparison_sources').set_value(['XrayDB','Henke','Henke_Br_swap','NIST','NIST_Br_swap']).run()
    at.selectbox(key='element_b').set_value('C').run()
    at.button(key='calculate').click().run()
    assert not at.exception and not at.error
    comparisons=at.session_state['result']['comparisons']
    values={meta['source']['source']:frame.delta_W_e2.iloc[0] for frame,meta in comparisons}
    assert set(values)=={'XrayDB','Henke','Henke_Br_swap','NIST','NIST_Br_swap'}
    assert values['Henke']!=values['Henke_Br_swap']
    assert values['NIST']!=values['NIST_Br_swap']
