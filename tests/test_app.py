"""The four v1 GUI behaviors retained under the deliberate v1.1 control layout."""
from pathlib import Path
from streamlit.testing.v1 import AppTest
APP=Path(__file__).resolve().parents[1]/'app.py'
def test_calculator_and_invalid_inputs():
    at=AppTest.from_file(str(APP),default_timeout=90).run()
    assert not at.exception
    at.button(key='calculate').click().run()
    assert not at.exception and not at.error
    assert '-155.2758466592' in at.metric[2].value
    at.number_input(key='q').set_value(-1)
    at.button(key='calculate').click().run()
    assert at.error and not at.exception

def test_navigation_and_comparison():
    at=AppTest.from_file(str(APP),default_timeout=90).run()
    at.sidebar.radio(key='page').set_value('Source Comparison').run()
    at.checkbox(key='enable_energy').set_value(True).run()
    at.checkbox(key='enable_q').set_value(True).run()
    at.button(key='calculate').click().run()
    assert not at.exception and not at.error
    assert len(at.get('plotly_chart'))==3  # ΔW(E), ΔW(q), selected-source ΔW(q)
    at.sidebar.radio(key='page').set_value('About / Scientific Notes').run()
    assert not at.exception
    at.sidebar.radio(key='page').set_value('Home').run()
    assert not at.exception

def test_extended_warning_and_source_domain_failure():
    at=AppTest.from_file(str(APP),default_timeout=90).run()
    at.selectbox(key='domain_mode').set_value('Extended')
    at.number_input(key='q').set_value(7.)
    at.button(key='calculate').click().run()
    assert not at.error and not at.exception
    assert any('Outside Phase 1 validated q-range' in w.value for w in at.warning)
    at.selectbox(key='data_source').set_value('NIST').run()
    at.number_input(key='near').set_value(1000.)
    at.button(key='calculate').click().run()
    assert at.error and not at.exception

def test_single_point_w_export_controls():
    at=AppTest.from_file(str(APP),default_timeout=90).run()
    at.radio(key='calculation_mode').set_value('W').run()
    at.button(key='calculate').click().run()
    assert not at.error and not at.exception
    assert len(at.get('download_button'))==3  # only single-point CSV/JSON/metadata, no disabled scans
    assert len(at.metric)==1
