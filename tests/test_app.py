from pathlib import Path
from streamlit.testing.v1 import AppTest
APP=Path(__file__).resolve().parents[1]/'app.py'
def test_calculator_and_invalid_inputs():
    at=AppTest.from_file(str(APP),default_timeout=90).run()
    assert not at.exception
    at.button[0].click().run()
    assert not at.exception and not at.error
    assert '-155.2758466592' in at.metric[2].value
    at.number_input[0].set_value(-1)
    at.button[0].click().run()
    assert at.error and not at.exception

def test_navigation_and_comparison():
    at=AppTest.from_file(str(APP),default_timeout=90).run()
    at.sidebar.radio[0].set_value('Source Comparison').run()
    at.button[0].click().run()
    assert not at.exception and not at.error
    assert len(at.get('plotly_chart'))==4
    at.sidebar.radio[0].set_value('About / Scientific Notes').run()
    assert not at.exception
    at.sidebar.radio[0].set_value('Home').run()
    assert not at.exception

def test_extended_warning_and_source_domain_failure():
    at=AppTest.from_file(str(APP),default_timeout=90).run()
    at.selectbox[2].set_value('Extended')
    at.number_input[0].set_value(7.)
    at.button[0].click().run()
    assert not at.error and not at.exception
    assert any('Outside Phase 1 validated q-range' in w.value for w in at.warning)
    at.selectbox[1].set_value('NIST')
    at.number_input[2].set_value(13.47)
    at.button[0].click().run()
    assert at.error and not at.exception

def test_single_point_w_export_controls():
    at=AppTest.from_file(str(APP),default_timeout=90).run()
    at.radio[0].set_value('W')
    at.checkbox[0].set_value(False)
    at.button[0].click().run()
    assert not at.error and not at.exception
    assert len(at.get('download_button'))>=10
