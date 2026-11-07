"""Optional real-browser validation. Requires playwright and an installed Chromium.
Run a local Streamlit server first. No source data are downloaded by this script.
"""
import json,os
from pathlib import Path
from playwright.sync_api import sync_playwright,expect
root=Path(__file__).resolve().parents[1]/'results/v1.1'
(root/'figures').mkdir(parents=True,exist_ok=True)
with sync_playwright() as p:
    browser=p.chromium.launch(executable_path=os.environ.get('ATRXL_CHROMIUM','/usr/bin/chromium'),headless=True)
    page=browser.new_page(viewport={'width':1600,'height':1100},accept_downloads=True)
    errors=[];page.on('pageerror',lambda error:errors.append(str(error)))
    page.goto('http://127.0.0.1:8502/',wait_until='networkidle')
    expect(page.get_by_text('1 — Single-point Calculation',exact=True)).to_be_visible()
    expect(page.get_by_role('checkbox',name='Enable Energy Scan',exact=True)).not_to_be_checked()
    expect(page.get_by_label('Energy min (keV)',exact=True)).to_be_disabled()
    page.get_by_role('button',name='Calculate',exact=True).click()
    expect(page.locator('[data-testid="stMetricValue"]').last).to_contain_text('-155.2758466592',timeout=30000)
    page.locator('[data-testid="stMain"]').evaluate('(el) => el.scrollTop = 0')
    page.screenshot(path=str(root/'figures/gui_delta.png'),full_page=True)
    page.get_by_role('radiogroup',name='Calculation mode',exact=True).get_by_text('W',exact=True).click()
    expect(page.get_by_label('E (keV)',exact=True)).to_be_visible()
    expect(page.get_by_label('Efar (keV)',exact=True)).to_have_count(0)
    page.get_by_role('button',name='Calculate',exact=True).click()
    expect(page.locator('[data-testid="stMetricValue"]')).to_have_count(1)
    with page.expect_download() as info:
        page.get_by_text('single_point JSON',exact=True).click()
    payload_path=root/'browser_W_download.json';info.value.save_as(payload_path)
    payload=json.loads(payload_path.read_text())
    assert payload['metadata']['calculation_mode']=='W'
    assert 'W_e2' in payload['results'][0] and 'delta_W_e2' not in payload['results'][0]
    page.get_by_text('Enable Energy Scan',exact=True).click()
    expect(page.get_by_label('Energy min (keV)',exact=True)).to_be_enabled()
    expect(page.get_by_label('qmin (Å⁻¹)',exact=True)).to_be_disabled()
    page.get_by_role('button',name='Calculate',exact=True).click()
    expect(page.locator('[data-testid="stPlotlyChart"]')).to_have_count(1,timeout=30000)
    page.screenshot(path=str(root/'figures/gui_W_energy.png'),full_page=True)
    assert not errors,errors
    result={'browser':browser.version,'primary_delta_benchmark_display':'-155.2758466592',
            'W_download_calculation_mode':payload['metadata']['calculation_mode'],
            'W_energy_scan_plots':page.locator('[data-testid="stPlotlyChart"]').count(),
            'q_scan_disabled':page.get_by_label('qmin (Å⁻¹)',exact=True).is_disabled(),
            'console_page_errors':errors}
    print(json.dumps(result,indent=2))
    browser.close()
