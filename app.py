"""Local app. Controls rerun immediately; numerical work only starts on Calculate."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
import numpy as np
import pandas as pd
import streamlit as st
from atrxl.validation import CalculationError,q_grid,numeric,validate
from atrxl.sources import SOURCES,LEGACY_SOURCES,supported_elements,source_availability,pair_domain,source_metadata
from atrxl.export import calculate,calculate_w,energy_scan,csv_bytes,json_bytes,results_json
from atrxl.plotting import interactive,png

st.set_page_config(page_title='A-TRXL Calculator v1.1',page_icon='⚛',layout='wide')
st.title('A-TRXL Scattering Factor & Contrast Calculator v1.1')
st.caption('Neutral atoms • Physical q (Å⁻¹) • Photon energy (keV) • Pair coefficients (electrons²)')
page=st.sidebar.radio('Workspace',['Home','Calculator','Source Comparison','About / Scientific Notes'],index=1,key='page')
LABELS={'XrayDB':'XrayDB/Chantler 4.5.8','Henke':'CXRO/Henke — both elements','NIST':'NIST/Chantler archived — both elements',
        'Henke_Br_swap':'Phase 1 legacy: Henke Br-only swap','NIST_Br_swap':'Phase 1 legacy: NIST Br-only swap'}

if page=='Home':
    st.header('Atomic pair weights and anomalous energy contrast')
    st.latex(r'f_a=f_{0,a}(q)+f_a^{\prime}(E)+if_a^{\prime\prime}(E),\quad W_{ab}=\mathrm{Re}(f_af_b^*),\quad \Delta W=W_{near}-W_{far}')
    st.write('Choose two independently supported neutral elements. W uses one photon energy. ΔW uses near-minus-far energies. Energy and q scans are independent and optional.')
    st.info('Primary Phase 1 benchmark: Br–Br, q=4 Å⁻¹, 13.45−12 keV, ΔW=−155.27584665917914 e².')
elif page=='About / Scientific Notes':
    st.warning('This calculator evaluates atomic pair-weight coefficients, not molecular scattering intensity, pump–probe double-difference signals, or experimental SNR.')
    st.write('Common f0: XrayDB/Waasmaier–Kirfel neutral-atom parameterization evaluated at physical q/(4π). Both selected atoms use the selected source anomalous factors. Henke/NIST full f1 is converted to f′ by subtracting each atomic Z. Imaginary terms and complex conjugation are retained; no Debye factor 2 or structural sinc kernel is included.')
    st.write('Default numerical input domain: q=0.2–6 Å⁻¹ (plus q=0 diagnostic), E=11.83–13.46 keV. Phase 1 research validation covers Br/C/O/H only. Expanded species have source-level tests, not Phase 1 physical validation. Near-edge neutral-atom predictions do not model chemical state or experimental uncertainty.')
    st.write('Full-source names Henke/NIST now mean both elements. Explicit Henke_Br_swap/NIST_Br_swap APIs preserve the historical Br-only comparisons. The Source Comparison workspace can select these labeled legacy channels.')
    st.write('Native XrayDB 4.5.8 f1 depends on spline extent. Historical Br/C/O/H within 11.83–13.46 keV retain Phase 1 endpoint anchors. Other energies/species use scalar native API. See validation report for numerical boundaries and provenance.')
    for source in SOURCES:
        with st.expander(LABELS[source]+' supported species'):
            try: st.write(', '.join(supported_elements(source)))
            except CalculationError as exc: st.error(str(exc))
    st.caption('Original Phase 1 reference files are read-only. All full source tables are stored locally with URL/version/retrieval timestamp/SHA-256. No table download occurs during a calculation.')
else:
    st.header('1 — Single-point Calculation')
    a,b,c,d=st.columns(4)
    source=c.selectbox('Data source',SOURCES,format_func=LABELS.get,key='data_source')
    try:
        availability=source_availability(source)
        elements=tuple(availability['supported'])
        if availability['unavailable']:
            st.warning('Unavailable factors: '+str(availability['unavailable']))
        if not elements: raise CalculationError('No verified source factors available locally')
    except CalculationError as exc:
        st.error(str(exc));st.stop()
    for key in ('element_a','element_b'):
        if st.session_state.get(key) not in elements:
            st.session_state[key]='Br' if 'Br' in elements else elements[0]
    atom_a=a.selectbox('Element A',elements,key='element_a')
    atom_b=b.selectbox('Element B',elements,key='element_b')
    mode=d.selectbox('Domain mode',['Default','Extended'],key='domain_mode')
    pair=(atom_a,atom_b)
    st.caption(f'{LABELS[source]}: {len(elements)} verified neutral species. f0 is common native XrayDB/Waasmaier–Kirfel for both atoms.')
    output=st.radio('Calculation mode',['ΔW','W'],horizontal=True,key='calculation_mode')
    domain=None
    try:
        domain=pair_domain(pair,source)
        if mode=='Extended':
            st.info(f"Actual common energy-domain envelope: {domain['min_keV']:.10g}–{domain['max_keV']:.10g} keV ({atom_a}/{atom_b}, {source})")
            st.warning('Outside Phase 1 input domain is numerical exploration. Expanded species are not Phase 1 validated. '+domain['warning'])
            with st.expander('Element domains, usable energy bands, edges and interpolation boundaries'):
                st.json(domain)
    except CalculationError as exc: st.error(str(exc))
    a,b,c=st.columns(3)
    q=a.number_input('Physical q (Å⁻¹)',value=4.,format='%.6f',key='q')
    if output=='W':
        energy=b.number_input('E (keV)',value=13.45,format='%.6f',key='energy')
        near=far=None
    else:
        far=b.number_input('Efar (keV)',value=12.,format='%.6f',key='far')
        near=c.number_input('Enear (keV)',value=13.45,format='%.6f',key='near')
        energy=None
    st.divider()
    st.header('2 — Energy Scan')
    enable_energy=st.checkbox('Enable Energy Scan',value=False,key='enable_energy')
    a,b=st.columns(2)
    emin=a.number_input('Energy min (keV)',value=11.83,format='%.6f',disabled=not enable_energy,key='emin')
    emax=b.number_input('Energy max (keV)',value=13.46,format='%.6f',disabled=not enable_energy,key='emax')
    st.caption('W: W(E). ΔW: ΔW(Enear; fixed Efar). The scan must remain within a single supported energy band.')
    st.divider()
    st.header('3 — q Scan')
    enable_q=st.checkbox('Enable q Scan',value=False,key='enable_q')
    a,b,c=st.columns(3)
    qlo=a.number_input('qmin (Å⁻¹)',value=.2,disabled=not enable_q,key='qlo')
    qhi=b.number_input('qmax (Å⁻¹)',value=6.,disabled=not enable_q,key='qhi')
    qstep=c.number_input('q step (Å⁻¹)',value=.02,format='%.6f',disabled=not enable_q,key='qstep')
    selected=[]
    if page=='Source Comparison':
        selected=st.multiselect('Comparison sources',list(SOURCES)+list(LEGACY_SOURCES[1:]),default=list(SOURCES),format_func=LABELS.get,key='comparison_sources')
        st.caption('Full Henke/NIST use both source atoms. Labeled legacy channels replace Br only and preserve the original Phase 1 comparison.')
    signature=(page,source,pair,mode,output,q,energy,near,far,enable_energy,emin,emax,enable_q,qlo,qhi,qstep,tuple(selected))
    clicked=st.button('Calculate',type='primary',key='calculate')
    if clicked:
        st.session_state.pop('result',None)
        try:
            if domain is None: raise CalculationError('Selected pair/source domain unavailable')
            def run(values):
                if output=='W': return calculate_w(pair,values,energy,source,mode)
                return calculate(pair,values,near,far,source,mode)
            point,meta=run(q)
            result={'point':(point,meta),'signature':signature,'output':output,'source':source,'pair':pair}
            if enable_q:
                grid=q_grid(qlo,qhi,qstep)
                result['q']=run(grid)
                result['q'][1]['curve_kind']='q_scan'
            if enable_energy:
                low,high=float(numeric(emin,'Energy min')),float(numeric(emax,'Energy max'))
                if high<=low: raise CalculationError('Energy max must exceed Energy min')
                es=np.linspace(low,high,160)
                result['energy']=energy_scan(pair,q,es,source,mode,calculation_mode=output,e_far=far)
            if selected:
                values=result['q'][0].q_Ainv.to_numpy() if 'q' in result else np.array([q])
                comparisons=[]
                for other in selected:
                    try:
                        if output=='W': frame,m=calculate_w(pair,values,energy,other,mode)
                        else: frame,m=calculate(pair,values,near,far,other,mode)
                        frame['source']=other;comparisons.append((frame,m))
                    except CalculationError as exc: st.warning(f'{LABELS[other]} unavailable for this request: {exc}')
                if not comparisons: raise CalculationError('No comparison source supports these inputs')
                result['comparisons']=comparisons
            elif page=='Source Comparison': raise CalculationError('Select at least one comparison source')
            st.session_state['result']=result
        except (CalculationError,ValueError,OverflowError) as exc: st.error(str(exc))
    result=st.session_state.get('result')
    if result and result['signature']!=signature:
        st.info('Inputs changed. Click Calculate to generate results for these settings.')
    if result and result['signature']==signature:
        point,meta=result['point'];output=result['output']
        st.divider();st.subheader('Single-point result · '+output)
        for warning in meta['warnings']: st.warning(warning)
        if output=='W': st.metric('W(q,E) (e²)',f'{point.W_e2.iloc[0]:.10f}')
        else:
            for col,key,label in zip(st.columns(3),['W_far_e2','W_near_e2','delta_W_e2'],['W far','W near','ΔW near − far']):
                col.metric(label+' (e²)',f'{point[key].iloc[0]:.10f}')
        st.dataframe(point,hide_index=True)
        def downloads(frame,metadata,prefix):
            st.download_button(prefix+' CSV',csv_bytes(frame),prefix+'.csv','text/csv',key=prefix+'_csv')
            st.download_button(prefix+' JSON',results_json(frame,metadata),prefix+'.json','application/json',key=prefix+'_json')
            st.download_button(prefix+' metadata',json_bytes(metadata),prefix+'_metadata.json','application/json',key=prefix+'_metadata')
        downloads(point,meta,'single_point')
        def plot(x,series,xlabel,title,key):
            st.plotly_chart(interactive(x,series,xlabel,title),width='stretch',key=key)
            st.download_button(title+' PNG',png(x,series,xlabel,title),key+'.png','image/png',key=key+'_png')
        ykey='W_e2' if output=='W' else 'delta_W_e2'
        for key,title in [('energy','Energy Scan'),('q','q Scan')]:
            if key in result:
                frame,metadata=result[key];st.subheader(title+' · '+output)
                x=frame['energy_keV' if output=='W' else 'energy_near_keV'] if key=='energy' else frame.q_Ainv
                label=output+(' (Enear − fixed Efar)' if key=='energy' and output=='ΔW' else '')
                plot(x,{label:frame[ykey]},'Photon energy (keV)' if key=='energy' else 'Physical q (Å⁻¹)',f"{title}: {output} · {meta['pair']} · {result['source']}",key+'_plot')
                downloads(frame,metadata,key+'_scan')
        if 'comparisons' in result:
            st.subheader('Source Comparison · '+output)
            comparisons=result['comparisons'];frame=pd.concat([f for f,_ in comparisons],ignore_index=True)
            metadata={'calculation_mode':output,'calculations':[m for _,m in comparisons]}
            st.dataframe(frame,hide_index=True)
            if 'q' in result:
                plot(comparisons[0][0].q_Ainv,{LABELS[m['source']['source']]:f[ykey] for f,m in comparisons},'Physical q (Å⁻¹)','Source Comparison · '+output,'comparison_plot')
            downloads(frame,metadata,'source_comparison')
        with st.expander('Full-precision calculation metadata'): st.json(meta)
