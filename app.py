"""Run with: streamlit run app.py (after pip install -e .)."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
import numpy as np
import pandas as pd
import streamlit as st
from atrxl.validation import PAIRS, CalculationError, q_grid, guard_edge
from atrxl.sources import SOURCES, source_metadata
from atrxl.weights import pair_weight, delta_pair_weight
from atrxl.export import calculate, csv_bytes, json_bytes, results_json
from atrxl.plotting import interactive, png

st.set_page_config(page_title='A-TRXL Calculator',page_icon='⚛',layout='wide')
st.title('A-TRXL Scattering Factor & Contrast Calculator v1.0')
st.caption('Neutral-atom coefficients • Physical q in Å⁻¹ • Native XrayDB 4.5.8 / Chantler')
page=st.sidebar.radio('Workspace',['Home','Calculator','Source Comparison','About / Scientific Notes'],index=1)
if page=='Home':
    st.header('Atomic contrast, with explicit scientific conventions')
    st.write('Choose Calculator for full-precision atomic factors, pair weights, q scans and energy curves. Source Comparison evaluates Br-only tabulation sensitivity.')
    st.latex(r'f_a=f_{0,a}(q)+f_a^{\prime}(E)+if_a^{\prime\prime}(E),\quad W_{ab}=\mathrm{Re}(f_af_b^*),\quad \Delta W=W_{near}-W_{far}')
    st.info('Default benchmark: Br–Br, q=4 Å⁻¹, 12.00 → 13.45 keV; ΔW = −155.27584665917914 electrons².')
elif page=='About / Scientific Notes':
    st.header('Scientific scope and provenance')
    st.warning('This calculator evaluates atomic pair-weight coefficients, not molecular scattering intensity, pump–probe double-difference signals, or experimental SNR.')
    st.write('Neutral atoms only. Each coefficient is for one unordered pair: no Debye multiplicity 2, structural sinc kernel, molecular geometry, or intensity normalization. XrayDB f0 receives q/(4π); photon energies are converted from keV to eV. Imaginary factors and complex conjugation are retained.')
    st.write('Phase 1 validated domain: q=0.2–6.0 Å⁻¹, E=11.83–13.46 keV. q=0 is a forward-limit diagnostic. Extended evaluations are numerical predictions outside that research scope, not validated physical accuracy. Br K-edge crossing is blocked. Near-edge chemical state and finite bandwidth are not modeled.')
    st.write('Henke and NIST are Br-only source swaps. Non-Br factors and all f0 values remain XrayDB. Archived adapters use linear full f1 minus Z and log-log f2, within verified pre-edge domains; no extrapolation.')
    for source in SOURCES:
        with st.expander(source+' provenance'):
            try: st.json(source_metadata(source))
            except CalculationError as exc: st.error(str(exc))
    st.write('Available Phase 1 inputs: five reports, contrast_weights.nc, calculation_metrics.json, xraydb_runtime_comparison.csv. A-TRXL Python Plan.pdf and the original sparse source excerpts were not supplied; raw tables were independently retrieved from reported URLs. See validation report for independent agreement checks.')
else:
    comparison=page=='Source Comparison'
    with st.form('inputs'):
        a,b,c=st.columns(3)
        pair=a.selectbox('Atomic pair',PAIRS)
        source=b.selectbox('Data source',SOURCES,format_func=lambda s: 'XrayDB/Chantler' if s=='XrayDB' else s+' Br-only source swap')
        mode=c.selectbox('Domain mode',['Default','Extended'])
        a,b,c=st.columns(3)
        q=a.number_input('Physical q (Å⁻¹)',value=4.0,format='%.6f')
        far=b.number_input('Efar (keV)',value=12.0,format='%.6f')
        near=c.number_input('Enear (keV)',value=13.45,format='%.6f')
        output=st.radio('Output',['ΔW','W'],horizontal=True)
        scan=st.checkbox('Scan q range',value=True)
        a,b,c=st.columns(3)
        lo=a.number_input('qmin (Å⁻¹)',value=0.2)
        hi=b.number_input('qmax (Å⁻¹)',value=6.0)
        step=c.number_input('q step (Å⁻¹)',value=0.02,format='%.6f')
        st.caption('Energy curve uses this requested interval; it must remain inside the selected source domain and on one side of Br K-edge.')
        a,b=st.columns(2)
        emin=a.number_input('Energy curve min (keV)',value=11.83,format='%.6f')
        emax=b.number_input('Energy curve max (keV)',value=13.46,format='%.6f')
        selected=st.multiselect('Comparison sources',SOURCES,default=list(SOURCES)) if comparison else [source]
        submitted=st.form_submit_button('Calculate',type='primary')
    if submitted:
        try:
            point,meta=calculate(pair,q,near,far,source,mode)
            grid=q_grid(lo,hi,step) if scan else np.array([q])
            curve,curve_meta=calculate(pair,grid,near,far,source,mode)
            if emax<=emin: raise CalculationError('Energy curve max must exceed min')
            energies=np.linspace(emin,emax,160);guard_edge(energies)
            energy_values=pair_weight(pair,q,energies,source,mode=mode)
            series={}
            comparison_metadata={}
            if comparison:
                if not selected: raise CalculationError('Select at least one comparison source')
                for s in selected:
                    try:
                        series[s+(' Br-only swap' if s!='XrayDB' else '')]=delta_pair_weight(pair,grid,near,far,s,mode=mode) if output=='ΔW' else pair_weight(pair,grid,near,s,mode=mode)
                        comparison_metadata[s]=source_metadata(s)
                    except CalculationError as exc: st.warning(str(exc))
                if not series: raise CalculationError('No comparison source supports these inputs')
            st.session_state['calculation']=(point,meta,curve,curve_meta,energies,energy_values,series,output,comparison_metadata)
        except (CalculationError,ValueError,OverflowError) as exc:
            st.session_state.pop('calculation',None);st.error(str(exc))
    if 'calculation' in st.session_state:
        point,meta,curve,curve_meta,energies,energy_values,series,output,comparison_metadata=st.session_state['calculation']
        st.subheader('Single-point results · selected output: '+output)
        st.caption('Displayed results belong to the last successful Calculate action; submit again after changing inputs.')
        for w in meta['warnings']: st.warning(w)
        cols=st.columns(3)
        for col,key,title in zip(cols,['W_far_e2','W_near_e2','delta_W_e2'],['W far','W near','ΔW near − far']):
            col.metric(title+' (e²)',f'{point[key].iloc[0]:.10f}')
        st.dataframe(point,hide_index=True)
        st.caption('f_real + i f_imag columns specify complex atomic factors. Relative change = 100 ΔW/Wfar; undefined denominators export as null/empty.')
        st.download_button('Single-point CSV',csv_bytes(point),'atrxl_point.csv','text/csv')
        st.download_button('Single-point JSON',results_json(point,meta),'atrxl_point.json','application/json')
        st.download_button('Calculation metadata',json_bytes(meta),'atrxl_metadata.json','application/json')
        grid=curve.q_Ainv.to_numpy()
        def show_plot(x,ys,xlabel,title,key):
            st.plotly_chart(interactive(x,ys,xlabel,title),width='stretch',key=key)
            st.download_button('Download '+title+' PNG',png(x,ys,xlabel,title),key+'.png','image/png',key=key+'_png')
        show_plot(grid,{f"W far ({meta['energy_far_keV']:g} keV)":curve.W_far_e2,f"W near ({meta['energy_near_keV']:g} keV)":curve.W_near_e2},'Physical q (Å⁻¹)',f"W · {meta['pair']} · {meta['source']['source']}",'weights')
        show_plot(grid,{'ΔW':curve.delta_W_e2},'Physical q (Å⁻¹)','Near − far contrast · '+meta['source']['source'],'contrast')
        st.download_button('Curve CSV',csv_bytes(curve),'atrxl_curve.csv','text/csv')
        st.download_button('Curve metadata',json_bytes(curve_meta),'atrxl_curve_metadata.json','application/json')
        show_plot(energies,{'W':energy_values},'Photon energy (keV)',f"Energy-dependent W · q={meta['q_Ainv']} Å⁻¹ · {meta['source']['source']}",'energy')
        energy_meta={**meta,'energy_curve_keV':energies.tolist(),'curve_kind':'energy-dependent W'}
        energy_frame=pd.DataFrame({'calculation_id':meta['calculation_id'],'energy_keV':energies,'W_e2':energy_values})
        st.download_button('Energy curve CSV',csv_bytes(energy_frame),'atrxl_energy.csv','text/csv')
        st.download_button('Energy metadata',json_bytes(energy_meta),'atrxl_energy_metadata.json','application/json')
        if series:
            show_plot(grid,series,'Physical q (Å⁻¹)','Source comparison · '+output,'sources')
            compare_frame=pd.DataFrame({'q_Ainv':grid,'calculation_id':curve_meta['calculation_id'],**series})
            st.download_button('Comparison CSV',csv_bytes(compare_frame),'atrxl_comparison.csv','text/csv')
            st.download_button('Comparison metadata',json_bytes({**curve_meta,'comparison_sources':comparison_metadata,'output':output}),'atrxl_comparison_metadata.json','application/json')
        with st.expander('Full reproducibility metadata'): st.json(meta)
