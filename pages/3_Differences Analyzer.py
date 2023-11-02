import streamlit as st
from streamlit_javascript import st_javascript
import time

# Import relevant libraries
import nidap_dashboard_lib as ndl   # Useful functions for dashboards connected to NIDAP
import basic_phenotyper_lib as bpl  # Useful functions for phenotyping collections of cells

def main():
    '''
    Main function for running the page
    '''
    # Use the whole page width
    st.set_page_config(
        page_title="Differences Analyzer",
        layout="wide"
    )

    for key, val in st.session_state.items():
        if (not key.endswith('__do_not_persist')) and (not key.startswith('FormSubmitter:')):
            st.session_state[key] = val

    st.header('Differences Analyzer\nNCATS-NCI-DMAP')

    # Toggles for different figures
    figToggle1, figToggle2, figToggle3 = st.columns([1, 1, 2])
    with figToggle1:
        st.radio("Show UMAP Density or Clusters?",
                        ['Density', 'Clusters'],
                        key = 'UMAPFigType',
                        horizontal = True)

    with figToggle2:
        st.radio("Filter by Phenotypes or Markers?",
                        ['Phenotypes', 'Markers'],
                        key = 'lineageDisplayToggle',
                        horizontal = True)

    if st.session_state.lineageDisplayToggle == 'Phenotypes':
        st.session_state.umaplineages = st.session_state.umapPheno
    elif st.session_state.lineageDisplayToggle == 'Markers':
        st.session_state.umaplineages = st.session_state.umapMarks

    if st.session_state.umapCompleted:
        st.session_state = ndl.setFigureObjs_UMAPDifferences(st.session_state)

    # Large UMAP Columns
    umapViz1, umapViz2 = st.columns(2)

    # FULL UMAP
    with umapViz1:
        st.header('Full Spatial UMAP')
        if st.session_state.umapCompleted:
            st.pyplot(st.session_state.UMAPFig)

    # Inspection UMAP
    with umapViz2:
        umapInsCol1, umapInsCol2 = st.columns(2)

        with umapInsCol1:
            inciOutSel   = st.selectbox('Experimental Outcomes', options = st.session_state.umapOutcomes, key = 'umapInspect_Feat')
        with umapInsCol2:
            inciPhenoSel = st.selectbox(st.session_state.lineageDisplayToggle, options = st.session_state.umaplineages, key = 'umapInspect_Ver')

        if st.session_state.umapCompleted:
            st.pyplot(st.session_state.UMAPFigInsp)

    # Difference Measures
    st.header('Difference Measures')

    diffUMAP1, diffUMAP2, diffUMAP3 = st.columns(3)

    with diffUMAP1:
        diffUMAPFeat = st.selectbox('Feature', options = st.session_state.umapOutcomes, key = 'diffUMAPSel_Feat')
        st.header('UMAP A')
        if st.session_state.umapCompleted:
            st.pyplot(st.session_state.UMAPFigDiff1)
    with diffUMAP2:
        diffUMAPVer  = st.selectbox(st.session_state.lineageDisplayToggle, options = st.session_state.umaplineages, key = 'diffUMAPSel_Ver')
        st.header('UMAP B')
        if st.session_state.umapCompleted:
            st.pyplot(st.session_state.UMAPFigDiff2)
    with diffUMAP3:
        st.write('#')
        st.write('###')
        st.write('###')
        st.header('UMAP A - UMAP B')
        if st.session_state.umapCompleted:
            st.pyplot(st.session_state.UMAPFigDiff3)

    # Clustering Columns
    HeatIncCol1, HeatIncCol2 = st.columns(2)
    # Heatmap
    with HeatIncCol1:
        st.header('Phenotype/Cluster Heatmap')
        curHeatRadio = st.radio("Normalize along features?",
                    ['No Norm', 'Norm within Clusters', 'Norm within Phenotypes'],
                    key = 'NormHeatRadio',
                    horizontal = True)
        if st.session_state.umapCompleted:
            st.pyplot(st.session_state.heatmapfig)

    # Incidence Plot        
    with HeatIncCol2:
        st.header('Incidence Lineplot')
        inciSel1, inciSel2 = st.columns(2)
        with inciSel1:
            inciOutSel   = st.selectbox('Experimental Outcomes', options = st.session_state.inciOutcomes, key = 'inciOutcomeSel')
        with inciSel2:
            inciPhenoSel = st.selectbox(st.session_state.lineageDisplayToggle, options = st.session_state.umaplineages, key = 'inciPhenoSel')

        if st.session_state.umapCompleted:
            st.pyplot(st.session_state.inciFig)


if __name__ == '__main__':
    main()
