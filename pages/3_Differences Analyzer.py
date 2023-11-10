import streamlit as st
from streamlit_javascript import st_javascript
import time

# Import relevant libraries
import nidap_dashboard_lib as ndl   # Useful functions for dashboards connected to NIDAP
import basic_phenotyper_lib as bpl  # Useful functions for phenotyping collections of cells

def reset_phenotype_selection():
    st.session_state.umapInspect_Ver = st.session_state.defLineageOpt
    st.session_state.diffUMAPSel_Ver = st.session_state.defLineageOpt

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

    ### SIDE BAR ORGANIZATION ###
    with st.sidebar:
        url = st_javascript("await fetch('').then(r => window.parent.location.href)")
        st.write(f'''[Open app in new Tab]({url})\n (MS Edge/ Google Chrome)''')

    # Toggles for different figures
    figToggle1, figToggle2, figToggle3 = st.columns([1, 1, 2])
    with figToggle2:
        st.radio("Show UMAP Density or Clusters?",
                        ['Density', 'Clusters'],
                        key = 'UMAPFigType',
                        horizontal = True)

    with figToggle1:
        st.radio("Filter by Phenotypes or Markers?",
                        ['Phenotypes', 'Markers'],
                        key = 'lineageDisplayToggle',
                        horizontal = True,
                        on_change = reset_phenotype_selection)

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
            st.selectbox('Feature', options = st.session_state.umapOutcomes, key = 'umapInspect_Feat')
        with umapInsCol2:
            st.selectbox(st.session_state.lineageDisplayToggle, options = st.session_state.umaplineages, key = 'umapInspect_Ver')

        if st.session_state.umapCompleted:
            st.pyplot(st.session_state.UMAPFigInsp)

    # Difference Measures
    st.header('Difference Measures')

    diff_umap_col = st.columns(3)

    with diff_umap_col[0]:
        st.selectbox('Feature', options = st.session_state.umapOutcomes, key = 'diffUMAPSel_Feat')
        st.header('UMAP A')
        if st.session_state.umapCompleted:
            st.pyplot(st.session_state.UMAPFigDiff0_Dens)
            st.pyplot(st.session_state.UMAPFigDiff0_Clus)
    with diff_umap_col[1]:
        st.selectbox(st.session_state.lineageDisplayToggle, options = st.session_state.umaplineages, key = 'diffUMAPSel_Ver')
        st.header('UMAP B')
        if st.session_state.umapCompleted:
            st.pyplot(st.session_state.UMAPFigDiff1_Dens)
            st.pyplot(st.session_state.UMAPFigDiff1_Clus)
    with diff_umap_col[2]:
        st.write('#')
        st.write('###')
        st.write('###')
        st.header('UMAP A - UMAP B')
        if st.session_state.umapCompleted:
            st.pyplot(st.session_state.UMAPFigDiff2_Dens)

if __name__ == '__main__':
    main()
