'''
Streamlit page for showing UMAP difference figures
'''
import streamlit as st

# Import relevant libraries
import nidap_dashboard_lib as ndl   # Useful functions for dashboards connected to NIDAP
import app_top_of_page as top
import streamlit_dataframe_editor as sde

def reset_phenotype_selection():
    '''
    Quick function to reset the UMAP differences visualizations
    '''
    st.session_state.umapInspect_Ver = st.session_state.defLineageOpt
    st.session_state.diffUMAPSel_Ver = st.session_state.defLineageOpt

def main():
    '''
    Main function for running the page
    '''
    # Use the whole page width
    st.set_page_config(
        page_title="UMAP Differences Analyzer",
        layout="wide"
    )

    # Run streamlit-dataframe-editor library initialization tasks at the top of the page
    st.session_state = sde.initialize_session_state(st.session_state)

    # Run Top of Page (TOP) functions
    st.session_state = top.top_of_page_reqs(st.session_state)

    st.header('UMAP Differences Analyzer\nNCATS-NCI-DMAP')

    # Toggles for different figures
    fig_toggle = st.columns([1, 1, 2])
    with fig_toggle[0]:
        st.radio("Show UMAP Density or Clusters?",
                 ['Density', 'Clusters'],
                 key = 'UMAPFigType',
                 horizontal = True)

    with fig_toggle[1]:
        st.radio("Filter by Phenotypes or Markers?",
                 ['Phenotypes', 'Markers'],
                 key = 'lineageDisplayToggle',
                 horizontal = True,
                 on_change = reset_phenotype_selection)

    if st.session_state.lineageDisplayToggle == 'Phenotypes':
        st.session_state.umaplineages = st.session_state.umapPheno
    elif st.session_state.lineageDisplayToggle == 'Markers':
        st.session_state.umaplineages = st.session_state.umapMarks

    if st.session_state.umap_completed:
        st.session_state = ndl.setFigureObjs_UMAPDifferences(st.session_state)
    else:
        st.warning('No spatial UMAP analysis detected. Please complete Neighborhood Profiles')

    # Large UMAP Columns
    umap_viz = st.columns(2)

    # FULL UMAP
    with umap_viz[0]:
        st.header('Full Spatial UMAP')
        if st.session_state.umap_completed:
            st.pyplot(st.session_state.UMAPFig)

    # Inspection UMAP
    with umap_viz[1]:
        umap_insp_col = st.columns(2)

        with umap_insp_col[0]:
            st.selectbox('Feature', options = st.session_state.umapOutcomes, key = 'umapInspect_Feat')
        with umap_insp_col[1]:
            st.selectbox(st.session_state.lineageDisplayToggle, options = st.session_state.umaplineages, key = 'umapInspect_Ver')

        if st.session_state.umap_completed:
            st.pyplot(st.session_state.UMAPFigInsp)

    # Difference Measures
    st.header('Difference Measures')

    diff_umap_col = st.columns(3)

    with diff_umap_col[0]:
        st.selectbox('Feature', options = st.session_state.umapOutcomes, key = 'diffUMAPSel_Feat')
        st.header('UMAP A')
        if st.session_state.umap_completed:
            st.pyplot(st.session_state.UMAPFigDiff0_Dens)
            st.pyplot(st.session_state.UMAPFigDiff0_Clus)
    with diff_umap_col[1]:
        st.selectbox(st.session_state.lineageDisplayToggle, options = st.session_state.umaplineages, key = 'diffUMAPSel_Ver')
        st.header('UMAP B')
        if st.session_state.umap_completed:
            st.pyplot(st.session_state.UMAPFigDiff1_Dens)
            st.pyplot(st.session_state.UMAPFigDiff1_Clus)
    with diff_umap_col[2]:
        st.write('#')
        st.write('###')
        st.write('###')
        st.header('UMAP A - UMAP B')
        if st.session_state.umap_completed:
            st.pyplot(st.session_state.UMAPFigDiff2_Dens)

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)

if __name__ == '__main__':
    main()
