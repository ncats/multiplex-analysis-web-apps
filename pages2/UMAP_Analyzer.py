'''
This is script which creates the UMAP Differences Analyzer page (MAWA).
'''
import streamlit as st
from streamlit_extras.add_vertical_space import add_vertical_space
import nidap_dashboard_lib as ndl   # Useful functions for dashboards connected to NIDAP

def reset_phenotype_selection():
    '''
    Quick function to reset the UMAP differences visualizations
    '''
    st.session_state.umapInspect_Ver = st.session_state.def_lineage_opt
    st.session_state.diffUMAPSel_Ver = st.session_state.def_lineage_opt

def main():
    '''
    Main function for running the page
    '''

    # Make a generic check to avoid neeeding to hunt down individual checks
    rdy_to_plot = st.session_state.cluster_completed

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

    if rdy_to_plot:
        st.session_state = ndl.setFigureObjs_UMAPDifferences(st.session_state)
    else:
        st.warning('No spatial UMAP analysis detected. Please complete Neighborhood Profiles')

    # Full UMAP Settings
    umap_sett_cols = st.columns(2)

    with umap_sett_cols[0]:
        st.header('Full Spatial UMAP')
    with umap_sett_cols[1]:
        umap_insp_col = st.columns(2)

        with umap_insp_col[0]:
            st.selectbox('Feature',
                         options = st.session_state.umapOutcomes, key = 'umapInspect_Feat')
        with umap_insp_col[1]:
            st.selectbox(st.session_state.lineageDisplayToggle,
                         options = st.session_state.umaplineages, key = 'umapInspect_Ver')

        if st.session_state.umap_ins_msg is not None:
            st.error(st.session_state.umap_ins_msg)
        else:
            add_vertical_space(2)

    # Large UMAP Columns
    umap_viz = st.columns(2)

    # FULL UMAP
    with umap_viz[0]:
        if rdy_to_plot:
            st.pyplot(st.session_state.UMAPFig)

    # Inspection UMAP
    with umap_viz[1]:
        if rdy_to_plot:
            st.pyplot(st.session_state.UMAPFigInsp)

    # Difference Measures
    st.header('Difference Measures')

    umap_diff_sett_cols = st.columns([2, 1])
    with umap_diff_sett_cols[0]:
        umap_diff_sett_subcols = st.columns(2)
        with umap_diff_sett_subcols[0]:
            st.selectbox('Feature',
                         options = st.session_state.umapOutcomes, key = 'diffUMAPSel_Feat')
        with umap_diff_sett_subcols[1]:
            st.selectbox(st.session_state.lineageDisplayToggle,
                         options = st.session_state.umaplineages, key = 'diffUMAPSel_Ver')

        if st.session_state.umap_diff_msg is not None:
            st.error(st.session_state.umap_diff_msg)
        else:
            add_vertical_space(2)

    diff_umap_col = st.columns(3)
    with diff_umap_col[0]:
        st.header('UMAP A')
        if rdy_to_plot:
            st.pyplot(st.session_state.UMAPFigDiff0_Dens)
            st.pyplot(st.session_state.UMAPFigDiff0_Clus)
    with diff_umap_col[1]:
        st.header('UMAP B')
        if rdy_to_plot:
            st.pyplot(st.session_state.UMAPFigDiff1_Dens)
            st.pyplot(st.session_state.UMAPFigDiff1_Clus)
    with diff_umap_col[2]:
        st.write('#')
        st.write('###')
        st.write('###')
        st.header('UMAP A - UMAP B')
        if rdy_to_plot:
            st.pyplot(st.session_state.UMAPFigDiff2_Dens)

if __name__ == '__main__':
    main()
