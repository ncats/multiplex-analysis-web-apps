import streamlit as st
from st_pages import show_pages_from_config, add_indentation
from streamlit_extras.app_logo import add_logo

# Import relevant libraries
import nidap_dashboard_lib as ndl   # Useful functions for dashboards connected to NIDAP
import basic_phenotyper_lib as bpl  # Useful functions for phenotyping collections of cells
import app_top_of_page as top

def reset_phenotype_selection():
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

    # Remove key values from session_state that should not persist
    for key, val in st.session_state.items():
        if (not key.endswith('__do_not_persist')) and (not key.startswith('FormSubmitter:')):
            st.session_state[key] = val

    # Apply pages order and indentation
    add_indentation()
    show_pages_from_config()

    # Sidebar organization
    with st.sidebar:
        st.write('**:book: [Documentation](https://ncats.github.io/multiplex-analysis-web-apps)**')

    # Add logo to page
    add_logo('app_images/mawa_logo-width315.png', height=150)

    # Run Top of Page (TOP) functions
    st.session_state = top.check_for_platform(st.session_state)

    if 'init' not in st.session_state:
        settings_yaml_file = 'config_files/OMAL_REEC.yml'
        # Initialize session_state values for streamlit processing
        st.session_state = ndl.init_session_state(st.session_state, settings_yaml_file)

    st.header('UMAP Differences Analyzer\nNCATS-NCI-DMAP')

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
    else:
        st.warning('No spatial UMAP analysis detected. Please complete Neighborhood Profiles')

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