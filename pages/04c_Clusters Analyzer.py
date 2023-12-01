import streamlit as st
from st_pages import show_pages_from_config, add_indentation
from streamlit_extras.app_logo import add_logo
import time

# Import relevant libraries
import nidap_dashboard_lib as ndl   # Useful functions for dashboards connected to NIDAP
import basic_phenotyper_lib as bpl  # Useful functions for phenotyping collections of cells

def reset_phenotype_selection():
    st.session_state.inciPhenoSel = st.session_state.defLineageOpt

def main():
    '''
    Main function for running the page
    '''
    # Use the whole page width
    st.set_page_config(
        page_title="Clusters Analyzer",
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

    if 'init' not in st.session_state:
        settings_yaml_file = 'config_files/OMAL_REEC.yml'
        # Initialize session_state values for streamlit processing
        st.session_state = ndl.init_session_state(st.session_state, settings_yaml_file)

    st.header('Clusters Analyzer\nNCATS-NCI-DMAP')

    st.radio("Filter by Phenotypes or Markers?",
            ['Phenotypes', 'Markers'],
            key = 'lineageDisplayToggle_clus',
            horizontal = True, 
            on_change = reset_phenotype_selection)
    
    if st.session_state.lineageDisplayToggle_clus == 'Phenotypes':
        st.session_state.cluslineages = st.session_state.umapPheno
    elif st.session_state.lineageDisplayToggle_clus == 'Markers':
        st.session_state.cluslineages = st.session_state.umapMarks

    if st.session_state.umapCompleted:
        st.session_state = ndl.setFigureObjs_UMAPDifferences(st.session_state)
    else:
        st.warning('No spatial UMAP analysis detected. Please complete Neighborhood Profiles')

    # Clustering Columns
    clusterfigs = st.columns(2)
    
    ### HEATMAP ###
    with clusterfigs[0]:
        st.header('Phenotype/Cluster Heatmap')
        st.radio("Normalize Heatmap",
                 options = ['No Norm', 'Norm within Clusters', 'Norm within Phenotypes'],
                 key = 'NormHeatRadio',
                 horizontal = True)
        if st.session_state.umapCompleted:
            st.pyplot(st.session_state.heatmapfig)

    ### INCIDENCE PLOT ###       
    with clusterfigs[1]:
        st.header('Incidence Lineplot')

        inciSel1, inciSel2 = st.columns(2)
        # Feature Select Box
        with inciSel1:
            st.selectbox('Feature', options = st.session_state.inciOutcomes, key = 'inciOutcomeSel')
        # Phenotype Select Box
        with inciSel2:
            st.selectbox(st.session_state.lineageDisplayToggle_clus, options = st.session_state.cluslineages, key = 'inciPhenoSel')

        if st.session_state.inciOutcomeSel == st.session_state.definciOutcomes:
            inci_radio_disabled = True
        else:
            inci_radio_disabled = False
        st.radio('Display As:', options = ('Count Differences', 'Percentages', 'Ratios'), 
                 key = 'Inci_Value_display', horizontal=True, disabled = inci_radio_disabled)

        if st.session_state.umapCompleted:
            st.pyplot(st.session_state.inciFig)

if __name__ == '__main__':
    main()
