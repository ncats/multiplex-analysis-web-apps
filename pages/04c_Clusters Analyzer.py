import streamlit as st

# Import relevant libraries
import nidap_dashboard_lib as ndl   # Useful functions for dashboards connected to NIDAP
import basic_phenotyper_lib as bpl  # Useful functions for phenotyping collections of cells
import app_top_of_page as top
import streamlit_dataframe_editor as sde

def reset_phenotype_selection():
    st.session_state.inciPhenoSel = st.session_state.defLineageOpt

def main():
    '''
    Main function for running the page
    '''

    st.radio("Filter by Phenotypes or Markers?",
             ['Phenotypes', 'Markers'],
             key = 'lineageDisplayToggle_clus',
             horizontal = True,
             on_change = reset_phenotype_selection)

    if st.session_state.lineageDisplayToggle_clus == 'Phenotypes':
        st.session_state.cluslineages = st.session_state.umapPheno
    elif st.session_state.lineageDisplayToggle_clus == 'Markers':
        st.session_state.cluslineages = st.session_state.umapMarks

    if st.session_state.umap_completed:
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
        if st.session_state.umap_completed:
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

        if st.session_state.umap_completed:
            st.pyplot(st.session_state.inciFig)

if __name__ == '__main__':

    # Set a wide layout
    st.set_page_config(page_title="Clusters Analyzer",
                       layout="wide"
    )
    st.title('Clusters Analyzer')

    # Run streamlit-dataframe-editor library initialization tasks at the top of the page
    st.session_state = sde.initialize_session_state(st.session_state)

    # Run Top of Page (TOP) functions
    st.session_state = top.top_of_page_reqs(st.session_state)

    main()

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)
