'''
This is script which creates the Clusters Analyzer page (MAWA).
'''

import streamlit as st
import nidap_dashboard_lib as ndl   # Useful functions for NIDAP dashboards

def reset_phenotype_selection():
    '''
    Quick function to reset the selections
    '''
    st.session_state.inciPhenoSel = st.session_state.def_lineage_opt

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

    # Make the figures for the Heatmap and Incidence Plot
    if st.session_state.umap_completed:
        st.session_state = ndl.set_figure_objs_clusters_analyzer(st.session_state)
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
        st.toggle('Subset Heatmap by Feature', value = False, key = 'toggle_heatmap_filter_feat')
        if st.session_state['toggle_heatmap_filter_feat']:
            nei_feat_filt_col = st.columns([2,2])
            with nei_feat_filt_col[0]:
                st.selectbox('Feature', options = st.session_state.umapOutcomes, key='heatmap_filter_feat')
            with nei_feat_filt_col[1]:
                selected_feat = st.session_state['heatmap_filter_feat']
                unique_values = st.session_state.spatial_umap.df_umap[selected_feat].unique()
                st.selectbox('Value', options = unique_values, key='heatmap_filter_value')
        if st.session_state.umap_completed:
            st.pyplot(st.session_state.heatmapfig)

    ### INCIDENCE PLOT ###
    with clusterfigs[1]:
        st.header('Incidence Lineplot')

        inci_sel_col = st.columns(2)
        # Feature Select Box
        with inci_sel_col[0]:
            st.selectbox('Feature', options = st.session_state.inciOutcomes, key = 'inciOutcomeSel')
        # Phenotype Select Box
        with inci_sel_col[1]:
            st.selectbox(st.session_state.lineageDisplayToggle_clus,
                         options = st.session_state.cluslineages, key = 'inciPhenoSel')

        if st.session_state.inciOutcomeSel == st.session_state.def_inci_feature:
            inci_radio_disabled = True
        else:
            inci_radio_disabled = False
        st.radio('Display As:', options = ('Count Differences', 'Percentages', 'Ratios'),
                 key = 'Inci_Value_display', horizontal=True, disabled = inci_radio_disabled)

        if st.session_state.umap_completed:
            st.pyplot(st.session_state.inci_fig)

if __name__ == '__main__':
    main()
