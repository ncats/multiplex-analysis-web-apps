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

def reset_heatmap_feature_values():
    '''
    Quick callback function to reset the heatmap feature values
    '''
    selected_feat = st.session_state['heatmap_filter_feat']
    unique_values = st.session_state.spatial_umap.df_umap[selected_feat].unique()
    st.session_state.heatmap_filter_value = unique_values[0]

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

    # Clustering Columns
    clust_sett_col = st.columns(2)

    ### HEATMAP SETTINGS ###
    with clust_sett_col[0]:
        st.header('Phenotype/Cluster Heatmap')
        st.radio("Normalize Heatmap",
                 options = ['No Norm', 'Norm within Clusters', 'Norm within Phenotypes'],
                 key = 'NormHeatRadio',
                 horizontal = True)
        st.session_state.toggle_heatmap_filter_feat = False
        # st.toggle('Subset Heatmap by Feature', value = False,
        #           key = 'toggle_heatmap_filter_feat', disabled= not st.session_state.umap_completed)
        # if st.session_state['toggle_heatmap_filter_feat']:
        #     nei_feat_filt_col = st.columns([2,2])
        #     with nei_feat_filt_col[0]:
        #         st.selectbox('Feature', options = st.session_state.umapOutcomes,
        #                      key='heatmap_filter_feat', on_change=reset_heatmap_feature_values)
        #     with nei_feat_filt_col[1]:
        #         selected_feat = st.session_state['heatmap_filter_feat']
        #         unique_values = st.session_state.spatial_umap.df_umap[selected_feat].unique()
        #         st.selectbox('Value', options = unique_values, key='heatmap_filter_value')

    ### INCIDENCE FIGURE SETTINGS ###
    with clust_sett_col[1]:
        st.header('Incidence Figure')

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

        inci_sel_col = st.columns(2)
        with inci_sel_col[0]:
            st.radio('Display As:', options = ('Count Differences', 'Percentages', 'Ratios'),
                     key = 'Inci_Value_display', horizontal=True, disabled = inci_radio_disabled)
        with inci_sel_col[1]:
            st.toggle('Show Raw Counts', key='inci_fig_show_raw_counts',
                      value=False, disabled=inci_radio_disabled)

    # Make the figures for the Heatmap and Incidence Figure
    if st.session_state.umap_completed:
        st.session_state = ndl.set_figure_objs_clusters_analyzer(st.session_state)

        # Clustering Columns
        cluster_figs = st.columns(2)
        with cluster_figs[0]:
            st.pyplot(st.session_state.heatmapfig)
        with cluster_figs[1]:
            st.plotly_chart(st.session_state.inci_fig)
    else:
        st.warning('No spatial UMAP analysis detected. Please complete Neighborhood Profiles')

if __name__ == '__main__':
    main()
