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

    if 'def_lineage_opt' not in st.session_state:
        st.session_state.def_lineage_opt = 'All Phenotypes'
    if 'def_umap_feature' not in st.session_state:
        st.session_state.def_umap_feature = 'phenotype'
    if 'def_inci_feature' not in st.session_state:
        st.session_state.def_inci_feature = 'Cell Counts'

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

            st.subheader('Incidence Figure Interpretation')
            if st.session_state.inci_appro_feat:
                if st.session_state.inciOutcomeSel == st.session_state.def_inci_feature:
                    st.write('''In this default state of the Incidence Figure, we display a
                            bar chart showing the number of cells that were included in each cluster.
                            Filtering by phenotype will instead show the number of cells of a given
                            phenotype that is assigned to a given cluster.
                            ''')
                elif st.session_state.Inci_Value_display == 'Count Differences':
                    st.write(f'''The Counts Differences mode shows the difference in the number of
                            cells that were identified as {st.session_state.inci_true_msg} to
                            the number identified as {st.session_state.inci_fals_msg} in each
                            cluster. The blue line represents the difference in the number of cells
                            for a given cluster, that is the number of cells for {st.session_state.inci_true_msg},
                            minus the number of cells that were identified as {st.session_state.inci_fals_msg}.
                            By clicking on the 'Show Raw Counts' toggle, you can view the actual counts for each cluster,
                            in each condition of {st.session_state.inciOutcomeSel}.
                            \nAlthough negative count differences are shown it is of course impossible to
                            have a negative quantity of cells. Negative numbers are used here to indicate a values towards
                            a given Feature condition.
                            \nAdditionally, for Features that are a pure binary, data was split directly between
                            the True and False conditions. For features that are a continuous variable, the data
                            are split along the median value. Any numerical annotations or notes associated with a greater than or less than 
                            symbol is comparing the median value of that data feature. For features that are a categorical variable,
                            the Feature is identified as inappropriate for splitting.
                            ''')
                elif st.session_state.Inci_Value_display == 'Percentages':
                    st.write(f'''The Percentages mode shows the percentage of cells that match a condition, which are placed in
                            each cluster. In this way, the sum of each of the points on the line sums to 100%. For example, 
                            one way to interpret the Percentages is to say, {st.session_state.inci_df.iloc[0, 4]:.2f}% of the {st.session_state.inci_true_msg} cells
                            can be found in {st.session_state.inci_df.index[0]}. At the moment, we are only displaying the percentages for
                            the more positive condition of the selected condition and not the inverse. We will soon have an
                            option to flip the view to the more negative condition. An example of a value that is not displayed
                            presently (but is calculated) is {st.session_state.inci_df.iloc[0, 5]:.2f}% of the 
                            {st.session_state.inci_fals_msg} cells can be found in {st.session_state.inci_df.index[0]}.
                            \nAdditionally, for Features that are a pure binary, data was split directly between
                            the True and False conditions. For features that are a continuous variable, the data
                            are split along the median value. Any numerical annotations or notes associated with a greater than or less than 
                            symbol is comparing the median value of that data feature. For features that are a categorical variable,
                            the Feature is identified as inappropriate for splitting.
                            ''')
                elif st.session_state.Inci_Value_display == 'Ratios':
                    st.write(f'''The Ratios mode measures the log10 ratio between the Feature-Condition percentages across all clusters.
                             That is to say that if {st.session_state.inci_df.iloc[0, 4]:.2f}% of the {st.session_state.inci_true_msg} cells
                             can be found in {st.session_state.inci_df.index[0]} and {st.session_state.inci_df.iloc[0, 5]:.2f}% of the
                            {st.session_state.inci_fals_msg} cells also be found in {st.session_state.inci_df.index[0]}, the calculated
                            ratio of {st.session_state.inci_df.index[0]} is {st.session_state.inci_df.iloc[0, 8]:.2f}. You may find that
                            this ratio calculation is not a pure reflection of these two percentages. That is because we needed to account for
                            conditions where one of the percentages was 0%. Therefore the true ratio is the log10 division between adjusted 
                            percentages. Please refer to the documentation for the full equation.
                            \nAll ratios are shown in relation to the more *positive* outome. Therefore all ratios larger than
                            1 show a higher incidence of cells for the upper conditions. Ratios less than 1 show a higher incidence
                            for cells in the lower conditions.
                            \nAdditionally, for Features that are a pure binary, data was split directly between
                            the True and False conditions. For features that are a continuous variable, the data
                            are split along the median value. Any numerical annotations or notes associated with a greater than or less than 
                            symbol is comparing the median value of that data feature. For features that are a categorical variable,
                            the Feature is identified as inappropriate for splitting.
                            ''')
                else:
                    st.write('Invalid Display Option')
            else:
                st.write('''The Incidence Figure is not appropriate for the selected Feature. Please select a different Feature.''')
    else:
        st.warning('No spatial UMAP analysis detected. Please complete Neighborhood Profiles')

if __name__ == '__main__':
    main()
