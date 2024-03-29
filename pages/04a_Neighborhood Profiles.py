'''
This is the python script which produces the NEIGHBORHOOD PROFILES PAGE
'''
import streamlit as st
import numpy as np
from streamlit_extras.add_vertical_space import add_vertical_space

# Import relevant libraries
import nidap_dashboard_lib as ndl   # Useful functions for dashboards connected to NIDAP
import basic_phenotyper_lib as bpl  # Useful functions for phenotyping collections of cells
import app_top_of_page as top
import streamlit_dataframe_editor as sde
import PlottingTools as umPT

def init_spatial_umap():
    '''
    Initalizing the spatial_umap object
    '''

    # Reset the settings required for Neighborhood Analysis
    st.session_state = ndl.reset_neigh_profile_settings(st.session_state)

    st.session_state.bc.startTimer()
    with st.spinner('Calculating Cell Counts and Areas'):
        st.session_state.spatial_umap = bpl.setup_Spatial_UMAP(st.session_state.df,
                                                               st.session_state.marker_multi_sel,
                                                               st.session_state.phenoOrder)

        st.session_state.spatial_umap = bpl.perform_density_calc(st.session_state.spatial_umap,
                                                                 st.session_state.bc,
                                                                 st.session_state.cpu_pool_size)
    st.write('Done Calculating Cell Counts and Areas')

    # Record time elapsed
    st.session_state.bc.set_value_df('time_to_run_counts', st.session_state.bc.elapsedTime())

    st.session_state.cell_counts_completed = True

def apply_umap(umap_style):
    '''
    Call back function for applying the UMAP functions
    '''
    clust_minmax = [1, 40]
    st.session_state.bc.startTimer()
    with st.spinner('Calculating UMAP'):
        st.session_state.spatial_umap = bpl.perform_spatialUMAP(st.session_state.spatial_umap,
                                                                st.session_state.bc,
                                                                umap_style)
    st.write('Done Calculating Spatial UMAP')

    # Record time elapsed
    st.session_state.bc.printElapsedTime(msg = 'Performing UMAP')
    st.session_state.bc.set_value_df('time_to_run_UMAP', st.session_state.bc.elapsedTime())

    # List of possible UMAP Lineages as defined by the completed UMAP
    st.session_state.umapPheno = [st.session_state.defLineageOpt]
    st.session_state.umapPheno.extend(st.session_state.pheno_summ['phenotype'])
    st.session_state.umapMarks = [st.session_state.defLineageOpt]
    st.session_state.umapMarks.extend(st.session_state.spatial_umap.markers)
    st.session_state.umapMarks.extend(['Other'])

    # Identify all of the features in the dataframe
    st.session_state.outcomes = st.session_state.spatial_umap.cells.columns

    # List of possible outcome variables as defined by the config yaml files
    st.session_state.umapOutcomes = [st.session_state.defumapOutcomes]
    st.session_state.umapOutcomes.extend(st.session_state.outcomes)
    st.session_state.inciOutcomes = [st.session_state.definciOutcomes]
    st.session_state.inciOutcomes.extend(st.session_state.outcomes)

    st.session_state.spatial_umap.prepare_df_umap_plotting(st.session_state.outcomes)

    st.session_state.df_umap = st.session_state.spatial_umap.df_umap

    # Perform possible cluster variations with the completed UMAP
    # st.session_state.bc.startTimer()
    # with st.spinner('Calculating Possible Clusters'):
    #     st.session_state.clust_range, st.session_state.wcss = bpl.measure_possible_clust(st.session_state.spatial_umap, clust_minmax)
    # st.session_state.bc.printElapsedTime(msg = 'Calculating possible clusters')

    st.session_state.wcss_calc_completed = True
    st.session_state.umapCompleted = True

    filter_and_plot()

def set_clusters():
    '''
    Callback function for setting the number of clusters
    and applying them to the UMAP/dataset
    '''
    st.session_state.bc.startTimer()
    st.session_state.spatial_umap = bpl.perform_clusteringUMAP(st.session_state.spatial_umap,
                                                               st.session_state.slider_clus_val)
    st.session_state.selected_nClus = st.session_state.slider_clus_val
    st.write('Done Calculating Clusters')

    # Record time elapsed
    st.session_state.bc.printElapsedTime(msg = 'Setting Clusters')
    st.session_state.bc.set_value_df('time_to_run_cluster', st.session_state.bc.elapsedTime())

    st.session_state.clustering_completed = True

    st.session_state.df_umap = st.session_state.spatial_umap.df_umap

    filter_and_plot()

def slide_id_prog_left_callback():
    '''
    callback function when the left Cell_ID progression button is clicked
    '''
    if st.session_state['idxSlide ID'] > 0:
        st.session_state['idxSlide ID'] -=1
        st.session_state['selSlide ID'] = st.session_state['uniSlide ID'][st.session_state['idxSlide ID']]
        st.session_state['selSlide ID_short'] = st.session_state['uniSlide ID_short'][st.session_state['idxSlide ID']]
        filter_and_plot()

def slide_id_prog_right_callback():
    '''
    callback function when the right Cell_ID progression button is clicked
    '''
    if st.session_state['idxSlide ID'] < st.session_state['numSlide ID']-1:
        st.session_state['idxSlide ID'] +=1
        st.session_state['selSlide ID'] = st.session_state['uniSlide ID'][st.session_state['idxSlide ID']]
        st.session_state['selSlide ID_short'] = st.session_state['uniSlide ID_short'][st.session_state['idxSlide ID']]
        filter_and_plot()

def slide_id_callback():
    '''
    Call back for setting the current Slide_ID index
    '''
    # st.session_state['idxSlide ID'] = st.session_state['uniSlide ID_short'].index(st.session_state['selSlide ID_short'])
    idx = st.session_state['idxSlide ID'] = st.session_state['uniSlide ID_short'].index(st.session_state['selSlide ID_short'])
    st.session_state['selSlide ID'] = st.session_state['uniSlide ID'][idx]

    # After correct index is selected, reapply filters, redraw plots
    filter_and_plot()

def filter_and_plot():
    '''
    function to update the filtering and the figure plotting
    '''
    st.session_state.prog_left_disabeled  = False
    st.session_state.prog_right_disabeled = False

    if st.session_state['idxSlide ID'] == 0:
        st.session_state.prog_left_disabeled = True

    if st.session_state['idxSlide ID'] == st.session_state['numSlide ID']-1:
        st.session_state.prog_right_disabeled = True

    if st.session_state.umapCompleted:
        st.session_state.df_umap_filt = st.session_state.df_umap.loc[st.session_state.df_umap['Slide ID'] == st.session_state['selSlide ID'], :]
        st.session_state = ndl.setFigureObjs_UMAP(st.session_state)

def main():
    '''
    Main function for running the page
    '''
    # Use the whole page width
    st.set_page_config(
        page_title="Neighborhood Profiles",
        layout="wide"
    )

    # Run streamlit-dataframe-editor library initialization tasks at the top of the page
    st.session_state = sde.initialize_session_state(st.session_state)

    # Run Top of Page (TOP) functions
    st.session_state = top.top_of_page_reqs(st.session_state)

    st.header('Neighborhood Profiles\nNCATS-NCI-DMAP')

    clust_minmax = [1, 40]
    npf_cols = st.columns([1, 1, 2])
    with npf_cols[0]:
        cellCountsButt = st.button('Perform Cell Counts/Areas Analysis')
        umap_butt       = st.button('Perform UMAP Analysis')
        st.toggle('Perform Clustering on UMAP Density Difference', value = False, key = 'toggle_clust_diff')
        st.slider('Number of K-means clusters', min_value=clust_minmax[0], max_value=clust_minmax[1], key = 'slider_clus_val')
        clustButt      = st.button('Perform Clustering Analysis')

    with npf_cols[1]:
        if st.session_state['toggle_clust_diff']:
            st.selectbox('Feature', options = ['Outcome', 'Survival_5yr'], key = 'dens_diff_feat_sel')
        if cellCountsButt:
            if st.session_state.phenotyping_completed:
                init_spatial_umap()
        if umap_butt:
            if st.session_state.cell_counts_completed:
                apply_umap(umap_style = 'Densities')
        if clustButt:
            if st.session_state.umapCompleted:
                set_clusters()

    with npf_cols[2]:
        if st.session_state.umapCompleted:
            vlim = .97
            n_bins = 200
            xx = np.linspace(np.min(st.session_state.df_umap['X']), np.max(st.session_state.df_umap['X']), n_bins + 1)
            yy = np.linspace(np.min(st.session_state.df_umap['Y']), np.max(st.session_state.df_umap['Y']), n_bins + 1)
            n_pad = 40

            w = None
            st.session_state.d_full = umPT.plot_2d_density(st.session_state.df_umap['X'],
                                                           st.session_state.df_umap['Y'],
                                                           bins=n_bins, w=w, return_matrix=True)

            st.session_state.UMAPFig = bpl.UMAPdraw_density(st.session_state.d_full, bins = [xx, yy], w=None, n_pad=n_pad, vlim=vlim)

            if st.session_state['toggle_clust_diff']:
                feat_comp1 = '= 1'
                feat_comp2 = '= 0'

                feat_label0 = f'{st.session_state.dens_diff_feat_sel} {feat_comp1} '
                feat_label1 = f'{st.session_state.dens_diff_feat_sel} {feat_comp2} '
                feat_labeld = f'{st.session_state.dens_diff_feat_sel} Difference '

                w = None

                st.session_state.df_umap_A = st.session_state.df_umap.loc[st.session_state.df_umap[st.session_state.dens_diff_feat_sel] == 1, :]
                st.session_state.df_umap_D = st.session_state.df_umap.loc[st.session_state.df_umap[st.session_state.dens_diff_feat_sel] == 0, :]

                st.session_state.d_A = umPT.plot_2d_density(st.session_state.df_umap_A['X'],
                                                            st.session_state.df_umap_A['Y'],
                                                            bins=n_bins, w=w, return_matrix=True)
                
                st.session_state.d_D = umPT.plot_2d_density(st.session_state.df_umap_D['X'],
                                                            st.session_state.df_umap_D['Y'],
                                                            bins=n_bins, w=w, return_matrix=True)
                
                st.session_state.d_diff = st.session_state.d_A - st.session_state.d_D

                st.session_state.UMAPFigDiff0_Dens = bpl.UMAPdraw_density(st.session_state.d_A, bins = [xx, yy], w=w, n_pad=n_pad, vlim=vlim, feat = feat_label0)
                st.session_state.UMAPFigDiff1_Dens = bpl.UMAPdraw_density(st.session_state.d_D, bins = [xx, yy], w=w, n_pad=n_pad, vlim=vlim, feat = feat_label1)
                st.session_state.UMAPFigDiff2_Dens = bpl.UMAPdraw_density(st.session_state.d_diff, bins = [xx, yy], w=w, n_pad=n_pad, vlim=vlim, feat = feat_labeld, diff = True)

                exp_cols = st.columns(3)
                with exp_cols[1]:
                    st.pyplot(fig=st.session_state.UMAPFig)
                diff_cols = st.columns(3)
                with diff_cols[0]:
                    st.pyplot(fig=st.session_state.UMAPFigDiff0_Dens)
                with diff_cols[1]:
                    st.pyplot(fig=st.session_state.UMAPFigDiff1_Dens)
                with diff_cols[2]:
                    st.pyplot(fig=st.session_state.UMAPFigDiff2_Dens)

            ### Clustering Meta Analysis and Description ###
            # with st.expander('Cluster Meta-Analysis', ):
            #     wcss_cols = st.columns(2)
            #     with wcss_cols[0]:
            #         st.markdown('''The within-cluster sum of squares (WCSS) is a measure of the
            #                     variability of the observations within each cluster. In general,
            #                     a cluster that has a small sum of squares is more compact than a
            #                     cluster that has a large sum of squares. Clusters that have higher
            #                     values exhibit greater variability of the observations within the
            #                     cluster.''')
            #     with wcss_cols[1]:
            #         if st.session_state.umapCompleted:
            #             elbowFig = bpl.draw_wcss_elbow_plot(st.session_state.clust_range, st.session_state.wcss, st.session_state.selected_nClus)
            #             st.pyplot(elbowFig)

    if not st.session_state.phenotyping_completed:
        st.warning('Step 0: Please complete phentoyping analysis (See Phenotyping Page)', icon="⚠️")
    elif not st.session_state.cell_counts_completed:
        st.warning('Step 1: Please complete Cell Counts and Areas analysis', icon="⚠️")
    elif not st.session_state.umapCompleted:
        st.warning('Step 2: Please run UMAP analysis', icon="⚠️")
    elif not st.session_state.clustering_completed:
        st.warning('Step 3: Please run clustering analysis', icon="⚠️")
    else:
        add_vertical_space(2)

    ### Visualizations ###
    uScatCol, uNeighPCol = st.columns(2)

    # Scatterplot Figure Column
    with uScatCol:
        # Print a column header
        st.header('Clusters Plot')

        clustOPheno = st.radio('Plot Colors by: ',
                               ('Clusters', 'Phenotype'),
                               horizontal = True, index = 0, key = 'clustOPheno')

        imageProgCol = st.columns([3, 1, 1, 2])
        with imageProgCol[0]:
            st.selectbox('Slide ID',
                         (st.session_state['uniSlide ID_short']),
                         key = 'selSlide ID_short',
                         on_change=slide_id_callback)
        with imageProgCol[1]:
            add_vertical_space(2)
            st.button('←', on_click=slide_id_prog_left_callback, disabled=st.session_state.prog_left_disabeled)
        with imageProgCol[2]:
            add_vertical_space(2)
            st.button('→', on_click=slide_id_prog_right_callback, disabled=st.session_state.prog_right_disabeled)
        with imageProgCol[3]:
            add_vertical_space(2)
            st.write(f'Image {st.session_state["idxSlide ID"]+1} of {st.session_state["numSlide ID"]}')

        if st.session_state.umapCompleted:
            if clustOPheno == 'Clusters':
                st.pyplot(st.session_state.seabornFig_clust)
            else:
                st.pyplot(st.session_state.phenoFig)

            clust_fig_col = st.columns([2, 1])
            with clust_fig_col[0]:
                st.text_input('.png file suffix (Optional)', key = 'cluster_scatter_suffix')
            with clust_fig_col[1]:
                add_vertical_space(2)
                if st.button('Append Export List', key = 'appendexportbutton_clusterscatter__do_not_persist'):
                    ndl.save_png(st.session_state.seabornFig_clust, 'Cluster Scatterplot', st.session_state.cluster_scatter_suffix)
                    st.toast(f'Added {st.session_state.cluster_scatter_suffix} to export list')

    # Neighborhood Profiles Figure Column
    with uNeighPCol:
        st.header('Neighborhood Profiles')
        if 'spatial_umap' in st.session_state:
            sel_npf_fig = st.selectbox('Select a cluster to view',
                                       list(range(st.session_state.selected_nClus)))
            if st.session_state.clustering_completed:

                npf_fig = bpl.neighProfileDraw(st.session_state.spatial_umap, sel_npf_fig)
                st.pyplot(fig=npf_fig)

                neigh_prof_col = st.columns([2, 1])
                with neigh_prof_col[0]:
                    st.text_input('.png file suffix (Optional)', key = 'neigh_prof_line_suffix')
                with neigh_prof_col[1]:
                    add_vertical_space(2)
                    if st.button('Append Export List', key = 'appendexportbutton_neighproline__do_not_persist'):

                        ndl.save_png(npf_fig, 'Neighborhood Profiles', st.session_state.neigh_prof_line_suffix)
                        st.toast(f'Added {st.session_state.neigh_prof_line_suffix} to export list')

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)

if __name__ == '__main__':
    main()
