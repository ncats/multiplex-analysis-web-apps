'''
This is the python script which produces the NEIGHBORHOOD PROFILES PAGE
'''
import os
import datetime
from copy import copy
import streamlit as st
import numpy as np
import dill
from streamlit_extras.add_vertical_space import add_vertical_space

# Import relevant libraries
import nidap_dashboard_lib as ndl   # Useful functions for dashboards connected to NIDAP
import basic_phenotyper_lib as bpl  # Useful functions for phenotyping collections of cells
import app_top_of_page as top
import streamlit_dataframe_editor as sde
from neighborhood_profiles import NeighborhoodProfiles, UMAPDensityProcessing

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

    st.session_state.density_completed = True

    # Save checkpoint for Neighborhood Profile structure
    save_neipro_struct()

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

    # Perform possible cluster variations with the completed UMAP
    # st.session_state.bc.startTimer()
    # with st.spinner('Calculating Possible Clusters'):
    #     st.session_state.clust_range, st.session_state.wcss = bpl.measure_possible_clust(st.session_state.spatial_umap, clust_minmax)
    # st.session_state.bc.printElapsedTime(msg = 'Calculating possible clusters')

    st.session_state.wcss_calc_completed = True
    st.session_state.umap_completed = True

    # Create Neighborhood Profiles Object
    st.session_state.npf = NeighborhoodProfiles(bc = st.session_state.bc)

    # Create Full UMAP example
    st.session_state.udp_full = UMAPDensityProcessing(st.session_state.npf, st.session_state.spatial_umap.df_umap)
    st.session_state.UMAPFig = st.session_state.udp_full.UMAPdraw_density()

    # Plot results
    filter_and_plot()

    # Save checkpoint for Neighborhood Profile structure
    save_neipro_struct()

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

    st.session_state.cluster_completed = True

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

    if st.session_state.umap_completed:
        st.session_state.spatial_umap.df_umap_filt = st.session_state.spatial_umap.df_umap.loc[st.session_state.spatial_umap.df_umap['Slide ID'] == st.session_state['selSlide ID'], :]
        if st.session_state['toggle_clust_diff']:
            palette = 'bwr'
        else:
            palette = 'tab20'
        st.session_state = ndl.setFigureObjs_UMAP(st.session_state, palette = palette)

def load_neipro_struct(file_name):
    '''
    Function to load the neighborhood profile structure
    '''

    # Load the Neighborhood Profile structure
    with open(f'{st.session_state.checkpoint_dir}/{file_name}', "rb") as dill_file:
        print(f'Loading Neighborhood Profiles Checkpoint-{file_name}')
        st.session_state.spatial_umap = dill.load(dill_file)

    st.session_state.phenotyping_completed = st.session_state.spatial_umap.phenotyping_completed
    st.session_state.density_completed     = st.session_state.spatial_umap.density_completed
    st.session_state.umap_completed        = st.session_state.spatial_umap.umap_completed
    st.session_state.cluster_completed     = st.session_state.spatial_umap.cluster_completed

    if st.session_state.umap_completed:
        # Create Neighborhood Profiles Object
        st.session_state.npf = NeighborhoodProfiles(bc = st.session_state.bc)

        # Create Full UMAP example
        st.session_state.udp_full = UMAPDensityProcessing(npf = st.session_state.npf, df = st.session_state.spatial_umap.df_umap)
        st.session_state.UMAPFig = st.session_state.udp_full.UMAPdraw_density()

        filter_and_plot()

def save_neipro_struct():
    '''
    Function to save the neighborhood profile structure
    '''

    file_name = f'neighborhood_profiles_checkpoint_{datetime.datetime.now().strftime("%d-%m-%Y %H-%M-%S")}'
    # Save the Neighborhood Profile structure
    with open(f'{st.session_state.checkpoint_dir}/{file_name}', "wb") as dill_file:
        print(f'Pickling Neighborhood Profiles Checkpoint-{file_name}')
        dill.dump(st.session_state.spatial_umap, dill_file)

def main():
    '''
    Main function for running the page
    '''

    clust_minmax = [1, 40]
    npf_cols = st.columns([1, 1, 2])
    with npf_cols[0]:
        nei_pro_tabs = st.tabs(['Analyze from Phenotyping', 'Load Previous Analysis'])
        with nei_pro_tabs[0]:
            dens_butt = st.button('Perform Cell Counts/Areas Analysis')
            umap_butt = st.button('Perform UMAP Analysis')
            st.toggle('Perform Clustering on UMAP Density Difference', value = False, key = 'toggle_clust_diff')
            if st.session_state['toggle_clust_diff'] is False: # Run Clustering Normally
                st.slider('Number of K-means clusters', min_value=clust_minmax[0], max_value=clust_minmax[1], key = 'slider_clus_val')
                clust_butt_disabled = False
            else:
                sep_clust_cols = st.columns(2)
                with sep_clust_cols[0]:
                    st.number_input('Number of Clusters for False Condition', min_value = 1, max_value = 10, value = 3, step = 1, key = 'num_clus_0')
                with sep_clust_cols[1]:
                    st.number_input('Number of Clusters for True Condition', min_value = 1, max_value = 10, value = 3, step = 1, key = 'num_clus_1')
                clust_butt_disabled = True
            clust_butt = st.button('Perform Clustering Analysis', disabled=clust_butt_disabled)
        with nei_pro_tabs[1]:
            neipro_checkpoint_files = os.listdir(st.session_state.checkpoint_dir)
            st.selectbox('Select Previous UMAP Results', options = neipro_checkpoint_files, key = 'sel_prev_umap')
            st.button('Load Selected UMAP Results', on_click=load_neipro_struct, args = (st.session_state.sel_prev_umap,))
            add_vertical_space(12)
    with npf_cols[1]:
        if dens_butt:
            if st.session_state.phenotyping_completed:
                init_spatial_umap()
        if umap_butt:
            if st.session_state.density_completed:
                apply_umap(umap_style = 'Densities')
        if clust_butt:
            if st.session_state.umap_completed:
                set_clusters()
        if st.session_state['toggle_clust_diff']:
            st.selectbox('Feature', options = st.session_state.outcomes, key = 'dens_diff_feat_sel')
            st.number_input('Cutoff Percentage', min_value = 0.01, max_value = 0.99, value = 0.2, step = 0.01, key = 'dens_diff_cutoff')

    with npf_cols[2]:
        if st.session_state.umap_completed:

            # Showing off the full UMAP
            with st.columns(3)[1]:
                st.pyplot(fig=st.session_state.UMAPFig)

            if st.session_state['toggle_clust_diff']:

                col_type = ndl.identify_col_type(st.session_state.spatial_umap.df_umap[st.session_state.dens_diff_feat_sel])

                if col_type == 'not_bool':
                    # Identify UMAP by Condition
                    median = np.round(st.session_state.spatial_umap.df_umap[st.session_state.dens_diff_feat_sel].median(), 2)
                    st.session_state.df_umap_fals = st.session_state.spatial_umap.df_umap.loc[st.session_state.spatial_umap.df_umap[st.session_state.dens_diff_feat_sel] <= median, :]
                    st.session_state.df_umap_true = st.session_state.spatial_umap.df_umap.loc[st.session_state.spatial_umap.df_umap[st.session_state.dens_diff_feat_sel] > median, :]
                    fals_msg = f'<= {median}'
                    true_msg = f'> {median}'
                    appro_feat = True
                elif col_type == 'bool':
                    # Identify UMAP by Condition
                    values = st.session_state.spatial_umap.df_umap[st.session_state.dens_diff_feat_sel].unique()
                    st.session_state.df_umap_fals = st.session_state.spatial_umap.df_umap.loc[st.session_state.spatial_umap.df_umap[st.session_state.dens_diff_feat_sel] == values[0], :]
                    st.session_state.df_umap_true = st.session_state.spatial_umap.df_umap.loc[st.session_state.spatial_umap.df_umap[st.session_state.dens_diff_feat_sel] == values[1], :]
                    fals_msg = f'= {values[0]}'
                    true_msg = f'= {values[1]}'
                    appro_feat = True
                else:
                    appro_feat = False
                    st.write('Feature must be boolean or numeric to perform density difference analysis')

                if appro_feat:
                    # Perform Density Calculations for each Condition
                    udp_fals = UMAPDensityProcessing(st.session_state.npf, st.session_state.df_umap_fals, xx=st.session_state.udp_full.xx, yy=st.session_state.udp_full.yy)
                    udp_true = UMAPDensityProcessing(st.session_state.npf, st.session_state.df_umap_true, xx=st.session_state.udp_full.xx, yy=st.session_state.udp_full.yy)

                    ## Copy over
                    udp_diff = copy(udp_fals)
                    ## Perform difference calculation
                    udp_diff.dens_mat = udp_true.dens_mat - udp_fals.dens_mat
                    ## Rerun the min/max calcs
                    udp_diff.umap_summary_stats()
                    ## Set Feature Labels
                    udp_fals.set_feature_label(st.session_state.dens_diff_feat_sel, fals_msg)
                    udp_true.set_feature_label(st.session_state.dens_diff_feat_sel, true_msg)
                    udp_diff.set_feature_label(st.session_state.dens_diff_feat_sel, 'Difference')

                    # Draw UMAPS
                    st.session_state.UMAPFig_fals = udp_fals.UMAPdraw_density()
                    st.session_state.UMAPFig_true = udp_true.UMAPdraw_density()
                    st.session_state.UMAPFig_diff = udp_diff.UMAPdraw_density(diff=True)

                    # Assign Masking and plot
                    udp_mask = copy(udp_diff)
                    udp_mask.filter_density_matrix(st.session_state.dens_diff_cutoff)
                    udp_mask.set_feature_label(st.session_state.dens_diff_feat_sel, f'Difference- Masked, cutoff = {st.session_state.dens_diff_cutoff}')
                    st.session_state.UMAPFig_mask = udp_mask.UMAPdraw_density(diff=True)

                    # Perform Clustering
                    udp_clus = copy(udp_mask)
                    udp_clus.perform_clustering(dens_mat_cmp=udp_mask.dens_mat,
                                                num_clus_0=st.session_state.num_clus_0,
                                                num_clus_1=st.session_state.num_clus_1)
                    udp_clus.set_feature_label(st.session_state.dens_diff_feat_sel, f'Clusters, False-{st.session_state.num_clus_0}, True-{st.session_state.num_clus_1}')
                    st.session_state.UMAPFig_clus = udp_clus.UMAPdraw_density(diff=True, legendtype='legend')
                    st.session_state.cluster_dict = udp_clus.cluster_dict

                    # Add cluster label column to cells dataframe
                    st.session_state.spatial_umap.df_umap.loc[:, 'clust_label'] = 'No Cluster'
                    st.session_state.spatial_umap.df_umap.loc[:, 'cluster'] = 'No Cluster'
                    st.session_state.spatial_umap.df_umap.loc[:, 'Cluster'] = 'No Cluster'

                    for key, val in st.session_state.cluster_dict.items():
                        if key != 0:
                            bin_clust = np.argwhere(udp_clus.dens_mat == key)
                            bin_clust = bin_clust[:, [1, 0]] # Swapping columns to by y, x
                            bin_clust = [tuple(x) for x in bin_clust]

                            significant_groups = st.session_state.udp_full.bin_indices_df_group[st.session_state.udp_full.bin_indices_df_group.set_index(['indx', 'indy']).index.isin(bin_clust)]

                            umap_ind = significant_groups.index.values
                            st.session_state.spatial_umap.df_umap.loc[umap_ind, 'clust_label'] = val
                            st.session_state.spatial_umap.df_umap.loc[umap_ind, 'cluster'] = val
                            st.session_state.spatial_umap.df_umap.loc[umap_ind, 'Cluster'] = val

                    # After assigning cluster labels, perform mean calculations
                    st.session_state.spatial_umap.mean_measures()

                    # Create the Cluster Scatterplot
                    filter_and_plot()

                    diff_cols = st.columns(3)
                    with diff_cols[0]:
                        st.pyplot(fig=st.session_state.UMAPFig_fals)
                    with diff_cols[1]:
                        st.pyplot(fig=st.session_state.UMAPFig_diff)
                    with diff_cols[2]:
                        st.pyplot(fig=st.session_state.UMAPFig_true)
                    mor_cols = st.columns(2)
                    with mor_cols[0]:
                        st.pyplot(fig=st.session_state.UMAPFig_mask)
                    with mor_cols[1]:
                        st.pyplot(fig=st.session_state.UMAPFig_clus)
                    st.session_state.cluster_completed = True

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
                #         if st.session_state.umap_completed:
                #             elbowFig = bpl.draw_wcss_elbow_plot(st.session_state.clust_range, st.session_state.wcss, st.session_state.selected_nClus)
                #             st.pyplot(elbowFig)

    if not st.session_state.phenotyping_completed:
        st.warning('Step 0: Please complete phentoyping analysis (See Phenotyping Page)', icon="⚠️")
    elif not st.session_state.density_completed:
        st.warning('Step 1: Please complete Cell Counts and Areas analysis', icon="⚠️")
    elif not st.session_state.umap_completed:
        st.warning('Step 2: Please run UMAP analysis', icon="⚠️")
    elif not st.session_state.cluster_completed:
        st.warning('Step 3: Please run clustering analysis', icon="⚠️")
    else:
        add_vertical_space(2)

    ### Visualizations ###
    uScatCol, uNeighPCol = st.columns(2)

    # Scatterplot Figure Column
    with uScatCol:
        # Print a column header
        st.header('Clusters Plot')

        # Plot Colors by Clusters or Phenotype
        clust_or_pheno = st.radio('Plot Colors by: ',
                                  ('Clusters', 'Phenotype'),
                                  horizontal = True, index = 0, key = 'clust_or_pheno')

        image_prog_col = st.columns([3, 1, 1, 2])
        with image_prog_col[0]:
            st.selectbox('Slide ID',
                         (st.session_state['uniSlide ID_short']),
                         key = 'selSlide ID_short',
                         on_change=slide_id_callback)
        with image_prog_col[1]:
            add_vertical_space(2)
            st.button('←', on_click=slide_id_prog_left_callback, disabled=st.session_state.prog_left_disabeled)
        with image_prog_col[2]:
            add_vertical_space(2)
            st.button('→', on_click=slide_id_prog_right_callback, disabled=st.session_state.prog_right_disabeled)
        with image_prog_col[3]:
            add_vertical_space(2)
            st.write(f'Image {st.session_state["idxSlide ID"]+1} of {st.session_state["numSlide ID"]}')

        if st.session_state.umap_completed:
            if clust_or_pheno == 'Clusters':
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
        with st.expander('Neighborhood Profile Options'):
            nei_sett_col = st.columns([1, 2, 1])
            with nei_sett_col[0]:
                st.toggle('Hide "Other" Phenotype', value = False, key = 'toggle_hide_other')
            with nei_sett_col[1]:
                st.toggle('Hide "No Cluster" Neighborhood Profile', value = False, key = 'toggle_hide_no_cluster')

        # If the spatial-umap is completed...
        if 'spatial_umap' in st.session_state:
            # List of Clusters to display
            if st.session_state['toggle_clust_diff']:
                list_clusters = list(st.session_state.cluster_dict.values())
            else:
                list_clusters = list(range(st.session_state.selected_nClus))
            if st.session_state['toggle_hide_no_cluster']:
                list_clusters.remove('No Cluster')

            cluster_sel_col = st.columns([3, 1])
            # Compare Clusters Toggle
            with cluster_sel_col[1]:
                add_vertical_space(2)
                st.toggle('Compare Cluster Neighborhoods', value = False, key = 'toggle_compare_clusters')
                if st.session_state['toggle_compare_clusters']:
                    st.radio('Compare as:', ('Difference', 'Ratio'), index = 0, key = 'compare_clusters_as', horizontal=True)

            # Cluster Select Widgets
            with cluster_sel_col[0]:
                sel_npf_fig  = st.selectbox('Select a cluster to view', list_clusters)
                sel_npf_fig2 = None
                if st.session_state['toggle_compare_clusters']:
                    sel_npf_fig2 = st.selectbox('Select a cluster to compare', list_clusters)

            if st.session_state.cluster_completed:
                # Draw the Neighborhood Profile
                npf_fig = bpl.neighProfileDraw(st.session_state.spatial_umap,
                                               sel_clus = sel_npf_fig,
                                               cmp_clus = sel_npf_fig2,
                                               cmp_style=st.session_state['compare_clusters_as'],
                                               hide_other = st.session_state['toggle_hide_other'],
                                               hide_no_cluster = st.session_state['toggle_hide_no_cluster'])
                st.pyplot(fig=npf_fig)

                # Create widgets for exporting the Neighborhood Profile images
                neigh_prof_col = st.columns([2, 1])
                with neigh_prof_col[0]:
                    st.text_input('.png file suffix (Optional)', key = 'neigh_prof_line_suffix')
                with neigh_prof_col[1]:
                    add_vertical_space(2)
                    if st.button('Append Export List', key = 'appendexportbutton_neighproline__do_not_persist'):

                        ndl.save_png(npf_fig, 'Neighborhood Profiles', st.session_state.neigh_prof_line_suffix)
                        st.toast(f'Added {st.session_state.neigh_prof_line_suffix} to export list')


if __name__ == '__main__':

    # Set a wide layout
    st.set_page_config(page_title="Neighborhood Profiles",
                       layout="wide")
    st.title('Neighborhood Profiles')

    # Run streamlit-dataframe-editor library initialization tasks at the top of the page
    st.session_state = sde.initialize_session_state(st.session_state)

    # Run Top of Page (TOP) functions
    st.session_state = top.top_of_page_reqs(st.session_state)

    main()

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)
