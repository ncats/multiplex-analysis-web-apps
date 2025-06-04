'''
This is the python script which produces the NEIGHBORHOOD PROFILES PAGE
'''

from copy import copy
import streamlit as st
import numpy as np
import dill
from streamlit_extras.add_vertical_space import add_vertical_space
import matplotlib.pyplot as plt
import pandas as pd
from natsort import natsorted
import nidap_dashboard_lib as ndl   # Useful functions for dashboards connected to NIDAP
import basic_phenotyper_lib as bpl  # Useful functions for phenotyping collections of cells
from neighborhood_profiles import NeighborhoodProfiles, UMAPDensityProcessing

def get_spatialUMAP(spatial_umap, bc, umap_subset_per_fit, umap_subset_toggle, umap_subset_per):
    '''
    Extract precomputed UMAP from the file

    Args:
        spatial_umap (spatial_umap): spatial_umap object
        bc (benchmark_collector): Benchmark Collector object
        UMAPStyle (str): Style of UMAP to use
    
    Returns:
        spatial_umap: spatial_umap object with the UMAP analysis performed
    '''

    min_image_size = spatial_umap.smallest_image_size
    n_fit = int(min_image_size*umap_subset_per_fit/100)
    n_tra = n_fit + int(min_image_size*umap_subset_per/100)

    # set training and "test" cells for umap training and embedding, respectively
    print('Setting Train/Test Split')
    spatial_umap.set_train_test(n_fit=n_fit, n_tra = n_tra, groupby_label = 'TMA_core_id', seed=54321, umap_subset_toggle = umap_subset_toggle)

    # fit umap on training cells
    bc.startTimer()
    print('Fitting Model')
    spatial_umap.umap_fit = spatial_umap.cells.loc[spatial_umap.cells['umap_train'].values, ['UMAP_1_20230327_152849', 'UMAP_2_20230327_152849']].values.reshape((spatial_umap.cells['umap_train'].sum(), -1))
    bc.printElapsedTime(f'      Fitting {np.sum(spatial_umap.cells["umap_train"] == 1)} points to a model')

    # Transform test cells based on fitted model
    bc.startTimer()
    print('Transforming Data')
    spatial_umap.umap_test = spatial_umap.cells.loc[spatial_umap.cells['umap_test'].values, ['UMAP_1_20230327_152849', 'UMAP_2_20230327_152849']].values.reshape((spatial_umap.cells['umap_test'].sum(), -1))
    bc.printElapsedTime(f'      Transforming {np.sum(spatial_umap.cells["umap_test"] == 1)} points with the model')

    spatial_umap.umap_completed = True
    
    # import pickle
    # with open('../Edits/spatial_umap_original_precomp.pkl', 'wb') as f:
    #     pickle.dump(spatial_umap, f)

    return spatial_umap

def init_spatial_umap():
    '''
    Initalizing the spatial_umap object
    '''

    # Reset the settings required for Neighborhood Analysis
    st.session_state = ndl.reset_neigh_profile_settings(st.session_state)

    if not st.session_state['calc_unique_areas_toggle']:
        area_filter = 0
    else:
        area_filter = st.session_state['area_filter_per']

    st.session_state.bc.startTimer()
    with st.spinner('Calculating Cell Counts and Areas', show_time=True):
        st.session_state.spatial_umap = bpl.setup_Spatial_UMAP(st.session_state.df,
                                                               st.session_state.marker_multi_sel,
                                                               st.session_state.phenoOrder,
                                                               st.session_state.datafile_min_img_size)

        st.session_state.spatial_umap = bpl.perform_density_calc(st.session_state.spatial_umap,
                                                                 st.session_state.bc,
                                                                 st.session_state.calc_unique_areas_toggle,
                                                                 st.session_state.cpu_pool_size,
                                                                 area_filter)

        # Record time elapsed
        st.session_state.bc.set_value_df('time_to_run_counts', st.session_state.bc.elapsedTime())

        st.session_state.density_completed = True

        # Save checkpoint for Neighborhood Profile structure
        # save_neipro_struct()

def apply_umap(umap_style):
    '''
    Call back function for applying the UMAP functions
    '''

    st.session_state.bc.startTimer()
    # if toggle for loading pre-generated UMAP is selected extract UMAP from file, works only with a specific dataset
    if st.session_state['load_generated_umap_toggle']:
        st.session_state.spatial_umap = get_spatialUMAP(st.session_state.spatial_umap,
                                                        st.session_state.bc,
                                                        st.session_state.umap_subset_per_fit,
                                                        st.session_state.umap_subset_toggle,
                                                        st.session_state.umap_subset_per)
    else:
        with st.spinner('Calculating UMAP', show_time=True):
            st.session_state.spatial_umap = bpl.perform_spatialUMAP(st.session_state.spatial_umap,
                                                                    st.session_state.bc,
                                                                    st.session_state.umap_subset_per_fit,
                                                                    st.session_state.umap_subset_toggle,
                                                                    st.session_state.umap_subset_per)

    # Record time elapsed
    st.session_state.bc.printElapsedTime(msg = 'Performing UMAP')
    st.session_state.bc.set_value_df('time_to_run_UMAP', st.session_state.bc.elapsedTime())

    # List of possible UMAP Lineages as defined by the completed UMAP
    st.session_state.umapPheno = [st.session_state.def_lineage_opt]
    st.session_state.umapPheno.extend(st.session_state.pheno_summ['phenotype'])
    st.session_state.umapMarks = [st.session_state.def_lineage_opt]
    st.session_state.umapMarks.extend(st.session_state.spatial_umap.markers)
    st.session_state.umapMarks.extend(['Other'])

    # Identify all of the features in the dataframe
    st.session_state.outcomes = st.session_state.spatial_umap.cells.columns
    st.session_state.spatial_umap.outcomes = st.session_state.spatial_umap.cells.columns

    st.session_state.dens_diff_feat_sel = st.session_state.outcomes[0]

    # List of possible outcome variables as defined by the config yaml files
    st.session_state.umapOutcomes = [st.session_state.def_umap_feature]
    st.session_state.umapOutcomes.extend(st.session_state.outcomes)
    st.session_state.inciOutcomes = [st.session_state.def_inci_feature]
    st.session_state.inciOutcomes.extend(st.session_state.outcomes)

    # creates the df_umap dataframe for plotting
    # creates the df_umap dataframe for plotting
    st.session_state.spatial_umap.prepare_df_umap_plotting(st.session_state.outcomes)

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
    # save_neipro_struct()

def set_clusters():
    '''
    Callback function for setting the number of clusters
    and applying them to the UMAP/dataset
    '''
    st.session_state.bc.startTimer()

    with st.spinner('Calculating Clusters', show_time=True):
        # If clustering is to be performed on the UMAP density difference
        if st.session_state['toggle_clust_diff'] and st.session_state['appro_feat']:

            # Split the UMAP by the selected values of the feature
            split_dict_full = st.session_state.udp_full.split_df_by_feature(st.session_state.dens_diff_feat_sel,
                                                                            st.session_state.feature_value_fals,
                                                                            st.session_state.feature_value_true,
                                                                            st.session_state.clust_diff_vals_code)

            # Perform Density Calculations for each Condition
            udp_fals = UMAPDensityProcessing(st.session_state.npf, split_dict_full['df_umap_fals'], xx=st.session_state.udp_full.xx, yy=st.session_state.udp_full.yy)
            udp_true = UMAPDensityProcessing(st.session_state.npf, split_dict_full['df_umap_true'], xx=st.session_state.udp_full.xx, yy=st.session_state.udp_full.yy)

            ## Copy over
            udp_diff = copy(udp_fals)
            ## Perform difference calculation
            udp_diff.dens_mat = np.log10(udp_fals.dens_mat) - np.log10(udp_true.dens_mat)
            ## Rerun the min/max calcs
            udp_diff.umap_summary_stats()
            ## Set Feature Labels
            udp_fals.set_feature_label(st.session_state.dens_diff_feat_sel, split_dict_full['fals_msg'])
            udp_true.set_feature_label(st.session_state.dens_diff_feat_sel, split_dict_full['true_msg'])
            udp_diff.set_feature_label(st.session_state.dens_diff_feat_sel, 'Difference')

            # Draw UMAPS
            st.session_state.UMAPFig_fals = udp_fals.UMAPdraw_density()
            st.session_state.UMAPFig_true = udp_true.UMAPdraw_density()
            st.session_state.UMAPFig_diff = udp_diff.UMAPdraw_density(diff=True)

            # Assign Masking and plot
            udp_mask = copy(udp_diff)
            udp_mask.filter_density_matrix(st.session_state.dens_diff_cutoff, st.session_state.udp_full.empty_bin_ind)
            udp_mask.set_feature_label(st.session_state.dens_diff_feat_sel, f'Difference- Masked, \ncutoff = {st.session_state.dens_diff_cutoff}')
            st.session_state.UMAPFig_mask = udp_mask.UMAPdraw_density(diff=True)

            # Perform Clustering
            udp_clus = copy(udp_mask)
            udp_clus.perform_clustering(dens_mat_cmp=udp_mask.dens_mat,
                                        num_clus_0=st.session_state.num_clus_0,
                                        num_clus_1=st.session_state.num_clus_1,
                                        clust_minmax=st.session_state.clust_minmax,
                                        cpu_pool_size=3)
            udp_clus.set_feature_label(st.session_state.dens_diff_feat_sel, f'Clusters, False-{st.session_state.num_clus_0}, True-{st.session_state.num_clus_1}')
            st.session_state.UMAPFig_clus = udp_clus.UMAPdraw_density(diff=True, legendtype='legend')
            st.session_state.cluster_dict = udp_clus.cluster_dict
            st.session_state.palette_dict = udp_clus.palette_dict
            st.session_state.elbow_fig_0 = udp_clus.elbow_fig_0
            st.session_state.elbow_fig_1 = udp_clus.elbow_fig_1

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

            # Benchmark how long it took to untangle indicies
            st.session_state.bc.printElapsedTime('Untangling bin indicies with UMAP indicies', split = True)

            # After assigning cluster labels, perform mean calculations
            st.session_state.spatial_umap.mean_measures()
            st.session_state.bc.printElapsedTime('Performing Mean Measures', split = True)

            # Average False condition and Average True Condition
            dens_df_fals = st.session_state.spatial_umap.dens_df_mean.loc[st.session_state.spatial_umap.dens_df_mean['clust_label'].str.contains('False'), :]
            dens_df_true = st.session_state.spatial_umap.dens_df_mean.loc[st.session_state.spatial_umap.dens_df_mean['clust_label'].str.contains('True'), :]

            dens_df_fals['clust_label'] = 'Average False'
            dens_df_mean_fals = dens_df_fals.groupby(['clust_label', 'phenotype', 'dist_bin'], as_index=False).mean()

            dens_df_true['clust_label'] = 'Average True'
            dens_df_mean_true = dens_df_true.groupby(['clust_label', 'phenotype', 'dist_bin'], as_index=False).mean()

            st.session_state.spatial_umap.dens_df_mean = pd.concat([st.session_state.spatial_umap.dens_df_mean, dens_df_mean_fals, dens_df_mean_true], axis=0)

            st.session_state.cluster_completed_diff = True

        else:
            st.session_state.spatial_umap = bpl.umap_clustering(spatial_umap = st.session_state.spatial_umap,
                                                                n_clusters = st.session_state.slider_clus_val,
                                                                clust_minmax = st.session_state.clust_minmax,
                                                                cpu_pool_size = 3)

            st.session_state.appro_feat = True
            st.session_state.cluster_completed_diff = False
            st.session_state.cluster_dict = st.session_state.spatial_umap.cluster_dict
            st.session_state.palette_dict = st.session_state.spatial_umap.palette_dict
            st.session_state.selected_nClus = st.session_state.slider_clus_val

        # Draw the 2D histogram UMAP colored by the clusters
        st.session_state.udp_full.cluster_dict = st.session_state.cluster_dict
        st.session_state.udp_full.palette_dict = st.session_state.palette_dict
        st.session_state.diff_clust_Fig = st.session_state.udp_full.umap_draw_clusters()

    st.session_state.bc.printElapsedTime(msg = 'Setting Clusters')
    st.session_state.bc.set_value_df('time_to_run_cluster', st.session_state.bc.elapsedTime())

    st.session_state.cluster_completed = True

    filter_and_plot()

def check_feature_approval_callback():
    '''
    Simple callback to test the current value of 
    st.session_state.dens_diff_feat_sel
    '''

    if not st.session_state['toggle_clust_diff']:
        st.session_state.appro_feat = True
    else:

        # Check feature values
        st.session_state.clust_diff_vals_code = st.session_state.udp_full.check_feature_values(st.session_state.dens_diff_feat_sel)

        if st.session_state.clust_diff_vals_code == 0:
            st.session_state.appro_feat = False
        else:
            st.session_state.appro_feat = True

        feature_vals = natsorted(st.session_state.udp_full.df[st.session_state.dens_diff_feat_sel].unique())
        if st.session_state.clust_diff_vals_code == 2:
            options_fals = [feature_vals[0]]
            options_true = [feature_vals[1]]
        elif st.session_state.clust_diff_vals_code > 2 and st.session_state.clust_diff_vals_code <= 15:
            options_fals = feature_vals
            options_true = feature_vals
        elif st.session_state.clust_diff_vals_code == 100:
            median = np.round(st.session_state.udp_full.df[st.session_state.dens_diff_feat_sel].median(), decimals = 2)

            options_fals = [median]
            options_true = [median]
        else:
            options_fals = ['None']
            options_true = ['None']

        st.session_state.clus_diff_vals_fals = options_fals
        st.session_state.clus_diff_vals_true = options_true
        st.session_state.feature_value_fals = options_fals[0]
        st.session_state.feature_value_true = options_true[0]

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
    callback function to update the filtering and the 
    figure plotting
    '''

    st.session_state.prog_left_disabled  = False
    st.session_state.prog_right_disabled = False

    if st.session_state['idxSlide ID'] == 0:
        st.session_state.prog_left_disabled = True

    if st.session_state['idxSlide ID'] == st.session_state['numSlide ID']-1:
        st.session_state.prog_right_disabled = True

    if st.session_state.umap_completed:
        st.session_state.spatial_umap.df_umap_filt = st.session_state.spatial_umap.df_umap.loc[st.session_state.spatial_umap.df_umap['Slide ID'] == st.session_state['selSlide ID'], :]
        st.session_state = ndl.setFigureObjs_UMAP(st.session_state, palette = st.session_state.palette_dict)

def load_neipro_struct():
    '''
    Function to load the neighborhood profile structure
    '''

    # Load the Neighborhood Profile structure
    with open(f'{st.session_state.checkpoint_dir}/neighborhood_profiles_checkpoint', "rb") as dill_file:
        print('Loading Neighborhood Profiles Checkpoint-neighborhood_profiles_checkpoint.pkl')
        st.session_state.spatial_umap = dill.load(dill_file)

    st.session_state.phenotyping_completed = st.session_state.spatial_umap.phenotyping_completed
    st.session_state.density_completed     = st.session_state.spatial_umap.density_completed
    st.session_state.umap_completed        = st.session_state.spatial_umap.umap_completed
    st.session_state.cluster_completed     = st.session_state.spatial_umap.cluster_completed

    # Slide ID Progression Initializeion
    st.session_state['uniSlide ID'] = st.session_state.spatial_umap.df_umap['Slide ID'].unique()
    st.session_state['selSlide ID'] = st.session_state['uniSlide ID'][0]
    st.session_state['idxSlide ID'] = 0
    st.session_state['numSlide ID'] = len(st.session_state['uniSlide ID'])
    st.session_state['uniSlide ID_short'] = st.session_state['uniSlide ID']
    st.session_state['selSlide ID_short'] = st.session_state['uniSlide ID_short'][0]

    st.session_state.prog_left_disabled = True
    st.session_state.prog_right_disabled = False
    if st.session_state['numSlide ID'] == 1:
        st.session_state.prog_right_disabled = True

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

    file_name = 'neighborhood_profiles_checkpoint'
    # Save the Neighborhood Profile structure
    with open(f'{st.session_state.checkpoint_dir}/{file_name}', "wb") as dill_file:
        print(f'Pickling Neighborhood Profiles Checkpoint-{file_name}.pkl')
        dill.dump(st.session_state.spatial_umap, dill_file)
    
    st.toast('Neighborhood Profiles Analysis Checkpoint Saved')

def diff_density_analysis():
    '''
    Function to perform the density difference analysis
    '''
    
    # If clustering is to be performed on the UMAP density difference
    if st.session_state['toggle_clust_diff']:
        
        # Idenfify the column type that is splitting the UMAP
        col_type = ndl.identify_col_type(st.session_state.spatial_umap.df_umap[st.session_state.dens_diff_feat_sel])

        if col_type == 'not_bool' or col_type == 'bool':
            st.session_state.feature_appro_message = None

            if col_type == 'not_bool':
                # Identify UMAP by Condition
                median = np.round(st.session_state.spatial_umap.df_umap[st.session_state.dens_diff_feat_sel].median(), 2)
                st.session_state.df_umap_fals = st.session_state.spatial_umap.df_umap.loc[st.session_state.spatial_umap.df_umap[st.session_state.dens_diff_feat_sel] <= median, :]
                st.session_state.df_umap_true = st.session_state.spatial_umap.df_umap.loc[st.session_state.spatial_umap.df_umap[st.session_state.dens_diff_feat_sel] > median, :]
                fals_msg = f'<= {median}'
                true_msg = f'> {median}'
            elif col_type == 'bool':
                # Identify UMAP by Condition
                values = st.session_state.spatial_umap.df_umap[st.session_state.dens_diff_feat_sel].unique()
                st.session_state.df_umap_fals = st.session_state.spatial_umap.df_umap.loc[st.session_state.spatial_umap.df_umap[st.session_state.dens_diff_feat_sel] == values[0], :]
                st.session_state.df_umap_true = st.session_state.spatial_umap.df_umap.loc[st.session_state.spatial_umap.df_umap[st.session_state.dens_diff_feat_sel] == values[1], :]
                fals_msg = f'= {values[0]}'
                true_msg = f'= {values[1]}'

            # Perform Density Calculations for each Condition
            st.session_state.udp_fals = UMAPDensityProcessing(st.session_state.npf, st.session_state.df_umap_fals, xx=st.session_state.udp_full.xx, yy=st.session_state.udp_full.yy)
            st.session_state.udp_true = UMAPDensityProcessing(st.session_state.npf, st.session_state.df_umap_true, xx=st.session_state.udp_full.xx, yy=st.session_state.udp_full.yy)

            ## Copy over
            st.session_state.udp_diff = copy(st.session_state.udp_fals)
            ## Perform difference calculation
            st.session_state.udp_diff.dens_mat = np.log10(st.session_state.udp_fals.dens_mat) - np.log10(st.session_state.udp_true.dens_mat)
            ## Rerun the min/max calcs
            st.session_state.udp_diff.umap_summary_stats()
            ## Set Feature Labels
            st.session_state.udp_fals.set_feature_label(st.session_state.dens_diff_feat_sel, fals_msg)
            st.session_state.udp_true.set_feature_label(st.session_state.dens_diff_feat_sel, true_msg)
            st.session_state.udp_diff.set_feature_label(st.session_state.dens_diff_feat_sel, 'Difference')

            # Draw UMAPS
            st.session_state.UMAPFig_fals = st.session_state.udp_fals.UMAPdraw_density()
            st.session_state.UMAPFig_true = st.session_state.udp_true.UMAPdraw_density()
            st.session_state.UMAPFig_diff = st.session_state.udp_diff.UMAPdraw_density(diff=True)

            # Assign Masking and plot
            st.session_state.udp_mask = copy(st.session_state.udp_diff)
            st.session_state.udp_mask.filter_density_matrix(st.session_state.dens_diff_cutoff)
            st.session_state.udp_mask.set_feature_label(st.session_state.dens_diff_feat_sel, f'Difference- Masked, cutoff = {st.session_state.dens_diff_cutoff}')
            st.session_state.UMAPFig_mask = st.session_state.udp_mask.UMAPdraw_density(diff=True)

            diff_density_apply_masking()
        else:
            st.session_state.feature_appro_message = 'Feature must be boolean or numeric to perform density difference analysis'

def diff_density_apply_masking():
    '''
    Function to apply masking to the density difference analysis
    '''
    # Assign Masking and plot
    st.session_state.udp_mask = copy(st.session_state.udp_diff)
    st.session_state.udp_mask.filter_density_matrix(st.session_state.dens_diff_cutoff)
    st.session_state.udp_mask.set_feature_label(st.session_state.dens_diff_feat_sel, f'Difference- Masked, cutoff = {st.session_state.dens_diff_cutoff}')
    st.session_state.UMAPFig_mask = st.session_state.udp_mask.UMAPdraw_density(diff=True)

    # Perform Clustering
    diff_density_perform_clustering()

def diff_density_perform_clustering():
    '''
    Function to set the clusters for the density difference analysis
    '''
    # Perform Clustering
    st.session_state.udp_clus = copy(st.session_state.udp_mask)
    st.session_state.udp_clus.perform_clustering(dens_mat_cmp=st.session_state.udp_mask.dens_mat,
                                                    num_clus_0=st.session_state.num_clus_0,
                                                    num_clus_1=st.session_state.num_clus_1)
    st.session_state.udp_clus.set_feature_label(st.session_state.dens_diff_feat_sel, f'Clusters, False-{st.session_state.num_clus_0}, True-{st.session_state.num_clus_1}')
    st.session_state.UMAPFig_clus = st.session_state.udp_clus.UMAPdraw_density(diff=True, legendtype='legend')
    st.session_state.cluster_dict = st.session_state.udp_clus.cluster_dict

    # Add cluster label column to cells dataframe
    st.session_state.spatial_umap.df_umap.loc[:, 'clust_label'] = 'No Cluster'
    st.session_state.spatial_umap.df_umap.loc[:, 'cluster'] = 'No Cluster'
    st.session_state.spatial_umap.df_umap.loc[:, 'Cluster'] = 'No Cluster'

    for key, val in st.session_state.cluster_dict.items():
        if key != 0:
            bin_clust = np.argwhere(st.session_state.udp_clus.dens_mat == key)
            bin_clust = bin_clust[:, [1, 0]] # Swapping columns to by y, x
            bin_clust = [tuple(x) for x in bin_clust]

            significant_groups = st.session_state.udp_full.bin_indices_df_group[st.session_state.udp_full.bin_indices_df_group.set_index(['indx', 'indy']).index.isin(bin_clust)]

            umap_ind = significant_groups.index.values
            st.session_state.spatial_umap.df_umap.loc[umap_ind, 'clust_label'] = val
            st.session_state.spatial_umap.df_umap.loc[umap_ind, 'cluster'] = val
            st.session_state.spatial_umap.df_umap.loc[umap_ind, 'Cluster'] = val

    # After assigning cluster labels, perform mean calculations
    # print(st.session_state.spatial_umap.df_umap['clust_label'].value_counts())
    st.session_state.spatial_umap.mean_measures()

    # Organize the mean density data
    dens_df_fals = st.session_state.spatial_umap.dens_df_mean.loc[st.session_state.spatial_umap.dens_df_mean['clust_label'].str.contains('False'), :]
    dens_df_true = st.session_state.spatial_umap.dens_df_mean.loc[st.session_state.spatial_umap.dens_df_mean['clust_label'].str.contains('True'), :]

    dens_df_fals['clust_label'] = 'Average Deceased'
    dens_df_mean_fals = dens_df_fals.groupby(['clust_label', 'phenotype', 'dist_bin'], as_index=False).mean()

    dens_df_true['clust_label'] = 'Average Alive'
    dens_df_mean_true = dens_df_true.groupby(['clust_label', 'phenotype', 'dist_bin'], as_index=False).mean()

    st.session_state.spatial_umap.dens_df_mean = pd.concat([st.session_state.spatial_umap.dens_df_mean, dens_df_mean_fals, dens_df_mean_true], axis=0)

    # Create the Cluster Scatterplot
    filter_and_plot()

def reset_neipro_feature_values():
    '''
    Quick callback function to reset the feature values used
    filter the Neighborhood Profiles figure
    '''

    selected_feat = st.session_state['NeiPro_filter_feat']
    unique_values = st.session_state.spatial_umap.df_umap[selected_feat].unique()
    st.session_state.NeiPro_filter_value = unique_values[0]

def main():
    '''
    Main function for running the page
    '''

    nei_pro_tabs = st.tabs(['Analyze from Phenotyping', 'Load Previous Analysis'])
    with nei_pro_tabs[0]:

        # Neighborhood Profiles Settings
        with st.expander('Neighborhood Profiles Settings', expanded = False):
            neipro_settings = st.columns([1, 1, 1, 1])
            with neipro_settings[0]:
                st.number_input('Number of CPUs', min_value = 1, max_value= 8, step = 1,
                                key = 'cpu_pool_size',
                                help = '''Number of CPUs to use for parallel processing.
                                This effects the speed of the Cell Density Analysis''')
                st.toggle('Calculate unique areas', key = 'calc_unique_areas_toggle',
                          help = '''For Multiplex studies measuring the density of cells
                          surrounding each cell, specific radii are drawn to measure how
                          density changes over the distance from the cell. If the same radii
                          is used for each cell, than the areas encompassed by the radii are
                          also the same. However, the actual areas surrounding a given cell may
                          be different due to edge effects of the tissue. When this toggle is
                          false, all areas are assumed to be the same. When this toggle is true,
                          the areas are calculated for each cell based on how the radii
                          interact with the edges of the tissue. Turning this on will make 
                          the density calculation run longer.
                          ''')
                st.number_input('Area Filter Ratio', min_value = 0.001, max_value = 1.0, step = 0.001,
                                format="%.3f", key = 'area_filter_per',
                                disabled=not st.session_state['calc_unique_areas_toggle'],
                                help = '''The area filter ratio identifies how much of an area surrounding
                                a cell can be to be considered large enough to be included in the density calculations.
                                Small values of the ratio (close to 0) include more cells, and large values
                                of the ratio (close to 1) include fewer cells. This can be useful for removing
                                cells that are on the edge of the tissue. This affects the number of cells included
                                in the Spatial UMAP processing and the speed of that step.''')
            with neipro_settings[1]:
                st.markdown(f'''Smallest image in dataset is {st.session_state.datafile_min_img_size} cells.
                            What percentage from each image should be used for the UMAP fitting step?''')
                st.number_input('Percentage of cells to Subset for Fitting Step',
                                min_value = 20, max_value = 80, step = 10,
                                key = 'umap_subset_per_fit')

            with neipro_settings[2]:
                st.toggle('Subset data transformed by UMAP',
                          key = 'umap_subset_toggle',
                          help = '''The UMAP model is always trained on a percentage of data included
                          in the smallest image. You can choose to transform the entire dataset using 
                          this trained model, or only transform a percentage of the data. This can be 
                          useful for large datasets. If a percentage is chosen for transformation, it 
                          is always a different sample than what the model was trained on.''')
                add_vertical_space(2)
                st.number_input('Percentage of cells to Subset for Transforming Step',
                                min_value = 20, max_value = 80, step = 10,
                                key = 'umap_subset_per', disabled = not st.session_state.umap_subset_toggle)

            with neipro_settings[3]:
                st.toggle('Load pre-generated UMAP',
                          value = False, key = 'load_generated_umap_toggle',)

        # Main Neighborhood Profiles Columns for the buttons and UMAP preview
        npf_cols = st.columns([2, 3])

        # Analysis Buttons
        with npf_cols[0]:
            butt_cols = st.columns(2)
            with butt_cols[0]:

                dens_butt  = st.button('Perform Cell Density Analysis')
                umap_butt  = st.button('Perform UMAP Analysis')
                clust_butt = st.button('Perform Clustering Analysis')

            # Button results and difference settings
            with butt_cols[1]:
                if st.session_state.phenotyping_completed:
                    if dens_butt:
                        init_spatial_umap()
                    if not st.session_state.density_completed:
                        st.write(':x: Step 1: Perform Cell Density')
                    else:
                        st.write(':white_check_mark: Density Analysis Complete')
                    add_vertical_space(1)

                if st.session_state.phenotyping_completed:
                    if st.session_state.density_completed:
                        if umap_butt:
                            apply_umap(umap_style = 'Densities')
                        if not st.session_state.umap_completed:
                            st.write(':x: Step 2: Perform UMAP')
                        else:
                            st.write(':white_check_mark: UMAP Analysis Completed')
                        add_vertical_space(1)

                if st.session_state.phenotyping_completed:
                    if st.session_state.umap_completed:

                        if clust_butt:
                            set_clusters()
                        if not st.session_state.cluster_completed:
                            st.write(':x: Step 3: Perform Clustering')
                        else:
                            st.write(':white_check_mark: Clustering Analysis Completed')

            # If UMAP is completed, display the clustering settings
            if st.session_state.umap_completed:
                with st.expander('Clustering Settings', expanded = True):
                    st.toggle('Perform Clustering on UMAP Density Difference',
                              value = False, key = 'toggle_clust_diff',
                              help = '''Perform clustering on the density difference between
                                        two levels of a dataset feature.''',
                                        on_change=check_feature_approval_callback)

                    clust_exp_col = st.columns(2)
                    with clust_exp_col[0]:

                        # Perform clustering on split features
                        if st.session_state['toggle_clust_diff'] is True:
                            st.selectbox('Feature', options = st.session_state.spatial_umap.outcomes,
                                         key = 'dens_diff_feat_sel',
                                         help = '''Select the feature to split the UMAP by.''',
                                         on_change=check_feature_approval_callback)

                            st.selectbox('Values for False Condition', key = 'feature_value_fals',
                                         options = st.session_state.clus_diff_vals_fals)
                            st.number_input('Number of Clusters for False Condition', min_value = 1, max_value = 10, value = 3, step = 1, key = 'num_clus_0')
                            if st.session_state.elbow_fig_0 is not None:
                                st.pyplot(st.session_state.elbow_fig_0)
                        # Perform clustering normally
                        else:
                            st.number_input('Number of K-means clusters',
                                    min_value=st.session_state.clust_minmax[0],
                                    max_value=st.session_state.clust_minmax[1],
                                    key = 'slider_clus_val')
                            if st.session_state.spatial_umap.elbow_fig is not None:
                                st.pyplot(st.session_state.spatial_umap.elbow_fig)

                    with clust_exp_col[1]:
                        if st.session_state['toggle_clust_diff'] is True:
                            st.number_input('Cutoff Percentage', min_value = 0.01, max_value = 0.99, value = 0.01, step = 0.01, key = 'dens_diff_cutoff')
                            st.selectbox('Values for True Condition', key = 'feature_value_true',
                                         options = st.session_state.clus_diff_vals_true)
                            st.number_input('Number of Clusters for True Condition', min_value = 1, max_value = 10, value = 3, step = 1, key = 'num_clus_1')
                            if st.session_state.elbow_fig_1 is not None:
                                st.pyplot(st.session_state.elbow_fig_1)
                    if st.session_state.cluster_completed:
                        st.markdown('''The within-cluster sum of squares (WCSS) is a measure of the
                                        variability of the observations within each cluster. In general,
                                        a cluster that has a small sum of squares is more compact than a
                                        cluster that has a large sum of squares. Clusters that have higher
                                        values exhibit greater variability of the observations within the
                                        cluster.''')

        # UMAP Density Preview
        with npf_cols[1]:

            # As long as the UMAP is completed, perform the density difference analysis
            if st.session_state.umap_completed:

                # Display visualization off the full UMAP
                with st.columns(3)[1]:
                    st.pyplot(fig=st.session_state.UMAPFig)

                if st.session_state.cluster_completed_diff and st.session_state['toggle_clust_diff']:
                    if st.session_state['appro_feat']:

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
                            st.pyplot(fig=st.session_state.diff_clust_Fig)
                if st.session_state['appro_feat'] is False:
                    st.write('Feature must be boolean or numeric to perform density difference analysis')

    # Tab for Loading Previous UMAP Results
    with nei_pro_tabs[1]:
        st.write('Feature coming soon!')
        # st.write('Checkpoint file: neighborhood_profiles_checkpoint.pkl')
        # st.button('Load checkpointed UMAP results', on_click=load_neipro_struct)
        add_vertical_space(19)

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
    viz_cols = st.columns(2)

    # Scatterplot Figure Column
    with viz_cols[0]:
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
            st.button('←', on_click=slide_id_prog_left_callback, disabled=st.session_state.prog_left_disabled)
        with image_prog_col[2]:
            add_vertical_space(2)
            st.button('→', on_click=slide_id_prog_right_callback, disabled=st.session_state.prog_right_disabled)
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
    with viz_cols[1]:
        st.header('Neighborhood Profiles')
        with st.expander('Neighborhood Profile Options'):
            nei_sett_col = st.columns([1, 1, 1, 1])
            with nei_sett_col[0]:
                st.toggle('Manual Y-axis scaling',
                          key = 'toggle_manual_y_axis_scaling_main')
                st.toggle('Hide "Other" Phenotype', value = False, key = 'toggle_hide_other')
            with nei_sett_col[1]:
                st.number_input('Y-axis Min', key = 'y_axis_min_main',
                                value = 0.1, step = 0.01,)
                st.toggle('Hide "No Cluster" Neighborhood Profile', value = False, key = 'toggle_hide_no_cluster')
            with nei_sett_col[2]:
                st.number_input('Y-axis Max', key = 'y_axis_max_main',
                                value = 10000, step = 10,)
            with nei_sett_col[3]:
                st.checkbox('Log Scale', key = 'nei_pro_toggle_log_scale', value = True)
            st.toggle('Subset Neighbourhood Profiles by Feature', value = False,
                      key = 'toggle_NeiPro_filter_feat', disabled=not st.session_state.umap_completed)
            if st.session_state['toggle_NeiPro_filter_feat']:
                nei_feat_filt_col = st.columns([2,2])
                with nei_feat_filt_col[0]:
                    st.selectbox('Feature', options = st.session_state.umapOutcomes,
                                 key='NeiPro_filter_feat', on_change=reset_neipro_feature_values)
                with nei_feat_filt_col[1]:
                    selected_feat = st.session_state['NeiPro_filter_feat']
                    unique_values = st.session_state.spatial_umap.df_umap[selected_feat].unique()
                    st.selectbox('Value', options = unique_values, key='NeiPro_filter_value')


        # If the spatial-umap is completed...
        if 'spatial_umap' in st.session_state:
            # List of Clusters to display
            list_clusters = list(st.session_state.spatial_umap.dens_df_mean['clust_label'].unique())
            if st.session_state['toggle_hide_no_cluster']:
                list_clusters.remove('No Cluster')

            cluster_sel_col = st.columns([3, 1])
            # Compare Clusters Toggle
            with cluster_sel_col[1]:
                add_vertical_space(2)
                st.toggle('Compare Cluster Neighborhoods', value = False, key = 'toggle_compare_clusters')
                if st.session_state['toggle_compare_clusters']:
                    st.radio('Compare as:', ('Ratio', 'Difference'), index = 0, key = 'compare_clusters_as', horizontal=True)

            # Cluster Select Widgets
            with cluster_sel_col[0]:
                st.selectbox('Select a cluster to view', list_clusters, key='sel_npf_fig')
                if st.session_state['toggle_compare_clusters']:
                    st.selectbox('Select a cluster to compare', list_clusters, key='sel_npf_fig2')
                else:
                    st.session_state['sel_npf_fig2'] = None

            if st.session_state.cluster_completed and st.session_state.appro_feat:

                # Draw the Neighborhood Profile
                npf_fig, ax = bpl.draw_scatter_fig(figsize=(14, 16))

                bpl.draw_neigh_profile_fig(st.session_state.spatial_umap,
                                           ax = ax,
                                           sel_clus = st.session_state['sel_npf_fig'],
                                           cmp_clus = st.session_state['sel_npf_fig2'],
                                           cmp_style=st.session_state['compare_clusters_as'],
                                           hide_other = st.session_state['toggle_hide_other'],
                                           hide_no_cluster = st.session_state['toggle_hide_no_cluster'])

                if st.session_state['nei_pro_toggle_log_scale']:
                    ax.set_yscale('log')

                if st.session_state['toggle_manual_y_axis_scaling_main']:
                    ax.set_ylim(st.session_state['y_axis_min_main'], st.session_state['y_axis_max_main'])

                if st.session_state['sel_npf_fig'] == st.session_state['sel_npf_fig2']:
                    st.markdown('## Please choose two different clusters to compare')
                else:
                    # Display the Neighborhood Profile
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

    # Drawing the subplots of Neighborhood Profiles per cluster combinations
    if st.session_state['appro_feat'] and st.session_state.cluster_completed_diff:

        supp_neipro_col = st.columns([4, 2])
        with supp_neipro_col[0]:
            with st.expander('Neighborhood Profile Subplots Settings', expanded = False):
                st.toggle('Manual Y-axis scaling',
                          key = 'toggle_manual_y_axis_scaling_supplemental')

                st.write('Settings for Individual Cluster Plots')
                neipro_exp_full = st.columns([3, 3, 1])
                with neipro_exp_full[0]:
                    st.number_input('Y-axis Min', key = 'y_axis_min_supplemental',
                                    value = 0.1, step = 0.001,)
                with neipro_exp_full[1]:
                    st.number_input('Y-axis Max', key = 'y_axis_max_supplemental',
                                    value = 10000, step = 10,)
                with neipro_exp_full[2]:
                    st.checkbox('Log Scale', key = 'log_scale_supplemental', value = True)

                st.write('Settings for Individual Cluster Ratios')
                neipro_exp_full = st.columns([3, 3, 1])
                with neipro_exp_full[0]:
                    st.number_input('Y-axis Min', key = 'y_axis_min_indiratio_supplemental',
                                    value = 0.1, step = 0.001,)
                with neipro_exp_full[1]:
                    st.number_input('Y-axis Max', key = 'y_axis_max_indiratio_supplemental',
                                    value = 10, step = 1,)
                with neipro_exp_full[2]:
                    st.checkbox('Log Scale', key = 'log_scale_indiratio_supplemental', value = True)

                st.write('Settings for Aggregated Cluster Ratios')
                neipro_exp_full = st.columns([3, 3, 1])
                with neipro_exp_full[0]:
                    st.number_input('Y-axis Min', key = 'y_axis_min_aggratio_supplemental',
                                    value = 0.1, step = 0.001,)
                with neipro_exp_full[1]:
                    st.number_input('Y-axis Max', key = 'y_axis_max_aggratio_supplemental',
                                    value = 10, step = 1,)
                with neipro_exp_full[2]:
                    st.checkbox('Log Scale', key = 'log_scale_aggratio_supplemental', value = True)

        npf_fig_big = plt.figure(figsize=(16, 45), facecolor = '#0E1117')

        title_supp = [f'DATASET: {st.session_state.datafile}',
                 f'FEATURE: {st.session_state.dens_diff_feat_sel}',
                 f'FALSE Val: {st.session_state.feature_value_fals}, TRUE Val: {st.session_state.feature_value_true}',]

        list_figures = [['Average False', None, 'Individual Cluster Plots'],
                        ['Average True', None, 'Individual Cluster Plots'],
                        ['Average False', 'Average True', 'Aggregate Cluster Ratios'],
                        ['False Cluster 1', None, 'Individual Cluster Plots'],
                        ['False Cluster 2', None, 'Individual Cluster Plots'],
                        ['False Cluster 3', None, 'Individual Cluster Plots'],
                        ['True Cluster 1', None, 'Individual Cluster Plots'],
                        ['True Cluster 2', None, 'Individual Cluster Plots'],
                        ['True Cluster 3', None, 'Individual Cluster Plots'],
                        ['False Cluster 1', 'True Cluster 1', 'Individual Cluster Ratios'],
                        ['False Cluster 2', 'True Cluster 1', 'Individual Cluster Ratios'],
                        ['False Cluster 3', 'True Cluster 1', 'Individual Cluster Ratios'],
                        ['False Cluster 1', 'True Cluster 2', 'Individual Cluster Ratios'],
                        ['False Cluster 2', 'True Cluster 2', 'Individual Cluster Ratios'],
                        ['False Cluster 3', 'True Cluster 2', 'Individual Cluster Ratios'],
                        ['False Cluster 1', 'Average True', 'Aggregate Cluster Ratios'],
                        ['False Cluster 2', 'Average True', 'Aggregate Cluster Ratios'],
                        ['False Cluster 3', 'Average True', 'Aggregate Cluster Ratios'],
                        ['True Cluster 1', 'Average False', 'Aggregate Cluster Ratios'],
                        ['True Cluster 2', 'Average False', 'Aggregate Cluster Ratios'],
                        ['True Cluster 3', 'Average False', 'Aggregate Cluster Ratios'],
                        ]

        num_figs = len(list_figures)
        num_cols = 3
        num_rows = np.ceil(num_figs/3).astype(int)
        for ii, cluster in enumerate(list_figures):
            axii = npf_fig_big.add_subplot(num_rows, 3, ii+1, facecolor = '#0E1117')

            if ii == ((num_rows*num_cols)-3):
                legend_flag = True
            else:
                legend_flag = False

            bpl.draw_neigh_profile_fig(st.session_state.spatial_umap,
                                       ax = axii,
                                       sel_clus = cluster[0],
                                       cmp_clus = cluster[1],
                                       cmp_style = 'Ratio',
                                       hide_other = st.session_state['toggle_hide_other'],
                                       hide_no_cluster = st.session_state['toggle_hide_no_cluster'],
                                       legend_flag = legend_flag)

            if st.session_state['toggle_manual_y_axis_scaling_supplemental']:
                if cluster[2] == 'Individual Cluster Plots':
                    axii.set_ylim([st.session_state['y_axis_min_supplemental'],
                                   st.session_state['y_axis_max_supplemental']])
                    if st.session_state['log_scale_supplemental']:
                        axii.set_yscale('log')
                elif cluster[2] == 'Individual Cluster Ratios':
                    axii.set_ylim([st.session_state['y_axis_min_indiratio_supplemental'],
                                   st.session_state['y_axis_max_indiratio_supplemental']])
                    if st.session_state['log_scale_indiratio_supplemental']:
                        axii.set_yscale('log')
                elif cluster[2] == 'Aggregate Cluster Ratios':
                    axii.set_ylim([st.session_state['y_axis_min_aggratio_supplemental'],
                                   st.session_state['y_axis_max_aggratio_supplemental']])
                    if st.session_state['log_scale_aggratio_supplemental']:
                        axii.set_yscale('log')
            else:
                axii.set_yscale('log')

        plot_title = ''
        for i in title_supp:
            plot_title = plot_title + i + '\n'

        # Super title with more information
        npf_fig_big.suptitle(t = plot_title, color = '#FAFAFA', x = 0.1, y = 0.9,
                             horizontalalignment = 'left',
                             verticalalignment = 'top',)

        st.pyplot(fig=npf_fig_big)

if __name__ == '__main__':
    main()
