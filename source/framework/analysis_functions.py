import time
import math
import os
import basic_phenotyper_lib as bpl  # Useful functions for phenotyping collections of cells
from copy import copy
import numpy as np
import pandas as pd

def run_analysis_job(function_name, inputs, job_dir):
    try:
        outputs_dir = os.path.join(job_dir, "outputs")  # This demonstrates that for a potentially asynchronous job that generates files, you should place the results in /tmp/multiplex_analysis_web_apps/job_data/<JOB_ID>/outputs specifically so the results are stored together with the worker output results in memory.
        if function_name == "find_primes_up_to":
            function_to_run = find_primes_up_to
        elif function_name == "init_spatial_umap":
            function_to_run = init_spatial_umap
        elif function_name == "apply_umap":
            function_to_run = apply_umap
        elif function_name == "set_clusters":
            function_to_run = set_clusters
        elif function_name == "clust_umap_dens_diff":
            function_to_run = clust_umap_dens_diff
        elif function_name == "run_sit_workflow":
            function_to_run = run_sit_workflow
        outputs = function_to_run(**inputs, results_topdir=outputs_dir)
        return outputs
    except Exception as e:
        print(f"Error occurred while running analysis job {function_name}: {e}")
        return None


def find_primes_up_to(limit, results_subdir, results_topdir):
    """
    Find all prime numbers up to a given limit using trial division.
    Returns the list of primes and timing information.

    On my laptop this takes about 6-8 seconds: primes, duration = find_primes_up_to(4000000).
    """
    start_time = time.time()

    if limit < 2:
        return [], 0

    primes = []

    for num in range(2, limit + 1):
        is_prime = True

        # Check if num is prime by testing divisibility
        for i in range(2, int(math.sqrt(num)) + 1):
            if num % i == 0:
                is_prime = False
                break

        if is_prime:
            primes.append(num)

    end_time = time.time()
    duration = end_time - start_time

    results_dir = os.path.join(results_topdir, results_subdir)

    # Create results directory if it doesn't exist.
    os.makedirs(results_dir, exist_ok=True)

    # Save results to text files
    with open(os.path.join(results_dir, "primes.txt"), "w") as f:
        f.write(f"Found {len(primes)} primes up to {limit} in {duration:.2f} seconds\n")
        f.write(f"First 10 primes: {primes[:10]}\n")
        f.write(f"Last 10 primes: {primes[-10:]}\n")

    return {"primes": primes, "duration": duration}


def init_spatial_umap(calc_unique_areas_toggle, area_filter_per, df, marker_multi_sel, phenoOrder, datafile_min_img_size, cpu_pool_size, results_topdir):
    '''
    Initalizing the spatial_umap object
    '''

    # Andrew moved this to its own button in Neighborhood_Profiles.py; search for "Initialize Neighborhood Profiles" there.
    # Reset the settings required for Neighborhood Analysis
    # st.session_state = ndl.reset_neigh_profile_settings(st.session_state)

    if not calc_unique_areas_toggle:
        area_filter = 0
    else:
        area_filter = area_filter_per

    # bc.startTimer()
    # with st.spinner('Calculating Cell Counts and Areas', show_time=True):
    spatial_umap = bpl.setup_Spatial_UMAP(df = df,
                                                            marker_names = marker_multi_sel,
                                                            pheno_order = phenoOrder,
                                                            smallest_image_size = datafile_min_img_size)

    spatial_umap = bpl.perform_density_calc(spatial_umap,
                                                            #  bc,
                                                                calc_unique_areas_toggle,
                                                                cpu_pool_size,
                                                                area_threshold = area_filter)

    # Record time elapsed
    # bc.set_value_df('time_to_run_counts', bc.elapsedTime())

    return {"spatial_umap": spatial_umap, "density_completed": True}

    # Save checkpoint for Neighborhood Profile structure
    # save_neipro_struct()


def get_spatialUMAP(spatial_umap, umap_subset_per_fit, umap_subset_toggle, umap_subset_per):
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
    # bc.startTimer()
    # print('Fitting Model')
    spatial_umap.umap_fit = spatial_umap.cells.loc[spatial_umap.cells['umap_train'].values, ['UMAP_1_20230327_152849', 'UMAP_2_20230327_152849']].values.reshape((spatial_umap.cells['umap_train'].sum(), -1))
    # bc.printElapsedTime(f'      Fitting {np.sum(spatial_umap.cells["umap_train"] == 1)} points to a model')

    # Transform test cells based on fitted model
    # bc.startTimer()
    # print('Transforming Data')
    spatial_umap.umap_test = spatial_umap.cells.loc[spatial_umap.cells['umap_test'].values, ['UMAP_1_20230327_152849', 'UMAP_2_20230327_152849']].values.reshape((spatial_umap.cells['umap_test'].sum(), -1))
    # bc.printElapsedTime(f'      Transforming {np.sum(spatial_umap.cells["umap_test"] == 1)} points with the model')

    spatial_umap.umap_completed = True
    
    # import pickle
    # with open('../Edits/spatial_umap_original_precomp.pkl', 'wb') as f:
    #     pickle.dump(spatial_umap, f)

    return spatial_umap

def apply_umap(spatial_umap, umap_subset_per_fit, umap_subset_toggle, 
               umap_subset_per, load_generated_umap_toggle, results_topdir):
    '''
    Call back function for applying the UMAP functions
    '''

    #bc.startTimer()
    # if toggle for loading pre-generated UMAP is selected extract UMAP from file, works only with a specific dataset
    if load_generated_umap_toggle:
        spatial_umap = get_spatialUMAP(spatial_umap,
                                                        #bc,
                                                        umap_subset_per_fit,
                                                        umap_subset_toggle,
                                                        umap_subset_per)
    else:
        spatial_umap = bpl.perform_spatialUMAP(spatial_umap,
                                                                #bc,
                                                                umap_subset_per_fit,
                                                                umap_subset_toggle,
                                                                umap_subset_per)

    return {"spatial_umap": spatial_umap, "umap_completed": True}


def set_clusters(spatial_umap, slider_clus_val, clust_minmax, results_topdir):
    spatial_umap = bpl.umap_clustering(spatial_umap = spatial_umap,
                                                                n_clusters = slider_clus_val,
                                                                clust_minmax = clust_minmax,
                                                                cpu_pool_size = 3)
    spatial_umap.mean_measures()
    
    return {"spatial_umap": spatial_umap, "cluster_completed": True, 
            "appro_feat": True, "cluster_completed_diff": False}

def clust_umap_dens_diff(udp_full, dens_diff_feat_sel, feature_value_fals, 
                         feature_value_true, clust_diff_vals_code, npf,
                         dens_diff_cutoff, 
                         num_clus_0, num_clus_1, clust_minmax, spatial_umap, results_topdir):
    # Import here to avoid circular import
    from neighborhood_profiles import UMAPDensityProcessing

    print(udp_full, flush=True)
    print(dens_diff_feat_sel, flush=True)
    print(feature_value_fals, flush=True)
    print(feature_value_true, flush=True)
    print(clust_diff_vals_code, flush=True)
    print(npf, flush=True)
    print(dens_diff_cutoff, flush=True)
    print(num_clus_0, flush=True)
    print(num_clus_1, flush=True)
    print(clust_minmax, flush=True)
    print(spatial_umap, flush=True)
    
    # Split the UMAP by the selected values of the feature
    split_dict_full = udp_full.split_df_by_feature(dens_diff_feat_sel,
                                                    feature_value_fals,
                                                    feature_value_true,
                                                    clust_diff_vals_code)

    # Perform Density Calculations for each Condition
    udp_fals = UMAPDensityProcessing(npf, split_dict_full['df_umap_fals'], xx=udp_full.xx, yy=udp_full.yy)
    udp_true = UMAPDensityProcessing(npf, split_dict_full['df_umap_true'], xx=udp_full.xx, yy=udp_full.yy)

    ## Copy over
    udp_diff = copy(udp_fals)
    ## Perform difference calculation
    udp_diff.dens_mat = np.log10(udp_fals.dens_mat) - np.log10(udp_true.dens_mat)
    ## Rerun the min/max calcs
    udp_diff.umap_summary_stats()
    ## Set Feature Labels
    udp_fals.set_feature_label(dens_diff_feat_sel, split_dict_full['fals_msg'])
    udp_true.set_feature_label(dens_diff_feat_sel, split_dict_full['true_msg'])
    udp_diff.set_feature_label(dens_diff_feat_sel, 'Difference')

    # Draw UMAPS
    UMAPFig_fals = udp_fals.UMAPdraw_density()
    UMAPFig_true = udp_true.UMAPdraw_density()
    UMAPFig_diff = udp_diff.UMAPdraw_density(diff=True)

        # Assign Masking and plot
    # Assign Masking and plot
    udp_mask = copy(udp_diff)
    udp_mask.filter_density_matrix(dens_diff_cutoff, udp_full.empty_bin_ind)
    udp_mask.set_feature_label(dens_diff_feat_sel, f'Difference- Masked, \ncutoff = {dens_diff_cutoff}')
    UMAPFig_mask = udp_mask.UMAPdraw_density(diff=True)

    # Perform Clustering
    udp_clus = copy(udp_mask)
    udp_clus.perform_clustering(dens_mat_cmp=udp_mask.dens_mat,
                                num_clus_0=num_clus_0,
                                num_clus_1=num_clus_1,
                                clust_minmax=clust_minmax,
                                cpu_pool_size=3)
    udp_clus.set_feature_label(dens_diff_feat_sel, f'Clusters, False-{num_clus_0}, True-{num_clus_1}')
    UMAPFig_clus = udp_clus.UMAPdraw_density(diff=True, legendtype='legend')
    cluster_dict = udp_clus.cluster_dict
    palette_dict = udp_clus.palette_dict
    elbow_fig_0 = udp_clus.elbow_fig_0
    elbow_fig_1 = udp_clus.elbow_fig_1

    # Add cluster label column to cells dataframe
    spatial_umap.df_umap.loc[:, 'clust_label'] = 'No Cluster'
    spatial_umap.df_umap.loc[:, 'cluster'] = 'No Cluster'
    spatial_umap.df_umap.loc[:, 'Cluster'] = 'No Cluster'

    for key, val in cluster_dict.items():
        if key != 0:
            bin_clust = np.argwhere(udp_clus.dens_mat == key)
            bin_clust = bin_clust[:, [1, 0]] # Swapping columns to by y, x
            bin_clust = [tuple(x) for x in bin_clust]

            significant_groups = udp_full.bin_indices_df_group[udp_full.bin_indices_df_group.set_index(['indx', 'indy']).index.isin(bin_clust)]

            umap_ind = significant_groups.index.values
            spatial_umap.df_umap.loc[umap_ind, 'clust_label'] = val
            spatial_umap.df_umap.loc[umap_ind, 'cluster'] = val
            spatial_umap.df_umap.loc[umap_ind, 'Cluster'] = val

    # Benchmark how long it took to untangle indicies
    #bc.printElapsedTime('Untangling bin indicies with UMAP indicies', split = True)

    # After assigning cluster labels, perform mean calculations
    spatial_umap.mean_measures()
    #bc.printElapsedTime('Performing Mean Measures', split = True)

    # Average Left condition and Average Right Condition
    dens_df_fals = spatial_umap.dens_df_mean.loc[spatial_umap.dens_df_mean['clust_label'].str.contains('Left'), :]
    dens_df_true = spatial_umap.dens_df_mean.loc[spatial_umap.dens_df_mean['clust_label'].str.contains('Right'), :]

    dens_df_fals['clust_label'] = 'Average Left'
    dens_df_mean_fals = dens_df_fals.groupby(['clust_label', 'phenotype', 'dist_bin'], as_index=False).mean()

    dens_df_true['clust_label'] = 'Average Right'
    dens_df_mean_true = dens_df_true.groupby(['clust_label', 'phenotype', 'dist_bin'], as_index=False).mean()

    spatial_umap.dens_df_mean = pd.concat([spatial_umap.dens_df_mean, dens_df_mean_fals, dens_df_mean_true], axis=0)

    return {"spatial_umap": spatial_umap, "cluster_completed_diff": True,
            "UMAPFig_fals": UMAPFig_fals,
            "UMAPFig_true": UMAPFig_true,
            "UMAPFig_diff": UMAPFig_diff,
            "UMAPFig_mask": UMAPFig_mask,
            "cluster_dict": cluster_dict,
            "palette_dict": palette_dict,
            "elbow_fig_0": elbow_fig_0,
            "elbow_fig_1": elbow_fig_1,
            "cluster_completed": True,
            "udp_full": udp_full
            }

def run_sit_workflow(dataset_obj, project_dir, allow_compound_species,
                     thickness_new, use_analytical_significance, n_neighs,
                     radius_instead_of_knn, workflow_bools, num_workers, 
                     block_names, use_multiprocessing,log_pval_range, num_valid_centers_minimum,
                     weight_rois_by_num_valid_centers, 
                     input_datafile, save_heatmap_data,
                     annotations_csv_files, phenotyping_method,
                     phenotype_identification_file,
                     annotation_coord_units_in_microns, annotation_microns_per_integer_unit,
                     settings__analysis__thickness, 
                     min_log_pval_for_plotting,
                     results_topdir):

    import time_cell_interaction_lib as tci
    slices = tci.TIMECellInteraction(
        dataset_obj,
        project_dir=results_topdir,
        allow_compound_species=allow_compound_species,
        thickness_new=thickness_new,
        use_analytical_significance=use_analytical_significance,
        n_neighs=n_neighs,
        radius_instead_of_knn=radius_instead_of_knn
    )

        # Plot every ROI using slices.df_roi_plotting_data in parallel
    iblock = 1
    if workflow_bools[iblock]:
        print('**** {}... ****'.format(block_names[iblock]))
        #st.write('**:sparkles: {}...**'.format(block_names[iblock]))
        start_time = time.time()
        slices.plot_rois(
            nworkers=num_workers,
            use_multiprocessing=use_multiprocessing
            )
        benchmarking_message = '...plot_rois() took {} seconds using {} CPU(s) {} hyperthreading'.format(int(np.round(time.time() - start_time)), (num_workers if use_multiprocessing else 1), ('WITH' if use_multiprocessing else 'WITHOUT'))
        print('')
        print('BENCHMARKING: {}'.format(benchmarking_message))
        print('')
        #st.write(benchmarking_message)

    # Calculate the P values from the coordinates of the species in every ROI in every slide. Note that this creates slices.df_density_pvals (the flattened metrics dataframe). This is a parallel function
    iblock = 2
    if workflow_bools[iblock]:
        print('**** {}... ****'.format(block_names[iblock]))
        #st.write('**:sparkles: {}...**'.format(block_names[iblock]))
        start_time = time.time()
        slices.calculate_metrics(
            nworkers=num_workers,
            use_multiprocessing=use_multiprocessing
            )
        benchmarking_message = '...calculate_metrics() took {} seconds using {} CPU(s) {} hyperthreading'.format(int(np.round(time.time() - start_time)), (num_workers if use_multiprocessing else 1), ('WITH' if use_multiprocessing else 'WITHOUT'))
        print('')
        print('BENCHMARKING: {}'.format(benchmarking_message))
        print('')
        #st.write(benchmarking_message)

    # Check the metrics, impose plotting requests, and convert arrays to numpy format; this creates slices.df_density_pvals_arrays and is not a parallel function
    iblock = 3
    if workflow_bools[iblock]:
        print('**** {}... ****'.format(block_names[iblock]))
        #st.write('**:sparkles: {}...**'.format(block_names[iblock]))
        slices.check_and_prepare_metrics_for_plotting(
            log_pval_range=log_pval_range,
            num_valid_centers_minimum=num_valid_centers_minimum
        )

    # Plot every density heatmap using slices.df_density_pvals in parallel
    iblock = 4
    if workflow_bools[iblock]:
        print('**** {}... ****'.format(block_names[iblock]))
        #st.write('**:sparkles: {}...**'.format(block_names[iblock]))
        start_time = time.time()
        slices.plot_dens_pvals_per_roi(
            nworkers=num_workers,
            use_multiprocessing=use_multiprocessing
            )
        benchmarking_message = '...plot_dens_pvals_per_roi() took {} seconds using {} CPU(s) {} hyperthreading'.format(int(np.round(time.time() - start_time)), (num_workers if use_multiprocessing else 1), ('WITH' if use_multiprocessing else 'WITHOUT'))
        print('')
        print('BENCHMARKING: {}'.format(benchmarking_message))
        print('')
        #st.write(benchmarking_message)

    # Plot ROI outlines individually for each ROI over the whole slides
    iblock = 5
    if workflow_bools[iblock]:
        print('**** {}... ****'.format(block_names[iblock]))
        #st.write('**:sparkles: {}...**'.format(block_names[iblock]))
        start_time = time.time()
        slices.plot_outline_for_single_roi_on_whole_slide(
            nworkers=num_workers,
            use_multiprocessing=use_multiprocessing
        )
        benchmarking_message = '...plot_outline_for_single_roi_on_whole_slide() took {} seconds using {} CPU(s) {} hyperthreading'.format(int(np.round(time.time() - start_time)), (num_workers if use_multiprocessing else 1), ('WITH' if use_multiprocessing else 'WITHOUT'))
        print('')
        print('BENCHMARKING: {}'.format(benchmarking_message))
        print('')
        #st.write(benchmarking_message)

    # Average the P values for each slide over the corresponding ROIs containing valid data. Note this creates slices.df_log_dens_pvals_arr_per_slide
    # Can run something like this afterward: slices.df_log_dens_pvals_arr_per_slide.drop('log_dens_pvals_arr', axis='columns').to_excel('/home/weismanal/transfer/slide_response_variables.xlsx')
    # Can also impute missing data in slices.df_log_dens_pvals_arr_per_slide, negate the arrays, and normalize to a 0-1 range as input for stats/ML
    iblock = 6
    if workflow_bools[iblock]:
        print('**** {}... ****'.format(block_names[iblock]))
        #st.write('**:sparkles: {}...**'.format(block_names[iblock]))
        slices.average_dens_pvals_over_rois_for_each_slide(
            weight_rois_by_num_valid_centers=weight_rois_by_num_valid_centers,
            input_datafile=input_datafile,  # this is just needed to get the input data filename to save to disk along with the df_log_dens_pvals_arr_per_slide for later read-in by the correlation analyzer
            save_heatmap_data=save_heatmap_data,
        )

    # Plot the ROIs on each slide
    iblock = 7
    if workflow_bools[iblock]:
        print('**** {}... ****'.format(block_names[iblock]))
        #st.write('**:sparkles: {}...**'.format(block_names[iblock]))
        slices.plot_whole_slide_patches()

    # Average the density P value data over all ROIs for each annotation region type and plot the final and intermediate results
    iblock = 8
    if workflow_bools[iblock]:
        print('**** {}... ****'.format(block_names[iblock]))
        #st.write('**:sparkles: {}...**'.format(block_names[iblock]))
        slices.average_over_rois_per_annotation_region(
            annotations_csv_files=annotations_csv_files,
            phenotyping_method=phenotyping_method,
            phenotype_identification_file=phenotype_identification_file,
            annotation_coord_units_in_microns=annotation_coord_units_in_microns,
            annotation_microns_per_integer_unit=annotation_microns_per_integer_unit,
            settings__analysis__thickness=settings__analysis__thickness,
            min_log_pval_for_plotting=min_log_pval_for_plotting
            )

    # Plot the density P values for each ROI over spatial plots of the slides; this probably overwrites existing plots, but it doesn't take long to regenerate them
    iblock = 9
    if workflow_bools[iblock]:
        print('**** {}... ****'.format(block_names[iblock]))
        #st.write('**:sparkles: {}...**'.format(block_names[iblock]))
        start_time = time.time()
        slices.plot_density_pvals_over_slides(
            nworkers=num_workers,
            use_multiprocessing=use_multiprocessing
        )
        benchmarking_message = '...plot_density_pvals_over_slides() took {} seconds using {} CPU(s) {} hyperthreading'.format(int(np.round(time.time() - start_time)), (num_workers if use_multiprocessing else 1), ('WITH' if use_multiprocessing else 'WITHOUT'))
        print('')
        print('BENCHMARKING: {}'.format(benchmarking_message))
        print('')
        #st.write(benchmarking_message)


    return {"slices": slices, "sit_workflow_completed": True}
