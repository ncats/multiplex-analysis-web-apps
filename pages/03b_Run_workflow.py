# Import relevant libraries
import streamlit as st
import os
import numpy as np
import time_cell_interaction_lib as tci  # import the TIME library stored in time_cell_interaction_lib.py
import time
import streamlit_utils
import app_top_of_page as top
import streamlit_dataframe_editor as sde

def main():
    '''
    Main function for running the page
    '''

    # If 'input_dataset' isn't in the session state, print an error message and return
    if 'input_dataset' not in st.session_state:
        st.error('An input dataset has not yet been opened. Please do so using the "Open File" page in the sidebar.')
        return

    # Set default widget values
    streamlit_utils.assign_default_values_in_session_state('num_workers', 7)
    streamlit_utils.assign_default_values_in_session_state('use_multiprocessing', True)
    block_names = ['Instantiate TIME class', 'Plot ROIs', 'Calculate P values', 'Check metrics, impose plotting settings, and convert to numpy format', 'Plot density heatmaps per ROI', 'Plot ROI outlines individually on the whole slides', 'Average density P values over ROIs for each slide', 'Plot all ROI outlines on the whole slides', 'Average density P values over ROIs for each annotation region type', 'Plot density P values for each ROI over slide spatial plot']
    component_bool_defaults = [True, True, True, True, True, True, True, True, False, False]
    component_checkbox_disabled = [False, False, False, False, False, False, False, False, False, False]
    component_checkbox_disabled[-2] = st.session_state['annotation_coordinate_units_is_disabled']
    component_help = [None, None, None, None, None, None, None, None, None, None]
    if component_checkbox_disabled[-2]:
        st.session_state['Average density P values over ROIs for each annotation region type'] = False
        component_help[-2] = 'Note at least one annotation file on the previous page must be selected'
    for iblock_name, block_name in enumerate(block_names):
        streamlit_utils.assign_default_values_in_session_state(block_name, component_bool_defaults[iblock_name])

    # Use only a single column
    col_settings, col_output = st.columns(2)
    with col_settings:

        # Section title
        st.subheader('Workflow settings')

        # Choose which components of the tool to run
        st.write('**Select components of the tool to run:**')
        workflow_bools = []
        for iblock, block_name in enumerate(block_names):
            st.checkbox(block_name, key=block_name, disabled=component_checkbox_disabled[iblock], help=component_help[iblock])
            workflow_bools.append(st.session_state[block_name])
        print('Workflow blocks that will be run:')
        for iblock_to_run in np.where(workflow_bools)[0]:
            print('  ({}): {}'.format(2 ** iblock_to_run, block_names[iblock_to_run]))

        # Separate section for job execution parameters
        st.write('**Select job execution parameters:**')

        # Determine whether we should employ threading
        use_multiprocessing = st.checkbox('Should we use multiple logical CPUs to speed up the calculations?', key='use_multiprocessing')

        # Get the number of threads to use for the calculations
        num_workers = st.number_input('Select number of threads for calculations:', min_value=1, max_value=os.cpu_count(), step=1, key='num_workers', disabled=(not use_multiprocessing))

    with col_output:

        # Section title
        st.subheader('Workflow execution')

        if st.button('Run workflow'):

            # Read in the dataset object
            if 'dataset_obj' in st.session_state:
                dataset_obj = st.session_state['dataset_obj']
            else:
                st.error('Dataset not yet loaded; please click on the "Load dataset" button', icon="🚨")
                st.stop()  # may only be available in the latest version of Streamlit!! (as of 5/11/23, version 1.22.0)

            # Instantiate the TIME class, which mainly loads the main datafile into a Pandas dataframe (or reads in a simulated one in the case of simulated data) and performs some preprocessing on it. Note that this creates slices.df_roi_plotting_data. This is *not* a parallel function
            iblock = 0
            if workflow_bools[iblock]:
                print('**** {}... ****'.format(block_names[iblock]))
                st.write('**:sparkles: {}...**'.format(block_names[iblock]))
                slices = tci.TIMECellInteraction(
                    dataset_obj,
                    project_dir=os.path.realpath(os.path.join(os.getcwd(), '..')),
                    allow_compound_species=st.session_state['sit__used_settings']['analysis']['allow_compound_species'],
                    thickness_new=st.session_state['sit__used_settings']['analysis']['thickness'],
                    use_analytical_significance=st.session_state['sit__used_settings']['analysis']['use_analytical_significance'],
                    n_neighs=st.session_state['sit__used_settings']['analysis']['n_neighs'],
                    radius_instead_of_knn=st.session_state['sit__used_settings']['analysis']['radius_instead_of_knn']
                )

            # Plot every ROI using slices.df_roi_plotting_data in parallel
            iblock = 1
            if workflow_bools[iblock]:
                print('**** {}... ****'.format(block_names[iblock]))
                st.write('**:sparkles: {}...**'.format(block_names[iblock]))
                start_time = time.time()
                slices.plot_rois(
                    nworkers=num_workers,
                    use_multiprocessing=use_multiprocessing
                    )
                benchmarking_message = '...plot_rois() took {} seconds using {} CPU(s) {} hyperthreading'.format(int(np.round(time.time() - start_time)), (num_workers if use_multiprocessing else 1), ('WITH' if use_multiprocessing else 'WITHOUT'))
                print('')
                print('BENCHMARKING: {}'.format(benchmarking_message))
                print('')
                st.write(benchmarking_message)

            # Calculate the P values from the coordinates of the species in every ROI in every slide. Note that this creates slices.df_density_pvals (the flattened metrics dataframe). This is a parallel function
            iblock = 2
            if workflow_bools[iblock]:
                print('**** {}... ****'.format(block_names[iblock]))
                st.write('**:sparkles: {}...**'.format(block_names[iblock]))
                start_time = time.time()
                slices.calculate_metrics(
                    nworkers=num_workers,
                    use_multiprocessing=use_multiprocessing
                    )
                benchmarking_message = '...calculate_metrics() took {} seconds using {} CPU(s) {} hyperthreading'.format(int(np.round(time.time() - start_time)), (num_workers if use_multiprocessing else 1), ('WITH' if use_multiprocessing else 'WITHOUT'))
                print('')
                print('BENCHMARKING: {}'.format(benchmarking_message))
                print('')
                st.write(benchmarking_message)

            # Check the metrics, impose plotting requests, and convert arrays to numpy format; this creates slices.df_density_pvals_arrays and is not a parallel function
            iblock = 3
            if workflow_bools[iblock]:
                print('**** {}... ****'.format(block_names[iblock]))
                st.write('**:sparkles: {}...**'.format(block_names[iblock]))
                slices.check_and_prepare_metrics_for_plotting(
                    log_pval_range=st.session_state['sit__used_settings']['plotting']['log_pval_range'],
                    num_valid_centers_minimum=st.session_state['sit__used_settings']['plotting']['num_valid_centers_minimum']
                )

            # Plot every density heatmap using slices.df_density_pvals in parallel
            iblock = 4
            if workflow_bools[iblock]:
                print('**** {}... ****'.format(block_names[iblock]))
                st.write('**:sparkles: {}...**'.format(block_names[iblock]))
                start_time = time.time()
                slices.plot_dens_pvals_per_roi(
                    nworkers=num_workers,
                    use_multiprocessing=use_multiprocessing
                    )
                benchmarking_message = '...plot_dens_pvals_per_roi() took {} seconds using {} CPU(s) {} hyperthreading'.format(int(np.round(time.time() - start_time)), (num_workers if use_multiprocessing else 1), ('WITH' if use_multiprocessing else 'WITHOUT'))
                print('')
                print('BENCHMARKING: {}'.format(benchmarking_message))
                print('')
                st.write(benchmarking_message)

            # Plot ROI outlines individually for each ROI over the whole slides
            iblock = 5
            if workflow_bools[iblock]:
                print('**** {}... ****'.format(block_names[iblock]))
                st.write('**:sparkles: {}...**'.format(block_names[iblock]))
                start_time = time.time()
                slices.plot_outline_for_single_roi_on_whole_slide(
                    nworkers=num_workers,
                    use_multiprocessing=use_multiprocessing
                )
                benchmarking_message = '...plot_outline_for_single_roi_on_whole_slide() took {} seconds using {} CPU(s) {} hyperthreading'.format(int(np.round(time.time() - start_time)), (num_workers if use_multiprocessing else 1), ('WITH' if use_multiprocessing else 'WITHOUT'))
                print('')
                print('BENCHMARKING: {}'.format(benchmarking_message))
                print('')
                st.write(benchmarking_message)

            # Average the P values for each slide over the corresponding ROIs containing valid data. Note this creates slices.df_log_dens_pvals_arr_per_slide
            # Can run something like this afterward: slices.df_log_dens_pvals_arr_per_slide.drop('log_dens_pvals_arr', axis='columns').to_excel('/home/weismanal/transfer/slide_response_variables.xlsx')
            # Can also impute missing data in slices.df_log_dens_pvals_arr_per_slide, negate the arrays, and normalize to a 0-1 range as input for stats/ML
            iblock = 6
            if workflow_bools[iblock]:
                print('**** {}... ****'.format(block_names[iblock]))
                st.write('**:sparkles: {}...**'.format(block_names[iblock]))
                slices.average_dens_pvals_over_rois_for_each_slide(
                    weight_rois_by_num_valid_centers=st.session_state['sit__used_settings']['plotting']['weight_rois_by_num_valid_centers'],
                    input_datafile=st.session_state['input_metadata']['datafile_path']  # this is just needed to get the input data filename to save to disk along with the df_log_dens_pvals_arr_per_slide for later read-in by the correlation analyzer
                )

            # Plot the ROIs on each slide
            iblock = 7
            if workflow_bools[iblock]:
                print('**** {}... ****'.format(block_names[iblock]))
                st.write('**:sparkles: {}...**'.format(block_names[iblock]))
                slices.plot_whole_slide_patches()

            # Average the density P value data over all ROIs for each annotation region type and plot the final and intermediate results
            iblock = 8
            if workflow_bools[iblock]:
                print('**** {}... ****'.format(block_names[iblock]))
                st.write('**:sparkles: {}...**'.format(block_names[iblock]))
                slices.average_over_rois_per_annotation_region(
                    annotations_csv_files=st.session_state['sit__used_settings']['annotation']['csv_files'],
                    phenotyping_method=st.session_state['sit__used_settings']['phenotyping']['method'],
                    phenotype_identification_file=st.session_state['sit__used_settings']['dataset']['phenotype_identification_tsv_file'],
                    annotation_coord_units_in_microns=st.session_state['sit__used_settings']['annotation']['annotation_coord_units_in_microns'],
                    annotation_microns_per_integer_unit=st.session_state['sit__used_settings']['annotation']['annotation_microns_per_integer_unit'],
                    settings__analysis__thickness=st.session_state['sit__used_settings']['analysis']['thickness'],
                    min_log_pval_for_plotting=st.session_state['sit__used_settings']['plotting']['min_log_pval']
                    )

            # Plot the density P values for each ROI over spatial plots of the slides; this probably overwrites existing plots, but it doesn't take long to regenerate them
            iblock = 9
            if workflow_bools[iblock]:
                print('**** {}... ****'.format(block_names[iblock]))
                st.write('**:sparkles: {}...**'.format(block_names[iblock]))
                start_time = time.time()
                slices.plot_density_pvals_over_slides(
                    nworkers=num_workers,
                    use_multiprocessing=use_multiprocessing
                )
                benchmarking_message = '...plot_density_pvals_over_slides() took {} seconds using {} CPU(s) {} hyperthreading'.format(int(np.round(time.time() - start_time)), (num_workers if use_multiprocessing else 1), ('WITH' if use_multiprocessing else 'WITHOUT'))
                print('')
                print('BENCHMARKING: {}'.format(benchmarking_message))
                print('')
                st.write(benchmarking_message)

# Call the main function
if __name__ == '__main__':

    # Set a wide layout and display the page heading
    st.set_page_config(layout="wide")
    st.title('Run workflow')

    # Run streamlit-dataframe-editor library initialization tasks at the top of the page
    st.session_state = sde.initialize_session_state(st.session_state)

    # Run Top of Page (TOP) functions
    st.session_state = top.top_of_page_reqs(st.session_state)

    main()

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)
