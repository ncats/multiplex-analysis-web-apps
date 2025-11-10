# Import relevant libraries
import streamlit as st
import os
import numpy as np
import time_cell_interaction_lib as tci  # import the TIME library stored in time_cell_interaction_lib.py
import time
import streamlit_utils
import framework.utils as framework_utils
import framework.analysis_framework as analysis_framework


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
        # This isn't actually a good fix because it's only the Squidpy enrichment that shouldn't have multiprocessing, not the entire workflow, but we need to implement that in the future
        # use_multiprocessing = st.checkbox('Should we use multiple logical CPUs to speed up the calculations?', key='use_multiprocessing', disabled=(st.session_state['settings__analysis__significance_calculation_method'] != 'Poisson (radius)'))

        # Get the number of threads to use for the calculations
        num_workers = st.number_input('Select number of threads for calculations:', min_value=1, max_value=os.cpu_count(), step=1, key='num_workers', disabled=(not use_multiprocessing))

    with col_output:

        # Section title
        st.subheader('Workflow execution')
    
    with col_output:
         # Section title
        st.subheader('Async Workflow execution')
        # Read in the dataset object
        if 'dataset_obj' in st.session_state:
            dataset_obj = st.session_state['dataset_obj']
            project_dir = os.path.realpath(os.path.join(os.getcwd(), '..'))
            # Submit the job to run the workflow asynchronously
            analysis_framework.job_submission(
                job_name="run_sit_workflow",
                inputs = dict(
                    dataset_obj=dataset_obj, 
                    project_dir=project_dir, 
                    allow_compound_species=st.session_state['sit__used_settings']['analysis']['allow_compound_species'],
                    thickness_new=st.session_state['sit__used_settings']['analysis']['thickness'],
                    use_analytical_significance=st.session_state['sit__used_settings']['analysis']['use_analytical_significance'],
                    n_neighs=st.session_state['sit__used_settings']['analysis']['n_neighs'],
                    radius_instead_of_knn=st.session_state['sit__used_settings']['analysis']['radius_instead_of_knn'],
                    workflow_bools=workflow_bools, 
                    num_workers=num_workers, 
                    block_names=block_names, 
                    use_multiprocessing=use_multiprocessing,
                    log_pval_range=st.session_state['sit__used_settings']['plotting']['log_pval_range'], 
                    num_valid_centers_minimum=st.session_state['sit__used_settings']['plotting']['num_valid_centers_minimum'],
                    weight_rois_by_num_valid_centers=st.session_state['sit__used_settings']['plotting']['weight_rois_by_num_valid_centers'], 
                    input_datafile=st.session_state['input_metadata']['datafile_path'], 
                    save_heatmap_data=st.session_state['sit__used_settings']['plotting']['save_heatmap_data'],
                    annotations_csv_files=st.session_state['sit__used_settings']['annotation']['csv_files'], 
                    phenotyping_method=st.session_state['sit__used_settings']['phenotyping']['method'],
                    phenotype_identification_file=st.session_state['sit__used_settings']['dataset']['phenotype_identification_tsv_file'],
                    annotation_coord_units_in_microns=st.session_state['sit__used_settings']['annotation']['annotation_coord_units_in_microns'], 
                    annotation_microns_per_integer_unit=st.session_state['sit__used_settings']['annotation']['annotation_microns_per_integer_unit'],
                    settings__analysis__thickness=st.session_state['sit__used_settings']['analysis']['thickness'], 
                    min_log_pval_for_plotting=st.session_state['sit__used_settings']['plotting']['min_log_pval']
                    ),
                    analysis_purpose = "run sit workflow",
                    st_key_prefix = "",
            )
            diff_clust_key = 'run_sit_workflow_results'
            if diff_clust_key not in st.session_state:
                st.warning("Density difference clustering analysis results are not yet available.")
                return
            # Set shortcuts to the job results.
            slices = st.session_state[diff_clust_key]['slices']


            
        else:
            st.error('Dataset not yet loaded; please click on the "Load dataset" button', icon="🚨")
            st.stop()  # may only be available in the latest version of Streamlit!! (as of 5/11/23, version 1.22.0)
        
        



# Call the main function
if __name__ == '__main__':
    main()
