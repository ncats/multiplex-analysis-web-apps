# Import relevant libraries
import streamlit as st
import os
import numpy as np
import pprint
import time_cell_interaction_lib as tci  # import the TIME library stored in time_cell_interaction_lib.py
import time
from pathlib import Path
from streamlit_extras.app_logo import add_logo
import streamlit_utils
from st_pages import show_pages_from_config, add_indentation
import app_top_of_page as top
import streamlit_dataframe_editor as sde

def main():

    # Set a wide layout
    st.set_page_config(layout="wide")

    # Run streamlit-dataframe-editor library initialization tasks at the top of the page
    st.session_state = sde.initialize_session_state(st.session_state)

    # Apply pages order and indentation
    add_indentation()
    show_pages_from_config()

    # Sidebar organization
    with st.sidebar:
        st.write('**:book: [Documentation](https://ncats.github.io/multiplex-analysis-web-apps)**')

    # Add logo to page
    add_logo('app_images/mawa_logo-width315.png', height=150)

    # Run Top of Page (TOP) functions
    st.session_state = top.check_for_platform(st.session_state)

    # Constant
    input_directory = os.path.join('.', 'input')
    
    # Display page heading
    st.title('Run workflow')

    if ('settings__input_datafile__filename' in st.session_state) and (st.session_state['settings__input_datafile__filename'] is not None):

        # Set default widget values
        # streamlit_utils.assign_default_values_in_session_state('settings_yaml_file', 'gmb.yml')
        streamlit_utils.assign_default_values_in_session_state('num_workers', 8)
        streamlit_utils.assign_default_values_in_session_state('use_multiprocessing', True)
        # streamlit_utils.assign_default_values_in_session_state('platform_set_up', False)  # likely get this out of here and put in a startup script!
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

        # # First thing, optionally run any platform-specific startup files; only do this once
        # if not st.session_state['platform_set_up']:
        #     # from pathlib import Path
        #     # import os
        #     sit_startup_file = os.path.join(Path.home(), ".dmap-dashboards", "sit_startup.py")
        #     if os.path.exists(sit_startup_file):
        #         with open(sit_startup_file) as f:
        #             exec(f.read())
        #     st.session_state['platform_set_up'] = True

        # Use only a single column
        col_settings, col_output = st.columns(2)
        with col_settings:

            # ---- Path checks and dataset loading --------------------------------------------------------------

            # Section title
            st.subheader('Dataset selection')

            # # Obtain a preset set of settings for the tool
            # config_file_options = [x for x in os.listdir('config/project') if x.endswith('.yml')]
            # settings_yaml_file = st.selectbox('Choose file containing project configuration settings:', options=config_file_options, key='settings_yaml_file')
            # settings_yaml_file = os.path.join('config/project', settings_yaml_file)
            # with open(settings_yaml_file, mode='rt') as file:
            #     settings = yaml.load(file, yaml.UnsafeLoader)
            # print('Imported settings for the workflow:')
            # pprint.pprint(settings, sort_dicts=False)

            # # Determine the project directory from the location of this Jupyter notebook
            project_dir = os.path.realpath(os.path.join(os.getcwd(), '..'))  # define the project directory
            # print('Project directory: {}'.format(project_dir))

            # # Initialize the boolean determining whether all the checks passed
            # all_checks_passed = True

            # # Check that the GitHub repository is correctly cloned to the local filesystem
            # tci_library_pathname = os.path.join(project_dir, 'repo', 'time_cell_interaction_lib.py')
            # repo_pathname = os.path.dirname(tci_library_pathname)
            # if os.path.exists(tci_library_pathname):  # check for existence of just one of the multiple files in the clone of the GitHub repository
            #     print('GitHub repository seems to correctly exist in the {} directory'.format(repo_pathname))
            # else:
            #     print('ERROR: GitHub repository does not seem to be cloned to the {} directory'.format(repo_pathname))
            #     st.error('GitHub repository does not seem to be cloned to the {} directory'.format(repo_pathname), icon="🚨")
            #     all_checks_passed = False

            # # Check that the results directory has been created
            # results_pathname = os.path.join(project_dir, 'results')
            # if os.path.exists(results_pathname):
            #     print('Results directory found: {}'.format(results_pathname))
            # else:
            #     print('ERROR: Results directory {} not found'.format(results_pathname))
            #     st.error('Results directory {} not found'.format(results_pathname), icon="🚨")
            #     all_checks_passed = False

            # # Check that the input datafile exists
            # input_datafile = os.path.join('.', 'input', settings['dataset']['input_datafile'])
            # if os.path.exists(input_datafile):
            #     print('Input datafile has been found: {}'.format(input_datafile))
            # else:
            #     print('ERROR: Input datafile {} not found'.format(input_datafile))
            #     st.error('Input datafile {} not found'.format(input_datafile), icon="🚨")
            #     all_checks_passed = False

            # # Print whether all checks passed
            # if all_checks_passed:
            #     print('\nEverything seems to be set up correctly')
            #     load_data_button_disabled = False
            # else:
            #     print('\nStop! At least one of the checks above failed. See the line(s) above beginning with "ERROR:" to determine what needs to be fixed')
            #     st.error('Stop! At least one of the checks above failed. See the error(s) above to learn what to address', icon="🚨")
            #     load_data_button_disabled = True

            # Map the new settings values to the old settings formats which are used in the workflow options below
            orig_settings = dict()
            orig_settings['dataset'], orig_settings['analysis'], orig_settings['plotting'], orig_settings['annotation'] = dict(), dict(), dict(), dict()
            orig_settings['dataset']['input_datafile'] = os.path.join('.', 'input', st.session_state['settings__input_datafile__filename'])
            orig_settings['dataset']['format'] = dict(zip(['HALO', 'Native', 'GMBSecondGeneration', 'QuPath', 'Steinbock'], ['OMAL', 'Native', 'GMBSecondGeneration', 'QuPath', 'Steinbock']))[st.session_state['settings__input_datafile__format']]
            orig_settings['dataset']['coord_units_in_microns'] = st.session_state['settings__input_datafile__coordinate_units']
            orig_settings['dataset']['sep'] = (',' if orig_settings['dataset']['input_datafile'].endswith('.csv') else '\t')
            orig_settings['analysis']['allow_compound_species'] = (False if st.session_state['settings__phenotyping__method'] == 'Marker' else True)
            orig_settings['dataset']['phenotype_identification_tsv_file'] = (os.path.join(input_directory, 'phenotypes', st.session_state['settings__phenotyping__phenotype_identification_file']) if st.session_state['settings__phenotyping__method'] == 'Custom' else None)
            orig_settings['dataset']['roi_width'] = (st.session_state['settings__analysis__roi_width'] if st.session_state['settings__analysis__partition_slides_into_rois'] else None)
            orig_settings['dataset']['overlap'] = (st.session_state['settings__analysis__roi_overlap'] if st.session_state['settings__analysis__partition_slides_into_rois'] else 0)
            orig_settings['analysis']['thickness'] = st.session_state['settings__analysis__neighbor_radius']
            orig_settings['analysis']['use_analytical_significance'] = (True if st.session_state['settings__analysis__significance_calculation_method'] == 'Poisson (radius)' else False)
            orig_settings['analysis']['n_neighs'] = (st.session_state['settings__analysis__n_neighs'] if st.session_state['settings__analysis__significance_calculation_method'] == 'Permutation (k-nearest neighbors)' else -1)
            orig_settings['analysis']['radius_instead_of_knn'] = (True if st.session_state['settings__analysis__significance_calculation_method'] == 'Permutation (radius)' else False)
            orig_settings['plotting']['num_valid_centers_minimum'] = st.session_state['settings__analysis__min_num_valid_centers']
            orig_settings['plotting']['weight_rois_by_num_valid_centers'] = st.session_state['settings__analysis__weight_by_num_valid_centers']
            orig_settings['plotting']['log_pval_range'] = (st.session_state['settings__analysis__log_pval_minimum'], st.session_state['settings__analysis__log_pval_maximum'])  # used for analysis and plotting both, except for the annotation-specific plots
            orig_settings['annotation']['csv_files'] = st.session_state['settings__annotation__used_annotation_files']
            orig_settings['annotation']['annotation_coord_units_in_microns'] = st.session_state['settings__annotation__coordinate_units']
            orig_settings['annotation']['annotation_microns_per_integer_unit'] = st.session_state['settings__annotation__microns_per_integer_unit']
            # orig_settings['annotation']['marker_column_names_list'] = st.session_state['settings__annotation__markers_designating_valid_objects']
            orig_settings['analysis']['images_to_analyze'] = st.session_state['settings__analysis__images_to_analyze']
            orig_settings['plotting']['min_log_pval'] = st.session_state['settings__plotting__min_log_pval']

            # Print the settings to screen in both the new and old formats
            print('Selected settings for the workflow (new format):')
            pprint.pprint(streamlit_utils.get_current_settings(), sort_dicts=False)
            print('Selected settings for the workflow (old format):')
            pprint.pprint(orig_settings, sort_dicts=False)

            # Get dataset loading columns for the button and an output message
            dataset_loading_col1, dataset_loading_col2 = st.columns(2)

            # Create a dataset loading button
            ready_to_preprocess_data = True
            for x in ['input_datafile', 'format', 'coord_units_in_microns', 'sep', 'phenotype_identification_tsv_file', 'roi_width', 'overlap']:
                if x not in orig_settings['dataset']:
                    ready_to_preprocess_data = False
                    break
            with dataset_loading_col1:
                dataset_just_loaded = st.button('Load dataset', disabled=(not ready_to_preprocess_data))

            # Load the dataset into a dataset object defined in dataset_formats.py if the button was pressed
            if dataset_just_loaded:
                if 'sit__using_gated_phenotypes' in st.session_state:
                    use_gated_phenotypes = st.session_state['sit__using_gated_phenotypes']
                else:
                    if orig_settings['dataset']['input_datafile'].split('/')[-1].startswith('orig_datafile_plus_gated_phenotypes'):
                        use_gated_phenotypes = True
                    else:
                        use_gated_phenotypes = False
                dataset_obj = tci.preprocess_dataset(
                    format=orig_settings['dataset']['format'],
                    input_datafile=orig_settings['dataset']['input_datafile'],
                    coord_units_in_microns=orig_settings['dataset']['coord_units_in_microns'],
                    phenotype_identification_tsv_file=orig_settings['dataset']['phenotype_identification_tsv_file'],
                    sep=orig_settings['dataset']['sep'],
                    roi_width=orig_settings['dataset']['roi_width'],
                    overlap=orig_settings['dataset']['overlap'],
                    images_to_analyze=orig_settings['analysis']['images_to_analyze'],
                    use_gated_phenotypes=use_gated_phenotypes
                    )
                st.session_state['dataset_obj'] = dataset_obj
                
            # If the button was pressed and therefore the data was loaded, then say so
            with dataset_loading_col2:
                if dataset_just_loaded:
                    st.info('Dataset loaded', icon="ℹ️")

            # ---- Main workflow --------------------------------------------------------------------------------

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

            # # Get job execution columns
            # job_exec_col1, job_exec_col2 = st.columns(2)

            # Determine whether we should employ threading
            # with job_exec_col1:
            use_multiprocessing = st.checkbox('Should we use multiple logical CPUs to speed up the calculations?', key='use_multiprocessing')

            # Get the number of threads to use for the calculations
            # with job_exec_col2:
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

                # Instantiate the TIME class, which mainly loads the main datafile into a Pandas dataframe (or reads in a simulated one in the case of simulated data) and performs some preprocessing on it. Note that this creates slices.df_roi_plotting_data. This is *n
                # iot* a parallel function
                # Took 40 min on 3625 ROIs on 7/5/22
                iblock = 0
                if workflow_bools[iblock]:
                    print('**** {}... ****'.format(block_names[iblock]))
                    st.write('**:sparkles: {}...**'.format(block_names[iblock]))
                    slices = tci.TIMECellInteraction(
                        dataset_obj,
                        project_dir=project_dir,
                        allow_compound_species=orig_settings['analysis']['allow_compound_species'],
                        thickness_new=orig_settings['analysis']['thickness'],
                        use_analytical_significance=orig_settings['analysis']['use_analytical_significance'],
                        n_neighs=orig_settings['analysis']['n_neighs'],
                        radius_instead_of_knn=orig_settings['analysis']['radius_instead_of_knn']
                    )
                    # st.session_state['main_analysis_input_data'] = slices.data
                    st.session_state['df_data_by_roi'] = slices.df_data_by_roi
                    st.session_state['thickness'] = orig_settings['analysis']['thickness']

                # Plot every ROI using slices.df_roi_plotting_data in parallel
                # Took <2 min on 3625 ROIs on 7/5/22
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
                # Took ~50 min on 3625 ROIs on 7/5/22
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
                        log_pval_range=orig_settings['plotting']['log_pval_range'],
                        num_valid_centers_minimum=orig_settings['plotting']['num_valid_centers_minimum']
                    )
                    st.session_state['df_density_pvals_arrays'] = slices.df_density_pvals_arrays

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
                        weight_rois_by_num_valid_centers=orig_settings['plotting']['weight_rois_by_num_valid_centers'],
                        input_datafile=orig_settings['dataset']['input_datafile']  # this is just needed to get the input data filename to save to disk along with the df_log_dens_pvals_arr_per_slide for later read-in by the correlation analyzer
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
                        annotations_csv_files=orig_settings['annotation']['csv_files'],
                        phenotyping_method=st.session_state['settings__phenotyping__method'],
                        phenotype_identification_file=orig_settings['dataset']['phenotype_identification_tsv_file'],
                        # marker_column_names_list=orig_settings['annotation']['marker_column_names_list'],
                        annotation_coord_units_in_microns=orig_settings['annotation']['annotation_coord_units_in_microns'],
                        annotation_microns_per_integer_unit=orig_settings['annotation']['annotation_microns_per_integer_unit'],
                        settings__analysis__thickness=orig_settings['analysis']['thickness'],
                        min_log_pval_for_plotting=orig_settings['plotting']['min_log_pval']
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

        # # Optionally save df_data_by_roi to disk
        # if st.button('Save df_data_by_roi to disk'):
        #     # Can be read back in using:
        #     # import pickle
        #     # with open('/home/weismanal/projects/spatial-interaction-tool/app-dev/repo/df_data_by_roi.pkl', 'rb') as f:
        #     #     df_data_by_roi = pickle.load(f)
        #     import pickle
        #     with open('/home/weismanal/projects/spatial-interaction-tool/app-dev/repo/two_dfs.pkl', 'wb') as f:
        #         pickle.dump((st.session_state['df_data_by_roi'], st.session_state['df_density_pvals_arrays']), f)

    else:
        st.warning('An input datafile must be selected on the previous page', icon='⚠️')

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)

if __name__ == '__main__':
    main()
