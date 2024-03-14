# Import relevant libraries
import os
import yaml
import pandas as pd
import streamlit as st
import streamlit_utils
import pprint
import platform_io
import time_cell_interaction_lib as tci  # import the TIME library stored in time_cell_interaction_lib.py
import utils
import app_top_of_page as top
import streamlit_dataframe_editor as sde
import dataset_formats

# Input/output directory initializations
input_directory = os.path.join('.', 'input')
output_directory = os.path.join('.', 'output')

def load_phenotyping_settings_from_thresholded_phenotyper():
    """
    Load the settings from the thresholded Phenotyper and map them to the corresponding ones in the session state.
    """

    # If the Phenotyper's phenotyping method is "Custom", create a phenotype assignments file (and assign the corresponding setting) from its phenotype assignment table
    if st.session_state['phenoMeth'] == 'Custom':
        df_pheno_assignments = st.session_state['pheno__de_phenotype_assignments'].reconstruct_edited_dataframe()
        df_pheno_assignments = df_pheno_assignments[df_pheno_assignments['phenotype'] != 'unassigned']
        st.session_state['settings__phenotyping__phenotype_identification_file'] = create_phenotype_assignments_file_from_phenotyper(df_pheno_assignments)

    # Set the phenotyping method from the Phenotyper and update its dependencies
    st.session_state['settings__phenotyping__method'] = st.session_state['phenoMeth']
    update_dependencies_of_phenotyping_method()

# Define the function to copy an input file from the output directory to the input directory
def copy_input_file_from_output_dir_to_input_dir(input_filename, input_subdir=None):
    import shutil
    if input_subdir is not None:
        input_directory2 = os.path.join(input_directory, input_subdir)
        if not os.path.exists(input_directory2):
            os.makedirs(input_directory2)
    else:
        input_directory2 = input_directory
    shutil.copy(os.path.join(output_directory, input_filename), os.path.join(input_directory2, input_filename))

def update_dependencies_of_input_datafile_filename():

    # This should no longer be used as of 3/13/24

    if st.session_state['settings__input_datafile__filename'] is not None:

        # Set full pathname to input datafile
        input_datafile_path = os.path.join(input_directory, st.session_state['settings__input_datafile__filename'])

        # Update input_datafile__format value
        _, _, _, _, file_format, _ = dataset_formats.extract_datafile_metadata(input_datafile_path)
        if file_format is not None:
            st.session_state['settings__input_datafile__format'] = file_format

        # Update the image options
        st.session_state['options_for_images'] = list(dataset_formats.get_image_series_in_datafile(input_datafile_path).unique())

        # Update analysis__images_to_analyze value
        st.session_state['settings__analysis__images_to_analyze'] = st.session_state['options_for_images']
        update_dependencies_of_analysis_images_to_analyze()

def update_dependencies_of_phenotyping_method():
    st.session_state['phenotyping_phenotype_identification_file_is_disabled'] = (True if (st.session_state['settings__phenotyping__method'] in ['Species', 'Marker']) else False)

def update_dependencies_of_analysis_images_to_analyze():

    # annotation__used_annotation_files options
    annotations_dir_listing = ([x for x in os.listdir(os.path.join(input_directory, 'annotations')) if x.endswith('.csv')] if os.path.exists(os.path.join(input_directory, 'annotations')) else [])
    st.session_state['options_for_annotation_files'] = [x for x in annotations_dir_listing if x.split('__')[0] in st.session_state['settings__analysis__images_to_analyze']]

    # annotation__used_annotation_files values
    st.session_state['settings__annotation__used_annotation_files'] = st.session_state['options_for_annotation_files']
    update_dependencies_of_annotation_used_annotation_files()

def update_dependencies_of_analysis_partition_slides_into_rois():
    if st.session_state['settings__analysis__partition_slides_into_rois']:
        st.session_state['analysis_roi_width_is_disabled'] = False
        st.session_state['analysis_roi_overlap_is_disabled'] = False
    else:
        st.session_state['analysis_roi_width_is_disabled'] = True
        st.session_state['analysis_roi_overlap_is_disabled'] = True

def update_dependencies_of_annotation_used_annotation_files():
    if len(st.session_state['settings__annotation__used_annotation_files']) == 0:
        st.session_state['annotation_coordinate_units_is_disabled'] = True
        st.session_state['annotation_coord_units_are_pixels_is_disabled'] = True
        st.session_state['annotation_microns_per_integer_unit_is_disabled1'] = True
        st.session_state['plotting_min_log_pval_is_disabled'] = True
    else:
        st.session_state['annotation_coordinate_units_is_disabled'] = False
        st.session_state['annotation_coord_units_are_pixels_is_disabled'] = False
        st.session_state['annotation_microns_per_integer_unit_is_disabled1'] = False
        st.session_state['plotting_min_log_pval_is_disabled'] = False

def update_dependencies_of_annotation_coordinate_units():
    if st.session_state['settings__annotation__coord_units_are_pixels']:
        st.session_state['settings__annotation__microns_per_integer_unit'] = st.session_state['settings__annotation__coordinate_units']

def update_dependencies_of_annotation_coord_units_are_pixels():
    if st.session_state['settings__annotation__coord_units_are_pixels']:
        st.session_state['settings__annotation__microns_per_integer_unit'] = st.session_state['settings__annotation__coordinate_units']
        st.session_state['annotation_microns_per_integer_unit_is_disabled2'] = True
    else:
        st.session_state['annotation_microns_per_integer_unit_is_disabled2'] = False

def update_dependencies_of_analysis_significance_calculation_method():
    if st.session_state['settings__analysis__significance_calculation_method'] == 'Permutation (k-nearest neighbors)':
        st.session_state['analysis_neighbor_radius_is_disabled'] = True
        st.session_state['analysis_n_neighs_is_disabled'] = False
    else:
        st.session_state['analysis_neighbor_radius_is_disabled'] = False
        st.session_state['analysis_n_neighs_is_disabled'] = True

def set_session_state_key(settings, str1, str2):
    if str2 in settings[str1]:
        st.session_state['settings__{}__{}'.format(str1, str2)] = settings[str1][str2]

def create_phenotype_assignments_file_from_phenotyper(df_phenotype_assignments):

    # Import relevant libraries
    from datetime import datetime
    import numpy as np

    # Set the full directory path to the phenotypes files
    phenotypes_path = os.path.join(input_directory, 'phenotypes')

    # Create this path if it doesn't already exist
    if not os.path.exists(phenotypes_path):
        os.makedirs(phenotypes_path)

    # Set the filename of the phenotype assignments file to write
    filename = 'phenotype_assignments_from_phenotyper-{}.tsv'.format(datetime.now().strftime("date%Y_%m_%d_time%H_%M_%S"))

    # Assign a new dataframe as a subset of the one containing the phenotype assignments
    df_phenotype_assignments_to_write = df_phenotype_assignments[['species_count', 'species_percent', 'species_name_short', 'phenotype', 'species_name_long']]

    # Set a dummy index; hopefully this should have no effect on anything
    df_phenotype_assignments_to_write.index = [-1] * len(df_phenotype_assignments_to_write)

    # Modify the columns to spec
    df_phenotype_assignments_to_write = df_phenotype_assignments_to_write[df_phenotype_assignments_to_write['species_name_long'].apply(lambda species_name_long: sum([0 if x[-1] == '-' else 1 for x in species_name_long.split(' ')]) != 0)]
    df_phenotype_assignments_to_write = df_phenotype_assignments_to_write.drop('species_name_long', axis='columns')
    df_phenotype_assignments_to_write['species_name_short'] = df_phenotype_assignments_to_write['species_name_short'].apply(lambda x: sorted(x.rstrip('+').split('+ ')))
    df_phenotype_assignments_to_write['species_percent'] = (df_phenotype_assignments_to_write['species_count'] / df_phenotype_assignments_to_write['species_count'].sum() * 100).apply(lambda x: np.round(x, decimals=8))

    # If any values in the "phenotype" column contain hyphens, output a warning that hyphens indicate compound phenotypes, which isn't (yet!) implemented in the Phenotyper and will be removed, though compound phenotypes can still be specified by modifying the TSV file manually, and note that the hyphens are being replaced with "(dash)". Perform this replacement.
    if df_phenotype_assignments_to_write['phenotype'].apply(lambda x: '-' in x).sum() > 0:
        st.warning('Hyphens in the assigned custom phenotype names in the Phenotyper app indicate compound phenotypes, which aren\'t (yet!) implemented in the Phenotyper and will be replaced with "(dash)". Compound phenotypes can still be specified by modifying the phenotype assignments TSV file manually. So if you don\'t want your hyphens replaced with "(dash)", please modify your custom phenotype names in the Phenotyper app and click the "Load relevant settings from Phenotyper" button again.')
        df_phenotype_assignments_to_write['phenotype'] = df_phenotype_assignments_to_write['phenotype'].apply(lambda x: x.replace('-', '(dash)'))

    # Write the dataframe to disk, also making a copy in the output directory so it can optionally be saved for good
    df_phenotype_assignments_to_write.to_csv(path_or_buf=os.path.join(phenotypes_path, filename), sep='\t', header=False)
    df_phenotype_assignments_to_write.to_csv(path_or_buf=os.path.join(output_directory, filename), sep='\t', header=False)

    # Return the filename of the written file
    return filename

def write_dataframe_to_disk(df, prefix='phenotyped_datafile_from_gater'):

    # Import relevant library
    from datetime import datetime
    import shutil

    # Set the filename of the phenotype assignments file to write
    filename = '{}-{}.csv'.format(prefix, datetime.now().strftime("date%Y_%m_%d_time%H_%M_%S"))

    # Save the dataframe to disk in both the output and input directories (the former for posterity, the latter so that it can be read in later)
    filepath_to_write = os.path.join(output_directory, filename)
    df.to_csv(path_or_buf=filepath_to_write, index=False)
    shutil.copy(filepath_to_write, input_directory)

    # Return the filename of the written file
    return filename

def load_relevant_settings_from_phenotyper():

    # As of 3/13/24 this should no longer be used, thankfully.

    # Get the datafile from the phenotyper, which may have old phenotype columns as "Phenotype_orig " and gated phenotypes as "Phenotype " if the Gater were used first
    new_df = st.session_state.df
    new_df_columns = new_df.columns

    # If the main data dataframe contains both "Phenotype_orig " and "Phenotype " columns, then the Phenotyper must have loaded the dataframe from the Gater, which means there is no datafile, so we must create one, and set its filename as the input datafile...
    phenotype_orig_columns_exist = len([column for column in new_df_columns if column.startswith('Phenotype_orig ')]) > 0
    phenotype_columns_exist = len([column for column in new_df_columns if column.startswith('Phenotype ')]) > 0
    if phenotype_orig_columns_exist and phenotype_columns_exist:

        # Grab the original datafile on which the multiaxial gating was based from the multiaxial gater
        orig_filename = st.session_state['mg__input_datafile_filename']

        # Get the type of text file that is the original datafile
        sep = (',' if orig_filename.endswith('.csv') else '\t')

        # If the original datafile is not present in the input directory, then don't do anything
        if not os.path.exists(os.path.join(input_directory, orig_filename)):
            st.error(f'The input file {orig_filename} does not appear to be present in the input directory. Please ensure it is there and try again.')
            return

        # Read in this original datafile from disk
        orig_df = pd.read_csv(os.path.join(input_directory, orig_filename), sep=sep)

        # Obtain the columns from the new datafile to paste on to the end of the original one
        new_df_to_add = new_df[[column for column in new_df_columns if column.startswith('Phenotype ')]]

        # Change the new phenotype columns to something uniquely identifiable
        columns_to_add = new_df_to_add.columns
        transform = dict(zip(columns_to_add, [column.replace('Phenotype ', 'Phenotype-from-multiaxial-gater ', 1) for column in columns_to_add]))
        new_df_to_add = new_df_to_add.rename(columns=transform)

        # Append these columns to the original dataframe and write the result to disk, storing the filename in the settings for the SIT
        st.session_state['settings__input_datafile__filename'] = write_dataframe_to_disk(pd.concat([orig_df, new_df_to_add], axis='columns'), prefix='orig_datafile_plus_gated_phenotypes')

        # Save the input datafile coordinate units in microns from the Gater (since Dante hasn't ported it yet [phenotyping_micron_coordinate_units key], load from the Gater for now instead of the Phenotyper)
        st.session_state['settings__input_datafile__coordinate_units'] = st.session_state['mg__input_datafile_coordinate_units']

    # Otherwise, the datafile was likely read in from disk (as opposed to from memory via Streamlit), so set that filename as the input datafile
    else:

        # If the Phenotyper has not been run, then don't do anything
        if 'datafileU' not in st.session_state:
            st.error('The Phenotyper does not appear to have been run. Please run it and try again.')
            return

        # Set the filename for the SIT as that in the Phenotyper's widget
        st.session_state['settings__input_datafile__filename'] = st.session_state['datafileU']

        # Save the input datafile coordinate units in microns from the Phenotyper
        st.session_state['settings__input_datafile__coordinate_units'] = st.session_state['phenotyping_micron_coordinate_units']

    # Update the dependencies of the input datafile filename since it has likely changed
    update_dependencies_of_input_datafile_filename()

    # If the Phenotyper's phenotyping method is "Custom", create a phenotype assignments file (and assign the corresponding setting) from its phenotype assignment table
    if st.session_state['phenoMeth'] == 'Custom':
        df_pheno_assignments = st.session_state['pheno__de_phenotype_assignments'].reconstruct_edited_dataframe()
        df_pheno_assignments = df_pheno_assignments[df_pheno_assignments['phenotype'] != 'unassigned']
        st.session_state['settings__phenotyping__phenotype_identification_file'] = create_phenotype_assignments_file_from_phenotyper(df_pheno_assignments)

    # Set the phenotyping method from the Phenotyper and update its dependencies
    st.session_state['settings__phenotyping__method'] = st.session_state['phenoMeth']
    update_dependencies_of_phenotyping_method()

def load_dataset_and_settings(checkpoints_exist, existing_dirs_to_delete, orig_settings):
    '''
    Callback for when the "Load dataset and settings" button is pressed
    '''

    # Delete any existing checkpoints so that both the preprocessing and the rest of the workflow will run from scratch
    if checkpoints_exist:
        platform_io.delete_selected_files_and_dirs('output', existing_dirs_to_delete)

    # Save to memory all settings that are actually being used to run the SIT
    st.session_state['sit__used_settings'] = orig_settings.copy()

    # Make a modified copy of the dataset object that suits the SIT
    dataset_obj = st.session_state['input_dataset'].copy()  # making a copy since we might modify the dataset_obj.data attribute by adding rows for patching potentially, filtering to images, and trimming the columns

    # Filter to just the necessary columns
    dataset_obj.data = dataset_formats.trim_dataframe_basic(dataset_obj.data)  # this returns a copy of the input dataframe

    # Filter to just the necessary rows
    dataset_obj.data = dataset_obj.data.loc[dataset_obj.data['Slide ID'].apply(lambda curr_image: curr_image in st.session_state['sit__used_settings']['analysis']['images_to_analyze'])]  # this returns a copy that will overwrite dataset_obj.data

    # Set a potentially important attribute
    dataset_obj.phenotype_identification_tsv_file = st.session_state['sit__used_settings']['dataset']['phenotype_identification_tsv_file']

    # Potentially patch the dataset into ROIs (generally increasing the number of cells)
    if st.session_state['sit__used_settings']['dataset']['roi_width'] is not None:
        dataset_obj.data = dataset_formats.potentially_apply_patching(dataset_obj.data, ['Cell X Position', 'Cell Y Position'], st.session_state['sit__used_settings']['dataset']['roi_width'], st.session_state['sit__used_settings']['dataset']['overlap'], lambda df_tmp: (df_tmp / 0.2).astype(int), lambda x: int(x / 0.2))  # see dataset_formats.py for 0.2 justification

    # Do any extra processing, namely, deleting ROIs with just a single coordinate
    dataset_obj.extra_processing()

    # Save the new dataset to the session state
    st.session_state['dataset_obj'] = dataset_obj

    # Reset any image path extraction in subsequent tabs
    if 'df_paths_per_roi' in st.session_state:
        del st.session_state['df_paths_per_roi']
    if 'df_paths_per_slide' in st.session_state:
        del st.session_state['df_paths_per_slide']
    if 'overlay_info' in st.session_state:
        del st.session_state['overlay_info']

    # If the button was pressed and therefore the data was loaded, then say so
    st.toast('Dataset loaded')

# Define the main function
def main():

    # Data needed for widget options
    options_for_parameter_files =                [os.path.join(input_directory, x) for x in os.listdir(input_directory) if x.endswith('.yml')] + \
                                                 [os.path.join(output_directory, x) for x in os.listdir(output_directory) if (x.endswith('.yml') and ('environment_as_of_' not in x))]
    options_for_phenotype_identification_files = ([x for x in os.listdir(os.path.join(input_directory, 'phenotypes')) if x.endswith('.tsv')] if os.path.exists(os.path.join(input_directory, 'phenotypes')) else [])
    options_for_phenotyping_methods = ['Species', 'Marker', 'Custom']
    options_for_significance_calculation_methods = ['Poisson (radius)', 'Permutation (radius)', 'Permutation (k-nearest neighbors)']

    # Set page settings
    st.set_page_config(layout="wide", page_title='Tool parameter selection')
    st.title('Tool parameter selection')

    # Run streamlit-dataframe-editor library initialization tasks at the top of the page
    st.session_state = sde.initialize_session_state(st.session_state)

    # Run Top of Page (TOP) functions
    st.session_state = top.top_of_page_reqs(st.session_state)

    # If 'input_dataset' isn't in the session state, print an error message and return
    if 'input_dataset' not in st.session_state:
        st.error('An input dataset has not yet been opened. Please do so using the "Open file(s)" page in the sidebar.')
        return

    # Assign default keys
    streamlit_utils.assign_default_values_in_session_state('preset_parameter_file', utils.get_first_element_or_none(options_for_parameter_files))

    # Load or reload default settings
    if st.button('Reload default settings'):
        st.session_state['settings__phenotyping__method'] = utils.get_first_element_or_none(options_for_phenotyping_methods)
        update_dependencies_of_phenotyping_method()
        st.session_state['settings__phenotyping__phenotype_identification_file'] = utils.get_first_element_or_none(options_for_phenotype_identification_files)
        if 'settings__analysis__images_to_analyze' not in st.session_state:
            st.session_state['options_for_images'] = []
            st.session_state['settings__analysis__images_to_analyze'] = st.session_state['options_for_images']
        update_dependencies_of_analysis_images_to_analyze()
        st.session_state['settings__analysis__partition_slides_into_rois'] = True
        update_dependencies_of_analysis_partition_slides_into_rois()
        st.session_state['settings__analysis__roi_width'] = 400
        st.session_state['settings__analysis__roi_overlap'] = 80
        st.session_state['settings__analysis__significance_calculation_method'] = utils.get_first_element_or_none(options_for_significance_calculation_methods)
        update_dependencies_of_analysis_significance_calculation_method()
        st.session_state['settings__analysis__neighbor_radius'] = 40
        st.session_state['settings__analysis__n_neighs'] = 6
        st.session_state['settings__analysis__min_num_valid_centers'] = 10
        st.session_state['settings__analysis__weight_by_num_valid_centers'] = False
        st.session_state['settings__analysis__log_pval_minimum'] = -50
        st.session_state['settings__analysis__log_pval_maximum'] = 0
        if 'settings__annotation__used_annotation_files' not in st.session_state:
            st.session_state['options_for_annotation_files'] = []
            st.session_state['settings__annotation__used_annotation_files'] = st.session_state['options_for_annotation_files']
        update_dependencies_of_annotation_used_annotation_files()
        st.session_state['settings__annotation__coordinate_units'] = 0.25
        # update_dependencies_of_annotation_coordinate_units() --> unnecessary to run *here*
        st.session_state['settings__annotation__coord_units_are_pixels'] = True
        update_dependencies_of_annotation_coord_units_are_pixels()
        if 'settings__annotation__microns_per_integer_unit' not in st.session_state:
            st.session_state['settings__annotation__microns_per_integer_unit'] = 0.25
        st.session_state['settings__plotting__min_log_pval'] = -50

    # Allow the user to select an existing parameter file
    columns = st.columns(2)
    with columns[0]:
        if st.session_state['preset_parameter_file'] not in options_for_parameter_files:
            if len(options_for_parameter_files) > 0:
                st.session_state['preset_parameter_file'] = options_for_parameter_files[0]
            else:
                st.session_state['preset_parameter_file'] = None
        st.selectbox('Optionally, instead load settings from available parameter file:', options_for_parameter_files, key='preset_parameter_file')

    # Optionally load the selected parameter file and map the loaded settings to the corresponding ones in the session state
    if st.button(':arrow_down: Load selected parameter file', disabled=(st.session_state['preset_parameter_file'] is None)):
        with open(st.session_state['preset_parameter_file'], mode='rt') as file:
            preset_settings = yaml.load(file, yaml.UnsafeLoader)
        settings = utils.validate_presets_and_map_to_settings(preset_settings, options_for_phenotype_identification_files, st.session_state['options_for_annotation_files'], options_for_phenotyping_methods, options_for_significance_calculation_methods, st.session_state['options_for_images'], message_function=streamlit_utils.streamlit_write_function)
        set_session_state_key(settings, 'phenotyping', 'method')
        update_dependencies_of_phenotyping_method()
        set_session_state_key(settings, 'phenotyping', 'phenotype_identification_file')
        set_session_state_key(settings, 'analysis', 'images_to_analyze')
        update_dependencies_of_analysis_images_to_analyze()
        set_session_state_key(settings, 'analysis', 'partition_slides_into_rois')
        update_dependencies_of_analysis_partition_slides_into_rois()
        set_session_state_key(settings, 'analysis', 'roi_width')
        set_session_state_key(settings, 'analysis', 'roi_overlap')
        set_session_state_key(settings, 'analysis', 'significance_calculation_method')
        update_dependencies_of_analysis_significance_calculation_method()
        set_session_state_key(settings, 'analysis', 'neighbor_radius')
        set_session_state_key(settings, 'analysis', 'n_neighs')
        set_session_state_key(settings, 'analysis', 'min_num_valid_centers')
        set_session_state_key(settings, 'analysis', 'weight_by_num_valid_centers')
        set_session_state_key(settings, 'analysis', 'log_pval_minimum')
        set_session_state_key(settings, 'analysis', 'log_pval_maximum')
        set_session_state_key(settings, 'annotation', 'used_annotation_files')
        update_dependencies_of_annotation_used_annotation_files()
        set_session_state_key(settings, 'annotation', 'coordinate_units')
        update_dependencies_of_annotation_coordinate_units()
        set_session_state_key(settings, 'annotation', 'coord_units_are_pixels')
        update_dependencies_of_annotation_coord_units_are_pixels()
        set_session_state_key(settings, 'annotation', 'microns_per_integer_unit')
        set_session_state_key(settings, 'plotting', 'min_log_pval')

    # Logical divider
    st.divider()

    # Create columns
    parameter_columns = st.columns(3)

    # Phenotyping settings widgets
    with parameter_columns[0]:
        st.header('Phenotyping settings')
        st.button('Load phenotyping settings from the thresholded phenotyper', on_click=load_phenotyping_settings_from_thresholded_phenotyper, disabled=('phenoMeth' not in st.session_state))
        st.selectbox('Method:', options_for_phenotyping_methods, key='settings__phenotyping__method', help='Species: phenotypes defined by the unique combinations of markers present in the input file. Marker: each marker is its own phenotype. Custom: custom phenotyping using a text file.', on_change=update_dependencies_of_phenotyping_method)
        if st.session_state['settings__phenotyping__phenotype_identification_file'] is not None:
            if (not os.path.exists(os.path.join(input_directory, 'phenotypes', st.session_state['settings__phenotyping__phenotype_identification_file']))) and (os.path.exists(os.path.join(output_directory, st.session_state['settings__phenotyping__phenotype_identification_file']))):
                copy_input_file_from_output_dir_to_input_dir(st.session_state['settings__phenotyping__phenotype_identification_file'], input_subdir='phenotypes')
                st.rerun()
        st.selectbox('Phenotype identification file:', options_for_phenotype_identification_files, key='settings__phenotyping__phenotype_identification_file', help='See [here](https://github.com/ncats/spatial-interaction-tool/blob/4e1240fba45cb3bc2290be18903af03f9d3fdf6a/config/phenotype_identifications/gmb_phenotype_ids_as_of_2022-07-05.tsv) for a sample phenotype identification file, which must be present in the input/phenotypes directory and have a .tsv extension.', disabled=st.session_state['phenotyping_phenotype_identification_file_is_disabled'])

    # Analysis settings widgets
    with parameter_columns[1]:
        st.header('Analysis settings')
        st.multiselect('Images to analyze:', st.session_state['options_for_images'], placeholder='', key='settings__analysis__images_to_analyze', on_change=update_dependencies_of_analysis_images_to_analyze)
        st.checkbox('Partition slides into regions of interest (ROIs)', key='settings__analysis__partition_slides_into_rois', on_change=update_dependencies_of_analysis_partition_slides_into_rois)
        st.number_input('ROI width (microns):', min_value=0.0, key='settings__analysis__roi_width', disabled=st.session_state['analysis_roi_width_is_disabled'])
        st.number_input('ROI overlap (microns):', min_value=0.0, key='settings__analysis__roi_overlap', disabled=st.session_state['analysis_roi_overlap_is_disabled'], help='In order to analyze the neighbors around every cell exactly once (as you should), this should be set to twice the neighbor radius (adjusted below).')
        st.selectbox('Significance calculation method:', options_for_significance_calculation_methods, key='settings__analysis__significance_calculation_method', on_change=update_dependencies_of_analysis_significance_calculation_method)
        st.number_input('Neighbor radius (microns):', min_value=0.0, key='settings__analysis__neighbor_radius', help='Radius around each center cell within which to count neighboring cells.', disabled=st.session_state['analysis_neighbor_radius_is_disabled'])
        st.number_input('Number of nearest neighbors:', min_value=0, key='settings__analysis__n_neighs', disabled=st.session_state['analysis_n_neighs_is_disabled'])
        st.number_input('Minimum number of valid centers:', min_value=1, key='settings__analysis__min_num_valid_centers', help='Valid centers are defined as cells inside an analysis radius from the ROI edges. If a particular phenotype in a ROI has fewer than the minimum number of valid cells, it is excluded from downstream analysis as a center species (not as a neighbor species). I.e., the user may want to define a minimum sample size.')
        st.checkbox('For just the "Display average heatmaps" page, weight by the number of valid centers', key='settings__analysis__weight_by_num_valid_centers', help='This refers to the averaging performed on the "Display average heatmaps" page, in which the rows in the density P value heatmaps for each ROI are weighted by the number of valid centers in each row. (This *is not* center-neighbor symmetric because in a given analysis, the number valid neighbors is the total number of neighbors in the ROI whereas the number of valid centers is the number of centers within an analysis radius buffer from the ROI edge.) This has no effect on the averaging done per annotation region type, i.e., the averaging performed on the "Display average heatmaps per annotation" page, in which the ROIs are uniformly weighted depending on the weighting method selected on that page, for each annotation region type. (This *is* center-neighbor symmetric because the same weight is assigned to all the center-neighbor pairs for each ROI.) So the averaging done on these two pages is fundamentally different, even aside from annotations being involved.')

        # Input of clipping information
        st.write('Clipping range of the log (base 10) density P values:')
        subcolumns = st.columns(2)
        with subcolumns[0]:
            st.number_input('Minimum:', max_value=0.0, key='settings__analysis__log_pval_minimum')
        with subcolumns[1]:
            st.number_input('Maximum:', max_value=0.0, key='settings__analysis__log_pval_maximum', help='For density P value heatmaps for individual ROIs this does not matter and is just a plotting property. But whenever averaging is performed over the ROIs (i.e., for all ROIs in a slide or all ROIs in an annotation region type), these values matter because the clipping is deliberately performed prior to the averaging.')

    # Annotations settings widgets
    with parameter_columns[2]:
        st.header('Annotation settings', help='We recommend averaging over annotation type because slides are generally heterogeneous. Each annotation filename (the files should reside in input/annotations) should contain the image number at the very beginning, followed by a *double* underscore "__". The annotation region name should be at the very end (aside from the arbitrary extension), preceded by a *double* underscore. E.g.: "12345__panel_3_20230913__tumor.csv". The image numbers in the annotation filenames should match those in the "Image Location" column of the original input datafile for analysis, where, as a reminder, the image number is at the very end (aside from the arbitrary extension), preceded by a *single* underscore, e.g., "\\\\\\10.1.1.1\\path\\to\\multiplex_panel_3_12345.tif".')
        st.multiselect('Annotation files to use:', st.session_state['options_for_annotation_files'], key='settings__annotation__used_annotation_files', placeholder='', on_change=update_dependencies_of_annotation_used_annotation_files)
        st.number_input('x-y coordinate units (microns):', min_value=0.0, key='settings__annotation__coordinate_units', help='E.g., if the coordinates in the annotation datafiles were pixels, this number would be a conversion to microns in units of microns/pixel. Note: This does not have to be the same value as for the input datafile, though it likely is.', format='%.4f', disabled=st.session_state['annotation_coordinate_units_is_disabled'], on_change=update_dependencies_of_annotation_coordinate_units, step=0.0001)
        st.checkbox('Are pixels the units of the annotation file coordinates?', key='settings__annotation__coord_units_are_pixels', disabled=st.session_state['annotation_coord_units_are_pixels_is_disabled'], on_change=update_dependencies_of_annotation_coord_units_are_pixels)
        st.number_input('Number of microns per integer unit (microns):', min_value=0.0, key='settings__annotation__microns_per_integer_unit', format='%.4f', disabled=(st.session_state['annotation_microns_per_integer_unit_is_disabled1'] or st.session_state['annotation_microns_per_integer_unit_is_disabled2']), help='Conversion of the coordinates in microns to an integer (not necessarily pixels). If the input coordinates are pixels (which are obviously integers), then this value ("microns_per_integer_unit") could be the same as the number of microns per coordinate unit ("coordinate_units") set above--which is basically trivial--so the box above should be checked. This is the most common case. But if, say, the input coordinate units were microns, then coordinate_units=1 and microns_per_integer_unit is most likely a different number, e.g., the number of microns per pixel. The point is to be able to convert the annotation coordinates to integer units, which don\'t necessarily need to be (but usually are) pixels.', step=0.0001)
        st.number_input('Minimum density log P value for plotting:', max_value=0.0, key='settings__plotting__min_log_pval', disabled=st.session_state['plotting_min_log_pval_is_disabled'])

    # Separate out the data loading button from the rest of the widgets above
    st.divider()

    # Convert widget settings above to those compatible with the original SIT API
    orig_settings = dict()
    orig_settings['dataset'], orig_settings['analysis'], orig_settings['plotting'], orig_settings['annotation'], orig_settings['phenotyping'] = dict(), dict(), dict(), dict(), dict()
    orig_settings['phenotyping']['method'] = st.session_state['settings__phenotyping__method']
    orig_settings['analysis']['allow_compound_species'] = (False if orig_settings['phenotyping']['method'] == 'Marker' else True)
    orig_settings['dataset']['phenotype_identification_tsv_file'] = (os.path.join(input_directory, 'phenotypes', st.session_state['settings__phenotyping__phenotype_identification_file']) if orig_settings['phenotyping']['method'] == 'Custom' else None)
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
    orig_settings['analysis']['images_to_analyze'] = st.session_state['settings__analysis__images_to_analyze']
    orig_settings['plotting']['min_log_pval'] = st.session_state['settings__plotting__min_log_pval']

    # Save the currently selected settings, in the original API format, to memory
    st.session_state['sit__workflow_settings'] = orig_settings

    # Print the settings to screen in both the new and old formats
    print('Selected settings for the workflow (new format):')
    pprint.pprint(streamlit_utils.get_current_settings(), sort_dicts=False)
    print('Selected settings for the workflow (old format):')
    pprint.pprint(orig_settings, sort_dicts=False)

    # Get dataset loading columns for the button and an output message
    dataset_loading_cols = st.columns(3)

    # Assess whether tci.preprocess_dataset() is ready to be called
    ready_to_preprocess_data = True
    for x in ['phenotype_identification_tsv_file', 'roi_width', 'overlap']:
        if x not in orig_settings['dataset']:
            ready_to_preprocess_data = False
            break
    for x in ['images_to_analyze']:
        if x not in orig_settings['analysis']:
            ready_to_preprocess_data = False
            break

    # Determine if any checkpoints (which are directories of pickle files or images) exist
    output_dir_listing = os.listdir(output_directory)
    dirs_to_delete = ['checkpoints', 'images', 'logs']
    existing_dirs_to_delete = set(dirs_to_delete).intersection(set(output_dir_listing))
    if len(existing_dirs_to_delete) > 0:
        checkpoints_exist = True
        help_message = 'WARNING: Clicking this button will delete these directories in the local `output` directory: {}. If you wish, first back them up in the "Data Import and Export" tab at left.'.format(existing_dirs_to_delete)
        potential_icon_prefix = '⚠️'
    else:
        checkpoints_exist = False
        help_message = None
        potential_icon_prefix = ''

    # Create a dataset (and settings) loading button
    with dataset_loading_cols[1]:
        if 'sit__used_settings' in st.session_state:
            if st.session_state['sit__used_settings'] != orig_settings:
                st.warning('The current settings differ from those used when the tool was last run. Click the button below to reload the dataset and settings, and then click the "Run workflow" button on the next page to rerun the workflow using the updated data/settings.', icon="⚠️")
        st.button('{} Load dataset and settings'.format(potential_icon_prefix), disabled=(not ready_to_preprocess_data), help=help_message, on_click=load_dataset_and_settings, kwargs={'checkpoints_exist': checkpoints_exist, 'existing_dirs_to_delete': existing_dirs_to_delete, 'orig_settings': orig_settings}, use_container_width=True)

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)

# Call the main function
if __name__ == '__main__':
    main()
