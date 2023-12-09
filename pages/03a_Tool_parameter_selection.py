import os
import yaml
import streamlit as st
from streamlit_extras.app_logo import add_logo
from st_pages import show_pages_from_config, add_indentation
import streamlit_utils

# Import relevant libraries
import utils
import app_top_of_page as top
import streamlit_dataframe_editor as sde

def main():

    # Constant
    input_directory = os.path.join('.', 'input')
    output_directory = os.path.join('.', 'output')

    # Data needed for widget options
    options_for_parameter_files =                [os.path.join(input_directory, x) for x in os.listdir(input_directory) if x.endswith('.yml')] + \
                                                 [os.path.join(output_directory, x) for x in os.listdir(output_directory) if (x.endswith('.yml') and ('environment_as_of_' not in x))]
    options_for_input_datafiles =                [x for x in os.listdir(input_directory) if x.endswith(('.csv', '.tsv'))]
    options_for_phenotype_identification_files = ([x for x in os.listdir(os.path.join(input_directory, 'phenotypes')) if x.endswith('.tsv')] if os.path.exists(os.path.join(input_directory, 'phenotypes')) else [])
    # options_for_annotation_files = st.session_state['options_for_annotation_files']  <-- dynamically updated options
    options_for_input_datafile_formats = ['HALO', 'Native', 'GMBSecondGeneration', 'REEC', 'QuPath']
    options_for_phenotyping_methods = ['Species', 'Marker', 'Custom']
    options_for_significance_calculation_methods = ['Poisson (radius)', 'Permutation (radius)', 'Permutation (k-nearest neighbors)']
    # options_for_images = st.session_state['options_for_images']  <-- dynamically updated options

    def update_dependencies_of_input_datafile_filename():
    
        # Import relevant library
        import dataset_formats

        if st.session_state['settings__input_datafile__filename'] is not None:

            # Set full pathname to input datafile
            input_datafile_path = os.path.join(input_directory, st.session_state['settings__input_datafile__filename'])

            # Update input_datafile__format value
            _, _, _, _, file_format, _ = dataset_formats.extract_datafile_metadata(input_datafile_path)
            if file_format is not None:
                st.session_state['settings__input_datafile__format'] = file_format

            # Update analysis__images_to_analyze options
            # st.session_state['options_for_images'] = utils.get_unique_image_ids_from_datafile(input_datafile_path)
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
        filename = 'phenotype_assignments_from_phenotyper-{}.tsv'.format(datetime.now().strftime("%Y%m%d_%H%M%S"))

        # Assign a new dataframe as a subset of the one containing the phenotype assignments
        df_phenotype_assignments_to_write = df_phenotype_assignments[['species_count', 'species_percent', 'species_name_short', 'phenotype', 'species_name_long']]

        # Set the index to the "Species int" value in time_cell_interaction_lib.py
        full_marker_list = [x.removeprefix('Phenotype ') for x in st.session_state['df'].columns if x.startswith('Phenotype ')]
        num_all_markers = len(full_marker_list)
        df_phenotype_assignments_to_write.index = df_phenotype_assignments_to_write['species_name_short'].apply(lambda species_name_short: (sum([2 ** (num_all_markers - full_marker_list.index(positive_marker) - 1) for positive_marker in species_name_short[:-1].split('+ ')]) if species_name_short != 'Other' else 0))

        # Modify the columns to spec
        df_phenotype_assignments_to_write = df_phenotype_assignments_to_write[df_phenotype_assignments_to_write['species_name_long'].apply(lambda species_name_long: sum([0 if x[-1] == '-' else 1 for x in species_name_long.split(' ')]) != 0)]
        df_phenotype_assignments_to_write = df_phenotype_assignments_to_write.drop('species_name_long', axis='columns')
        df_phenotype_assignments_to_write['species_name_short'] = df_phenotype_assignments_to_write['species_name_short'].apply(lambda x: sorted(x.rstrip('+').split('+ ')))
        df_phenotype_assignments_to_write['species_percent'] = (df_phenotype_assignments_to_write['species_count'] / df_phenotype_assignments_to_write['species_count'].sum() * 100).apply(lambda x: np.round(x, decimals=8))

        # Write the dataframe to disk
        df_phenotype_assignments_to_write.to_csv(path_or_buf=os.path.join(phenotypes_path, filename), sep='\t', header=False)

        # Return the filename of the written file
        return filename

    def load_relevant_settings_from_phenotyper():
        st.session_state['settings__input_datafile__filename'] = st.session_state['datafileU']
        update_dependencies_of_input_datafile_filename()
        st.session_state['settings__input_datafile__coordinate_units'] = st.session_state['phenotyping_micron_coordinate_units']
        if st.session_state['phenoMeth'] == 'Custom':
            phenotyper_assignments_filename = create_phenotype_assignments_file_from_phenotyper(st.session_state['spec_summ_dataeditor'])
            st.session_state['settings__phenotyping__phenotype_identification_file'] = phenotyper_assignments_filename
        st.session_state['settings__phenotyping__method'] = st.session_state['phenoMeth']
        update_dependencies_of_phenotyping_method()

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

    # Display page heading
    st.title('Tool parameter selection')

    # Assign default keys
    streamlit_utils.assign_default_values_in_session_state('preset_parameter_file', utils.get_first_element_or_none(options_for_parameter_files))

    # Load or reload default settings
    if st.button('Reload default settings') or ('settings__input_datafile__filename' not in st.session_state):
        st.session_state['settings__input_datafile__filename'] = utils.get_first_element_or_none(options_for_input_datafiles)
        update_dependencies_of_input_datafile_filename()
        if 'settings__input_datafile__format' not in st.session_state:
            st.session_state['settings__input_datafile__format'] = utils.get_first_element_or_none(options_for_input_datafile_formats)
        st.session_state['settings__input_datafile__coordinate_units'] = 0.25
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
        settings = utils.validate_presets_and_map_to_settings(preset_settings, options_for_input_datafiles, options_for_input_datafile_formats, options_for_phenotype_identification_files, st.session_state['options_for_annotation_files'], options_for_phenotyping_methods, options_for_significance_calculation_methods, st.session_state['options_for_images'], message_function=streamlit_utils.streamlit_write_function)
        set_session_state_key(settings, 'input_datafile', 'filename')
        update_dependencies_of_input_datafile_filename()
        set_session_state_key(settings, 'input_datafile', 'format')
        set_session_state_key(settings, 'input_datafile', 'coordinate_units')
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

    st.button(':arrow_right: Load relevant settings from Phenotyper', on_click=load_relevant_settings_from_phenotyper)

    # Logical divider
    st.divider()

    # Create columns
    parameter_columns = st.columns(3)

    # Input datafile settings widgets
    with parameter_columns[0]:
        st.header('Input datafile settings')
        st.selectbox('Filename:', options_for_input_datafiles, key='settings__input_datafile__filename', help='An input datafile must be present in the "input" directory and have a .csv or .tsv extension.', on_change=update_dependencies_of_input_datafile_filename)
        st.selectbox('Format:', options_for_input_datafile_formats, key='settings__input_datafile__format')
        st.number_input('x-y coordinate units (microns):', min_value=0.0, key='settings__input_datafile__coordinate_units', help='E.g., if the coordinates in the input datafile were pixels, this number would be a conversion to microns in units of microns/pixel.', format='%.4f', step=0.0001)

    # Phenotyping settings widgets
        st.header('Phenotyping settings')
        st.selectbox('Method:', options_for_phenotyping_methods, key='settings__phenotyping__method', help='Species: phenotypes defined by the unique combinations of markers present in the input file. Marker: each marker is its own phenotype. Custom: custom phenotyping using a text file.', on_change=update_dependencies_of_phenotyping_method)
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

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)

if __name__ == '__main__':
    main()
