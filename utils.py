def set_filename_corresp_to_roi(df_paths, roi_name, curr_colname, curr_dir, curr_dir_listing):
    """Update the path in a main paths-holding dataframe corresponding to a particular ROI in a particular directory.

    Args:
        df_paths (Pandas dataframe): Dataframe holding the absolute file paths.
        roi_name (str): Name of the particular ROI.
        curr_colname (str): Name of the corresponding column in the input dataframe df_paths.
        curr_dir (str): Directory in which to search for a file based on the current ROI.

    Returns:
        Pandas dataframe: Output df_paths filled in with the current file path (if it exists).
    """

    # Import relevant library
    import os

    # Obtain a list of the filenames corresponding to the input ROI
    filename_list = [x for x in curr_dir_listing if roi_name in x]

    # If no filenames were found, note this
    if len(filename_list) == 0:
        print('No filenames found for ROI {} in the directory {}'.format(roi_name, curr_dir))

    # If one filename was found, update the appropriate entry in df_paths
    elif len(filename_list) == 1:
        df_paths.loc[roi_name, curr_colname] = os.path.join(curr_dir, filename_list[0])

    # If more than one filename was found, then something may be wrong, and say so
    else:
        print('UH OH: More than one filename {} has been found for ROI {} in the directory {}'.format(filename_list, roi_name, curr_dir))

    # Return the updated paths dataframe
    return df_paths

def get_paths_for_rois():
    """Get the pathnames for all three types of ROI-based analyses.

    Args:
        radius_in_microns (int, optional): Distance in microns around each center cell. Defaults to 40.

    Returns:
        Pandas dataframe: Holder of absolute pathnames to all the ROI-based images.
    """

    # Import relevant libraries
    import os
    import pandas as pd

    # Obtain the directory holding the subdirectories containing various types of plots (in this case, three types)
    # plots_dir = os.path.join(os.getcwd(), '..', 'results', 'webpage', 'slices_1x{}'.format(radius_in_microns), 'real')
    plots_dir = os.path.join('.', 'output', 'images')

    # Obtain the paths to the subdirectories
    outlines_dir = os.path.join(plots_dir, 'single_roi_outlines_on_whole_slides')
    rois_dir = os.path.join(plots_dir, 'roi_plots')
    heatmaps_dir = os.path.join(plots_dir, 'dens_pvals_per_roi')

    # Get the directory listings
    outlines_dir_listing = os.listdir(outlines_dir)
    rois_dir_listing = os.listdir(rois_dir)
    heatmaps_dir_listing = os.listdir(heatmaps_dir)

    # Initialize an empty Pandas dataframe holding the full image pathnames where the index is the core ROI name
    df_paths = pd.DataFrame([x.rstrip('.png') for x in os.listdir(outlines_dir)], columns=['roi_name']).set_index('roi_name')

    # Determine the filenames in the various subdirectories corresponding to each ROI and store them in the df_paths dataframe
    df_paths['roi'] = ''
    df_paths['heatmap'] = ''
    df_paths['outline'] = ''
    for roi_name in df_paths.index:
        df_paths = set_filename_corresp_to_roi(df_paths=df_paths, roi_name=roi_name, curr_colname='roi', curr_dir=rois_dir, curr_dir_listing=rois_dir_listing)
        df_paths = set_filename_corresp_to_roi(df_paths=df_paths, roi_name=roi_name, curr_colname='heatmap', curr_dir=heatmaps_dir, curr_dir_listing=heatmaps_dir_listing)
        df_paths = set_filename_corresp_to_roi(df_paths=df_paths, roi_name=roi_name, curr_colname='outline', curr_dir=outlines_dir, curr_dir_listing=outlines_dir_listing)

    # Add columns containing the patient "case" ID and the slide "condition", in order to aid in sorting the data
    cases = []
    conditions = []
    for roi_name in df_paths.index:
        slide_id = roi_name.split('-')[0]
        cases.append(int(slide_id[:-1]))
        conditions.append(slide_id[-1])
    df_paths['case'] = cases
    df_paths['condition'] = conditions

    # Sort the data by case, then by condition, then by the ROI string
    df_paths = df_paths.sort_values(by=['case', 'condition', 'roi_name'])

    # Return the paths dataframe
    return df_paths

def get_paths_for_slides():
    """Get the pathnames for both types of slide-based analyses.

    5/13/23: Note I may want to implement the methodology in get_paths_for_rois() above, but for now the following, original way should suffice for the slides!

    Args:
        radius_in_microns (int, optional): Distance in microns around each center cell. Defaults to 40.

    Returns:
        Pandas dataframe: Holder of absolute pathnames to all the slide-based images.
    """

    # Import relevant libraries
    import os
    import pandas as pd

    # Obtain the directory holding the subdirectories containing various types of plots (in this case, two types)
    # plots_dir = os.path.join(os.getcwd(), '..', 'results', 'webpage', 'slices_1x{}'.format(radius_in_microns), 'real')
    plots_dir = os.path.join('.', 'output', 'images')

    # Obtain the paths to the subdirectories
    slides_dir = os.path.join(plots_dir, 'whole_slide_patches')
    heatmaps_dir = os.path.join(plots_dir, 'dens_pvals_per_slide')

    # List the contents of each directory
    slides_listing = os.listdir(slides_dir)
    slides_listing = [y for y in slides_listing if '-patched.png' not in y]
    heatmaps_listing = os.listdir(heatmaps_dir)

    # Initialize an empty Pandas dataframe holding the full image pathnames where the index is the core slide name
    df_paths = pd.DataFrame([x.rstrip('.png') for x in slides_listing], columns=['slide_name']).set_index('slide_name')

    # Determine the filenames of each of the image types corresponding to each slide name
    corresp_slide_filename = []
    corresp_slide_filename_patched = []
    corresp_heatmap_filename = []
    for slide_name in df_paths.index:
        for slide_filename, heatmap_filename in zip(slides_listing, heatmaps_listing):
            if slide_name in slide_filename:
                corresp_slide_filename.append(os.path.join(plots_dir, 'whole_slide_patches', slide_filename))
                corresp_slide_filename_patched.append(os.path.join(plots_dir, 'whole_slide_patches', '{}-patched.png'.format(slide_filename.rstrip('.png'))))
            if slide_name in heatmap_filename:
                corresp_heatmap_filename.append(os.path.join(plots_dir, 'dens_pvals_per_slide', heatmap_filename))

    # Add these paths to the main paths dataframe
    df_paths['slide'] = corresp_slide_filename
    df_paths['slide_patched'] = corresp_slide_filename_patched
    df_paths['heatmap'] = corresp_heatmap_filename

    # Add columns containing the patient "case" ID and the slide "condition", in order to aid in sorting the data
    cases = []
    conditions = []
    for slide_name in df_paths.index:
        slide_id = slide_name.split('-')[0]
        cases.append(int(slide_id[:-1]))
        conditions.append(slide_id[-1])
    df_paths['case'] = cases
    df_paths['condition'] = conditions

    # Sort the data by case, then by condition, then by the slide string
    df_paths = df_paths.sort_values(by=['case', 'condition', 'slide_name'])

    # Return the paths dataframe
    return df_paths

def get_overlay_info():
    """Get the items making up the overlay image filenames.

    Args:
        radius_in_microns (int, optional): Distance in microns around each center cell. Defaults to 40.

    Returns:
        dict: {str: Path to the overlays images, list: Correctly sorted slide names, list: Center species, list: Neighbor species, list: P value types}
    """

    # Import relevant libraries
    import os
    import pandas as pd

    # Obtain the directory holding the subdirectory of interest
    # plots_dir = os.path.join(os.getcwd(), '..', 'results', 'webpage', 'slices_1x{}'.format(radius_in_microns), 'real')
    plots_dir = os.path.join('.', 'output', 'images')

    # Obtain the path to the subdirectory
    overlays_dir = os.path.join(plots_dir, 'density_pvals_over_slide_spatial_plot')

    # List the contents of the directory
    overlays_listing = os.listdir(overlays_dir)

    # Create an empty dataframe containing the slide names as the index
    df_paths = pd.DataFrame(set([x.split('-with_log_dens_pvals_per_roi')[0] for x in overlays_listing]), columns=['slide_name']).set_index('slide_name')

    # Add columns containing the patient "case" ID and the slide "condition", in order to aid in sorting the data
    cases = []
    conditions = []
    for slide_name in df_paths.index:
        slide_id = slide_name.split('-')[0]
        cases.append(int(slide_id[:-1]))
        conditions.append(slide_id[-1])
    df_paths['case'] = cases
    df_paths['condition'] = conditions

    # Sort the data by case, then by condition, then by the slide string, returning the index as a list (i.e., correctly sorted slide names)
    slide_names = df_paths.sort_values(by=['case', 'condition', 'slide_name']).index.to_list()

    # Determine the rest of the unique identifiers making up the overlay image names
    center_species = list(set([x.split('__center_')[1].split('__neighbor_')[0] for x in overlays_listing]))
    center_species.sort()
    neighbor_species = list(set([x.split('__neighbor_')[1].split('__')[0] for x in overlays_listing]))
    neighbor_species.sort()
    pval_types = list(set([x.split('__')[3].split('_pvals.png')[0] for x in overlays_listing]))
    pval_types.sort()

    # Return the needed path and lists as a single dictionary
    return {'overlay_dir': overlays_dir, 'slide_names': slide_names, 'center_species': center_species, 'neighbor_species': neighbor_species, 'pval_types': pval_types}

def string_has_bad_values(string, bad_values=[' ', 'XMin', 'XMax', 'YMin', 'YMax']):
    is_bad = False
    for bad_value in bad_values:
        if bad_value in string:
            is_bad = True
    return is_bad

def detect_markers_in_annotation_files(selected_annotation_files):

    # Import relevant libraries
    import os
    import dataset_formats

    # Constant
    annotations_dir = os.path.join('input', 'annotations')

    # Initialize an empty list to hold all the annotation markers in the annotation files
    annotation_markers = []

    # For every file in the annotations directory...
    for annotation_file in selected_annotation_files:

        # Obtain the markers in the current annotation file
        _, _, _, _, _, current_annot_markers = dataset_formats.extract_datafile_metadata(os.path.join(annotations_dir, annotation_file))

        # For each of these column names...
        for current_annot_marker in current_annot_markers:

            # If it's not already in the running list, add it to the running list
            if current_annot_marker not in annotation_markers:
                annotation_markers.append(current_annot_marker)

    # Return the running list of markers in the annotation files
    return annotation_markers

# def validate_presets_and_map_to_settings(preset_settings, possible_input_datafiles, possible_input_datafile_formats, possible_phenotype_identification_files, markers_in_annotation_files, possible_annotation_files, possible_phenotyping_methods, possible_significance_calculation_methods, options_for_images, message_function=print):
def validate_presets_and_map_to_settings(preset_settings, possible_input_datafiles, possible_input_datafile_formats, possible_phenotype_identification_files, possible_annotation_files, possible_phenotyping_methods, possible_significance_calculation_methods, options_for_images, message_function=print):

    # Import relevant libraries
    import sys
    import os

    # Make it easy to do type checking
    def isnum(val):
        return isinstance(val, (int, float))
    def isstr(val):
        return isinstance(val, str)
    def isbool(val):
        return isinstance(val, bool)
    def isint(val):
        return isinstance(val, int)
    def istuple(val):
        return isinstance(val, tuple)
    def islist(val):
        return isinstance(val, list)
    
    # Validate a single preset setting based on a validation function, showing an error message and quitting if the check fails, and returning the value if it passes
    def validate_preset_setting(preset_settings, section, key, is_valid, error_message, mapper=None, transformation=None, message_function=message_function):
        val = preset_settings[section][key]
        if not is_valid(val):
            message_function('Error in preset setting [{}][{}]: {}'.format(section, key, error_message.replace('Value is not', 'Value ({}) is not'.format(val), 1)))
            sys.exit()
        else:
            if mapper is not None:  # allow for making an optional mapping
                val = mapper[val]
            if transformation is not None:  # allow for making an optional transformation
                val = transformation(val)
            return val
        
    # Constants
    orig_possible_dataset_formats = ['OMAL', 'Native', 'GMBSecondGeneration']  # current dataset types defined in dataset_formats.py
    input_directory = os.path.join('.', 'input')

    # Initialize the dictionary holding the actual settings dictionary to be used in the workflow
    settings = dict()
    settings['input_datafile'], settings['phenotyping'], settings['analysis'], settings['annotation'], settings['plotting'] = dict(), dict(), dict(), dict(), dict()

    # If the input settings file is the new format...
    if 'input_datafile' in preset_settings:

        # Validate the settings in the input_datafile (new format) section
        settings['input_datafile']['filename'] = validate_preset_setting(preset_settings, 'input_datafile', 'filename', lambda x: isstr(x) and (x in possible_input_datafiles), 'Value is not present in the "input" directory (present values are [{}])'.format(possible_input_datafiles))
        options_for_images, possible_annotation_files = get_updated_dynamic_options(input_directory, settings['input_datafile']['filename'])
        settings['input_datafile']['format'] = validate_preset_setting(preset_settings, 'input_datafile', 'format', lambda x: isstr(x) and (x in possible_input_datafile_formats), 'Value is not an acceptable value (options are {})'.format(possible_input_datafile_formats))
        settings['input_datafile']['coordinate_units'] = validate_preset_setting(preset_settings, 'input_datafile', 'coordinate_units', lambda x: isnum(x) and (x > 0), 'Value is not a positive number')

        # Validate the settings in the phenotyping (new format) section
        settings['phenotyping']['method'] = validate_preset_setting(preset_settings, 'phenotyping', 'method', lambda x: isstr(x) and (x in possible_phenotyping_methods), 'Value is not an acceptable value (options are {})'.format(possible_phenotyping_methods))
        if settings['phenotyping']['method'] == 'Custom':
            settings['phenotyping']['phenotype_identification_file'] = validate_preset_setting(preset_settings, 'phenotyping', 'phenotype_identification_file', lambda x: isstr(x) and (x in possible_phenotype_identification_files), 'Value is not present in the "input" directory (present values are [{}])'.format(possible_phenotype_identification_files))

        # Validate the settings in the analysis (new format) section
        settings['analysis']['images_to_analyze'] = validate_preset_setting(preset_settings, 'analysis', 'images_to_analyze', lambda x: islist(x) and (sum([(y in options_for_images) for y in x]) == len(x)), 'Value is not a list containing images that are present in the selected input datafile (options are {})'.format(options_for_images))
        settings['analysis']['partition_slides_into_rois'] = validate_preset_setting(preset_settings, 'analysis', 'partition_slides_into_rois', lambda x: isbool(x), 'Value is not a boolean')
        if settings['analysis']['partition_slides_into_rois']:
            settings['analysis']['roi_width'] = validate_preset_setting(preset_settings, 'analysis', 'roi_width', lambda x: isnum(x) and (x > 0), 'Value is not a positive number')
            settings['analysis']['roi_overlap'] = validate_preset_setting(preset_settings, 'analysis', 'roi_overlap', lambda x: isnum(x), 'Value is not a number')
        settings['analysis']['neighbor_radius'] = validate_preset_setting(preset_settings, 'analysis', 'neighbor_radius', lambda x: isnum(x) and (x > 0), 'Value is not a positive number')
        settings['analysis']['significance_calculation_method'] = validate_preset_setting(preset_settings, 'analysis', 'significance_calculation_method', lambda x: isstr(x) and (x in possible_significance_calculation_methods), 'Value is not an acceptable value (options are {})'.format(possible_significance_calculation_methods))
        settings['analysis']['n_neighs'] = validate_preset_setting(preset_settings, 'analysis', 'n_neighs', lambda x: isint(x) and (x >= 0), 'Value is not a non-negative integer')
        settings['analysis']['min_num_valid_centers'] = validate_preset_setting(preset_settings, 'analysis', 'min_num_valid_centers', lambda x: isint(x) and (x >= 1), 'Value is not a positive integer')
        settings['analysis']['weight_by_num_valid_centers'] = validate_preset_setting(preset_settings, 'analysis', 'weight_by_num_valid_centers', lambda x: isbool(x), 'Value is not a boolean')
        settings['analysis']['log_pval_minimum'] = validate_preset_setting(preset_settings, 'analysis', 'log_pval_minimum', lambda x: isnum(x) and (x < 0), 'Value is not a negative number')
        settings['analysis']['log_pval_maximum'] = validate_preset_setting(preset_settings, 'analysis', 'log_pval_maximum', lambda x: isnum(x) and (x <= 0), 'Value is not a non-positive number')

        # Validate the settings in the annotation (new format) section
        settings['annotation']['used_annotation_files'] = validate_preset_setting(preset_settings, 'annotation', 'used_annotation_files', lambda x: islist(x) and (sum([(y in possible_annotation_files) for y in x]) == len(x)), 'Value is not a list containing filenames that are present in the input/annotations directory (options are {})'.format(possible_annotation_files))
        if len(settings['annotation']['used_annotation_files']) > 0:
            settings['annotation']['coordinate_units'] = validate_preset_setting(preset_settings, 'annotation', 'coordinate_units', lambda x: isnum(x) and (x > 0), 'Value is not a positive number')
            settings['annotation']['coord_units_are_pixels'] = validate_preset_setting(preset_settings, 'annotation', 'coord_units_are_pixels', lambda x: isbool(x), 'Value is not a boolean')
            if settings['annotation']['coord_units_are_pixels']:
                settings['annotation']['microns_per_integer_unit'] = validate_preset_setting(preset_settings, 'annotation', 'microns_per_integer_unit', lambda x: isnum(x) and (x > 0), 'Value is not a positive number')
            # settings['annotation']['markers_designating_valid_objects'] = validate_preset_setting(preset_settings, 'annotation', 'markers_designating_valid_objects', lambda x: islist(x) and (sum([(y in markers_in_annotation_files) for y in x]) == len(x)), 'Value is not a list containing marker strings that are present in the annotation files (options are {})'.format(markers_in_annotation_files))

        settings['plotting']['min_log_pval'] = validate_preset_setting(preset_settings, 'plotting', 'min_log_pval', lambda x: isnum(x) and (x < 0), 'Value is not a negative number')

    # If the input settings file is in the original format...
    else:

        # Validate the settings in the input_datafile (new format) section
        settings['input_datafile']['filename'] = validate_preset_setting(preset_settings, 'dataset', 'input_datafile', lambda x: isstr(x) and (x in possible_input_datafiles), 'Value is not a string or is not present in the "input" directory (present values are [{}])'.format(possible_input_datafiles))
        options_for_images, possible_annotation_files = get_updated_dynamic_options(input_directory, settings['input_datafile']['filename'])
        settings['input_datafile']['format'] = validate_preset_setting(preset_settings, 'dataset', 'format', lambda x: isstr(x) and (x in orig_possible_dataset_formats), 'Value is not a string or is not in not an acceptable value (options are {})'.format(orig_possible_dataset_formats), mapper=dict(zip(orig_possible_dataset_formats, possible_input_datafile_formats)))
        settings['input_datafile']['coordinate_units'] = validate_preset_setting(preset_settings, 'dataset', 'coord_units_in_microns', lambda x: isnum(x) and (x > 0), 'Value is not a positive number')

        # Validate the settings in the phenotyping (new format) section
        preset__dataset__phenotype_identification_tsv_file = validate_preset_setting(preset_settings, 'dataset', 'phenotype_identification_tsv_file', lambda x: (x is None) or (isstr(x) and (x in possible_phenotype_identification_files)), 'Value is not None or is not a string that is present in the "input" directory (present values are [{}])'.format(possible_phenotype_identification_files))
        preset__analysis__allow_compound_species = validate_preset_setting(preset_settings, 'analysis', 'allow_compound_species', lambda x: isbool(x), 'Value is not a boolean')
        if (preset__dataset__phenotype_identification_tsv_file is None) and (preset__analysis__allow_compound_species is True):
            settings['phenotyping']['method'] = 'Species'
            # settings['phenotyping']['phenotype_identification_file'] = None
        elif (preset__dataset__phenotype_identification_tsv_file is None) and (preset__analysis__allow_compound_species is False):
            settings['phenotyping']['method'] = 'Marker'
            # settings['phenotyping']['phenotype_identification_file'] = None
        elif (preset__dataset__phenotype_identification_tsv_file is not None) and (preset__analysis__allow_compound_species is True):
            settings['phenotyping']['method'] = 'Custom'
            settings['phenotyping']['phenotype_identification_file'] = preset__dataset__phenotype_identification_tsv_file
        else:
            print('ERROR: Bad combination of original input parameters')
            sys.exit()

        # Validate the settings in the analysis (new format) section
        settings['analysis']['images_to_analyze'] = options_for_images
        preset__dataset__roi_width = validate_preset_setting(preset_settings, 'dataset', 'roi_width', lambda x: (x is None) or (isnum(x) and (x > 0)), 'Value is not None or is <= 0')
        preset__dataset__overlap = validate_preset_setting(preset_settings, 'dataset', 'overlap', lambda x: isnum(x) and (x >= 0), 'Value is not a non-negative number')
        if preset__dataset__roi_width is None:
            settings['analysis']['partition_slides_into_rois'] = False
        else:  # preset__dataset__roi_width is not None
            settings['analysis']['partition_slides_into_rois'] = True
            settings['analysis']['roi_width'] = preset__dataset__roi_width
            settings['analysis']['roi_overlap'] = preset__dataset__overlap
        settings['analysis']['neighbor_radius'] = validate_preset_setting(preset_settings, 'analysis', 'thickness', lambda x: isnum(x) and (x > 0), 'Value is not a positive number')
        settings['analysis']['significance_calculation_method'] = validate_preset_setting(preset_settings, 'analysis', 'use_analytical_significance', lambda x: isbool(x), 'Value is not a boolean', mapper={True: 'Poisson (radius)', False: 'Permutation (radius)'})
        settings['analysis']['n_neighs'] = validate_preset_setting(preset_settings, 'analysis', 'n_neighs', lambda x: isint(x) and (x >= 0), 'Value is not a non-negative integer')
        settings['analysis']['min_num_valid_centers'] = validate_preset_setting(preset_settings, 'plotting', 'num_valid_centers_minimum', lambda x: isint(x) and (x >= 1), 'Value is not a positive integer')
        settings['analysis']['weight_by_num_valid_centers'] = validate_preset_setting(preset_settings, 'plotting', 'weight_rois_by_num_valid_centers', lambda x: isbool(x), 'Value is not a boolean')
        preset__plotting__log_pval_range = validate_preset_setting(preset_settings, 'plotting', 'log_pval_range', lambda x: istuple(x) and (len(x) == 2) and isnum(x[0]) and (x[0] < 0) and isnum(x[1]) and (x[1] <= 0), 'Value is not a two-element tuple of numbers < 0 and <= 0, respectively')
        settings['analysis']['log_pval_minimum'] = preset__plotting__log_pval_range[0]
        settings['analysis']['log_pval_maximum'] = preset__plotting__log_pval_range[1]

        # Validate the settings in the annotation (new format) section
        settings['annotation']['used_annotation_files'] = possible_annotation_files
        if len(settings['annotation']['used_annotation_files']) != 0:
            settings['annotation']['coordinate_units'] = validate_preset_setting(preset_settings, 'annotation', 'annotation_coord_units_in_microns', lambda x: isnum(x) and (x > 0), 'Value is not a positive number')
            settings['annotation']['coord_units_are_pixels'] = False
            settings['annotation']['microns_per_integer_unit'] = validate_preset_setting(preset_settings, 'annotation', 'annotation_microns_per_integer_unit', lambda x: isnum(x) and (x > 0), 'Value is not a positive number')

        settings['plotting']['min_log_pval'] = preset__plotting__log_pval_range[0]

    # Return the nested dictionary of settings
    return settings

def get_first_element_or_none(options):
    if len(options) > 0:
        return options[0]
    else:
        return None

# def get_settings_defaults(options_for_input_datafiles, options_for_phenotype_identification_files, options_for_annotation_files, options_for_input_datafile_formats, options_for_phenotyping_methods, options_for_significance_calculation_methods, options_for_markers_in_annotation_files, options_for_images):
def get_settings_defaults(options_for_input_datafiles, options_for_phenotype_identification_files, options_for_annotation_files, options_for_input_datafile_formats, options_for_phenotyping_methods, options_for_significance_calculation_methods, options_for_images):

    # Initialize the dictionary holding the actual settings dictionary to be used in the workflow
    settings = dict()
    settings['input_datafile'], settings['phenotyping'], settings['analysis'], settings['annotation'] = dict(), dict(), dict(), dict()

    # Assign defaults
    settings['input_datafile']['filename'] = get_first_element_or_none(options_for_input_datafiles)

    import dataset_formats
    import os
    if settings['input_datafile']['filename'] is not None:
        settings['input_datafile']['format'] = dataset_formats.extract_datafile_metadata(os.path.join('.', 'input', settings['input_datafile']['filename']))[4]
    else:
        settings['input_datafile']['format'] = get_first_element_or_none(options_for_input_datafile_formats)

    settings['input_datafile']['coordinate_units'] = 0.25
    settings['phenotyping']['method'] = get_first_element_or_none(options_for_phenotyping_methods)
    settings['phenotyping']['phenotype_identification_file'] = get_first_element_or_none(options_for_phenotype_identification_files)
    settings['analysis']['images_to_analyze'] = options_for_images
    settings['analysis']['partition_slides_into_rois'] = True
    settings['analysis']['roi_width'] = 400
    settings['analysis']['roi_overlap'] = 80
    settings['analysis']['neighbor_radius'] = 40
    settings['analysis']['significance_calculation_method'] = get_first_element_or_none(options_for_significance_calculation_methods)
    settings['analysis']['min_num_valid_centers'] = 10
    settings['analysis']['weight_by_num_valid_centers'] = False
    settings['analysis']['log_pval_minimum'] = -50
    settings['analysis']['log_pval_maximum'] = 0
    settings['annotation']['used_annotation_files'] = options_for_annotation_files
    settings['annotation']['coordinate_units'] = 0.25
    settings['annotation']['coord_units_are_pixels'] = True
    settings['annotation']['microns_per_integer_unit'] = 0.25
    # settings['annotation']['markers_designating_valid_objects'] = options_for_markers_in_annotation_files

    # Return the nested dictionary of settings
    return settings

def get_unique_image_ids_from_datafile(datafile_path):

    # Import relevant libraries
    import pandas as pd
    import numpy as np
    import dataset_formats

    # Obtain the image number extraction parameters
    relevant_column_str, string_processing_func, _, _, _, _ = dataset_formats.extract_datafile_metadata(datafile_path)

    # If the datafile is some recognized file format...
    if relevant_column_str is not None:

        sep = (',' if datafile_path.split('.')[-1] == 'csv' else '\t')

        # Efficiently get just the relevant column of the datafile
        df = pd.read_csv(datafile_path, usecols=[relevant_column_str], sep=sep)

        # Extract from it the unique image IDs
        return [y.strip() for y in sorted([string_processing_func(x) for x in df[relevant_column_str].unique()])]

    # If the datafile is not a recognized file format...
    else:

        # Return None
        return None

def get_updated_dynamic_options(input_directory, input_datafile_filename):

    # options_for_images, options_for_annotation_files = get_updated_dynamic_options(input_directory, input_datafile_filename)
    
    # Import relevant libraries
    import os

    if input_datafile_filename is not None:

        # Set full pathname to input datafile
        input_datafile_path = os.path.join(input_directory, input_datafile_filename)

        # Update analysis__images_to_analyze options
        options_for_images = get_unique_image_ids_from_datafile(input_datafile_path)

        # annotation__used_annotation_files options
        annotations_dir_listing = ([x for x in os.listdir(os.path.join(input_directory, 'annotations')) if x.endswith('.csv')] if os.path.exists(os.path.join(input_directory, 'annotations')) else [])
        options_for_annotation_files = [x for x in annotations_dir_listing if x.split('__')[0] in options_for_images]

        return options_for_images, options_for_annotation_files
    
    else:

        return None, None

# Function to determine whether an n x n matrix is symmetric
def is_symmetric(matrix):
    import numpy as np
    return np.allclose(matrix, matrix.T, equal_nan=True)

# Function to determine whether all matrices corresponding to the first two dimensions of an arbitrarily large np.ndarray are symmetric
def is_symmetric_all(arr):
    
    # Import relevant library
    import numpy as np

    # Get the shape of the full ndarray
    arr_shape = arr.shape

    # Get the shape of the ndarray beyond the first two dimensions
    rest_of_dims = arr_shape[2:]

    # If there are dimensions beyond the first two (e.g., n x n x 1, n x n x 2 x 9, etc.)...
    if len(rest_of_dims) >= 1:

        # Assume all contained matrices are symmetric until proven otherwise
        all_are_symmetric = True

        # Get the total length of the other dimensions, calling it "m" below
        rest_of_dims_prod = np.prod(rest_of_dims)

        # Get a version of the original array reshaped to n x n x m
        arr_reshaped = arr.reshape(arr_shape[:2] + (rest_of_dims_prod,))

        # For every matrix in the array...
        for imatrix in range(rest_of_dims_prod):

            # Get the current (strictly n x n) matrix
            curr_matrix = arr_reshaped[:, :, imatrix]

            # If the current matrix is not symmetric, set the overall "are all symmetric" boolean to false, and get out of the loop
            if not is_symmetric(curr_matrix):
                all_are_symmetric = False
                break
    
    # If there are only two dimensions...
    else:

        # Set the overall "all are symmetric" boolean to whether the one single matrix is symmetric
        all_are_symmetric = is_symmetric(arr)

    # Return the boolean declaring whether all contained matrices are symmetric
    return all_are_symmetric

# def execute_data_parallelism_potentially(function=(lambda x: x), list_of_tuple_arguments=[(4444)], nworkers=0, task_description='', do_benchmarking=False, mp_start_method=None):
def execute_data_parallelism_potentially(function=(lambda x: x), list_of_tuple_arguments=[(4444,)], nworkers=0, task_description='', do_benchmarking=False, mp_start_method='forkserver'):  # spawn works with name=main block in Home.py
    # Note I forced mp_start_method = 'spawn' up until 4/27/23. Removing that and letting Python choose the default for the OS got parallelism working on NIDAP. I likely forced it to be spawn a long time ago maybe to get it working on Biowulf or my laptop or something like that. This worked in all scenarios including on my laptop (in WSL) though I get weird warnings I believe. I got confident about doing it this most basic way on 4/27/23 after reading Goyo's 2/7/23 example [here](https://discuss.streamlit.io/t/streamlit-session-state-with-multiprocesssing/29230/2) showing the same exact method I've been using except for forcing multiprocessing to use the "spawn" start method.

    # Import relevant library
    import multiprocessing as mp
    import time

    # Determine whether we're requesting to use the multiprocessing module in the first place
    if nworkers == 0:
        use_multiprocessing = False
    else:
        use_multiprocessing = True
        if mp_start_method is None:
            mp_start_method = mp.get_start_method()

    # Record the start time
    if do_benchmarking:
        start_time = time.time()

    # Farm out the function execution to multiple CPUs on different parts of the data
    if use_multiprocessing:
        print('Running {} function calls using {} workers for the {}'.format(len(list_of_tuple_arguments), nworkers, task_description))
        with mp.get_context(mp_start_method).Pool(nworkers) as pool:
            pool.map(function, list_of_tuple_arguments)

    # Execute fully in serial without use of the multiprocessing module
    else:
        print('NOT using multiprocessing for the {}'.format(task_description))
        for args_as_single_tuple in list_of_tuple_arguments:
            function(args_as_single_tuple)

    # Output how long the task took
    if do_benchmarking:
        elapsed_time = time.time() - start_time
        print('BENCHMARKING: The task took {} seconds using {} CPU(s) {} hyperthreading'.format(elapsed_time, (nworkers if use_multiprocessing else 1), ('WITH' if use_multiprocessing else 'WITHOUT')))