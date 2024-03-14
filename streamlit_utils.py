import streamlit as st
from streamlit_extras.app_logo import add_logo

# Define a function to easily assign default values to items in the session state
def assign_default_values_in_session_state(key, value, args=None):
    if key not in st.session_state:
        if callable(value):
            if args is None:
                st.session_state[key] = value()
            else:
                st.session_state[key] = value(args)
        else:
            st.session_state[key] = value

# Function to go to a previous option using a button
def previous_option(options, selected_option='selected_image_id', button_name='image_button'):
    curr_index = options.index(st.session_state[selected_option])
    if curr_index != 0:
        curr_index = curr_index - 1
    st.session_state[selected_option] = options[curr_index]
    update_button_activities(options, selected_option=selected_option, button_name=button_name)

# Function to go to a next option using a button
def next_option(options, selected_option='selected_image_id', button_name='image_button'):
    curr_index = options.index(st.session_state[selected_option])
    if curr_index != (len(options) - 1):
        curr_index = curr_index + 1
    st.session_state[selected_option] = options[curr_index]
    update_button_activities(options, selected_option=selected_option, button_name=button_name)

# Function to determine whether previous or next buttons should be disabled
def update_button_activities(options, selected_option='selected_image_id', button_name='image_button'):
    curr_index = options.index(st.session_state[selected_option])
    if curr_index == 0:
        st.session_state['disable_previous_{}'.format(button_name)] = True
    else:
        st.session_state['disable_previous_{}'.format(button_name)] = False
    if curr_index == (len(options) - 1):
        st.session_state['disable_next_{}'.format(button_name)] = True
    else:
        st.session_state['disable_next_{}'.format(button_name)] = False

def propagate_settings_to_session_state(settings):
    for key0 in settings:
        for key1 in settings[key0]:
            current_setting = settings[key0][key1]
            combined_key = 'settings__' + key0 + '__' + key1
            st.session_state[combined_key] = current_setting

def get_current_settings():

    # Initialize the dictionary holding the actual settings dictionary to be used in the workflow
    settings = dict()
    settings['phenotyping'], settings['analysis'], settings['annotation'], settings['plotting'] = dict(), dict(), dict(), dict()

    # Assign current values of the settings
    for key in st.session_state:
        if key.startswith('settings__'):
            key0 = key.split('__')[1]
            key1 = key.split('__')[2]
            settings[key0][key1] = st.session_state[key]

    # Return the nested dictionary of settings
    return settings

def streamlit_write_function(message):
    st.error(message, icon="🚨")

def update_session_state_keys(transformation_dict):
    for key in transformation_dict:
        st.session_state[key] = transformation_dict[key]

def load_input_dataset(datafile_path_or_df, coord_units_in_microns, input_dataset_key='input_dataset', input_metadata_key='input_metadata'):
    """
    Load and standardize the input datafile and save the settings in the session state.

    This should be called from an "Open file(s)" page in the sidebar.

    Args:
        datafile_path_or_df (str): The path to the input datafile.
        coord_units_in_microns (float): The conversion factor from the units of the input datafile to microns.
        input_dataset_key (str): The key to be used for the input dataset in the session state. Defaults to 'input_dataset'.
        input_metadata_key (str): The key to be used for the input metadata in the session state. Defaults to 'input_metadata'.

    Returns:
        None
    """

    # Import relevant libraries
    import utils
    import os

    # Define the metadata for the input data
    if isinstance(datafile_path_or_df, str):
        metadata = {
            'datafile_path': os.path.realpath(datafile_path_or_df),
            'coord_units_in_microns': coord_units_in_microns
        }
    else:
        metadata = {
            'datafile_path': None,
            'coord_units_in_microns': coord_units_in_microns
        }

    # Load and standardize the input datafile into the session state
    dataset_obj = utils.load_and_standardize_input_datafile(datafile_path_or_df, coord_units_in_microns)

    # Save the input data and metadata to the session state
    if dataset_obj is not None:
        st.session_state[input_dataset_key] = dataset_obj
        st.session_state[input_metadata_key] = metadata  # save the metadata to the session state as well

    # If the input data is in an unsupported format, assign null values and display an error message
    else:
        st.session_state[input_dataset_key] = None
        st.session_state[input_metadata_key] = None

def get_updated_dynamic_options(input_directory):

    # options_for_images, options_for_annotation_files = get_updated_dynamic_options(input_directory, input_datafile_filename)

    input_datafile_filename = st.session_state['input_metadata']['datafile_path']
    
    # Import relevant library
    import os

    if input_datafile_filename is not None:

        # Update analysis__images_to_analyze options
        # options_for_images = get_unique_image_ids_from_datafile(input_datafile_path)
        # options_for_images = list(dataset_formats.get_image_series_in_datafile(input_datafile_path).unique())
        options_for_images = st.session_state['input_dataset'].data['Slide ID'].unique()

        # annotation__used_annotation_files options
        annotations_dir_listing = ([x for x in os.listdir(os.path.join(input_directory, 'annotations')) if x.endswith('.csv')] if os.path.exists(os.path.join(input_directory, 'annotations')) else [])
        options_for_annotation_files = [x for x in annotations_dir_listing if x.split('__')[0] in options_for_images]

        return options_for_images, options_for_annotation_files
    
    else:

        return None, None
