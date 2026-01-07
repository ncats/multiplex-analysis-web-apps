'''
Python Script which produces the "Open File" page of the Multiplex Analysis Web Apps.
'''

import os
import streamlit as st
import streamlit_utils
import utils
import framework.utils as framework_utils
import framework.manage_sessions as manage_sessions

ST_KEY_PREFIX_APP = "app.py__"
ST_KEY_PREFIX_STARTUP = "startup.py__"


def reset_session_state(extra_keys_to_keep=[], delete_input_dir=True):
    """Load the session state from the session directory."""
    try:
        # Back up app session-specific (i.e., startup.py-defined) variables we ultimately don't want to overwrite.
        keys_to_keep = [ST_KEY_PREFIX_STARTUP + "app_session_id", "previous_page_name", "current_page_name", ST_KEY_PREFIX_APP + "app_initialized", "platform"]

        # Delete everything in the session state but the keys to keep.
        for key in list(st.session_state.keys()):
            if (key not in keys_to_keep) and (key not in extra_keys_to_keep):
                del st.session_state[key]

        # Delete everything from the input and output directories.
        if delete_input_dir:
            framework_utils.ensure_empty_directory(os.path.join(framework_utils.session_dir(), "input"))
        framework_utils.ensure_empty_directory(os.path.join(framework_utils.session_dir(), "output"))

        return True
    except Exception as e:
        st.error(f"Failed to reset session state: {e}")
        return False


def clear_session_state():
    '''
    Function to clear all session state variables 
    except for the ones that are needed for the app to run.
    '''
    keys_to_keep = ['session_selection']
    for key in st.session_state.keys():
        if key.startswith(('unifier__', 'opener__')):
            keys_to_keep.append(key)
    reset_session_state(extra_keys_to_keep=keys_to_keep, delete_input_dir=False)
    

def load_input_dataset():
    '''
    Function to run when the "Load the selected input dataset" button is clicked.
    '''
    clear_session_state()
    st.session_state['opener__load_input_dataset'] = True
    st.session_state['input_dataset'] = None


def toggle_changed():
    """
    Function to run when the "Load dataset from Datafile Unifier" toggle is changed.
    """
    if st.session_state['opener__load_from_datafile_unifier']:
        st.session_state['opener__opener__selected_input_file'] = None
        st.session_state['opener__microns_per_coordinate_unit'] = 1.0

def main():
    """
    Main function for the Open File page.
    """

    # Constant
    input_dir = os.path.join(framework_utils.session_dir(), 'input')
    num_rows_to_sample = 100

    # Initialization
    show_dataframe_updates = False

    # In the first third of the page, create the input file selection and loading widgets
    open_data_cols = st.columns([2, 3])
    with open_data_cols[0]:

        # Display a header
        st.header('Input options')

        # Create a toggle to load the dataset from Datafile Unifier
        if 'opener__load_from_datafile_unifier' not in st.session_state:
            if 'unifier__df' in st.session_state:
                st.session_state['opener__load_from_datafile_unifier'] = True
                toggle_changed()
            else:
                st.session_state['opener__load_from_datafile_unifier'] = False
        st.toggle('Load dataset from Datafile Unifier', key='opener__load_from_datafile_unifier', disabled=('unifier__df' not in st.session_state), on_change=toggle_changed)

        # Create a dropdown of the available files in the "input" directory
        available_input_files = [file for file in os.listdir(input_dir) if file.lower().endswith(('.csv', '.tsv'))]
        if 'opener__selected_input_file' not in st.session_state:
            if available_input_files:
                st.session_state['opener__selected_input_file'] = available_input_files[0]
        else:
            if st.session_state['opener__selected_input_file'] not in available_input_files:
                st.session_state['opener__selected_input_file'] = None
        st.selectbox('Select an available input file to load:', options=available_input_files, key='opener__selected_input_file', disabled=st.session_state['opener__load_from_datafile_unifier'])

        # Set the number of microns per coordinate unit
        if 'opener__microns_per_coordinate_unit' not in st.session_state:
            st.session_state['opener__microns_per_coordinate_unit'] = 1.0

        # Determine the input datafile or input dataframe
        if st.session_state['opener__load_from_datafile_unifier']:
            if 'unifier__df' in st.session_state:
                input_file_or_df = st.session_state['unifier__df']
            else:
                input_file_or_df = None
        else:
            if st.session_state['opener__selected_input_file'] is not None:
                input_file_or_df = os.path.join(input_dir, st.session_state['opener__selected_input_file'])
            else:
                input_file_or_df = None

        # Stop if no input file or dataframe is available
        if input_file_or_df is None:
            st.info('Please select a valid input source.')
            return
        
        # Read in the dataset either from memory or from a file
        if 'input_dataset' in st.session_state:
            help_message = 'WARNING: This will clear all downstream analyses. If you want to keep them, please save them first using the Data Import and Export tool at left.'
            button_text = '⚠️ Load the selected input dataset'
        else:
            help_message = None
            button_text = 'Load the selected input dataset'
        st.button(button_text, on_click=load_input_dataset, help=help_message)

        # Not using a callback because of st.spinner() and the return statement
        if 'opener__load_input_dataset' not in st.session_state:
            st.session_state['opener__load_input_dataset'] = False
        if st.session_state['opener__load_input_dataset']:
            st.session_state['opener__load_input_dataset'] = False
            with st.spinner('Loading the input dataset...'):
                streamlit_utils.load_input_dataset(input_file_or_df, st.session_state['opener__microns_per_coordinate_unit'])  # this assigns the input dataset to st.session_state['input_dataset'] and the metadata to st.session_state['input_metadata']
                if st.session_state['input_dataset'] is not None:
                    st.session_state['input_dataset'].data, st.session_state['input_dataframe_memory_usage_bytes'] = utils.downcast_dataframe_dtypes(st.session_state['input_dataset'].data, also_return_final_size=True)
                    # st.session_state['adata'] = utils.create_anndata_from_dataframe(st.session_state['input_dataset'].data, columns_for_data_matrix='float')
            if st.session_state['input_dataset'] is not None:
                st.info('The input data have been successfully loaded and validated.')
                show_dataframe_updates = True
            else:
                st.error('The input data are in an unsupported format. As long as you use the Datafile Unifier first (even for just a single file), this should work.')
                return

        # Stop if the input dataset is not yet loaded
        if ('input_dataset' not in st.session_state) or (st.session_state['input_dataset'] is None):
            st.info('To continue, please press the button above to load an input dataset.')
            return

    with open_data_cols[1]:
        # Now that we know the input dataset is assigned and is valid,
        # print out a sample of the main, input dataframe, plus some information about it
        st.header('Loaded dataset')

        # Assign shortcuts to the loaded dataset, metadata, and dataframe
        dataset_obj = st.session_state['input_dataset']
        metadata = st.session_state['input_metadata']

        # This dataframe (st.session_state['input_dataset'].data) 
        # is the dataframe that should be used throughout ehe entire MAWA suite
        df = dataset_obj.data

        # Get information about the loaded dataset
        if metadata['datafile_path'] is not None:
            datafile_path = os.path.basename(metadata['datafile_path'])
        else:
            datafile_path = 'loaded from Datafile Unifier'
        information = f'''
        Properties:

        :small_orange_diamond: Datafile: `{datafile_path}`  
        :small_orange_diamond: Coordinate units: `{metadata['coord_units_in_microns']} microns/coord`  
        :small_orange_diamond: Dataset format: `{type(dataset_obj)}`  
        :small_orange_diamond: Number of rows: `{df.shape[0]}`  
        :small_orange_diamond: Number of columns: `{df.shape[1]}`  
        :small_orange_diamond: Minimum coordinate spacing: `{dataset_obj.min_coord_spacing_:.4f} microns`  
        :small_orange_diamond: Loaded memory usage: `{st.session_state['input_dataframe_memory_usage_bytes'] / 1024 ** 2:.2f} MB`  
        '''

        # Display more information if the dataset has been preprocessed inside MAWA
        if metadata['preprocessing'] is not None:
            information += f':small_orange_diamond: Preprocessing: `{metadata["preprocessing"]}`'

        # Display the information and the sampled dataframe
        st.markdown(information)

    st.header('Dataframe sample')
    resample_dataframe = st.button('Refresh dataframe sample')
    if ('opener__sampled_df' not in st.session_state) or resample_dataframe or show_dataframe_updates:
        st.session_state['opener__sampled_df'] = utils.sample_df_without_replacement_by_number(df=df, n=num_rows_to_sample).sort_index()
    st.write(st.session_state['opener__sampled_df'])

# Run the main function
if __name__ == '__main__':
    main()
