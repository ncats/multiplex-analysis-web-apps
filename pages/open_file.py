# Import relevant libraries
import streamlit as st
import app_top_of_page as top
import streamlit_dataframe_editor as sde
import os
import streamlit_utils
import utils

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
    input_dir = os.path.join('.', 'input')
    num_rows_to_sample = 100

    # Initialization
    show_dataframe_updates = False

    # Display a header
    st.header('Input options')

    # In the first third of the page, create the input file selection and loading widgets
    with st.columns(3)[0]:

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

        # Create a number input for the number of microns per coordinate unit
        if 'opener__microns_per_coordinate_unit' not in st.session_state:
            st.session_state['opener__microns_per_coordinate_unit'] = 1.0
        if st.session_state['opener__load_from_datafile_unifier']:
            help_message = 'Remember that the dataset coordinates were converted to microns in the Datafile Unifier.'
        else:
            help_message = None
        st.number_input('Enter the number of microns per coordinate unit in the input:', min_value=0.0, key='opener__microns_per_coordinate_unit', format='%.4f', step=0.0001, disabled=st.session_state['opener__load_from_datafile_unifier'], help=help_message)

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
        if st.button('Load the selected input dataset'):
            with st.spinner('Loading the input dataset...'):
                streamlit_utils.load_input_dataset(input_file_or_df, st.session_state['opener__microns_per_coordinate_unit'])  # this assigns the input dataset to st.session_state['input_dataset'] and the metadata to st.session_state['input_metadata']
            if st.session_state['input_dataset'] is not None:
                st.session_state['input_dataset'].data, st.session_state['input_dataframe_memory_usage_bytes'] = utils.downcast_dataframe_dtypes(st.session_state['input_dataset'].data, also_return_final_size=True)
                st.info('The input data have been successfully loaded and validated.')
                show_dataframe_updates = True
            else:
                st.error('The input data are in an unsupported format.')
                return

        # Stop if the input dataset is not yet loaded
        if ('input_dataset' not in st.session_state) or (st.session_state['input_dataset'] is None):
            st.info('To continue, please press the button above to load an input dataset.')
            return

    # Now that we know the input dataset is assigned and is valid, print out a sample of the main, input dataframe, plus some information about it
    st.header('Loaded dataset')

    # Assign shortcuts to the loaded dataset, metadata, and dataframe
    dataset_obj = st.session_state['input_dataset']
    metadata = st.session_state['input_metadata']
    df = dataset_obj.data  # this (st.session_state['input_dataset'].data) is the dataframe that should be used in the entire app suite

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

    # Display the information and the sampled dataframe
    st.markdown(information)
    st.write('Sampled dataframe:')
    resample_dataframe = st.button('Refresh dataframe sample')
    if ('opener__sampled_df' not in st.session_state) or resample_dataframe or show_dataframe_updates:
        st.session_state['opener__sampled_df'] = df.sample(min(num_rows_to_sample, len(df))).sort_index()
    st.write(st.session_state['opener__sampled_df'])

# Run the main function
if __name__ == '__main__':

    # Set page settings
    page_name = 'Open File'
    st.set_page_config(layout='wide', page_title=page_name)
    st.title(page_name)

    # Run streamlit-dataframe-editor library initialization tasks at the top of the page
    st.session_state = sde.initialize_session_state(st.session_state)

    # Run Top of Page (TOP) functions
    st.session_state = top.top_of_page_reqs(st.session_state)

    main()

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)
