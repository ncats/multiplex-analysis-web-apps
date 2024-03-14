# Import relevant libraries
import streamlit as st
import app_top_of_page as top
import streamlit_dataframe_editor as sde
import os
import streamlit_utils
import pandas as pd

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

    # Set page settings
    page_name = 'Open File'
    st.set_page_config(layout='wide', page_title=page_name)
    st.title(page_name)

    # Run streamlit-dataframe-editor library initialization tasks at the top of the page
    st.session_state = sde.initialize_session_state(st.session_state)

    # Run Top of Page (TOP) functions
    st.session_state = top.top_of_page_reqs(st.session_state)

    # Constant
    input_dir = os.path.join('.', 'input')
    num_rows_to_sample = 100

    # Initialization
    show_dataframe_updates = False

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
        
        # Write out the selected input source, but if it is a dataframe, just say so
        if isinstance(input_file_or_df, str):
            st.write(f'Selected input: {os.path.basename(input_file_or_df)}')
        else:
            st.write(f'Selected input: Dataset from Datafile Unifier consisting of {input_file_or_df.shape[0]} rows and {input_file_or_df.shape[1]} columns')
        st.write(f'Coordinate units: {st.session_state["opener__microns_per_coordinate_unit"]} microns')

        # Read in the dataset either from memory or from a file
        if st.button('Load the selected input dataset'):
            streamlit_utils.load_input_dataset(input_file_or_df, st.session_state["opener__microns_per_coordinate_unit"])
            show_dataframe_updates = True

        # Stop if the input dataset is not yet loaded
        if 'input_dataset' not in st.session_state:
            st.info('To continue, please press the button above to load an input dataset.')
            return

        # Stop if the input dataset is in an unrecognized format
        if st.session_state['input_dataset'] is None:
            st.error('The input data are in an unsupported format.')
            return
    
    # Now that we know the input dataset st.session_state['input_dataset'] is assigned and is valid, print a sample of the main, input dataframe
    st.header('Sample of loaded dataframe')
    df = st.session_state['input_dataset'].data
    resample_dataframe = st.button('Refresh dataframe sample')
    if ('opener__sampled_df' not in st.session_state) or resample_dataframe or show_dataframe_updates:
        if not isinstance(df, pd.DataFrame):
            st.error('The input data is not a dataframe.')
            return
        else:
            st.session_state['opener__sampled_df'] = df.sample(min(num_rows_to_sample, len(df))).sort_index()
    st.write(st.session_state['opener__sampled_df'])
    st.write('The full loaded dataframe is part of the {} class and has {} rows and {} columns.'.format(type(st.session_state['input_dataset']), df.shape[0], df.shape[1]))

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)

# Run the main function
if __name__ == '__main__':
    main()
