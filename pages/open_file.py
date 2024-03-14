# Import relevant libraries
import streamlit as st
import app_top_of_page as top
import streamlit_dataframe_editor as sde
import os

def toggle_changed():
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
    if 'opener__selected_input_file' in st.session_state:
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
    st.number_input('Enter the number of microns per coordinate unit in the input:', key='opener__microns_per_coordinate_unit', disabled=st.session_state['opener__load_from_datafile_unifier'], help=help_message)

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
        st.write(f'Selected input: {input_file_or_df}')
    else:
        st.write(f'Selected input: Dataset from Datafile Unifier consisting of {input_file_or_df.shape[0]} rows and {input_file_or_df.shape[1]} columns')
    st.write(f'Coordinate units: {st.session_state["opener__microns_per_coordinate_unit"]} microns')

    pass

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)

# Run the main function
if __name__ == '__main__':
    main()
