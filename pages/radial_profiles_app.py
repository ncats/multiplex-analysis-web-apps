# Import relevant libraries
import streamlit as st
import app_top_of_page as top
import streamlit_dataframe_editor as sde
import radial_profiles
import time
import numpy as np


def main():
    """
    Main function for the page.
    """

    # Global variable
    st_key_prefix = 'radial_profiles__'

    # Ensure a dataset has been opened in the first place
    if 'input_dataset' not in st.session_state:
        st.warning('Please open a dataset from the Open File page at left.')
        return

    # Save a shortcut to the dataframe
    df = st.session_state['input_dataset'].data

    # On the first third vertical third of the page...
    with st.columns(3)[0]:
    
        # Checkbox for whether to run checks
        key = st_key_prefix + '__run_checks'
        if key not in st.session_state:
            st.session_state[key] = False
        st.checkbox('Run checks', key=key)

        # Number input for the threshold for the RawIntNorm check
        key = st_key_prefix + '__perc_thresh_rawintnorm_column_check'
        if key not in st.session_state:
            st.session_state[key] = 0.01
        if st.session_state[st_key_prefix + '__run_checks']:
            st.number_input('Threshold for the RawIntNorm column check (%):', min_value=0.0, max_value=100.0, key=key)

        # Number input to select the nuclear intensity channel
        key = st_key_prefix + '__nuclear_channel'
        if key not in st.session_state:
            st.session_state[key] = 2
        st.number_input('Nuclear channel:', min_value=1, key=key)

        # Checkbox for whether to apply the z-score filter
        key = st_key_prefix + '__do_z_score_filter'
        if key not in st.session_state:
            st.session_state[key] = True
        st.checkbox('Do z-score filter', key=key)

        # Number input for the z-score filter threshold
        key = st_key_prefix + '__z_score_filter_threshold'
        if key not in st.session_state:
            st.session_state[key] = 3
        if st.session_state[st_key_prefix + '__do_z_score_filter']:
            st.number_input('z-score filter threshold:', min_value=0.0, key=key)

        # If dataset preprocessing is desired...
        if st.button('Preprocess dataset'):

            # Record the start time
            start_time = time.time()

            # Preprocess the dataset
            df = radial_profiles.preprocess_dataset(
                df,
                perc_thresh_rawintnorm_column_check=st.session_state[st_key_prefix + '__perc_thresh_rawintnorm_column_check'],
                image_col='Slide ID',
                nuclear_channel=st.session_state[st_key_prefix + '__nuclear_channel'],
                do_z_score_filter=st.session_state[st_key_prefix + '__do_z_score_filter'],
                z_score_filter_threshold=st.session_state[st_key_prefix + '__z_score_filter_threshold'],
                run_checks=st.session_state[st_key_prefix + '__run_checks']
            )

            # Output the time taken
            st.write(f'Preprocessing took {int(np.round(time.time() - start_time))} seconds')

            # Calculate the memory usage of the transformed dataframe
            st.session_state['input_dataframe_memory_usage_bytes'] = df.memory_usage(deep=True).sum()

            # Update the preprocessing parameters
            st.session_state['input_metadata']['preprocessing'] = {
                'location': 'Radial Profiles app',
                'nuclear_channel': st.session_state[st_key_prefix + '__nuclear_channel'],
                'do_z_score_filter': st.session_state[st_key_prefix + '__do_z_score_filter'],
            }
            if st.session_state[st_key_prefix + '__do_z_score_filter']:
                st.session_state['input_metadata']['preprocessing']['z_score_filter_threshold'] = st.session_state[st_key_prefix + '__z_score_filter_threshold']

            # Display information about the new dataframe
            df.info()

            # In case df has been modified not-in-place in any way, reassign the input dataset as the modified df
            st.session_state['input_dataset'].data = df


# Run the main function
if __name__ == '__main__':

    # Set page settings
    page_name = 'Radial Profiles'
    st.set_page_config(layout='wide', page_title=page_name)
    st.title(page_name)

    # Run streamlit-dataframe-editor library initialization tasks at the top of the page
    st.session_state = sde.initialize_session_state(st.session_state)

    # Run Top of Page (TOP) functions
    st.session_state = top.top_of_page_reqs(st.session_state)

    # Call the main function
    main()

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)
