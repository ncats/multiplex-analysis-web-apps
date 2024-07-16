# Import relevant libraries
import streamlit as st
import radial_profiles
import time
import numpy as np

# Global variable
st_key_prefix = 'preprocessing__'


# Function to initialize the preprocessing section
def initialize_radial_profiles_preprocessing(df):

    # Checkbox for whether to run checks
    key = st_key_prefix + 'run_checks'
    if key not in st.session_state:
        st.session_state[key] = False
    run_checks = st.checkbox('Run checks', key=key)

    # Number input for the threshold for the RawIntNorm check
    key = st_key_prefix + 'perc_thresh_rawintnorm_column_check'
    if key not in st.session_state:
        st.session_state[key] = 0.01
    if run_checks:
        st.number_input('Threshold for the RawIntNorm column check (%):', min_value=0.0, max_value=100.0, key=key)
    perc_thresh_rawintnorm_column_check = st.session_state[key]

    # Number input to select the nuclear intensity channel
    key = st_key_prefix + 'nuclear_channel'
    if key not in st.session_state:
        st.session_state[key] = 1
    nuclear_channel = st.number_input('Nuclear channel:', min_value=1, key=key)

    # Checkbox for whether to apply the z-score filter
    key = st_key_prefix + 'do_z_score_filter'
    if key not in st.session_state:
        st.session_state[key] = True
    do_z_score_filter = st.checkbox('Do z-score filter', key=key)

    # Number input for the z-score filter threshold
    key = st_key_prefix + 'z_score_filter_threshold'
    if key not in st.session_state:
        st.session_state[key] = 3
    if do_z_score_filter:
        st.number_input('z-score filter threshold:', min_value=0.0, key=key)
    z_score_filter_threshold = st.session_state[key]

    # If dataset preprocessing is desired...
    if st.button('Preprocess dataset'):

        # Record the start time
        start_time = time.time()

        # Preprocess the dataset
        df = radial_profiles.preprocess_dataset(
            df,
            perc_thresh_rawintnorm_column_check=perc_thresh_rawintnorm_column_check,
            image_col='Slide ID',
            nuclear_channel=nuclear_channel,
            do_z_score_filter=do_z_score_filter,
            z_score_filter_threshold=z_score_filter_threshold,
            run_checks=run_checks
        )

        # If the dataset is None, it's likely preprocessing has already been performed, so display a warning and return
        if df is None:
            st.warning('It appears that the dataset has already been preprocessed because there is no "Label" column. If you would like to re-preprocess the dataset, please reload it from the Open File page.')
            return

        # Output the time taken
        st.write(f'Preprocessing took {int(np.round(time.time() - start_time))} seconds')

        # Calculate the memory usage of the transformed dataframe
        st.session_state['input_dataframe_memory_usage_bytes'] = df.memory_usage(deep=True).sum()

        # Update the preprocessing parameters
        st.session_state['input_metadata']['preprocessing'] = {
            'location': 'Radial Profiles app',
            'nuclear_channel': nuclear_channel,
            'do_z_score_filter': do_z_score_filter,
        }
        if do_z_score_filter:
            st.session_state['input_metadata']['preprocessing']['z_score_filter_threshold'] = z_score_filter_threshold

        # Display information about the new dataframe
        df.info()

        # In case df has been modified not-in-place in any way, reassign the input dataset as the modified df
        st.session_state['input_dataset'].data = df

    # Return the modified dataframe
    return df


def main():
    """
    Main function for the page.
    """

    st.header('ImageJ Output for Radial Profiles Analysis')

    # Ensure a dataset has been opened in the first place
    if 'input_dataset' not in st.session_state:
        st.warning('Please open a dataset from the Open File page at left.')
        return
    
    # Save a shortcut to the dataframe
    df = st.session_state['input_dataset'].data

    # Set up preprocessing
    with st.columns(3)[0]:
        df = initialize_radial_profiles_preprocessing(df)

    # Ensure the main dataframe is updated per the operations above
    st.session_state['input_dataset'].data = df


# Run the main function
if __name__ == '__main__':

    # Call the main function
    main()
