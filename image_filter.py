# Import relevant libraries
import streamlit as st
import pandas as pd
import numpy as np

# Global variable
st_key_prefix = 'test_filter__'


# This is an image filter that should behave somewhat like a Streamlit (macro) widget
def image_filter(df, selected_cols_for_filtering, st_key_prefix, key, color='red'):

    with st.expander(f'Image filter for :{color}[{key}] group:', expanded=False):

        # Build a widget for all selected filtering columns
        for col in selected_cols_for_filtering:
            build_multiselect(df, col, key)

        # Create a mask for each filter
        masks = [df[col].isin(st.session_state[st_key_prefix + 'filtering_multiselect_' + key + '_' + col]) for col in selected_cols_for_filtering if st.session_state[st_key_prefix + 'filtering_multiselect_' + key + '_' + col]]

        # Combine the masks
        combined_mask = np.logical_and.reduce(masks)

        # Apply the combined mask to the DataFrame
        if masks:
            df_masked = df[combined_mask]
        else:
            df_masked = df.copy()

        # Output the number of images that passed the filter
        st.write(f'Number of images filtered above: {len(df_masked)}')

        # Create an interactive dataframe to allow the user to customize the image selection
        df_selection = st.dataframe(df_masked, on_select='rerun', hide_index=True, key=st_key_prefix + key + '_df_selection')

        # Output the number of images that have been manually selected by the user
        st.write(f'Number of images selected above: {len(df_selection["selection"]["rows"])}')

        # Output the filenames of the selected images
        ser_selection = df_masked['input_filename'].iloc[df_selection['selection']['rows']]
        st.dataframe(ser_selection, hide_index=True)

        # Convert the list of selected images to a list
        selected_images = ser_selection.tolist()

        # Save it to the session state
        full_key = st_key_prefix + key + '_selected_images'
        st.session_state[full_key] = selected_images

    # Also return it
    return(selected_images)


# Reset the filtering columns
def reset_filtering_columns():
    for key in st.session_state.keys():
        if key.startswith(st_key_prefix + 'filtering_multiselect_'):
            st.session_state[key] = []
    if st_key_prefix + 'df_deduped' in st.session_state:
        del st.session_state[st_key_prefix + 'df_deduped']


# Build a multiselect widget for a given column
def build_multiselect(df, col, widget_key_prefix):
    unique_vals = df[col].unique()
    st.multiselect(f'Filter image on `{col}`:', unique_vals, key=st_key_prefix + 'filtering_multiselect_' + widget_key_prefix + '_' + col)


# Main function
def main():

    # Load the full dataframe from disk
    # Sample of how it's written to disk in the first place from preprocess_radial_profile_data.ipynb:
    #   import random
    #   input_filenames = df_transformed['input_filename'].unique()
    #   df_transformed[df_transformed['input_filename'].isin(random.sample(list(input_filenames), 10))].to_hdf(os.path.join('image_data.h5'), key='df_transformed_partial', mode='w', format='table', complevel=9)
    if st.button('Load data from disk'):
        st.session_state[st_key_prefix + 'df'] = pd.read_hdf('image_data.h5')
        st.info(f'Data of shape {st.session_state[st_key_prefix + "df"].shape} loaded successfully')

    # Ensure the full dataset has been loaded from disk
    if st_key_prefix + 'df' not in st.session_state:
        st.warning('Please load the data first')
        return
    
    # Get a shortcut to the full dataframe
    df = st.session_state[st_key_prefix + 'df']

    # Get the columns on which we will allow the user to create a filter
    df_cols = [col for col in df.columns.tolist() if col != 'input_filename']

    # Allow the user to select the columns on which they want to filter
    selected_cols_for_filtering = st.multiselect('Select columns on which to filter:', df_cols, key=st_key_prefix + 'selected_cols_for_filtering', on_change=reset_filtering_columns)

    # Simplify the dataframe to presumably just the essentially categorical columns
    if st.button('Prepare filtering data'):
        st.session_state[st_key_prefix + 'df_deduped'] = df[['input_filename'] + selected_cols_for_filtering].drop_duplicates().sort_values(selected_cols_for_filtering)

    # Ensure the deduplication based on the selected columns has been performed
    if st_key_prefix + 'df_deduped' not in st.session_state:
        st.warning('Please prepare the filtering data first')
        return
    
    # Get a shortcut to the deduplicated dataframe
    df_deduped = st.session_state[st_key_prefix + 'df_deduped']

    # Create two image filters
    selected_images_baseline = image_filter(df_deduped, selected_cols_for_filtering, st_key_prefix, key='baseline', color='blue')
    selected_images_signal = image_filter(df_deduped, selected_cols_for_filtering, st_key_prefix, key='signal', color='red')

    # Output the selected images in each group
    st.write('Selected images in the baseline group:')
    st.write(selected_images_baseline)
    st.write('Selected images in the signal group:')
    st.write(selected_images_signal)  # or could write, e.g., st.session_state[st_key_prefix + 'signal' + '_selected_images']


# Run the main function
if __name__ == '__main__':
    main()
