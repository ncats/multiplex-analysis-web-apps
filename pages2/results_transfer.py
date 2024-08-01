# Import relevant libraries
import streamlit as st
import os
import pandas as pd
from datetime import datetime
import pytz
import zipfile
import nidap_io

# Global variable
st_key_prefix = 'results_transfer__'


def zip_files(file_paths, output_zip):
    """
    Zips the given list of file paths into a single zip file.

    :param file_paths: List of file paths to include in the zip file.
    :param output_zip: The path to the output zip file.
    """
    with zipfile.ZipFile(output_zip, 'w') as zipf:
        for file in file_paths:
            # Add file to the zip file
            zipf.write(file, os.path.basename(file))


def main():
    """
    Main function for the page.
    """

    # Output directory
    output_dir = 'output'

    # Get the list of files in the directory
    dir_listing = os.listdir(output_dir)

    # Create a list to hold file information
    file_info = []

    # Define the EDT timezone
    edt = pytz.timezone('US/Eastern')

    # For each file in the directory...
    for file_name in dir_listing:

        # Get the full file path
        file_path = os.path.join(output_dir, file_name)

        # Get the last modification time
        mod_time = os.path.getmtime(file_path)

        # Convert the timestamp to a datetime object
        mod_time_dt = datetime.fromtimestamp(mod_time)

        # Localize the datetime object to UTC and then convert to EDT
        mod_time_edt = mod_time_dt.astimezone(edt)

        # Format the datetime object to a human-readable string
        mod_time_readable = mod_time_edt.strftime('%Y-%m-%d %H:%M:%S %Z')

        # Append the file information to the list
        file_info.append({'File Name': file_name, 'Last Modified': mod_time_readable})

    # Create a DataFrame from the file information
    st.write('Files in the output directory:')
    df = pd.DataFrame(file_info)

    # Display the DataFrame in Streamlit
    selected_files = st.dataframe(df, hide_index=True, on_select='rerun')

    # On the first of three columns...
    with st.columns(3)[0]:

        # Add a text input for the zip file name
        key = st_key_prefix + 'zipfile_name'
        if key not in st.session_state:
            st.session_state[key] = 'output.zip'
        st.text_input('Zip file name:', key=key)
        zipfile_path = os.path.join(output_dir, st.session_state[key])

        # Check if we're on NIDAP
        on_nidap = st.session_state['platform'] == 'nidap'

        # Add a button to zip (and transfer, if on NIDAP) the selected files
        if st.button('Zip and transfer to NIDAP' if on_nidap else 'Zip'):

            # Ensure selection exists
            if 'selection' in selected_files and 'rows' in selected_files['selection']:

                # Get the selected files
                selected_files_list = df.iloc[selected_files['selection']['rows']]['File Name'].tolist()

                # Check if any files were selected
                if selected_files_list:

                    # Zip the selected files
                    zip_files([os.path.join(output_dir, file) for file in selected_files_list], zipfile_path)

                    # Display a success message
                    st.success('Files zipped')

                    # If we're on NIDAP...
                    if on_nidap:

                        # Get the output dataset
                        dataset = nidap_io.get_foundry_dataset(alias='output')

                        # Upload the zip file to the dataset
                        nidap_io.upload_file_to_dataset(dataset, zipfile_path)

                        # Display a success message
                        st.success('Zip file transferred to NIDAP')

            # Otherwise, display an error message
                else:
                    st.warning('No files selected to zip!')
            else:
                st.warning('No files selected to zip!')


# Run the main function
if __name__ == '__main__':

    # Call the main function
    main()
