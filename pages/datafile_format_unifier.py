# Import relevant libraries
import os
import streamlit as st
import pandas as pd
import app_top_of_page as top
import streamlit_dataframe_editor as sde

def list_files(directory, extensions):
    """
    List files in a directory with given extension(s).

    Args:
        directory (str): The directory to search for files
        extensions (tuple): The file extensions for which to search

    Returns:
        list: A list of files in the directory with the given extension(s)
    """
    files = []
    for file in sorted(os.listdir(directory)):
        if file.endswith(extensions):
            files.append(file)
    return files

def main():

    # Set page settings
    st.set_page_config(layout='wide', page_title='Datafile Unifier')
    st.title('Datafile Unifier')

    # Run streamlit-dataframe-editor library initialization tasks at the top of the page
    st.session_state = sde.initialize_session_state(st.session_state)

    # Run Top of Page (TOP) functions
    st.session_state = top.top_of_page_reqs(st.session_state)

    # Constants
    directory = os.path.join('.', 'input')
    extensions = ('.csv', '.tsv')

    # Split the page into two columns
    main_columns = st.columns(2)

    # In the first column...
    with main_columns[0]:

        # Get a list of files with the given extensions in the requested directory
        files = list_files(directory, extensions)

        # Display a header for the datafile selection section
        st.header('Datafile selection')

        # If no files are found, write a message to the user
        if len(files) == 0:
            st.warning('No ".csv" or ".tsv" files found in the `input` directory.')

        # If files are found, display them in a dataframe editor
        else:
            
            # Write a message to the user
            num_files = len(files)
            if num_files == 1:
                st.write('Detected 1 ".csv" or ".tsv" file in the `input` directory:')
            else:
                st.write('Detected {} ".csv" and ".tsv" files in the `input` directory:'.format(num_files))

            # Create a dataframe from the list of files, including a column for the user to select files
            df_input_files = pd.DataFrame(files, columns=['Filename'])
            df_input_files.insert(0, 'Selected', False)

            # If the dataframe editor has been initialized, check to see if the detected files have changed from the last time the dataframe editor was initialized; if so, delete the dataframe editor so that it can be re-initialized
            if 'unifier__de_datafile_selection' in st.session_state:
                if set(st.session_state['unifier__de_datafile_selection'].reconstruct_edited_dataframe()['Filename']) != set(df_input_files['Filename']):
                    del st.session_state['unifier__de_datafile_selection']

            # If the dataframe editor has not been initialized, initialize it with the current file listing
            if 'unifier__de_datafile_selection' not in st.session_state:
                st.session_state['unifier__de_datafile_selection'] = sde.DataframeEditor(df_name='unifier__df_datafile_selection', default_df_contents=df_input_files)

            # Display the dataframe editor
            st.session_state['unifier__de_datafile_selection'].dataframe_editor(reset_data_editor_button_text='Reset file selections')

            # Extract the selected dataframe from the dataframe editor
            df_reconstructed = st.session_state['unifier__de_datafile_selection'].reconstruct_edited_dataframe()

            # Write a message to the user about the number of selected files
            selected_rows = df_reconstructed['Selected'] == True
            num_selected_rows = selected_rows.sum()
            if num_selected_rows == 1:
                st.write('There is 1 row selected.')
            else:
                st.write('There are {} rows selected.'.format(num_selected_rows))

            # Extract the selected files from the dataframe editor
            input_files = sorted(df_reconstructed[selected_rows]['Filename'].to_list())
            
            # Create a button to concatenate the selected files
            if st.button('Concatenate selected files'):

                # Efficiently check if the columns are equal for all input files
                columns_equal = True
                if len(input_files) > 1:
                    first_file_columns = pd.read_csv(os.path.join(directory, input_files[0]), nrows=0).columns
                    for input_file in input_files[1:]:
                        current_file_columns = pd.read_csv(os.path.join(directory, input_file), nrows=0).columns
                        if not first_file_columns.equals(current_file_columns):
                            st.error('Columns are not equal for files: {} and {}'.format(input_files[0], input_file))
                            columns_equal = False
                            break

                # If the columns are equal for all input files, concatenate all files into a single dataframe
                if columns_equal:
                    st.session_state['unifier__df'] = pd.concat([pd.read_csv(os.path.join(directory, input_file)) for input_file in input_files])
                    st.session_state['unifier__input_files'] = input_files

            # If the concatenated dataframe exists...
            if ('unifier__df' in st.session_state) and (st.session_state['unifier__input_files'] != input_files):
                st.warning('The input files have changed since the last time the dataframe was created. Please re-concatenate the files.')

    if 'unifier__df' in st.session_state:

        df = st.session_state['unifier__df']

        # In the second column...
        with main_columns[1]:
            st.write('second column')

        st.divider()
        st.header('Sample of concatenated dataframe')
        st.write('Click on the page and hit "r" to refresh the sample.')
        st.write(df.sample(100))
        st.write('The full concatenated dataframe has {} rows and {} columns.'.format(df.shape[0], df.shape[1]))

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)

# Run the main function
if __name__ == '__main__':
    main()
