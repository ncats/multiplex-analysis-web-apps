# Import relevant libraries
import os
import streamlit as st
import pandas as pd
import app_top_of_page as top
import streamlit_dataframe_editor as sde
import re
# from streamlit_dataframe_editor import DataframeEditor  # attempting Lens' fix from https://stackoverflow.com/questions/1412787/picklingerror-cant-pickle-class-decimal-decimal-its-not-the-same-object for the error PicklingError: Can't pickle <class 'streamlit_dataframe_editor.DataframeEditor'>: it's not the same object as streamlit_dataframe_editor.DataframeEditor --> does not work, same error

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
    main_columns = st.columns(3)

    # In the first column...
    with main_columns[0]:

        # Get a list of files with the given extensions in the requested directory
        files = list_files(directory, extensions)

        # Display a header for the datafile selection section
        st.header(':one: Select datafiles')

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
            if st.button(':star2: Concatenate selected files :star2:'):

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
                    st.session_state['unifier__df'] = pd.concat([pd.read_csv(os.path.join(directory, input_file)) for input_file in input_files], ignore_index=True)
                    st.session_state['unifier__input_files'] = input_files
                    # Delete some subsequent keys if present
                    for key_to_delete in ['unifier__columns_actually_used_to_drop_rows', 'unifier__columns_actually_used_to_define_slides']:
                        if key_to_delete in st.session_state:
                            del st.session_state[key_to_delete]
                    st.toast('Files concatenated successfully')

            # If the concatenated dataframe exists...
            if ('unifier__input_files' in st.session_state) and (st.session_state['unifier__input_files'] != input_files):
                st.warning('The input files have changed since the last time the dataframe was created. Please re-concatenate the files or adjust the input files to match the previous selection.')

    # If concatentation has been performed...
    if 'unifier__df' in st.session_state:

        # Get a shortcut to the concatenated dataframe
        df = st.session_state['unifier__df']

        # In the second column...
        with main_columns[1]:

            # Drop rows with `None` values in selected columns
            st.header(':two: (Optional) Drop null rows')
            st.write('Observe the concatenated dataframe at bottom and select columns by which to drop rows. Rows will be dropped if the selected columns have a value of `None`.')
            st.multiselect('Select columns by which to to drop rows:', df.columns, key='unifier__columns_to_drop_rows_by')
            if st.button(':star2: Drop rows :star2:'):
                st.session_state['unifier__df'] = df.dropna(subset=st.session_state['unifier__columns_to_drop_rows_by']).reset_index(drop=True).convert_dtypes()
                st.session_state['unifier__columns_actually_used_to_drop_rows'] = st.session_state['unifier__columns_to_drop_rows_by']
                # Delete some subsequent keys if present
                for key_to_delete in ['unifier__columns_actually_used_to_define_slides']:
                    if key_to_delete in st.session_state:
                        del st.session_state[key_to_delete]
                st.toast('Rows dropped successfully')
            if ('unifier__columns_actually_used_to_drop_rows' in st.session_state) and (set(st.session_state['unifier__columns_actually_used_to_drop_rows']) != set(st.session_state['unifier__columns_to_drop_rows_by'])):
                st.warning('The columns used to drop rows have changed since the last time rows were dropped. Please start from scratch (the dataframe has been overwritten) or adjust the columns to match the previous selection.')
            df = st.session_state['unifier__df']

            # Identify columns that combine to uniquely define slides
            st.header(':three: Combine columns to uniquely define slides')
            st.multiselect('Select string columns to combine to uniquely define slides:', df.select_dtypes(include=['object', 'string']).columns, key='unifier__columns_to_combine_to_uniquely_define_slides')
            if st.button(':star2: Combine columns :star2:'):
                subset_columns = st.session_state['unifier__columns_to_combine_to_uniquely_define_slides']
                unique_rows = df.drop_duplicates(subset=subset_columns)[subset_columns]
                df_from = unique_rows.apply(lambda x: '__'.join(x), axis='columns')
                df_to = unique_rows.apply(lambda x: '__'.join(x.apply(lambda y: re.split(r'[/\\]', y)[-1])).replace(' ', '_').replace('.', '_'), axis='columns')
                transformation = dict(zip(df_from, df_to))
                df.insert(0, 'Slide ID', df[subset_columns].apply(lambda x: transformation['__'.join(x)], axis='columns'))  # much faster than the commented line below
                # df.insert(0, 'Slide ID', df[subset_columns].apply(lambda x: '__'.join(x.apply(lambda y: re.split(r'[/\\]', y)[-1])).replace(' ', '_').replace('.', '_'), axis='columns'))  # takes much longer to do it directly as here
                st.session_state['unifier__df'] = df
                st.session_state['unifier__columns_actually_used_to_define_slides'] = subset_columns
                st.toast('Columns combined into a "Slide ID" column successfully')
            if ('unifier__columns_actually_used_to_define_slides' in st.session_state) and (set(st.session_state['unifier__columns_actually_used_to_define_slides']) != set(st.session_state['unifier__columns_to_combine_to_uniquely_define_slides'])):
                st.warning('The columns used to define slides have changed since the last time slides were defined. Please re-combine the columns or adjust them to match the previous selection.')
            df = st.session_state['unifier__df']

        with main_columns[2]:

            # Store the current columns in a variable
            df_columns = df.columns
            df_numeric_columns = df.select_dtypes(include='number').columns
            
            # Allow user to select a column as ROI identifier
            # st.header(':four: Select ROI identifier')
            st.selectbox('Select a column from to be used as a region of interest (ROI) identifier. If no particular identifier exists, use the "Slide ID" column:', df_columns, key='unifier__roi_identifier_column')

            # Allow user to select one or two columns for specifying coordinates
            # st.header(':four: Select coordinate columns')
            coordinate_options = ['One column (centroid)', 'Two columns (min and max)']
            st.radio('Select number of columns that specify one coordinate axis:', coordinate_options, key='unifier__number_of_coordinate_columns')
            if st.session_state['unifier__number_of_coordinate_columns'] == coordinate_options[0]:
                x_coordinate_column = st.selectbox('Select a column for the x-coordinate:', df_numeric_columns, key='unifier__x_coordinate_column')
                y_coordinate_column = st.selectbox('Select a column for the y-coordinate:', df_numeric_columns, key='unifier__y_coordinate_column')
                # TODO: Add code to handle one column selection
            else:
                x_min_coordinate_column = st.selectbox('Select a column for the minimum x-coordinate:', df_numeric_columns, key='unifier__x_min_coordinate_column')
                x_max_coordinate_column = st.selectbox('Select a column for the maximum x-coordinate:', df_numeric_columns, key='unifier__x_max_coordinate_column')
                y_min_coordinate_column = st.selectbox('Select a column for the minimum y-coordinate:', df_numeric_columns, key='unifier__y_min_coordinate_column')
                y_max_coordinate_column = st.selectbox('Select a column for the maximum y-coordinate:', df_numeric_columns, key='unifier__y_max_coordinate_column')
                # TODO: Add code to handle two columns selection

            # Allow user to select the coordinate units in microns
            # st.header(':four: Select coordinate units')
            st.number_input('Enter the number of microns per coordinate unit:', value=1.0, key='unifier__microns_per_coordinate_unit')


            # if st.button(':star2: Set ROI identifier :star2:'):
            #     st.session_state['unifier__roi_identifier'] = roi_identifier
            #     st.toast('ROI identifier set successfully')
            # if 'unifier__roi_identifier' in st.session_state:
            #     st.write('ROI identifier column: {}'.format(st.session_state['unifier__roi_identifier']))

        # Output a sample of the concatenated dataframe
        st.divider()
        st.header('Sample of concatenated dataframe')
        resample_dataframe = st.button('Resample dataframe')
        if ('sampled_df' not in st.session_state) or resample_dataframe:
            sampled_df = df.sample(100).sort_index()
            st.session_state['sampled_df'] = sampled_df
        sampled_df = st.session_state['sampled_df']
        st.write(sampled_df)
        st.write('The full concatenated dataframe has {} rows and {} columns.'.format(df.shape[0], df.shape[1]))

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)

# Run the main function
if __name__ == '__main__':
    main()
