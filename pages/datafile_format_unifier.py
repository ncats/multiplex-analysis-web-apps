# Import relevant libraries
import os
import streamlit as st
import pandas as pd
import app_top_of_page as top
import streamlit_dataframe_editor as sde
import re
import utils

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
    """
    Main function for the Datafile Unifier page.
    """

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

    # Write a message to the user
    st.write('After completing Section :one:, a sample of your dataset will be displayed at bottom. Use the sample to help you complete the rest of the sections on this page.')

    # Split the page into three main columns
    main_columns = st.columns(3)

    # In the first column...
    with main_columns[0]:

        # ---- 1. Select datafiles to combine --------------------------------------------------------------------------------------------------------------------------------

        # Display a header for the datafile selection section
        st.header(':one: Select datafiles to combine')

        # Get a list of files with the given extensions in the requested directory
        files = list_files(directory, extensions)

        # If no files are found, write a message to the user
        if len(files) == 0:
            st.warning('No ".csv" or ".tsv" files found in the `input` directory.')

        # If files are found, display them in a dataframe editor
        else:
            
            # Write messages to the user
            num_files = len(files)
            if num_files == 1:
                st.write('Detected 1 ".csv" or ".tsv" file in the `input` directory.')
            else:
                st.write('Detected {} ".csv" and ".tsv" files in the `input` directory.'.format(num_files))
            st.write('Note: You can double-click on any cell to see a full filename that is too long to fit in the cell.')
    
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

            # Preprocess the selected file selections
            selected_rows = df_reconstructed['Selected'] == True
            num_selected_rows = selected_rows.sum()
            if num_selected_rows == 1:
                st.write('There is 1 file selected.')
                button_text = ':star2: Combine selected files into single dataframe ⚠️ :star2:'
                button_help_message = 'Note: Even if just one file is selected, please click this button to proceed.'
            else:
                st.write('There are {} files selected.'.format(num_selected_rows))
                button_text = ':star2: Combine selected files into single dataframe :star2:'
                button_help_message = None
            with st.expander('Selected files:', expanded=False):
                st.dataframe(df_reconstructed.loc[selected_rows, 'Filename'], hide_index=True)

            # Extract the selected files from the dataframe editor
            input_files = sorted(df_reconstructed[selected_rows]['Filename'].to_list())
            st.session_state['unifier__input_files'] = input_files
            
            # Create a button to concatenate the selected files
            if st.button(button_text, help=button_help_message):

                # Render a progress spinner while the files are being combined
                with st.spinner('Combining files...'):

                    # Efficiently check if the columns are equal for all input files
                    columns_equal = True
                    if len(input_files) > 1:
                        sep = (',' if input_files[0].split('.')[-1] == 'csv' else '\t')
                        first_file_columns = pd.read_csv(os.path.join(directory, input_files[0]), nrows=0, sep=sep).columns
                        for input_file in input_files[1:]:
                            sep = (',' if input_file.split('.')[-1] == 'csv' else '\t')
                            current_file_columns = pd.read_csv(os.path.join(directory, input_file), nrows=0, sep=sep).columns
                            if not first_file_columns.equals(current_file_columns):
                                st.error('Columns are not equal for files: {} and {}'.format(input_files[0], input_file))
                                columns_equal = False
                                break

                    # If the columns are equal for all input files, concatenate all files into a single dataframe
                    if columns_equal:
                        sep = (',' if input_files[0].split('.')[-1] == 'csv' else '\t')
                        st.session_state['unifier__df'] = pd.concat([pd.read_csv(os.path.join(directory, input_file), sep=sep) for input_file in input_files], ignore_index=True)

                        # Save the setting used for this operation
                        st.session_state['unifier__input_files_actual'] = input_files

                # Display a success message
                if columns_equal:
                    st.toast('Files combined successfully')

            # If the selected input files have changed since the last time the combined dataframe was created, display a warning
            if ('unifier__input_files_actual' in st.session_state) and (st.session_state['unifier__input_files_actual'] != st.session_state['unifier__input_files']):
                st.warning('The selected input files have changed since the last time the combined dataframe was created. Please re-combine the files or adjust the input files to match the previous selection.')

    # If concatentation has been performed...
    if 'unifier__df' in st.session_state:

        # Get a shortcut to the concatenated dataframe
        df = st.session_state['unifier__df']

        # In the first column...
        with main_columns[0]:

            # ---- 2. (Optional) Delete null rows --------------------------------------------------------------------------------------------------------------------------------

            # Display a header for the null row deletion section
            st.header('(Optional) :two: Delete null rows')

            # Create an expander for the null row deletion section
            with st.expander('(Optional) Click to expand:', expanded=False):

                # Write a note to the user
                st.write('Observe the combined dataframe at bottom and select columns here by which to remove rows. Rows will be removed if the selected columns have a value of `None`.')

                # Allow the user to select the columns by which to delete rows
                if 'unifier__columns_to_drop_rows_by' not in st.session_state:
                    st.session_state['unifier__columns_to_drop_rows_by'] = []
                st.multiselect('Select columns by which to to delete rows:', df.columns, key='unifier__columns_to_drop_rows_by')

                # Create a button to delete rows from the dataframe
                if st.button(':star2: Delete rows from dataframe :star2:'):

                    # Render a progress spinner while the rows are being deleted
                    with st.spinner('Deleting rows...'):

                        row_count_before = len(df)
                        st.session_state['unifier__df'] = df.dropna(subset=st.session_state['unifier__columns_to_drop_rows_by']).reset_index(drop=True).convert_dtypes()
                        row_count_after = len(st.session_state['unifier__df'])

                        # Save the setting used for this operation
                        st.session_state['unifier__columns_to_drop_rows_by_actual'] = st.session_state['unifier__columns_to_drop_rows_by']

                    # Display a success message
                    st.toast(f'{row_count_after - row_count_before} rows deleted successfully')

                # If the selected columns to delete rows have changed since the last time rows were deleted, display a warning
                if ('unifier__columns_to_drop_rows_by_actual' in st.session_state) and (set(st.session_state['unifier__columns_to_drop_rows_by_actual']) != set(st.session_state['unifier__columns_to_drop_rows_by'])):
                    st.warning('The columns used to remove rows have changed since the last time rows were removed. Please start from scratch (the dataframe has been overwritten) or adjust the columns to match the previous selection.')

                # Get a shortcut to the concatenated dataframe
                df = st.session_state['unifier__df']

        # In the second column...
        with main_columns[1]:

            # ---- 3. Identify images --------------------------------------------------------------------------------------------------------------------------------

            # Display a header for the image identification section
            st.header(':three: Identify images')

            # Allow the user to select the columns by which to uniquely define images
            if 'unifier__columns_to_combine_to_uniquely_define_slides' not in st.session_state:
                st.session_state['unifier__columns_to_combine_to_uniquely_define_slides'] = []
            st.multiselect('Select columns to combine to uniquely define images:', df.columns, key='unifier__columns_to_combine_to_uniquely_define_slides')  # removing .select_dtypes(include=['object', 'string']) from df

            # Create a button to assign images to the dataframe
            if st.button(':star2: Assign images :star2:'):

                # Render a progress spinner while the images are being assigned
                with st.spinner('Assigning images...'):

                    # Perform the operation
                    subset_columns = st.session_state['unifier__columns_to_combine_to_uniquely_define_slides']
                    unique_rows = df.drop_duplicates(subset=subset_columns)[subset_columns]
                    df_from = unique_rows.apply(lambda x: '__'.join(x.apply(str)), axis='columns')
                    df_to = unique_rows.apply(lambda x: '__'.join(x.apply(str).apply(lambda y: re.split(r'[/\\]', y)[-1])).replace(' ', '_').replace('.', '_'), axis='columns')
                    transformation = dict(zip(df_from, df_to))
                    df_subset = df[subset_columns]
                    for column in subset_columns:
                        if df_subset[column].dtype not in ['string', 'object']:
                            df_subset[column] = df_subset[column].apply(str)
                    utils.dataframe_insert_possibly_existing_column(df, 0, 'Image ID (standardized)', df_subset.apply(lambda x: transformation['__'.join(x)], axis='columns'))

                    # Save this dataframe to memory
                    st.session_state['unifier__df'] = df

                    # Save the setting used for this operation
                    st.session_state['unifier__columns_to_combine_to_uniquely_define_slides_actual'] = subset_columns

                # Display a success message
                st.toast('Columns combined into a "Image ID (standardized)" column successfully')

            # If the selected columns to define images have changed since the last time images were defined, display a warning
            if ('unifier__columns_to_combine_to_uniquely_define_slides_actual' in st.session_state) and (set(st.session_state['unifier__columns_to_combine_to_uniquely_define_slides_actual']) != set(st.session_state['unifier__columns_to_combine_to_uniquely_define_slides'])):
                st.warning('The columns used to define slides have changed since the last time slides were defined. Please re-assign the image names or adjust the columns to match the previous selection.')

            # Get a shortcut to the concatenated dataframe
            df = st.session_state['unifier__df']

            # ---- 4. Identify regions of interest (ROIs) --------------------------------------------------------------------------------------------------------------------------------

            # Display a header for the region of interest (ROI) identification section
            st.header(':four: Identify regions of interest (ROIs)')

            # Allow the user to select the column by which to uniquely define ROIs
            if 'unifier__roi_explicitly_defined' not in st.session_state:
                st.session_state['unifier__roi_explicitly_defined'] = False
            if st.checkbox('Check here if the ROIs are explicitly defined by a column in the dataframe', key='unifier__roi_explicitly_defined', help='Note: This is not often the case, so this can usually be left unchecked.'):
                if 'unifier__roi_column' not in st.session_state:
                    st.session_state['unifier__roi_column'] = df.columns[0]
                st.selectbox('Select the column containing the ROI names:', df.columns, key='unifier__roi_column')

            # Assign some dynamic button-related variables
            if not st.session_state['unifier__roi_explicitly_defined']:
                button_text = ':star2: Assign ROIs ⚠️ :star2:'
                button_help_message = 'Note: Even if the ROIs are not explicitly defined by a column in the dataframe, please click this button to proceed.'
            else:
                button_text = ':star2: Assign ROIs :star2:'
                button_help_message = None

            # Create a button to assign ROIs to the dataframe
            if st.button(button_text, help=button_help_message):

                # Render a progress spinner while the ROIs are being assigned
                with st.spinner('Assigning ROIs...'):

                    # Perform the operation
                    if st.session_state['unifier__roi_explicitly_defined']:
                        utils.dataframe_insert_possibly_existing_column(df, 1, 'ROI ID (standardized)', df[st.session_state['unifier__roi_column']])
                        toast_message = '"ROI ID (standardized)" column created successfully'
                        num_roi_columns = 1
                    else:
                        if 'ROI ID (standardized)' in df.columns:
                            del df['ROI ID (standardized)']
                        toast_message = 'We have successfully ensured that the "ROI ID (standardized)" column is not present in the dataset'
                        num_roi_columns = 0

                    # Save this dataframe to memory
                    st.session_state['unifier__df'] = df

                    # Save the settings used for this operation
                    st.session_state['unifier__roi_explicitly_defined_actual'] = st.session_state['unifier__roi_explicitly_defined']
                    if st.session_state['unifier__roi_explicitly_defined']:
                        st.session_state['unifier__roi_column_actual'] = st.session_state['unifier__roi_column']
                    st.session_state['unifier__num_roi_columns_actual'] = num_roi_columns

                # Display a success message
                st.toast(toast_message)

            # If the selected columns to define ROIs have changed since the last time ROIs were defined, display a warning
            if st.session_state['unifier__roi_explicitly_defined']:
                if ('unifier__roi_explicitly_defined_actual' in st.session_state) and ((st.session_state['unifier__roi_explicitly_defined_actual'] != st.session_state['unifier__roi_explicitly_defined']) or (st.session_state['unifier__roi_column_actual'] != st.session_state['unifier__roi_column'])):
                    st.warning('The use of or column used to define ROIs has changed since the last time ROIs were defined. Please re-assign the ROI names or adjust the settings to match the previous ones.')
            else:
                if ('unifier__roi_explicitly_defined_actual' in st.session_state) and (st.session_state['unifier__roi_explicitly_defined_actual'] != st.session_state['unifier__roi_explicitly_defined']):
                    st.warning('The use of or column used to define ROIs has changed since the last time ROIs were defined. Please re-assign the ROI names or adjust the settings to match the previous ones.')

            # Get a shortcut to the concatenated dataframe
            df = st.session_state['unifier__df']

            # ---- 5. Identify coordinates --------------------------------------------------------------------------------------------------------------------------------

            # Display a header for the coordinate identification section
            st.header(':five: Identify coordinates')

            # Define some constants/variables
            df_numeric_columns = df.select_dtypes(include='number').columns
            coordinate_options = ['One column (centroid)', 'Two columns (min and max)']

            # Allow the user to select coordinate related settings
            smallest_columns = st.columns(2)
            with smallest_columns[0]:
                if 'unifier__number_of_coordinate_columns' not in st.session_state:
                    st.session_state['unifier__number_of_coordinate_columns'] = coordinate_options[0]
                st.radio('Select number of columns that specify one coordinate axis:', coordinate_options, key='unifier__number_of_coordinate_columns')
            with smallest_columns[1]:
                if 'unifier__microns_per_coordinate_unit' not in st.session_state:
                    st.session_state['unifier__microns_per_coordinate_unit'] = 1.0
                st.number_input('Enter the number of microns per coordinate unit in the columns below:', key='unifier__microns_per_coordinate_unit')
            if st.session_state['unifier__number_of_coordinate_columns'] == coordinate_options[0]:
                if 'unifier__x_coordinate_column' not in st.session_state:
                    st.session_state['unifier__x_coordinate_column'] = df_numeric_columns[0]
                if 'unifier__y_coordinate_column' not in st.session_state:
                    st.session_state['unifier__y_coordinate_column'] = df_numeric_columns[0]
                with smallest_columns[0]:
                    st.selectbox('Select a column for the x-coordinate:', df_numeric_columns, key='unifier__x_coordinate_column')
                    st.selectbox('Select a column for the y-coordinate:', df_numeric_columns, key='unifier__y_coordinate_column')
            else:
                if 'unifier__x_min_coordinate_column' not in st.session_state:
                    st.session_state['unifier__x_min_coordinate_column'] = df_numeric_columns[0]
                if 'unifier__y_min_coordinate_column' not in st.session_state:
                    st.session_state['unifier__y_min_coordinate_column'] = df_numeric_columns[0]
                if 'unifier__x_max_coordinate_column' not in st.session_state:
                    st.session_state['unifier__x_max_coordinate_column'] = df_numeric_columns[0]
                if 'unifier__y_max_coordinate_column' not in st.session_state:
                    st.session_state['unifier__y_max_coordinate_column'] = df_numeric_columns[0]
                with smallest_columns[0]:
                    st.selectbox('Select a column for the minimum x-coordinate:', df_numeric_columns, key='unifier__x_min_coordinate_column')
                with smallest_columns[1]:
                    st.selectbox('Select a column for the maximum x-coordinate:', df_numeric_columns, key='unifier__x_max_coordinate_column')
                with smallest_columns[0]:
                    st.selectbox('Select a column for the minimum y-coordinate:', df_numeric_columns, key='unifier__y_min_coordinate_column')
                with smallest_columns[1]:
                    st.selectbox('Select a column for the maximum y-coordinate:', df_numeric_columns, key='unifier__y_max_coordinate_column')

            # Create a button to assign coordinates to the dataframe, including converting them to microns and rounding them to the nearest 0.2 microns
            if st.button(':star2: Assign coordinates :star2:'):

                # Render a progress spinner while the coordinates are being assigned
                with st.spinner('Assigning coordinates...'):

                    # Perform the operation
                    if st.session_state['unifier__number_of_coordinate_columns'] == coordinate_options[0]:
                        centroid_x = df[st.session_state['unifier__x_coordinate_column']]
                        centroid_y = df[st.session_state['unifier__y_coordinate_column']]
                    else:
                        centroid_x = (df[st.session_state['unifier__x_min_coordinate_column']] + df[st.session_state['unifier__x_max_coordinate_column']]) / 2
                        centroid_y = (df[st.session_state['unifier__y_min_coordinate_column']] + df[st.session_state['unifier__y_max_coordinate_column']]) / 2
                    centroid_x = centroid_x * st.session_state['unifier__microns_per_coordinate_unit']
                    centroid_y = centroid_y * st.session_state['unifier__microns_per_coordinate_unit']
                    utils.dataframe_insert_possibly_existing_column(df, st.session_state['unifier__num_roi_columns_actual'] + 1, 'Centroid X (µm) (standardized)', (centroid_x / 0.2).round() * 0.2)
                    utils.dataframe_insert_possibly_existing_column(df, st.session_state['unifier__num_roi_columns_actual'] + 2, 'Centroid Y (µm) (standardized)', (centroid_y / 0.2).round() * 0.2)

                    # Save this dataframe to memory
                    st.session_state['unifier__df'] = df

                    # Save the settings used for this operation
                    st.session_state['unifier__number_of_coordinate_columns_actual'] = st.session_state['unifier__number_of_coordinate_columns']
                    st.session_state['unifier__microns_per_coordinate_unit_actual'] = st.session_state['unifier__microns_per_coordinate_unit']
                    if st.session_state['unifier__number_of_coordinate_columns'] == coordinate_options[0]:
                        st.session_state['unifier__x_coordinate_column_actual'] = st.session_state['unifier__x_coordinate_column']
                        st.session_state['unifier__y_coordinate_column_actual'] = st.session_state['unifier__y_coordinate_column']
                    else:
                        st.session_state['unifier__x_min_coordinate_column_actual'] = st.session_state['unifier__x_min_coordinate_column']
                        st.session_state['unifier__y_min_coordinate_column_actual'] = st.session_state['unifier__y_min_coordinate_column']
                        st.session_state['unifier__x_max_coordinate_column_actual'] = st.session_state['unifier__x_max_coordinate_column']
                        st.session_state['unifier__y_max_coordinate_column_actual'] = st.session_state['unifier__y_max_coordinate_column']

                # Display a success message
                st.toast('Coordinate columns created successfully')

            # If the selected columns to define coordinates have changed since the last time coordinates were defined, display a warning
            if st.session_state['unifier__number_of_coordinate_columns'] == coordinate_options[0]:
                if ('unifier__number_of_coordinate_columns_actual' in st.session_state) and any(st.session_state[key] != st.session_state[key + '_actual'] for key in ['unifier__number_of_coordinate_columns', 'unifier__microns_per_coordinate_unit', 'unifier__x_coordinate_column', 'unifier__y_coordinate_column']):
                    st.warning('The values of some coordinate settings have changed since the last time coordinates were assigned. Please re-assign the coordinates or adjust the settings to match the previous ones.')
            else:
                if ('unifier__number_of_coordinate_columns_actual' in st.session_state) and any(st.session_state[key] != st.session_state[key + '_actual'] for key in ['unifier__number_of_coordinate_columns', 'unifier__microns_per_coordinate_unit', 'unifier__x_min_coordinate_column', 'unifier__y_min_coordinate_column', 'unifier__x_max_coordinate_column', 'unifier__y_max_coordinate_column']):
                    st.warning('The values of some coordinate settings have changed since the last time coordinates were assigned. Please re-assign the coordinates or adjust the settings to match the previous ones.')

            # Get a shortcut to the concatenated dataframe
            df = st.session_state['unifier__df']

            # Allow the user to select the which columns correspond to the markers/phenotypes
            # st.header(':four: Select marker/phenotype columns')
            # st.multiselect('Optional: Select the categorical columns that correspond to the markers/phenotypes (i.e., thresholded intensities):', df_columns, key='unifier__marker_columns')
            # if st.button(':star2: Set ROI identifier :star2:'):
            #     st.session_state['unifier__roi_identifier'] = roi_identifier
            #     st.toast('ROI identifier set successfully')
            # if 'unifier__roi_identifier' in st.session_state:
            #     st.write('ROI identifier column: {}'.format(st.session_state['unifier__roi_identifier']))
            # cols_to_keep = ['Slide ID', 'tag', 'Cell X Position', 'Cell Y Position'] + df.loc[0, :].filter(regex='^Phenotype ').index.tolist()
            # cols_to_keep = ['Slide ID', 'tag', 'Cell X Position', 'Cell Y Position'] + df.loc[0, :].filter(regex='^Phenotype ').index.tolist() + extra_cols_to_keep
            # delete ROIs with single coordinates
            # optionally apply patching
            # do whatever is in HALO's pipeline
            # pick up going through dataset_formats.py line 218

        # Output a sample of the concatenated dataframe
        st.divider()
        st.header('Sample of unified dataframe')
        resample_dataframe = st.button('Refresh dataframe sample')
        if ('sampled_df' not in st.session_state) or resample_dataframe:
            sampled_df = df.sample(100).sort_index()
            st.session_state['sampled_df'] = sampled_df
        sampled_df = st.session_state['sampled_df']
        st.write(sampled_df)
        st.write('The full combined dataframe has {} rows and {} columns.'.format(df.shape[0], df.shape[1]))

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)

# Run the main function
if __name__ == '__main__':
    main()
