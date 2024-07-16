# Import relevant libraries
import pandas as pd
import os
import utils
import time
import scipy.stats


# Define a function to determine the common suffix in a series
def common_suffix(ser):
    strings = ser.tolist()
    reversed_strings = [s[::-1] for s in strings]
    reversed_lcs = os.path.commonprefix(reversed_strings)
    suffix = reversed_lcs[::-1]
    print(f'The common suffix is "{suffix}".')
    return suffix


# Create a function to check for duplicate columns
def has_duplicate_columns(df, sample_size=1000):
    print('NOTE: Sampling is being performed, so if this returns True, consider increasing the sample_size parameter (or manually check for duplicates). If it returns False, you know that there are no duplicate columns in the DataFrame and do not need to stress about the sampling.')
    df = df.sample(n=sample_size)
    df_transposed = df.T
    df_deduplicated = df_transposed.drop_duplicates().T
    return len(df.columns) != len(df_deduplicated.columns)  # if they're the same length (has_duplicate_columns() returns False), you know there are no duplicate columns


# Define a function to appropriately transform the dataframe
def transform_dataframe(df1, index, columns, values, repeated_columns, verbose=False):

    # Sample call
    # df2.equals(transform_dataframe(df1, index='a', columns='c', values=['d', 'e'], repeated_columns=['b', 'f'], verbose=True))

    # Pivot the specified values columns based on the specified index and columns
    if verbose:
        print(f'Pivoting the DataFrame based on the index "{index}", the columns "{columns}", and the values "{values}"...')
        start_time = time.time()
    pivot_df = df1.pivot_table(index=index, columns=columns, values=values, aggfunc='first', observed=False)
    if verbose:
        print(f'Pivoting the DataFrame took {time.time() - start_time:.2f} seconds.')

    # Flatten the multi-index columns
    pivot_df.columns = [f'{x[0]}{x[1]}' for x in pivot_df.columns]

    # Reset the index to make the index column a column again
    pivot_df.reset_index(inplace=True)

    # Merge with the repeated_columns from the original DataFrame
    if verbose:
        print(f'Merging the DataFrame with the repeated columns "{repeated_columns}"...')
        start_time = time.time()
    df2 = pd.merge(pivot_df, df1[[index] + repeated_columns].drop_duplicates(), on=index)
    if verbose:
        print(f'Merging the DataFrame took {time.time() - start_time:.2f} seconds.')

    # Return the transformed DataFrame
    return df2


# Transform the dataframes image-by-image to save both memory and time
def transform_dataframes_in_chunks(df, image_col, new_index_col, distinct_old_row_identifier_col, cols_with_unique_rows, cols_with_repeated_rows, verbose=False):
    for unique_image in df[image_col].unique():
        df_image = df[df[image_col] == unique_image]
        df_image2 = transform_dataframe(df_image, index=new_index_col, columns=distinct_old_row_identifier_col, values=cols_with_unique_rows, repeated_columns=cols_with_repeated_rows, verbose=verbose)
        yield df_image2


# Efficiently turn a series of strings (in the format of the ImageJ "Labe" column) into three series of the relevant parts
def process_label_column(ser_label):

    # Regular expression to capture the three parts
    pattern = r'^(.*?) - (T=.*?:c:.*/.*? t:.*/.*?) - (.*)$'

    # Extract parts into separate columns
    extracted_df = ser_label.str.extract(pattern)

    # Return the desired parts
    return extracted_df[0], extracted_df[1], extracted_df[2]


def preprocess_dataset(df, perc_thresh_rawintnorm_column_check=0.01, image_col='Slide ID', nuclear_channel=2, do_z_score_filter=True, z_score_filter_threshold=3, run_checks=False):

    # Note the input dataframe df should likely be that loaded from the Open File page
    # This function is essentially the same as the preprocess_radial_profile_data.ipynb Jupyter notebook

    # Constant
    run_repeated_rows_check = False

    # Output the initial dataframe shape
    print('Initial dataframe shape:', df.shape)

    # Determine whether preprocessing has already been performed
    if 'Label' not in df.columns:
        print('It appears that the dataset has already been preprocessed because there is no "Label" column. If you would like to re-preprocess the dataset, please reload it from the Open File page.')
        return

    # Efficiently save some data from the "Label" column
    tif_name, middle_data, last_name = process_label_column(df['Label'])

    # Show that the basename of the image name is always the same as the last name
    if run_checks:
        assert (tif_name.apply(lambda x: os.path.splitext(x)[0]) == last_name).all(), 'The basename of the image name is not always the same as the last name.'

    # Add the TIF image names to the dataframe
    df['tif_name'] = tif_name.astype(str)
    print('tif_name:', df['tif_name'])

    # Show that there are always five fields of the middle data when split using ":"
    if run_checks:
        assert list(middle_data.apply(lambda x: len(x.split(':'))).unique()) == [5], 'There are not always five fields of the middle data when split using ":".'

    # Add the time ("T") field to the observations dataframe
    df['T'] = middle_data.apply(lambda x: x.split(':')[0].removeprefix('T=')).astype(int)
    print('T:', df['T'])

    # Add the cell ID field to the observations dataframe
    df['cell_id'] = middle_data.apply(lambda x: x.split(':')[1]).astype(str)
    print('cell_id:', df['cell_id'])

    # Check that the "c" field is the exact same as the "Ch" column
    if run_checks:
        assert (df['Ch'] == middle_data.apply(lambda x: x.split(':')[3].split('/')[0]).astype(int)).all(), 'The "c" field is not the exact same as the "Ch" column.'

    # Check that the last field is exactly the same as one plus the "T" field
    if run_checks:
        assert (df['T'] == middle_data.apply(lambda x: x.split(':')[4].split('/')[0]).astype(int) - 1).all(), 'The last field is not exactly the same as one plus the "T" field.'

    # Check that the .tif basename is completely contained in the actual input filename
    df_small = df[['input_filename', 'tif_name']].drop_duplicates()
    if run_checks:
        assert df_small.apply(lambda x: os.path.splitext(x['tif_name'])[0].replace('REEEC', 'REEC') in x['input_filename'], axis='columns').all(), 'The .tif basename is not completely contained in the actual input filename.'

    # Determine (and remove from the small df) the common suffix in the input_filename field
    ser = df_small['input_filename']
    suffix1 = common_suffix(ser)
    df_small['input_filename'] = ser.str.removesuffix(suffix1)
    print('df_small:', df_small)

    # Ensure that the "T=X" part of the input filename is the same as the "T" field
    if run_checks:
        assert df['input_filename'].apply(lambda x: x.removesuffix(suffix1).split('=')[-1]).astype(int).equals(df['T']), 'The "T=X" part of the input filename is not the same as the "T" field.'

    # Determine the next common suffix aside from the T value
    ser = df_small['input_filename'].apply(lambda x: x.removesuffix(suffix1).split('=')[0])
    suffix2 = common_suffix(ser)
    df_small['input_filename'] = ser.str.removesuffix(suffix2)
    print('df_small:', df_small)

    # Get a series of the data to process from the "input_filename" field
    ser_remaining_data = df['input_filename'].apply(lambda x: x.removesuffix(suffix1).split('=')[0].removesuffix(suffix2))
    print('ser_remaining_data:', ser_remaining_data)

    # Add a column to identify if the cell was processed using the REEC
    df['REEC'] = ser_remaining_data.str.endswith('_REEC').astype(bool)
    print('REEC:', df['REEC'])

    # Remove the just-processed suffix
    ser_remaining_data = ser_remaining_data.str.removesuffix('_REEC')
    print('ser_remaining_data:', ser_remaining_data)

    # Add columns identifying the well and cell type
    df['well_id'] = ser_remaining_data.apply(lambda x: x.split('_')[0]).astype(str)
    df['cell_type'] = ser_remaining_data.apply(lambda x: x.split('_')[1]).astype(str)
    print('well_id and cell_type:', df[['well_id', 'cell_type']])

    # Normalize the total intensity using the cell areas
    df['RawIntNorm'] = df['RawIntDen'] / df['Area']
    print('RawIntNorm:', df['RawIntNorm'])

    # Check that this yields the same result already included in the datafile
    if run_checks:
        min_true_val = df['Mean'].min()
        thresh = min_true_val * perc_thresh_rawintnorm_column_check / 100
        mad = (df['RawIntNorm'] - df['Mean']).abs().max()
        if mad < thresh:
            print(f'The new and existing columns ARE equal to within {perc_thresh_rawintnorm_column_check}% of the minimum existing value (MAD: {mad}).')
        else:
            print(f'The new and existing columns are NOT equal to within {perc_thresh_rawintnorm_column_check}% of the minimum existing value (MAD: {mad}).')

    # Set a row ID as a combination of the input filename and the extracted cell ID
    # Note then that for each row ID there should be the same number of rows as there are channels, unless some cell IDs are duplicated, which is what the following cells test.
    # As may often be necessary, we are doing this on a per-image basis to save memory.
    unique_images = df[image_col].unique()
    for image in unique_images:
        image_loc = df[image_col] == image
        df.loc[image_loc, 'row_id'] = df.loc[image_loc, 'input_filename'].astype(str) + ' - ' + df.loc[image_loc, 'cell_id']

    # Display the rows corresponding to duplicated cell IDs
    num_channels = df['Ch'].nunique()
    vc = df['row_id'].value_counts()
    df_dupes = df[df['row_id'].isin(vc[vc > num_channels].index)]
    print('df_dupes:', df_dupes)

    # Create a dataframe to aid in transforming the duplicated cell IDs to something unique
    df_assist_deduping = df_dupes.groupby(by=['row_id', 'Ch'], observed=False)[' '].aggregate(sorted).to_frame().reset_index()
    df_assist_deduping['index'] = df_assist_deduping[' '].apply(lambda x: list(range(len(x))))
    df_assist_deduping_expanded = df_assist_deduping.explode(' ')
    df_assist_deduping_expanded['index'] = df_assist_deduping.explode('index')['index'].apply(lambda x: f'{x:04d}')
    print('df_assist_deduping_expanded:', df_assist_deduping_expanded)

    # Check that the row ID in combination with the "blank index" uniquely identifies all rows (the duplicated cell IDs are not used here)
    if run_checks:
        assert (df['row_id'] + ' - ' + df[' '].astype(str)).nunique() == len(df), 'The row ID in combination with the "blank index" does not uniquely identify all rows.'

    # Append string indices to the cell IDs in order to de-duplicate them
    for _, ser in df_assist_deduping_expanded.iterrows():
        row_id = ser['row_id']
        blank_index = ser[' ']
        index_to_append = ser['index']
        loc = (df['row_id'] == row_id) & (df[' '] == blank_index)
        assert loc.sum() == 1
        df.loc[loc, 'cell_id'] = df.loc[loc, 'cell_id'] + ':' + index_to_append

    # Output the de-duplicated versions of the previously duplicated rows (just `cell_id` is modified so far)
    if run_checks:
        print('df at the locations of the duplicates:', df.loc[df_dupes.index])

    # Now confirm that there are no more duplicated cell IDs
    for image in df[image_col].unique():
        image_loc = df[image_col] == image
        df.loc[image_loc, 'row_id'] = df.loc[image_loc, 'input_filename'].astype(str) + ' - ' + df.loc[image_loc, 'cell_id']
    vc = df['row_id'].value_counts()
    if run_checks:
        df_dupes = df[df['row_id'].isin(vc[vc > num_channels].index)]
        print('df_dupes:', df_dupes)

    # Really confirm there are no duplicated rows anymore
    if run_checks:
        assert list(vc.unique()) == [num_channels], 'There are still duplicated rows.'

    # Check for duplicate columns (ignoring MAWA-created columns)
    if run_checks:
        print('has_duplicate_columns:', has_duplicate_columns(df.drop(columns=['input_filename'])))

    # Check whether suspected columns are equal
    if run_checks:
        assert (df['IntDen'] == df['RawIntDen']).all(), 'The "IntDen" and "RawIntDen" columns are not equal.'

    # Drop the ostensibly more processed one
    df.drop(columns=['IntDen'], inplace=True)

    # Check again for duplicate columns
    # This should return False if just loading the data in a Jupyter notebook but not when dataset_formats.py is used (as in Open File) because it probably does add duplicate columns
    if run_checks:
        print('has_duplicate_columns:', has_duplicate_columns(df.drop(columns=['input_filename'])))

    # Show that in general, the cells are not strictly in sorted order
    if run_checks:
        is_sorted = df['cell_id'].iloc[::2].reset_index(drop=True).equals(df['cell_id'].iloc[1::2].reset_index(drop=True))
        print('is_sorted:', is_sorted)

    # Sort them appropriately and then check again
    df.sort_values(by=['well_id', 'cell_type', 'REEC', 'T', 'cell_id'], inplace=True)
    if run_checks:
        is_sorted = df['cell_id'].iloc[::2].reset_index(drop=True).equals(df['cell_id'].iloc[1::2].reset_index(drop=True))
        print('is_sorted:', is_sorted)

    # Sample usage of the transformation function
    if run_checks:
        df1 = pd.DataFrame(
            {
                'a': [1, 1, 2, 2, 3, 3],
                'c': [1, 2, 1, 2, 1, 2],
                'd': ['AA', 'BB', 'CC', 'DD', 'EE', 'FF'],
                'e': ['AAA', 'BBB', 'CCC', 'DDD', 'EEE', 'FFF'],
                'b': ['b1', 'b1', 'b2', 'b2', 'b3', 'b3'],
                'f': ['f1', 'f1', 'f2', 'f2', 'f3', 'f3'],
            }
        )
        df2 = pd.DataFrame(
            {
                'a': [1, 2, 3],
                'd1': ['AA', 'CC', 'EE'],
                'd2': ['BB', 'DD', 'FF'],
                'e1': ['AAA', 'CCC', 'EEE'],
                'e2': ['BBB', 'DDD', 'FFF'],
                'b': ['b1', 'b2', 'b3'],
                'f': ['f1', 'f2', 'f3'],
            }
        )
        assert df2.equals(transform_dataframe(df1, 'a', 'c', ['d', 'e'], ['b', 'f'], verbose=True)), 'The transformation function does not work as expected.'

    # Compress the dataframe prior to transforming it
    df = utils.downcast_dataframe_dtypes(df)
    print('Shape of dataframe prior to transforming it:', df.shape)

    # Programmatically get the columns with duplicated data in sequential rows and manually confirm the result
    distinct_old_row_identifier_col = 'Ch'
    new_index_col = 'row_id'
    cols_with_repeated_rows2 = []
    cols_with_unique_rows2 = []
    for column in df.columns:
        ser = df[column]
        if ser.iloc[::2].reset_index(drop=True).equals(ser.iloc[1::2].reset_index(drop=True)):
            cols_with_repeated_rows2.append(column)
        else:
            cols_with_unique_rows2.append(column)
    cols_with_repeated_rows2 = [col for col in cols_with_repeated_rows2 if col not in [distinct_old_row_identifier_col, new_index_col]]
    cols_with_unique_rows2 = [col for col in cols_with_unique_rows2 if col not in [distinct_old_row_identifier_col, new_index_col]]

    # These two definitions are part of a manual check that, if passing consistently, can likely be removed in the future
    if run_checks and run_repeated_rows_check:
        # These sample rows are expected a unified dataframe but not for one processed by dataset_formats.py (which adds more columns), so in general we'll force this check to be skipped
        cols_with_repeated_rows = ['Image ID_(standardized)', 'Centroid X (µm)_(standardized)', 'Centroid Y (µm)_(standardized)', 'Area', 'X', 'Y', 'Perim.', 'Circ.', 'AR', 'Round', 'Solidity', 'input_filename', 'tif_name', 'T', 'cell_id', 'REEC', 'well_id', 'cell_type']
        cols_with_unique_rows = [' ', 'Label', 'Mean', 'StdDev', 'RawIntDen', 'RawIntNorm']
        print('Repeated:', cols_with_repeated_rows2)
        print('Unique:', cols_with_unique_rows2)
        assert cols_with_repeated_rows2 == cols_with_repeated_rows, 'The columns with repeated rows are not as expected.'
        assert cols_with_unique_rows2 == cols_with_unique_rows, 'The columns with unique rows are not as expected.'

    # Create a generator that transforms the dataframes image-by-image
    df_generator = transform_dataframes_in_chunks(df, image_col, new_index_col=new_index_col, distinct_old_row_identifier_col=distinct_old_row_identifier_col, cols_with_unique_rows=cols_with_unique_rows2, cols_with_repeated_rows=cols_with_repeated_rows2, verbose=False)

    # Concatenate the DataFrames generated by the generator into a single DataFrame
    df_transformed = pd.concat(df_generator, ignore_index=True)

    # Print the memory usage of the new dataframe
    if run_checks:
        df_transformed.info(memory_usage='deep')
        print('Shape of df_transformed:', df_transformed.shape)

    # Get the signal intensity and signal Z score column names
    signal_channels = set(df['Ch'].unique()) - {nuclear_channel}
    signal_intensity_columns = [f'RawIntNorm{channel}' for channel in signal_channels]
    signal_zscore_columns = [f'{col}_zscore' for col in signal_intensity_columns]
    print('signal_intensity_columns:', signal_intensity_columns)
    print('signal_zscore_columns:', signal_zscore_columns)

    # Calculate the Z scores on each area-normalized signal channel for each image
    for curr_image in unique_images:
        for curr_column in signal_intensity_columns:
            curr_image_loc = df_transformed[image_col] == curr_image
            df_transformed.loc[curr_image_loc, curr_column + '_zscore'] = scipy.stats.zscore(df_transformed.loc[curr_image_loc, curr_column])
    print('Columns of df_transformed that end with "_zscore":', df_transformed.loc[:, df_transformed.columns.str.endswith('_zscore')])

    # Filter out cells for each image that are more than three standard deviations from the mean for any channel
    if do_z_score_filter:
        num_cells_before = df_transformed[image_col].value_counts()
        df_transformed = df_transformed[(df_transformed[signal_zscore_columns].abs() < z_score_filter_threshold).all(axis='columns')]
        num_cells_removed = num_cells_before - df_transformed[image_col].value_counts()
        print('Number of cells removed by Z score threshold:')
        print(num_cells_removed)

    # Return the transformed dataframe
    return df_transformed
