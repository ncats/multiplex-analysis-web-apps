# Import relevant libraries
import numpy as np
import pandas as pd
import time
import neighbors_counts_for_neighborhood_profiles_orig
import os
import streamlit as st
import streamlit_dataframe_editor as sde
import app_top_of_page as top


def calculate_densities(df, radius_edges=[0, 25, 50, 100, 150, 200], spatial_x_colname='Cell X Position', spatial_y_colname='Cell Y Position', image_colname='Slide ID', phenotype_colname='XXXX', debug_output=False, num_cpus_to_use=7, cast_to_float32=False, output_dir=os.path.join('.', 'output'), counts_matrix_csv_filename='counts_matrix.csv'):
    """
    Calculate the spatial density matrix for a given dataframe.

    Args:
        df (pandas.DataFrame): The input dataframe containing cell data.
        radius_edges (list): The edges of the radius bins for density calculation. Default is [0, 25, 50, 100, 150, 200].
        spatial_x_colname (str): The column name for the X position of cells. Default is 'Cell X Position'.
        spatial_y_colname (str): The column name for the Y position of cells. Default is 'Cell Y Position'.
        image_colname (str): The column name for the image ID. Default is 'Slide ID'.
        phenotype_colname (str): The column name for the phenotype. Default is 'XXXX'.
        debug_output (bool): Whether to enable debug output. Default is False.
        num_cpus_to_use (int): The number of CPUs to use for parallel processing. Default is 7.
        cast_to_float32 (bool): Whether to cast the density matrix to float32. Default is False.

    Returns:
        df_density_matrix (pandas.DataFrame): The spatial density matrix.
        timing_string_counts (str): The timing information for counting neighbors and creating the counts matrix.
    """
    
    # Variables
    image_names = df[image_colname].unique()
    phenotypes = df[phenotype_colname].unique()
    radii = np.array(radius_edges)
    num_ranges = len(radii) - 1
    range_strings = ['({}, {}]'.format(radii[iradius], radii[iradius + 1]) for iradius in range(num_ranges)]
    counts_matrix_csv_path = os.path.join(output_dir, counts_matrix_csv_filename)

    # Calculate the counts matrix for all images if it hasn't already been calculated
    if not os.path.exists(counts_matrix_csv_path):
        print(f'Counting neighbors around {df.shape[0]} cells...')
        start_time = time.time()
        df_counts_matrix = neighbors_counts_for_neighborhood_profiles_orig.calculate_density_matrix_for_all_images(image_names, df, phenotypes, phenotype_colname, image_colname, [spatial_x_colname, spatial_y_colname], radii, num_ranges, range_strings, debug_output=debug_output, num_cpus_to_use=num_cpus_to_use, swap_inequalities=True, cast_to_float32=cast_to_float32)
        timing_string_counts = f'Time to count neighbors resulting in a counts matrix of shape {df_counts_matrix.shape}: {int(time.time() - start_time)} seconds'
        print(timing_string_counts)
        df_counts_matrix.to_csv(counts_matrix_csv_path, index=False)
    else:
        df_counts_matrix = pd.read_csv(counts_matrix_csv_path)
        timing_string_counts = f'Counts matrix already exists at {counts_matrix_csv_path}; not re-calculated and simply read back in from disk'

    # Fill in any NaN values with 0 and convert to integers
    df_counts_matrix = df_counts_matrix.fillna(0).astype(int)

    # Get the areas of the 2D annuli/rings
    radius_range_holder = []
    area_holder = []
    for radius_range in list(set([column.split(' in range ')[1] for column in df_counts_matrix.columns])):
        r1 = float(radius_range.split(', ')[0][1:])
        r2 = float(radius_range.split(', ')[1][:-1])
        area = np.pi * (r2 ** 2 - r1 ** 2)
        radius_range_holder.append(radius_range)
        area_holder.append(area)
    ser_areas = pd.Series(area_holder, index=radius_range_holder).sort_values()

    # Get the areas for each column
    areas_for_columns = [ser_areas[column.split(' in range ')[1]] for column in df_counts_matrix.columns]

    # Divide the columns of df_counts_matrix by areas_for_columns to get the spatial density in units of counts per unit area
    df_density_matrix = df_counts_matrix.div(pd.Series(areas_for_columns, index=df_counts_matrix.columns), axis=1)

    # Rename the columns
    new_column_names = ['spatial_umap_density of ' + column.removeprefix('Number of neighbors of type ') for column in df_counts_matrix.columns]
    df_density_matrix = df_density_matrix.rename(columns=dict(zip(df_counts_matrix.columns, new_column_names)))

    # Print out the histogram of the final density matrix
    print(np.histogram(df_density_matrix))

    # Return the final density matrix
    return df_density_matrix, timing_string_counts


def sample_index(index, n_samples):
    """
    Randomly samples `n_samples` indices from the given `index` array without replacement, returning the sorted result.

    This should be applied per image.

    Args:
        index (array-like): The array of indices to sample from.
        n_samples (int): The number of indices to sample.

    Returns:
        numpy.array: An array of randomly sampled indices.
    """
    return np.sort(np.random.choice(index, n_samples, replace=False))


def get_train_subset_indices(df_train, smallest_image_size_frac, num_analysis_subsets, image_colname='Slide ID'):
    """
    Get the training subset indices for fitting and analyzing the UMAP model.

    Args:
        df_train (pandas.DataFrame): The training dataframe.
        smallest_image_size_frac (float): The fraction of the smallest image size to use as the number of samples per image.
        num_analysis_subsets (int): The number of analysis subsets to generate.

    Returns:
        tuple: A tuple containing the following:
            - df_subset_indices_train (pandas.DataFrame): A dataframe with the training and analysis subset indices.
            - num_samples_per_image (int): The number of samples per image.

    Raises:
        AssertionError: If the concatenation of the indices has added rows.

    """

    # Group the training data grouped by image
    df_train_by_image = df_train.groupby(image_colname)
    
    # Get the desired number of samples per image
    num_samples_per_image = int(df_train_by_image.size().min() * smallest_image_size_frac)

    # Get a sampling of df_train that samples rows from each image equally
    train_indices = df_train_by_image.apply(lambda df_image: sample_index(df_image.index, num_samples_per_image)).explode()  # goes into df_train. This is a series whose index is the image ID and the values are the indices of the samples. This should be used for fitting the UMAP model.

    # Do the same for num_analysis_subsets subsets but first ensuring that the training indices are not included in the possible analysis indices
    subset_indices_holder = [train_indices]
    df_train_only_analysis_indices_by_image = df_train.drop(train_indices).groupby(image_colname)
    for _ in range(num_analysis_subsets):
        subset_indices_holder.append(df_train_only_analysis_indices_by_image.apply(lambda df_image: sample_index(df_image.index, num_samples_per_image)).explode())

    # Get one big dataframe with all the indices, both "train" and analysis
    df_subset_indices_train = pd.concat(subset_indices_holder, axis='columns', columns=['train_subset'] + [f'analysis_subset_{i}' for i in range(num_analysis_subsets)])

    # Ensure that the concatenation hasn't added any rows
    assert len(df_subset_indices_train) == len(train_indices), f'The concatenation of the indices has added rows, going from {len(train_indices)} to {len(df_subset_indices_train)}'

    # Return the training indices and the number of samples per image
    return df_subset_indices_train, num_samples_per_image


def get_test_subset_indices(df_test, num_samples_per_image, num_test_subsets, image_colname='Slide ID'):
    """
    Get index subsets from a test DataFrame.

    Args:
        df_test (pandas.DataFrame): The DataFrame containing the test data.
        num_samples_per_image (int): The number of samples to be taken from each image.
        num_test_subsets (int): The number of test subsets to be generated.

    Returns:
        pandas.DataFrame: A DataFrame containing the test indices for each subset.

    Raises:
        AssertionError: If the concatenation of the indices has added rows.

    """

    # Get a sampling of df_test that samples rows from each image equally
    subset_indices_holder = []
    df_test_by_image = df_test.groupby(image_colname)
    for _ in range(num_test_subsets):
        subset_indices_holder.append(df_test_by_image.apply(lambda df_image: sample_index(df_image.index, num_samples_per_image)).explode())

    # Get one big dataframe with all the test indices
    df_subset_indices_test = pd.concat(subset_indices_holder, axis='columns', columns=[f'test_subset_{i}' for i in range(num_test_subsets)])

    # Ensure that the concatenation hasn't added any rows
    assert len(df_subset_indices_test) == len(subset_indices_holder[0]), f'The concatenation of the indices has added rows, going from {len(subset_indices_holder[0])} to {len(df_subset_indices_test)}'

    # Return the test indices
    return df_subset_indices_test


# Main function
def main():
    """
    Main function for the page.
    """

    # If 'input_dataset' isn't in the session state, print an error message and return
    if 'input_dataset' not in st.session_state:
        st.error('An input dataset has not yet been opened. Please do so using the "Open File" page in the sidebar.')
        return
    else:
        df_input = st.session_state['input_dataset'].data

    # Parameters
    # TODO: Make these widgets
    smallest_image_size_frac = 0.1
    num_analysis_subsets = 10
    num_test_subsets = 10

    # Split the input dataframe into training and testing data
    # TODO: Replace this with the actual split
    df_train = df_input[df_input['is_train']]
    df_test = df_input[~df_input['is_train']]

    # Get the subset indices for fitting and analyzing the UMAP model
    df_subset_indices_train, num_samples_per_image = get_train_subset_indices(df_train, smallest_image_size_frac, num_analysis_subsets)

    # Get the subset indices for testing the overall model
    df_subset_indices_test = get_test_subset_indices(df_test, num_samples_per_image, num_test_subsets)


# Call the main function
if __name__ == '__main__':

    # Set page settings
    page_name = 'Spatial UMAP Prediction'
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
