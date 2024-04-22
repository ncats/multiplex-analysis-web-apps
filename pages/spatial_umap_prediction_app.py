# Import relevant libraries
import numpy as np
import pandas as pd
import time
import neighbors_counts_for_neighborhood_profiles_orig
import os
import streamlit as st
import streamlit_dataframe_editor as sde
import app_top_of_page as top
import sit_03a_Tool_parameter_selection as sit
import new_phenotyping_lib
import utils
import anndata


def load_input_dataset(df_input):
    """Creates an AnnData object from a pandas DataFrame.

    Args:
        df_input (pandas.DataFrame): The input dataframe containing the dataset.

    Returns:
        anndata.AnnData: The loaded input dataset.
    """

    # Create an AnnData object from the input dataframe in the recommended way
    adata = anndata.AnnData(X=df_input[['Cell X Position', 'Cell Y Position']].values, obs=df_input.drop(columns=['Cell X Position', 'Cell Y Position']), var=pd.DataFrame(index=['Cell X Position', 'Cell Y Position']))

    # Return the AnnData object
    return adata


def load_input_dataset_interactive():
    """Loads the input dataset from the session state and returns the AnnData object.

    Returns:
        anndata.AnnData: The loaded input dataset.
    """

    # If 'input_dataset' isn't in the session state, print an error message and return
    if 'input_dataset' not in st.session_state:
        st.error('An input dataset has not yet been opened. Please do so using the "Open File" page in the sidebar.')
        return
    
    # Load the input dataset
    with st.button('Load input dataset'):
        with st.spinner('Loading input dataset...'):
            adata = load_input_dataset(st.session_state['input_dataset'].data)

    # Return the AnnData object
    return adata


def apply_phenotyping(phenotyping_method, df_input, df_pheno_assignments, remove_allneg_phenotypes=True):
    """
    Apply phenotyping to a given dataframe.

    Args:
        phenotyping_method (str): The phenotyping method to use.
        df_input (pandas.DataFrame): The input dataframe containing cell data.
        df_pheno_assignments (pandas.DataFrame): The dataframe containing phenotype assignments.
        remove_allneg_phenotypes (bool): Whether to remove all negative phenotypes. Default is True.

    Returns:
        df_input (pandas.DataFrame): The input dataframe with phenotyping applied.
        phenotype_colname (str): The column name for the phenotype.
    """

    # If the Phenotyper's phenotyping method is "Custom", create a phenotype assignments file (and assign the corresponding setting) from its phenotype assignment table
    if phenotyping_method == 'Custom':
        df_pheno_assignments = df_pheno_assignments[df_pheno_assignments['phenotype'] != 'unassigned']
        phenotype_identification_file = sit.create_phenotype_assignments_file_from_phenotyper(df_pheno_assignments)
    else:
        phenotype_identification_file = None

    # Apply phenotyping
    df_input, _, species_int_to_pheno_name = new_phenotyping_lib.apply_phenotyping(df_input, phenotyping_method, phenotype_identification_file, species_int_colname='Species int', remove_allneg_phenotypes=remove_allneg_phenotypes)

    # Convert the integers identifying the species to their corresponding useful names
    df_input['Species int'] = df_input['Species int'].map(species_int_to_pheno_name)
    phenotype_colname = 'pheno_' + utils.get_timestamp()
    df_input.rename(columns={'Species int': phenotype_colname}, inplace=True)

    # Return the dataframe and the name of the phenotype column
    return df_input, phenotype_colname


def apply_phenotyping_interactive(df_input, remove_allneg_phenotypes=True):
    """
    Apply phenotyping to a given dataframe using the phenotyping method (and potentially a phenotype assignments file) from the session state.

    Args:
        df_input (pandas.DataFrame): The input dataframe containing cell data.
        remove_allneg_phenotypes (bool): Whether to remove all negative phenotypes. Default is True.

    Returns:
        df_input (pandas.DataFrame): The input dataframe with phenotyping applied.
        phenotype_colname (str): The column name for the phenotype.
    """

    # If 'phenoMeth' isn't in the session state, print an error message and return
    if 'phenoMeth' not in st.session_state:
        st.error('A phenotyping method has not yet been selected. Please do so using the "Thresholded Intensities" page in the sidebar.')
        return

    # Create a button and spinner
    if st.button('Apply phenotyping'):
        with st.spinner('Applying phenotyping...'):

            # Get necessary variables from the session state
            phenotyping_method = st.session_state['phenoMeth']
            df_pheno_assignments = st.session_state['pheno__de_phenotype_assignments'].reconstruct_edited_dataframe() if phenotyping_method == 'Custom' else None

            # Apply the phenotyping
            df_input, phenotype_colname = apply_phenotyping(phenotyping_method, df_input, df_pheno_assignments, remove_allneg_phenotypes=remove_allneg_phenotypes)

    # Return the dataframe and the name of the phenotype column
    return df_input, phenotype_colname


def calculate_densities(df, radius_edges=[0, 25, 50, 100, 150, 200], spatial_x_colname='Cell X Position', spatial_y_colname='Cell Y Position', image_colname='Slide ID', phenotype_colname='pheno_20240327_152849', debug_output=False, num_cpus_to_use=7, cast_to_float32=False, output_dir=os.path.join('.', 'output'), counts_matrix_csv_filename='counts_matrix.csv'):
    """
    Calculate the spatial density matrix for a given dataframe.

    Args:
        df (pandas.DataFrame): The input dataframe containing cell data.
        radius_edges (list): The edges of the radius bins for density calculation. Default is [0, 25, 50, 100, 150, 200].
        spatial_x_colname (str): The column name for the X position of cells. Default is 'Cell X Position'.
        spatial_y_colname (str): The column name for the Y position of cells. Default is 'Cell Y Position'.
        image_colname (str): The column name for the image ID. Default is 'Slide ID'.
        phenotype_colname (str): The column name for the phenotype. Default is 'pheno_20240327_152849'.
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

    # Return the final density matrix and the timing string for counting neighbors and creating the counts matrix
    return df_density_matrix, timing_string_counts


def calculate_densities_interactive(df_input, radius_edges=[0, 25, 50, 100, 150, 200], spatial_x_colname='Cell X Position', spatial_y_colname='Cell Y Position', image_colname='Slide ID', phenotype_colname='pheno_20240327_152849', debug_output=False, num_cpus_to_use=7, cast_to_float32=False, output_dir=os.path.join('.', 'output'), counts_matrix_csv_filename='counts_matrix.csv'):
    """
    Calculate densities interactively.

    Args:
        df_input (pandas.DataFrame): Input dataframe containing spatial information.
        radius_edges (list, optional): List of radius edges for density calculation. Defaults to [0, 25, 50, 100, 150, 200].
        spatial_x_colname (str, optional): Column name for X position. Defaults to 'Cell X Position'.
        spatial_y_colname (str, optional): Column name for Y position. Defaults to 'Cell Y Position'.
        image_colname (str, optional): Column name for slide ID. Defaults to 'Slide ID'.
        phenotype_colname (str, optional): Column name for phenotype. Defaults to 'pheno_20240327_152849'.
        debug_output (bool, optional): Enable debug output. Defaults to False.
        num_cpus_to_use (int, optional): Number of CPUs to use for calculation. Defaults to 7.
        cast_to_float32 (bool, optional): Cast data to float32. Defaults to False.
        output_dir (str, optional): Output directory path. Defaults to './output'.
        counts_matrix_csv_filename (str, optional): Filename for counts matrix CSV. Defaults to 'counts_matrix.csv'.

    Returns:
        pandas.DataFrame: Final density matrix.
        str: Timing string for counting neighbors and creating the counts matrix.
    """

    # Create a button and spinner
    if st.button('Calculate densities'):
        with st.spinner('Calculating densities...'):

            # Calculate the densities
            df_density_matrix, timing_string_counts = calculate_densities(df_input, radius_edges=radius_edges, spatial_x_colname=spatial_x_colname, spatial_y_colname=spatial_y_colname, image_colname=image_colname, phenotype_colname=phenotype_colname, debug_output=debug_output, num_cpus_to_use=num_cpus_to_use, cast_to_float32=cast_to_float32, output_dir=output_dir, counts_matrix_csv_filename=counts_matrix_csv_filename)

    # Return the final density matrix and the timing string for counting neighbors and creating the counts matrix
    return df_density_matrix, timing_string_counts


def train_test_split(df_input):
    """
    Split the input dataframe into training and testing data.

    TODO: Replace this with the actual split.

    Args:
        df_input (pandas.DataFrame): The input dataframe.

    Returns:
        df_train (pandas.DataFrame): The training dataframe.
        df_test (pandas.DataFrame): The testing dataframe.
    """

    # TODO: Replace this with the actual split
    df_train = df_input
    df_test = df_input

    # Return the training and testing dataframes
    return df_train, df_test


def train_test_split_interactive(df_input):
    """
    Split the input dataframe into training and testing data interactively.

    Args:
        df_input (pandas.DataFrame): The input dataframe.

    Returns:
        df_train (pandas.DataFrame): The training dataframe.
        df_test (pandas.DataFrame): The testing dataframe.
    """

    # Create a button and spinner
    if st.button('Split input dataset into training and testing sets'):
        with st.spinner('Splitting input dataset into training and testing sets...'):

            # Split the input dataframe into training and testing data
            df_train, df_test = train_test_split(df_input)

    # Return the training and testing dataframes
    return df_train, df_test


def sample_index(index, n_samples):
    """
    Randomly sample `n_samples` indices from the given `index` array without replacement, returning the sorted result.

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


def get_train_subset_indices_interactive(df_train, smallest_image_size_frac, num_analysis_subsets, image_colname='Slide ID'):
    """
    Get the training subset indices for fitting and analyzing the UMAP model interactively.

    Args:
        df_train (pandas.DataFrame): The training dataframe.
        smallest_image_size_frac (float): The fraction of the smallest image size to use as the number of samples per image.
        num_analysis_subsets (int): The number of analysis subsets to generate.

    Returns:
        tuple: A tuple containing the following:
            - df_subset_indices_train (pandas.DataFrame): A dataframe with the training and analysis subset indices.
            - num_samples_per_image (int): The number of samples per image.
    """

    # Create a button and spinner
    if st.button('Get training subset indices'):
        with st.spinner('Getting training subset indices...'):

            # Get the training subset indices
            df_subset_indices_train, num_samples_per_image = get_train_subset_indices(df_train, smallest_image_size_frac, num_analysis_subsets, image_colname=image_colname)

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


def get_test_subset_indices_interactive(df_test, num_samples_per_image, num_test_subsets, image_colname='Slide ID'):
    """
    Get index subsets from a test DataFrame interactively.

    Args:
        df_test (pandas.DataFrame): The DataFrame containing the test data.
        num_samples_per_image (int): The number of samples to be taken from each image.
        num_test_subsets (int): The number of test subsets to be generated.

    Returns:
        pandas.DataFrame: A DataFrame containing the test indices for each subset.
    """

    # Create a button and spinner
    if st.button('Get test subset indices'):
        with st.spinner('Getting test subset indices...'):

            # Get the test subset indices
            df_subset_indices_test = get_test_subset_indices(df_test, num_samples_per_image, num_test_subsets, image_colname=image_colname)

    # Return the test indices
    return df_subset_indices_test


# Main function
def main():
    """
    Main function for the page.
    """

    # Parameters
    # TODO: Make these widgets when they're to first be used
    remove_allneg_phenotypes = False
    smallest_image_size_frac = 0.1
    num_analysis_subsets = 10
    num_test_subsets = 10

    # Load the input dataset. Since it's required to run Open File first, this will return None if it hasn't been run yet. This uses Python's "new" walrus operator :=
    if (df_input := load_input_dataset_interactive()) is None: return

    # Apply phenotyping to the input dataset. Since it's required to run the Phenotyper first, this will return None if it hasn't been run yet. Note the walrus operator := doesn't support tuple unpacking, which is why we can't make the following more concise
    df_input, phenotype_colname = apply_phenotyping_interactive(df_input, remove_allneg_phenotypes=remove_allneg_phenotypes) or (None, None)
    if df_input is None: return

    # Calculate the densities
    df_density_matrix, timing_string_counts = calculate_densities_interactive(df_input, phenotype_colname=phenotype_colname)

    # Split the input images into training and testing images
    # TODO: Replace this with the actual split
    df_train, df_test = train_test_split_interactive(df_input)

    # Get the subset indices for fitting and analyzing the UMAP model
    df_subset_indices_train, num_samples_per_image = get_train_subset_indices_interactive(df_train, smallest_image_size_frac, num_analysis_subsets)

    # Get the subset indices for testing the overall model
    df_subset_indices_test = get_test_subset_indices_interactive(df_test, num_samples_per_image, num_test_subsets)


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
