# Import relevant libraries
import numpy as np
import pandas as pd
import time
import neighbors_counts_for_neighborhood_profiles_orig


def calculate_densities(df, radius_edges=[0, 25, 50, 100, 150, 200], spatial_x_colname='Cell X Position', spatial_y_colname='Cell Y Position', image_colname='Slide ID', phenotype_colname='XXXX', debug_output=False, num_cpus_to_use=7, cast_to_float32=False):

    # Variables
    image_names = df[image_colname].unique()
    phenotypes = df[phenotype_colname].unique()
    radii = np.array(radius_edges)
    num_ranges = len(radii) - 1
    range_strings = ['({}, {}]'.format(radii[iradius], radii[iradius + 1]) for iradius in range(num_ranges)]

    # Calculate the counts matrix for all images
    print(f'Counting neighbors around {df.shape[0]} cells...')
    start_time = time.time()
    df_counts_matrix = neighbors_counts_for_neighborhood_profiles_orig.calculate_density_matrix_for_all_images(image_names, df, phenotypes, phenotype_colname, image_colname, [spatial_x_colname, spatial_y_colname], radii, num_ranges, range_strings, debug_output=debug_output, num_cpus_to_use=num_cpus_to_use, swap_inequalities=True, cast_to_float32=cast_to_float32)
    timing_string_counts = f'Time to count neighbors resulting in a counts matrix of shape {df_counts_matrix.shape}: {int(time.time() - start_time)} seconds'
    print(timing_string_counts)

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


# This should be applied per image
def sample_index(index, n_samples):
    """
    Randomly samples `n_samples` indices from the given `index` array without replacement, returning the sorted result.

    Parameters:
        index (array-like): The array of indices to sample from.
        n_samples (int): The number of indices to sample.

    Returns:
        numpy.array: An array of randomly sampled indices.
    """
    return np.sort(np.random.choice(index, n_samples, replace=False))


# Main function
def main():

    # Parameters
    smallest_image_size_frac = 0.1
    num_analysis_subsets = 10
    num_test_subsets = 10

    # Junk load the training and testing data
    df_train = pd.read_csv('data/train.csv')
    df_test = pd.read_csv('data/test.csv')

    # Group the training data grouped by image
    df_train_by_image = df_train.groupby('Slide ID')
    
    # Get the desired number of samples per image
    smallest_image_size = df_train_by_image.size().min()
    num_samples_per_image = int(smallest_image_size * smallest_image_size_frac)

    # Get a sampling of df_train that samples rows from each image equally
    train_indices = df_train_by_image.apply(lambda df_image: sample_index(df_image.index, num_samples_per_image)).explode()  # goes into df_train. This is a series whose index is the image ID and the values are the indices of the samples. This should be used for fitting the UMAP model.

    # Do the same for num_analysis_subsets subsets but first ensuring that the training indices are not included in the possible analysis indices
    subset_indices_holder = [train_indices]
    df_train_only_analysis_indices_by_image = df_train.drop(train_indices).groupby('Slide ID')
    for _ in range(num_analysis_subsets):
        subset_indices_holder.append(df_train_only_analysis_indices_by_image.apply(lambda df_image: sample_index(df_image.index, num_samples_per_image)).explode())

    # Get one big dataframe with all the indices, both "train" and analysis
    df_subset_indices_train = pd.concat(subset_indices_holder, axis='columns', columns=['train_subset'] + [f'analysis_subset_{i}' for i in range(num_analysis_subsets)])

    # Ensure that the concatenation hasn't added any rows
    assert len(df_subset_indices_train) == len(train_indices), f'The concatenation of the indices has added rows, going from {len(train_indices)} to {len(df_subset_indices_train)}'

    # vvvv TEST vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    # Do the same for num_test_subsets subsets
    subset_indices_holder = []
    df_test_by_image = df_test.groupby('Slide ID')
    for _ in range(num_test_subsets):
        subset_indices_holder.append(df_test_by_image.apply(lambda df_image: sample_index(df_image.index, num_samples_per_image)).explode())

    # Get one big dataframe with all the test indices
    df_subset_indices_test = pd.concat(subset_indices_holder, axis='columns', columns=[f'test_subset_{i}' for i in range(num_test_subsets)])

    # Ensure that the concatenation hasn't added any rows
    assert len(df_subset_indices_test) == len(subset_indices_holder[0]), f'The concatenation of the indices has added rows, going from {len(subset_indices_holder[0])} to {len(df_subset_indices_test)}'
    # ^^^^ TEST ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


# Call the main function
if __name__ == '__main__':
    main()
