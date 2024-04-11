# Import relevant libraries
import numpy as np
import pandas as pd
import PlottingTools_orig as PlottingTools
import matplotlib.pyplot as plt
import time
import neighbors_counts_for_neighborhood_profiles_orig
import umap
from sklearn.preprocessing import MinMaxScaler
import datetime
import utils


# Add columns to partition a dataframe into train and test sets for the UMAP, roughly similar to Giraldo et. al. 2021. Note frac_train and frac_test cannot be larger than 0.5
def get_umap_train_and_test_sets(df, frac_train=0.5, frac_test=0.5, image_colname='ShortName', num_umap_test_sets=1):

    # Calculate the size of the train and test sets based on the size of the smallest image
    min_num_cells_per_image = df.groupby(image_colname).size().min()
    num_train_cells_per_image = int(np.floor(min_num_cells_per_image * frac_train))
    num_test_cells_per_image = int(np.floor(min_num_cells_per_image * frac_test))

    # Get a set of train indices and assign them to 'umap_train'
    train_indices = df.groupby(image_colname).apply(lambda x: x.sample(n=num_train_cells_per_image, replace=False).index, include_groups=False).explode().values
    df['umap_train'] = False
    df.loc[train_indices, 'umap_train'] = True

    # Get the subset of df that is not in the training set
    df_not_train = df[~df['umap_train']]

    # Delete all currently existing umap test columns, if any
    df.drop(columns=[column for column in df.columns if column.startswith('umap_test_')], inplace=True)

    # Pre-calculate group indices
    group_indices = df_not_train.groupby(image_colname).indices
    group_indices = {key: df_not_train.index[value] for key, value in group_indices.items()}  # convert the indices to the actual indices in the dataframe

    # Add 'umap_test' columns to df for each umap test set and initialize them to False, creating a dictionary of new columns
    new_columns = {
        f'umap_test_{iumap_test_set}': pd.Series(False, index=df.index)
        for iumap_test_set in range(num_umap_test_sets)
    }

    # Create a new DataFrame that includes the new columns
    df = pd.concat([df, pd.DataFrame(new_columns)], axis=1)

    # Loop through each umap test set
    for iumap_test_set in range(num_umap_test_sets):

        # Get a set of test indices and assign them to 'umap_test_i'
        for _, indices in group_indices.items():
            test_indices = np.random.choice(indices, size=num_test_cells_per_image, replace=False)
            df.loc[test_indices, f'umap_test_{iumap_test_set}'] = True

    # Return the dataframe with the UMAP train and test columns added
    return df


# Calculate a dictionary of clusters as keys and list of bin tuples as values using the normalized histogram differences between two conditions for a single set of UMAP "test" data
def calculate_difference_clusters(df, diff_cutoff_frac, umap_x_colname='UMAP_1_20230327_152849', umap_y_colname='UMAP_2_20230327_152849', binary_colname='Survival_5yr', num_umap_bins=200, plot_manual_histogram_diff=False, plot_diff_matrix=True, umap_test_colname='umap_test_0'):

    # Get the x and y UMAP ranges for all UMAP test sets
    umap_test_colnames = [column for column in df.columns if column.startswith('umap_test_')]
    umap_x_range = [9999, -9999]
    umap_y_range = [9999, -9999]
    for colname in umap_test_colnames:
        umap_x_range[0] = min(umap_x_range[0], df.loc[df[colname], umap_x_colname].min())
        umap_x_range[1] = max(umap_x_range[1], df.loc[df[colname], umap_x_colname].max())
        umap_y_range[0] = min(umap_y_range[0], df.loc[df[colname], umap_y_colname].min())
        umap_y_range[1] = max(umap_y_range[1], df.loc[df[colname], umap_y_colname].max())

    # Get a universal set of edges for the UMAPs
    edges_x = np.linspace(umap_x_range[0], umap_x_range[1], num_umap_bins + 1)
    edges_y = np.linspace(umap_y_range[0], umap_y_range[1], num_umap_bins + 1)

    # Get subsets of the full data that are the entire test set and both binary subsets of the test set
    df_test = df[df[umap_test_colname]]
    binary_values = sorted(list(df[binary_colname].unique()))
    df_binary_subset_0 = df_test[df_test[binary_colname] == binary_values[0]]
    df_binary_subset_1 = df_test[df_test[binary_colname] == binary_values[1]]

    # Get the 2D histograms for each condition
    d_binary_subset_0 = PlottingTools.plot_2d_density(df_binary_subset_0[umap_x_colname], df_binary_subset_0[umap_y_colname], bins=[edges_x, edges_y], return_matrix=True)
    d_binary_subset_1 = PlottingTools.plot_2d_density(df_binary_subset_1[umap_x_colname], df_binary_subset_1[umap_y_colname], bins=[edges_x, edges_y], return_matrix=True)

    # Get the difference between the histograms for the binary subsets
    d_diff = d_binary_subset_1 - d_binary_subset_0

    # "Mask" the difference matrix based on a cutoff
    cutoff = np.abs(d_diff).max() * diff_cutoff_frac
    d_diff[d_diff > cutoff] = 1
    d_diff[d_diff < -cutoff] = -1
    d_diff[(d_diff >= -cutoff) & (d_diff <= cutoff)] = 0

    # TODO: Implement further clustering of the 1 and -1 classes!

    # Get the clusters based only on the cutoff, deliberately not performing any clustering
    # clusters = {0: [tuple([x[1], x[0]]) for x in np.transpose(np.where(np.isclose(d_diff, -1)))], 1: [tuple([x[1], x[0]]) for x in np.transpose(np.where(np.isclose(d_diff, 1)))]}
    clusters = {
        -1: [tuple([x[1], x[0]]) for x in np.transpose(np.where(np.isclose(d_diff, -1)))],
        0: [tuple([x[1], x[0]]) for x in np.transpose(np.where(np.isclose(d_diff, 0)))],
        1: [tuple([x[1], x[0]]) for x in np.transpose(np.where(np.isclose(d_diff, 1)))]
        }

    # Plot the difference matrix
    if plot_diff_matrix:
        fig_d_diff, ax_d_diff = plt.subplots()
        PlottingTools.plot_2d_density(d_diff, bins=[edges_x, edges_y], n_pad=30, circle_type='arch', cmap=plt.get_cmap('bwr'), ax=ax_d_diff)
    else:
        fig_d_diff = None

    # Optionally plot the difference matrix manually
    if plot_manual_histogram_diff:
        fig, ax = plt.subplots()
        c = ax.pcolormesh(edges_x, edges_y, d_diff, cmap='bwr')
        fig.colorbar(c, ax=ax)
        fig_d_diff_manual = fig
    else:
        fig_d_diff_manual = None

    # Return the calculated clusters, edges, and figures
    return clusters, edges_x, edges_y, fig_d_diff, fig_d_diff_manual


# Determine the bins of the UMAP x and y coordinates of an input dataframe
def perform_umap_binning(df, edges_x, edges_y, image_colname='ShortName', spatial_x_colname='spatial_x', spatial_y_colname='spatial_y', umap_x_colname='umap_x', umap_y_colname='umap_y', property_colnames=['property_a', 'property_b']):

    # Sample values of boolean_subset_on_cells:
    #   * pd.Series(True, index=cells.index)  # --> all cells
    #   * cells['survival'] == 'dead'         # --> subset of cells

    # Start by defining a test dataframe containing just the cells in the test set from the original dataframe, where the index has been reset to a range but the original indices are stored for reference
    df_test_with_bins = df[[spatial_x_colname, spatial_y_colname, umap_x_colname, umap_y_colname, image_colname] + property_colnames]
    df_test_with_bins.index.name = 'cells_dataframe_index'
    df_test_with_bins = df_test_with_bins.reset_index(drop=False)

    # Get the histogram and bins to which the cells in the test set correspond
    # counts, _, _ = np.histogram2d(df_test[umap_x_colname].values, df_test[umap_y_colname].values, bins=[edges_x, edges_y])  # not actually used, but useful to keep here to see how it relates to the following calls to np.digitize()
    bin_index_x = np.digitize(df_test_with_bins[umap_x_colname].values, edges_x[:-1]) - 1
    bin_index_y = np.digitize(df_test_with_bins[umap_y_colname].values, edges_y[:-1]) - 1

    # Add the bin indices to the test dataframe
    df_test_with_bins['bin_index_x'] = bin_index_x.flatten()
    df_test_with_bins['bin_index_y'] = bin_index_y.flatten()

    # Return the dataframe
    return df_test_with_bins


# Given the UMAP bins of each cell and the labels assigned to the bins, assign the labels to a per-cell dataframe and to a transformed dataframe grouped by bin, i.e., assign the label to each bin in that transformed dataframe
def assign_cluster_labels(df_test_with_bins, cluster_labels, image_colname='ShortName', min_cells_per_bin=1):

    # Say we perform clustering on the bins and get a dictionary of cluster labels (0, 1, 2, ..., k-1) as keys and the indices of the bins in each cluster as values
    # Note these bin indices must correspond to the ones coming out of np.digitize(), this is crucial!!
    # cluster_labels = {0: [(3, 4), (4, 5)], 1: [(0, 1), (1, 2), (2, 3)], 2: [(5, 1), (0, 7), (7, 8)]}

    # Loop through each cluster and add a column to df_test containing the cluster label for each cell, if a label has been applied to the cell
    for cluster_label, cluster_bins in cluster_labels.items():
        cells_dataframe_indices_for_curr_cluster = df_test_with_bins.loc[df_test_with_bins.set_index(['bin_index_x', 'bin_index_y']).index.isin(cluster_bins), 'cells_dataframe_index'].values
        df_test_with_bins.loc[df_test_with_bins['cells_dataframe_index'].isin(cells_dataframe_indices_for_curr_cluster), 'cluster_label'] = cluster_label

    # Invert cluster_labels to get a dictionary of bin indices as keys and cluster labels as values, and convert to a Series
    bin_cluster_labels = {bin_index: cluster_label for cluster_label, bins in cluster_labels.items() for bin_index in bins}
    bin_cluster_labels = pd.Series(bin_cluster_labels, name='cluster_label')

    # Group by bin and obtain the unique images in each group, in addition to the bin means, stds, and number of cells in each bin
    df_grouped = df_test_with_bins.rename(columns={'cluster_label': 'cluster_label_by_cell'}).groupby(['bin_index_x', 'bin_index_y'])
    bin_unique_images_within = df_grouped[image_colname].agg(set)
    bin_unique_images_within.name = "unique_images"
    df_grouped_columns_less_image_colname = [col for col in df_test_with_bins.columns if col != image_colname]  # get all the columns in df_grouped except for image_colname
    df_grouped_columns_less_image_colname = [col if col != 'cluster_label' else 'cluster_label_by_cell' for col in df_grouped_columns_less_image_colname]  # since we want to use the renamed version of "cluster_label"
    bin_means = df_grouped[df_grouped_columns_less_image_colname].mean()
    bin_stds = df_grouped[df_grouped_columns_less_image_colname].std()
    bin_counts = df_grouped[df_grouped_columns_less_image_colname].size()  # get the number of test cells in each bin. This could be useful if we want to e.g. only use a bin with a minimum number of cells
    bin_counts.name = 'num_cells'

    # Concatenate the grouped data together
    bin_means = pd.concat([bin_means, bin_counts, bin_cluster_labels, bin_unique_images_within], axis='columns')
    bin_stds = pd.concat([bin_stds, bin_counts, bin_cluster_labels, bin_unique_images_within], axis='columns')

    # Drop rows from bin_means and bin_stds where the number of cells in the bin is less than min_cells_per_bin
    bin_means = bin_means[bin_means['num_cells'] >= min_cells_per_bin]
    bin_stds = bin_stds[bin_stds['num_cells'] >= min_cells_per_bin]

    # Return the plotly figures
    return bin_means, df_test_with_bins


def run_density_calculation(df, radius_edges=[0, 25, 50, 100, 150, 200], spatial_x_colname='CentroidX', spatial_y_colname='CentroidY', image_colname='ShortName', phenotype_colname='pheno_20230327_152849', debug_output=False, num_cpus_to_use=7, cast_to_float32=False):

    # Variables
    image_names = df[image_colname].unique()
    phenotypes = df[phenotype_colname].unique()
    radii = np.array(radius_edges)
    num_ranges = len(radii) - 1
    range_strings = ['({}, {}]'.format(radii[iradius], radii[iradius + 1]) for iradius in range(num_ranges)]

    # Calculate the density matrix for all images
    print(f'Counting neighbors around {df.shape[0]} cells...')
    start_time = time.time()
    df_density_matrix = neighbors_counts_for_neighborhood_profiles_orig.calculate_density_matrix_for_all_images(image_names, df, phenotypes, phenotype_colname, image_colname, [spatial_x_colname, spatial_y_colname], radii, num_ranges, range_strings, debug_output=debug_output, num_cpus_to_use=num_cpus_to_use, swap_inequalities=True, cast_to_float32=cast_to_float32)
    timing_string_density = f'Time to count neighbors resulting in a counts matrix of shape {df_density_matrix.shape}: {int(time.time() - start_time)} seconds'
    print(timing_string_density)

    # Fill in any NaN values with 0 and convert to integers
    df_density_matrix = df_density_matrix.fillna(0).astype(int)

    # Get the areas of the 2D annuli/rings
    radius_range_holder = []
    area_holder = []
    for radius_range in list(set([column.split(' in range ')[1] for column in df_density_matrix.columns])):
        r1 = float(radius_range.split(', ')[0][1:])
        r2 = float(radius_range.split(', ')[1][:-1])
        area = np.pi * (r2 ** 2 - r1 ** 2)
        radius_range_holder.append(radius_range)
        area_holder.append(area)
    ser_areas = pd.Series(area_holder, index=radius_range_holder).sort_values()

    # Get the areas for each column
    areas_for_columns = [ser_areas[column.split(' in range ')[1]] for column in df_density_matrix.columns]

    # Divide the columns of df_density_matrix by areas_for_columns to get the spatial density in units of counts per unit area (despite it being called "density_matrix" above... it's really a counts matrix)
    df_density_matrix_divided = df_density_matrix.div(pd.Series(areas_for_columns, index=df_density_matrix.columns), axis=1)

    # Rename the columns
    new_column_names = ['spatial_umap_density of ' + column.removeprefix('Number of neighbors of type ') for column in df_density_matrix.columns]
    df_density_matrix_divided = df_density_matrix_divided.rename(columns=dict(zip(df_density_matrix.columns, new_column_names)))

    # # Normalize the density matrix --> do this before the UMAP instead
    # df_density_matrix_divided = df_density_matrix_divided / df_density_matrix_divided.max().max()

    # # Convert to float32 --> float64 seems to better represent the densities perhaps
    # df_density_matrix_divided = df_density_matrix_divided.astype(np.float32)

    # Print out the histogram of the final density matrix
    print(np.histogram(df_density_matrix_divided))

    # Return the final density matrix
    return df_density_matrix_divided, timing_string_density


# Calculate `umap_transform = umap_fit.transform(df_density)` in batches, since during an intermediate step, this tries to create an ndarray of size (len(df_density), num_rows_in_training_data) and therefore can easily run out of memory
def umap_transform_in_batches(umap_fit, df_density, nworkers=7, max_arr_size_in_mb=200):

    # Get a reasonable batch size that ensures the temporary matrix created during the transformation doesn't reach max_arr_size_in_mb
    max_arr_size_in_bytes = max_arr_size_in_mb * 1024 * 1024  # 1 MB = 1024 KB = 1024^2 bytes
    max_num_elements_in_arr = max_arr_size_in_bytes / 8  # 8 bytes per float64 element

    # Get the number of rows in the training data
    num_rows_in_training_data = umap_fit.n_samples_fit_
    print('Using as the number of rows in the training data:', num_rows_in_training_data)

    # Number of elements in a row of the array
    batch_size = int(max_num_elements_in_arr / num_rows_in_training_data)

    # Initialize the arguments to umap_fit.transform()... this is a list of tuples
    original_data = []

    # For each batch...
    # transformed_data = []  # --> if we didn't use parallelism
    for i in range(0, len(df_density), batch_size):

        # Store the corresponding slice of the dataframe to transform
        # print(f'Transforming batch {i // batch_size + 1} of {len(df_density) // batch_size + 1}...') = []  # --> if we didn't use parallelism
        batch = df_density.iloc[i:i+batch_size, :]
        original_data.append((batch,))
        # transformed_batch = umap_fit.transform(batch) = []  # --> if we didn't use parallelism
        # transformed_data.append(transformed_batch) = []  # --> if we didn't use parallelism

    # Perform the transformation in parallel
    transformed_data = utils.execute_data_parallelism_potentially(function=umap_fit.transform, list_of_tuple_arguments=original_data, nworkers=nworkers, task_description='UMAP transformation of a density matrix', do_benchmarking=True, mp_start_method=None, use_starmap=True)

    # Concatenate the results of all the batches together
    umap_transform = np.concatenate(transformed_data)

    # Put the UMAP results into a new dataframe
    current_datetime = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    df_umap = pd.DataFrame({
        f'UMAP_1_{current_datetime}': umap_transform[:, 0],
        f'UMAP_2_{current_datetime}': umap_transform[:, 1]}, index=df_density.index)

    # Return the UMAP dataframe
    return df_umap


# Perform UMAP on the density matrix of a dataframe, fitting on training cells and transforming on either all the testing cells or the entire dataset
# TODO: See options to implement below
def run_umap_calculation(df, nworkers=7, transform_full_dataset=False):

    # Options to implement:
    #   * Perform supervised UMAP
    #   * Convert to float32 before running UMAP
    #   * Understand why the transformations are taking so long in parallel and even in serial (I thought Dante said they were fast... maybe it's the area transformation that Giraldo does?)... I get roughly 33 min for whole dataset in serial and 16 min for whole data in parallel using 7 workers

    # Get a copy (since we're about to normalize) of the density matrix of just the density columns
    df_density = df.loc[:, [column for column in df.columns if column.startswith('spatial_umap_density of ')]].copy()

    # Get a boolean series identifying the cells to be used for training the UMAP
    training_set_locations = df['umap_train']

    # Normalize to [0, 1] range, taking care of the scaler returning a numpy array instead of a dataframe
    print(f'[min, max] before normalization: [{df_density.min().min()}, {df_density.max().max()}]')
    df_density = pd.DataFrame(MinMaxScaler().fit_transform(df_density), columns=df_density.columns, index=df_density.index)
    print(f'[min, max] after normalization: [{df_density.min().min()}, {df_density.max().max()}]')

    # Get the pre-defined (in get_umap_train_and_test_sets()) training set
    df_density_train = df_density.loc[training_set_locations, :]

    # Perform UMAP on the training set
    training_set_shape = df_density_train.shape
    print(f'Fitting UMAP to the training set of shape {training_set_shape}...')
    start_time = time.time()
    umap_fit = umap.UMAP().fit(df_density_train)
    timing_string_fit = f'Time to fit UMAP to the training set of shape {training_set_shape}: {int(time.time() - start_time)} seconds'
    print(timing_string_fit)

    # If we want to transform just the densities needed for any of the test sets...
    if not transform_full_dataset:

        # Get a boolean series identifying the cells to be used for "testing" the UMAP
        umap_test_columns = [column for column in df.columns if column.startswith('umap_test_')]
        testing_set_locations = df.loc[:, umap_test_columns].any(axis='columns')

        # Get the pre-defined (in get_umap_train_and_test_sets()) *combined* testing set
        df_density_test = df_density.loc[testing_set_locations, :]

        # Transform the testing set using the fitted UMAP
        testing_set_shape = df_density_test.shape
        print(f'Transforming the *combined* testing set of shape {testing_set_shape} using the fitted UMAP...')
        start_time = time.time()
        df_umap = umap_transform_in_batches(umap_fit, df_density_test, nworkers=nworkers)
        timing_string_transform = f'Time to transform the *combined* testing set of shape {testing_set_shape} using the fitted UMAP: {int(time.time() - start_time)} seconds'
        print(timing_string_transform)

    # If we want to transform the entire dataset...
    else:

        # Transform all cells using the fitted UMAP
        full_dataset_shape = df_density.shape
        print(f'Transforming all cells of shape {full_dataset_shape} using the fitted UMAP...')
        start_time = time.time()
        df_umap = umap_transform_in_batches(umap_fit, df_density, nworkers=nworkers)
        timing_string_transform = f'Time to transform all cells of shape {full_dataset_shape} using the fitted UMAP: {int(time.time() - start_time)} seconds'
        print(timing_string_transform)

    # Return the UMAP results
    return df_umap, timing_string_fit, timing_string_transform
