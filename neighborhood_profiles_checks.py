# Import relevant libraries
import numpy as np
import pandas as pd


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
