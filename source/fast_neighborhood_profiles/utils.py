import time
import numpy as np
import pandas as pd
import scipy.spatial


def fast_neighbors_counts_for_block2(df_image, image_name, coord_column_names, phenotypes, radii, phenotype_column_name, max_chunk_size_in_mb=200):
    # A block can be an image, ROI, etc. It's the entity over which it makes sense to calculate the neighbors of centers. Here, we're assuming it's an image, but in the SIT for e.g., we generally want it to refer to a ROI.

    # max_chunk_size_in_mb=200, for a 100K-cell dataset, will yield about 250-row chunks, which will yield about 400 chunks i.e. center KDTrees
    # max_chunk_size_in_mb=500, for a 100K-cell dataset, will yield about 650-row chunks, which will yield about 150 chunks i.e. center KDTrees
    # max_chunk_size_in_mb=1000, for a 100K-cell dataset, will yield about 1300-row chunks, which will yield about 80 chunks i.e. center KDTrees
    # max_chunk_size_in_mb=2000, for a 100K-cell dataset, will yield about 2600-row chunks, which will yield about 40 chunks i.e. center KDTrees
    # max_chunk_size_in_mb=5000, for a 100K-cell dataset, will yield about 6600-row chunks, which will yield about 15 chunks i.e. center KDTrees

    # Print the image name
    print(f'Calculating neighbor counts for image {image_name} ({len(df_image)} cells) using the new kdtree method...')

    # Record the start time
    start_time = time.time()

    # Store some properties of the input dataframe
    df_image_index = df_image.index
    num_cells = len(df_image)

    # Get the number of rows per chunk based on the maximum chunk size in MB and assuming every cell will be counted as a neighbor around every center (i.e., use upper limits)
    largest_sublist_size_in_mb = num_cells * 8  / 1024 ** 2  # 8 bytes per int64 --> units = mb / row
    num_rows_per_chunk = int(max_chunk_size_in_mb / largest_sublist_size_in_mb)

    # Get the corresponding integer indices to index things like the input image dataframe
    start_indices = np.arange(0, num_cells, num_rows_per_chunk)
    stop_indices = start_indices + num_rows_per_chunk

    # Initialize a list to hold the dataframes of neighbor counts for each radius (not each radius range)
    df_counts_holder = [pd.DataFrame(0, index=phenotypes, columns=df_image_index) for _ in radii]

    # Pre-calculate the neighbor tree for each phenotype
    neighbor_trees = []
    for neighbor_phenotype in phenotypes:

        # Get the boolean series identifying the current neighbor phenotype
        ser_curr_neighbor_phenotype = df_image[phenotype_column_name] == neighbor_phenotype

        # Construct the KDTree for the current phenotype in the entire current image. This represents the neighbors
        neighbor_trees.append(scipy.spatial.KDTree(df_image.loc[ser_curr_neighbor_phenotype, coord_column_names]))

    # For each chunk of centers...
    for start_index, stop_index in zip(start_indices, stop_indices):

        # Construct the KDTree for the current chunk. This represents the centers
        center_tree = scipy.spatial.KDTree(df_image.loc[df_image_index[start_index:stop_index], coord_column_names])

        # For each neighbor tree...
        for ineighbor_phenotype, curr_neighbor_tree in enumerate(neighbor_trees):

            # For each radius, which should be monotonically increasing and start with 0...
            for iradius, radius in enumerate(radii):

                # Get the list of lists containing the indices of the neighbors for each center
                neighbors_for_radius = center_tree.query_ball_tree(curr_neighbor_tree, radius)

                # In the correct dataframe (corresponding to the current radius), set the counts of neighbors (of the current phenotype) for each center
                df_counts_holder[iradius].iloc[ineighbor_phenotype, start_index:stop_index] = [len(neighbors_for_center) for neighbors_for_center in neighbors_for_radius]

    # For each annulus, i.e., each radius range...
    df_counts_holder_annulus = []
    for iradius in range(len(radii) - 1):

        # Get the counts of neighbors in the current annulus
        df_counts_curr_annulus = df_counts_holder[iradius + 1] - df_counts_holder[iradius]

        # Rename the index to reflect the current radius range
        radius_range_str = f'({radii[iradius]}, {radii[iradius + 1]}]'
        df_counts_curr_annulus.index = [f'{phenotype} in {radius_range_str}' for phenotype in phenotypes]

        # Add a transpose of this (so centers are in rows and phenotypes/radii are in columns) to the running list of annulus dataframes
        df_counts_holder_annulus.append(df_counts_curr_annulus.T)

    # Concatenate the annulus dataframes to get the final dataframe of neighbor counts for the current image
    df_curr_counts = pd.concat(df_counts_holder_annulus, axis='columns')

    # Convert the dataframe to a more efficient dtype (from int64)
    df_curr_counts = df_curr_counts.astype(np.int32)

    # Print the time taken to calculate the neighbor counts for the current image
    print(f'  ...finished calculating neighbor counts for image {image_name} ({len(df_image)} cells) in {time.time() - start_time:.2f} seconds')

    # Return the final dataframe of neighbor counts for the current image
    return df_curr_counts
