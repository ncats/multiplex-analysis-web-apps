import scipy.spatial
import pandas as pd
import time
import os
import numpy as np
import utils


# def fast_neighbors_counts_for_block(df_all_images, image_column_name, image_name, coord_column_names, phenotypes, radii, phenotype_column_name):
def fast_neighbors_counts_for_block(df_image, image_name, coord_column_names, phenotypes, radii, phenotype_column_name):
    # A block can be an image, ROI, etc. It's the entity over which it makes sense to calculate the neighbors of centers. Here, we're assuming it's an image, but in the SIT for e.g., we generally want it to refer to a ROI.

    # Print the image name
    print(f'Calculating neighbor counts for image {image_name} ({len(df_image)} cells)...')

    start_time = time.time()

    # Get the boolean series identifying the current image specified by image_name
    # ser_curr_image = df_all_images[image_column_name] == image_name

    # Get the number of cells in the current image
    # num_cells_in_image = ser_curr_image.sum()
    # num_cells_in_image = len(df_image)

    # Construct the KDTree for the entire current image. This represents the centers
    # center_tree = scipy.spatial.KDTree(df_all_images.loc[ser_curr_image, coord_column_names])
    center_tree = scipy.spatial.KDTree(df_image[coord_column_names])

    # Initialize a list to hold the dataframes of neighbor counts for each radius (not each radius range)
    # df_counts_holder = [pd.DataFrame(0, index=phenotypes, columns=range(num_cells_in_image)) for _ in radii]
    df_counts_holder = [pd.DataFrame(0, index=phenotypes, columns=df_image.index) for _ in radii]

    # For each phenotype...
    for neighbor_phenotype in phenotypes:

        # Get the boolean series identifying the current neighbor phenotype
        # ser_curr_neighbor_phenotype = df_all_images[phenotype_column_name] == neighbor_phenotype
        ser_curr_neighbor_phenotype = df_image[phenotype_column_name] == neighbor_phenotype

        # Construct the KDTree for the current phenotype in the entire current image. This represents the neighbors
        # curr_neighbor_tree = scipy.spatial.KDTree(df_all_images.loc[ser_curr_image & ser_curr_neighbor_phenotype, coord_column_names])
        curr_neighbor_tree = scipy.spatial.KDTree(df_image.loc[ser_curr_neighbor_phenotype, coord_column_names])

        # For each radius, which should be monotonically increasing and start with 0...
        for iradius, radius in enumerate(radii):

            # Get the list of lists containing the indices of the neighbors for each center
            neighbors_for_radius = center_tree.query_ball_tree(curr_neighbor_tree, radius)

            # In the correct dataframe (corresponding to the current radius), set the counts of neighbors (of the current phenotype) for each center
            df_counts_holder[iradius].loc[neighbor_phenotype, :] = [len(neighbors_for_center) for neighbors_for_center in neighbors_for_radius]

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

    # Set the index of the final dataframe to correspond to that of the original dataframe
    # df_curr_counts.index = df_all_images[ser_curr_image].index
    # df_curr_counts.index = df_image.index

    df_curr_counts = df_curr_counts.astype(np.int32)

    print(f'  ...finished calculating neighbor counts for image {image_name} ({len(df_image)} cells) in {time.time() - start_time:.2f} seconds')

    # Return the final dataframe of neighbor counts for the current image
    return df_curr_counts


def main():

    # Parameters
    input_file = os.path.join('.', 'input', 'Combo_CSVfiles_20230327_152849.csv')
    radii = np.array([0, 25, 50, 100, 150, 200])
    coord_column_names = ['CentroidX', 'CentroidY']
    image_column_name = 'ShortName'
    phenotype_column_name = 'pheno_20230327_152849'
    method = 'kdtree'

    # Read in the datafile
    df = utils.downcast_dataframe_dtypes(pd.read_csv(input_file))

    # To see mapping of phenotype names
    # print(df.iloc[:, 83:92].drop_duplicates())

    # Variables
    image_names = sorted(df[image_column_name].unique())
    phenotypes = df[phenotype_column_name].value_counts().index

    tuple_holder = []

    for image_name in image_names:
        # print(fast_neighbors_counts_for_block(df[df[image_column_name] == image_name], image_name, coord_column_names, phenotypes, radii, phenotype_column_name).shape)  # 9.7, 15.0, 7.8, 8.1, 9.5, 8.8, 8.0, 11.6, 7.7, 8.5, 8.0
        tuple_holder.append((df[df[image_column_name] == image_name], image_name, coord_column_names, phenotypes, radii, phenotype_column_name))

    nworkers = 7

    df_counts_holder = utils.execute_data_parallelism_potentially(function=fast_neighbors_counts_for_block, list_of_tuple_arguments=tuple_holder, nworkers=nworkers, task_description='calculation of the counts matrix for all images', do_benchmarking=True, mp_start_method=None, use_starmap=True)

    df_counts = pd.concat(df_counts_holder, axis='index')

    return df_counts


if __name__ == '__main__':
    main()
