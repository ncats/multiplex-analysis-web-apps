# Import relevant libraries
import numpy as np
import utils
import os
import pandas as pd
import multiprocessing
import time

def calculate_density_matrix_for_all_images(image_names, df, phenotypes, phenotype_column_name, image_column_name, coord_column_names, radii, num_ranges, range_strings, debug_output=False, num_cpus_to_use=1, swap_inequalities=False, cast_to_float32=False):
    """
    Calculate the density matrix for all images.

    Args:
        image_names (numpy.ndarray): The array of image names.
        df (pandas.DataFrame): The dataframe containing the data for all images.
        phenotypes (numpy.ndarray): The array of phenotypes.
        phenotype_column_name (str): The name of the column containing the phenotype information.
        image_column_name (str): The name of the column containing the image information.
        coord_column_names (list): The list of column names containing the coordinate information.
        radii (numpy.ndarray): The array of radii.
        num_ranges (int): The number of ranges.
        range_strings (list): The list of range strings.
        debug_output (bool, optional): Whether to print debug output.
        num_cpus_to_use (int, optional): The number of CPUs to use. Defaults to 1.

    Returns:
        pandas.DataFrame: The dataframe containing the density matrix for all images.
    """

    # Initialize a list of the keyword arguments
    kwargs_list = []

    # Loop through the images
    for image in image_names:

        # Create a dictionary for the variables
        kwargs_list.append(
            (
                df[df[image_column_name] == image][[phenotype_column_name] + coord_column_names].copy(),
                phenotypes,
                phenotype_column_name,
                image,
                coord_column_names,
                radii,
                num_ranges,
                range_strings,
                debug_output,
                swap_inequalities,
                cast_to_float32
            )
        )

    # Fan out the function to num_cpus_to_use CPUs
    results = utils.execute_data_parallelism_potentially(function=calculate_density_matrix_for_image, list_of_tuple_arguments=kwargs_list, nworkers=num_cpus_to_use, task_description='calculation of the counts matrix for neighborhood profiles checker', do_benchmarking=True, mp_start_method=None, use_starmap=True)

    # Concatenate the results into a single dataframe
    return pd.concat(results)

def calculate_density_matrix_for_image(df_image, phenotypes, phenotype_column_name, image, coord_column_names, radii, num_ranges, range_strings, debug_output=False, swap_inequalities=False, cast_to_float32=False):
    """
    Calculate the density matrix for a single image.

    Note that upon a cursory test, casting to float32 doesn't seem to speed it up at all but does produce (extremely minor) differences, so we should probably always keep it to False.

    Args:
        df_image (pandas.DataFrame): The dataframe containing the data for the current image.
        phenotypes (numpy.ndarray): The array of phenotypes.
        phenotype_column_name (str): The name of the column containing the phenotype information.
        image (str): The name of the current image.
        coord_column_names (list): The list of column names containing the coordinate information.
        radii (numpy.ndarray): The array of radii.
        num_ranges (int): The number of ranges.
        range_strings (list): The list of range strings.
        debug_output (bool, optional): Whether to print debug output.

    Returns:
        pandas.DataFrame: The dataframe containing the density matrix for the current image.
    """

    # Note that upon a cursory test, casting to float32 doesn't seem to speed it up at all but does produce (extremely minor) differences, so we should probably always keep it to False.
    print(f'The datatypes of the coordinate columns are: {df_image[coord_column_names].dtypes}')
    if cast_to_float32:
        df_image[coord_column_names] = df_image[coord_column_names].astype(np.float32)
        print(f'After casting, the datatypes of the coordinate columns are: {df_image[coord_column_names].dtypes}')

    # Initialize the start time
    start_time = time.time()

    # Initialize the dataframe to store the number of neighbors for the current image
    df_num_neighbors_image = pd.DataFrame(index=df_image.index)

    # Loop through the phenotypes as center phenotypes
    for center_phenotype in phenotypes:

        # Get the locations of the current image and center phenotype in the dataframe
        center_loc_for_image = df_image[phenotype_column_name] == center_phenotype

        # Get the total number of centers of the current type in the current image
        num_centers_in_image = center_loc_for_image.sum()

        # If there are no centers of the current type in the current image, print a message
        if num_centers_in_image == 0:
            if debug_output:
                print(f'No centers found for image {image} and phenotype {center_phenotype}')

        # Otherwise, calculate the number of neighbors of each type in the current image, for all neighbor phenotypes and all radii
        else:

            # Get the coordinates of the centers of the current type in the current image as a numpy array
            arr_image_center_phenotype = df_image.loc[center_loc_for_image, coord_column_names].to_numpy()

            # Loop through the phenotypes as neighbor phenotypes
            for neighbor_phenotype in phenotypes:

                # Get the locations of the current image and neighbor phenotype in the dataframe
                neighbor_loc_for_image = df_image[phenotype_column_name] == neighbor_phenotype

                # Get the total number of neighbors of the current type in the current image
                num_neighbors_in_image = neighbor_loc_for_image.sum()

                # If there are no neighbors of the current type in the current image, print a message
                if num_neighbors_in_image == 0:
                    if debug_output:
                        print(f'No neighbors found for image {image} and phenotype {neighbor_phenotype}')

                # Otherwise, calculate the number of neighbors of the current type in the current image, for all radii
                else:

                    # Print the number of centers and neighbors found for the current image and phenotypes
                    if debug_output:
                        print(f'Number of centers found for image {image} and phenotype {center_phenotype}: {num_centers_in_image}')
                        print(f'Number of neighbors found for image {image} and phenotype {neighbor_phenotype}: {num_neighbors_in_image}')

                    # Get the coordinates of the neighbors of the current type in the current image as a numpy array
                    arr_image_neighbor_phenotype = df_image.loc[neighbor_loc_for_image, coord_column_names].to_numpy()

                    # Calculate the number of neighbors around the centers of the current types in the current image, for all radii
                    nneighbors = utils.calculate_neighbor_counts_with_possible_chunking(center_coords=arr_image_center_phenotype, neighbor_coords=arr_image_neighbor_phenotype, radii=radii, single_dist_mat_cutoff_in_mb=200, verbose=False, test=False, swap_inequalities=swap_inequalities)  # (num_centers, num_ranges)

                    # Print the shape of the number of neighbors array
                    if debug_output:
                        print(nneighbors.shape)

                    # Add the number of neighbors to the dataframe
                    for iradius_range in range(num_ranges):
                        range_string = range_strings[iradius_range]
                        df_num_neighbors_image.loc[center_loc_for_image, f'Number of neighbors of type {neighbor_phenotype} in range {range_string}'] = nneighbors[:, iradius_range]  # note that since we are adding columns dynamically that the order of these columns may not be logical because sometimes there are no centers or no neighbors

    # Print the time taken to calculate the number of neighbors for the current image
    if debug_output:
        print(f'Time taken to calculate the number of neighbors for image {image} ({len(df_image)} rows) on a single CPU: {(time.time() - start_time) / 60:.2f} minutes')
    
    # Return the dataframe with the number of neighbors for the current image
    return df_num_neighbors_image

# Define the main function
def main():
    """
    This is a sample of how to calculate the density matrix for the entire dataset.
    """

    # Constants
    input_file = os.path.join('.', 'input', 'Combo_CSVfiles_20230327_152849.csv')
    radii_small_spacing = np.arange(0, 251, 25)
    radii_large_spacing = np.arange(0, 251, 50)
    coord_column_names = ['CentroidX', 'CentroidY']
    image_column_name = 'ShortName'
    phenotype_column_name = 'pheno_20230327_152849'

    # Read in the datafile
    df = pd.read_csv(input_file)

    # To see mapping of phenotype names
    # print(df.iloc[:, 83:92].drop_duplicates())

    # Variables
    image_names = df[image_column_name].unique()
    phenotypes = df[phenotype_column_name].unique()
    radii = radii_small_spacing
    # radii = radii_large_spacing
    debug_output = True
    num_ranges = len(radii) - 1
    # range_strings = ['[{}, {})'.format(radii[iradius], radii[iradius + 1]) for iradius in range(num_ranges)]
    range_strings = ['({}, {}]'.format(radii[iradius], radii[iradius + 1]) for iradius in range(num_ranges)]
    num_cpus_to_use = int(multiprocessing.cpu_count() / 2)

    # Calculate the density matrix for all images
    df_density_matrix = calculate_density_matrix_for_all_images(image_names, df, phenotypes, phenotype_column_name, image_column_name, coord_column_names, radii, num_ranges, range_strings, debug_output=debug_output, num_cpus_to_use=num_cpus_to_use, swap_inequalities=True)

    # Print the shape final density matrix dataframe, which can be concatenated with the original dataframe
    print(f'Shape of final density matrix: {df_density_matrix.shape}')

    # Fill in any NaN values with 0 and convert to integers
    df_density_matrix = df_density_matrix.fillna(0).astype(int)

    # To concatenate the density matrix with the original dataframe
    # pd.concat([df, df_density_matrix], axis='columns')

# Call the main function
if __name__ == '__main__':
    main()
