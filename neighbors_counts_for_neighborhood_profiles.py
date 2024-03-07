'''
Alternative method of performing the Cell density calculations
by Andrew
'''

# Import relevant libraries
import os
import time
import multiprocessing
import numpy as np
import pandas as pd
import utils

class dummySessionState:
    '''
    This is a simple class meant to mimic the SessionState class 
    from the streamlit library. It is used to store the state of the 
    app and its variables.
    '''
    def __init__(self):
        pass

    def __setitem__(self, key, value):
        setattr(self, key, value)

    def __getitem__(self, key):
        return getattr(self, key)

def calculate_density_matrix_for_all_images(struct, debug_output=False):
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
        range_strings (list): The list of range strings.
        debug_output (bool, optional): Whether to print debug output.
        num_cpus_to_use (int, optional): The number of CPUs to use. Defaults to 1.

    Returns:
        pandas.DataFrame: The dataframe containing the density matrix for all images.
    """

    df          = struct.cells
    phenotypes  = struct.species
    radii       = np.concatenate([[0], struct.dist_bin_px])

    num_cpus_to_use = int(multiprocessing.cpu_count() / 2)
    coord_column_names = ['Cell X Position', 'Cell Y Position']
    phenotype_column_name = 'Lineage'
    image_column_name     = 'Slide ID'
    image_names = df[image_column_name].unique()
    num_ranges = len(radii) - 1
    range_strings = [f'{radii[iradius]}, {radii[iradius + 1]})' for iradius in range(num_ranges)]

    # Initialize keyword arguments
    kwargs_list = []

    # Initialize the start time
    start_time = time.time()
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
                range_strings,
                debug_output
            )
        )

    # Get the number of CPUs to use
    print(f'Using {num_cpus_to_use} CPUs')

    # Create a pool of worker processes
    with multiprocessing.Pool(processes=num_cpus_to_use) as pool:

        # Apply the calculate_density_matrix_for_image function to each set of keyword arguments in kwargs_list
        # A single call would be something like: calculate_density_matrix_for_image(**kwargs_list[4])
        results = pool.starmap(calculate_density_matrix_for_image, kwargs_list)

    print(f'All images took {(time.time() - start_time) / 60:.2f} minutes to complete')

    df_density_matrix = pd.concat(results)
    full_array = None
    for ii, phenotype in enumerate(phenotypes):
        cols2Use = np.arange(0, num_ranges, 1) + (ii*(num_ranges))
        array_set = df_density_matrix.iloc[:, cols2Use].to_numpy()
        if full_array is None:
            full_array = array_set
        else:
            full_array = np.dstack((full_array, array_set))

    # Concatenate the results into a single dataframe
    return full_array

def calculate_density_matrix_for_image(df_image, phenotypes, phenotype_column_name, image, coord_column_names, radii, range_strings, debug_output=False):
    """
    Calculate the density matrix for a single image.

    Args:
        df_image (pandas.DataFrame): The dataframe containing the data for the current image.
        phenotypes (numpy.ndarray): The array of phenotypes.
        phenotype_column_name (str): The name of the column containing the phenotype information.
        image (str): The name of the current image.
        coord_column_names (list): The list of column names containing the coordinate information.
        radii (numpy.ndarray): The array of radii.
        range_strings (list): The list of range strings.
        debug_output (bool, optional): Whether to print debug output.

    Returns:
        pandas.DataFrame: The dataframe containing the density matrix for the current image.
    """

    # Initialize the start time
    start_time = time.time()

    # Get the number of range segments
    num_ranges = len(radii) - 1

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
                pass
                # print(f'No centers found for image {image} and phenotype {center_phenotype}')

        # Otherwise, calculate the number of neighbors of each type in the current image, for all neighbor phenotypes and all radii
        else:

            # Get the coordinates of the centers of the current type in the current image as a numpy array
            arr_image_center_phenotype = df_image[center_loc_for_image][coord_column_names].to_numpy()

            # Loop through the phenotypes as neighbor phenotypes
            for neighbor_phenotype in phenotypes:

                # Get the locations of the current image and neighbor phenotype in the dataframe
                neighbor_loc_for_image = df_image[phenotype_column_name] == neighbor_phenotype

                # Get the total number of neighbors of the current type in the current image
                num_neighbors_in_image = neighbor_loc_for_image.sum()

                # If there are no neighbors of the current type in the current image, print a message
                if num_neighbors_in_image == 0:
                    if debug_output:
                        pass
                        # print(f'No neighbors found for image {image} and phenotype {neighbor_phenotype}')

                # Otherwise, calculate the number of neighbors of the current type in the current image, for all radii
                else:

                    # Print the number of centers and neighbors found for the current image and phenotypes
                    if debug_output:
                        pass
                        # print(f'Number of centers found for image {image} and phenotype {center_phenotype}: {num_centers_in_image}')
                        # print(f'Number of neighbors found for image {image} and phenotype {neighbor_phenotype}: {num_neighbors_in_image}')

                    # Get the coordinates of the neighbors of the current type in the current image as a numpy array
                    arr_image_neighbor_phenotype = df_image[neighbor_loc_for_image][coord_column_names].to_numpy()

                    # Calculate the number of neighbors around the centers of the current types in the current image, in each radii range
                    nneighbors = utils.calculate_neighbor_counts_with_possible_chunking(center_coords = arr_image_center_phenotype,
                                                                                        neighbor_coords = arr_image_neighbor_phenotype,
                                                                                        radii = radii,
                                                                                        single_dist_mat_cutoff_in_mb = 200,
                                                                                        verbose = False,
                                                                                        test = False)  # (num_centers, num_ranges)

                    # Add the number of neighbors to the dataframe
                    for irange in range(num_ranges):
                        range_string = range_strings[irange]
                        # note that since we are adding columns dynamically that the order of these columns may not be logical because sometimes there are no centers or no neighbors
                        df_num_neighbors_image.loc[center_loc_for_image, f'Number of neighbors of type {neighbor_phenotype} in range {range_string}'] = nneighbors[:, irange]  

    # Print the time taken to calculate the number of neighbors for the current image
    if debug_output:
        print(f'Time to calculate neighbors for image {image} ({len(df_image)} rows) on a single CPU: {(time.time() - start_time) / 60:.2f} minutes')

    # Return the dataframe with the number of neighbors for the current image
    return df_num_neighbors_image

# Define the main function
def main():
    """
    This is a sample of how to calculate the density matrix for the entire dataset.
    """

    session_state = dummySessionState()

    # Constants
    num_cpus_to_use = int(multiprocessing.cpu_count() / 2)
    datafile = 'Combo_CSVfiles_20230327_152849.csv'
    input_file = os.path.join('.', 'input', datafile)

    # Read in the datafile
    df = pd.read_csv(input_file)
    radii = np.array([0, 25, 50, 100, 150, 200])

    image_column_name = 'ShortName'
    coord_column_names = ['CentroidX', 'CentroidY']
    phenotype_column_name = 'pheno_20230327_152849'

    # Variables
    image_names = df[image_column_name].unique()
    phenotypes = df[phenotype_column_name].unique()
    debug_output = True
    num_ranges = len(radii) - 1
    range_strings = [f'{radii[iradius]}, {radii[iradius + 1]})' for iradius in range(num_ranges)]

    # Calculate the density matrix for all images
    df_density_matrix = calculate_density_matrix_for_all_images(image_names, df, phenotypes, phenotype_column_name, image_column_name, coord_column_names, radii, range_strings, debug_output=debug_output, num_cpus_to_use=num_cpus_to_use)

    # Print shape of final density matrix dataframe
    print(f'Shape of final density matrix: {df_density_matrix.shape}')

    # Fill in any NaN values with 0 and convert to integers
    df_density_matrix = df_density_matrix.fillna(0).astype(int)

    df_density_matrix.to_csv('C:/Users/smithdaj/Desktop/AndrewMethod/Combo_CSVfiles_Out2.csv', index = False)

# Call the main function
if __name__ == '__main__':
    main()
