# Import relevant libraries
import pandas as pd
import dataset_formats
import numpy as np

def decompound_integer_field(df, integer_field_name, component_columns):
    """Modify a dataframe to decompound a single compound species.
    """

    # For every unique species in the dataset...
    for unique_int in df[integer_field_name].unique():

        # Log transform the species integer
        log2_transform = np.log2(unique_int)

        # If the log of the species integer is not an integer, then we're dealing with a compound species (i.e., more than one marker in the species)
        # I.e., if the current unique_int corresponds to a compound integer field...
        if log2_transform != np.round(log2_transform):

            # Get a copy of the slice of the original dataframe corresponding to the current compound species
            df_compound_group = df[df[integer_field_name] == unique_int].copy()

            # Get a series showing the makeup of the current compound species
            field_makeup = df_compound_group[component_columns].iloc[0, :]

            # Print this information
            print('{} {} corresponds to compound field:\n{}'.format(integer_field_name, unique_int, field_makeup))

            # Create a marker-zeroed-out version of the dataframe corresponding to the current group of compound species
            df_compound_group_zeroed = df_compound_group.copy()
            df_compound_group_zeroed.loc[:, component_columns] = 0

            # Initialize a running list of dataframes to concatenate with the main dataframe excluding the current compound species
            dfs_to_concat = [df[df[integer_field_name] != unique_int]]

            # For each positive marker in the current compound species...
            for positive_component in field_makeup[field_makeup == 1].index:

                # Create a copy of the zeroed dataframe
                df_compound_group = df_compound_group_zeroed.copy()

                # Turn on the bit corresponding to the current positive marker
                df_compound_group.loc[:, positive_component] = 1

                # Set the species integer for the group of the corresponding base-10 value
                df_compound_group.loc[:, integer_field_name] = int(''.join([str(x) for x in df_compound_group[component_columns].iloc[0, :]]), base=2)

                # Append the dataframe to the running list of dataframes to concatenate
                dfs_to_concat.append(df_compound_group)

            # Update the main dataframe as that decompounded for the current compound species
            df = pd.concat(dfs_to_concat)

    # Return the resulting dataframe
    return df

def map_species_to_possibly_compound_phenotypes(df, phenotype_identification_file, full_marker_list):

    # Import relevant library
    import ast

    # Load the phenotype assignments from the biologists
    df_phenotype_spec = pd.read_csv(phenotype_identification_file, sep='\t', header=None).iloc[:, -2:].rename({3: 'marker_list', 4: 'phenotypes'}, axis='columns')

    # From their possibly compounded phenotypes (where a compound phenotype consists of multiple phenotypes separated by a hyphen), determine the full list of possible phenotypes
    full_phenotype_list = []
    for spec_row in df_phenotype_spec.iterrows():
        curr_phenotypes = spec_row[1]['phenotypes']
        for phenotype in curr_phenotypes.split('-'):
            if phenotype not in full_phenotype_list:
                full_phenotype_list.append(phenotype)

    # Get prefixed column names for the individual phenotypes in order to avoid possible duplication of columns
    phenotype_colnames = ['phenotype ' + x for x in full_phenotype_list]

    # For each row in the biologists' phenotype specification file...
    num_updated_rows = 0
    for spec_row in df_phenotype_spec.iterrows():

        # Get a list of strings corresponding to the marker list and a string of the corresponding possibly compound phenotypes
        curr_marker_list = ast.literal_eval(spec_row[1]['marker_list'])
        curr_phenotypes = spec_row[1]['phenotypes']
        
        # Determine the species integer corresponding to the current marker list
        marker_str_list = ['0'] * len(full_marker_list)
        for marker in curr_marker_list:
            marker_str_list[full_marker_list.index(marker)] = '1'
        curr_species_int = int(''.join(marker_str_list), base=2)

        # Determine the phenotype integer corresponding to the current possibly compound phenotype
        phenotype_str_list = ['0'] * len(full_phenotype_list)
        for phenotype in curr_phenotypes.split('-'):
            phenotype_str_list[full_phenotype_list.index(phenotype)] = '1'
        curr_phenotype_int = int(''.join(phenotype_str_list), base=2)

        # Add individual component phenotype columns as well as the phenotype integer column to the dataframe
        curr_df_indexes = df[df['species_int'] == curr_species_int].index
        df.loc[curr_df_indexes, phenotype_colnames] = [int(x) for x in phenotype_str_list]  # new line
        df.loc[curr_df_indexes, 'phenotype_int'] = curr_phenotype_int
        num_updated_rows = num_updated_rows + len(curr_df_indexes)

    # Ensure the total number of rows modified equals the size of the dataframe itself
    assert num_updated_rows == len(df), 'ERROR: Not a one-to-one mapping of the rows'

    # Cast the new columns to integers since for some reason their default datatypes are floats
    cols_to_cast = phenotype_colnames + ['phenotype_int']
    for col_to_cast in cols_to_cast:
        df[col_to_cast] = df[col_to_cast].astype(int)

    # Return the final dataframe
    return df, phenotype_colnames

def apply_phenotyping(csv_file_path, method, phenotype_identification_file):
    """Load a datafile and apply one of three phenotyping methods: Species, Marker, or Custom.
    """

    # Apply some checks to the function parameters
    if phenotype_identification_file is not None:
        assert method == 'Custom', 'ERROR: The phenotype identification file is not None but it will not be used'
    if phenotype_identification_file is None:
        assert method in ['Species', 'Marker'], 'ERROR: The phenotype identification file is not specified but someone other than the Species or Marker phenotyping method (i.e., Custom) has been specified'

    # Determine the field seperator
    sep = (',' if csv_file_path.split('.')[-1] == 'csv' else '\t')

    # Read in the datafile
    df = pd.read_csv(csv_file_path, sep=sep)

    # From the detected datafile format, determine the coordinate columns and the marker information
    _, _, coord_cols, marker_prefix, _, markers_in_csv_file = dataset_formats.extract_datafile_metadata(csv_file_path)

    # Determine the marker columns from the datafile
    marker_cols = [marker_prefix + marker for marker in markers_in_csv_file]

    # Extract the relevant columns of the dataframe
    df = df[coord_cols + marker_cols]

    # If the marker columns are not 1s and 0s, map their values to 1s and 0s
    marker_cols_first_row = df[marker_cols].iloc[0, :].to_list()  # get just the first row of marker values
    if (0 not in marker_cols_first_row) and (1 not in marker_cols_first_row):
        df[marker_cols] = df[marker_cols].map(lambda x: ({'+': 1, '-': 0}[x[-1]] if isinstance(x, str) else 0))
        
    # Get the integer conversion of the markers lined up a base-two bits
    df['species_int'] = df[marker_cols].apply(lambda x: int(''.join([str(y) for y in list(x)]), base=2), axis='columns')

    # Drop objects that aren't positive for any markers of interest
    df = df[df['species_int'] != 0]

    # Print the initial species makeup of the dataframe prior to phenotyping
    value_counts = df['species_int'].value_counts()
    print('Dataframe makeup before {} phenotyping (species_int):'.format(method))
    print(value_counts)
    print('Total objects: {}'.format(value_counts.sum()))

    # If the user requests the Species phenotyping method...
    if method == 'Species':
        
        pass
        # # Set the best column names to those of the markers
        # phenotype_colnames = marker_cols
    
    # If the user requests the Marker phenotyping method...
    elif method == 'Marker':

        # Decompound the species
        df = decompound_integer_field(df, 'species_int', marker_cols)

        # # Set the best column names to those of the markers
        # phenotype_colnames = marker_cols
    
    # If the user requests the Custom phenotyping method...
    elif method == 'Custom':
        
        # Use the TSV file to assign possibly more than one phenotype to each species
        df, phenotype_colnames = map_species_to_possibly_compound_phenotypes(df, phenotype_identification_file, markers_in_csv_file)

        # Decompound the possibly compound phenotypes so there is a single phenotype per row, likely duplicating coordinates in the process
        df = decompound_integer_field(df, 'phenotype_int', phenotype_colnames)

    # Print the phenotype makeup of the dataframe after phenotyping
    if 'phenotype_int' in df.columns:
        col_to_count = 'phenotype_int'
    else:
        col_to_count = 'species_int'
    value_counts = df[col_to_count].value_counts()
    print('Dataframe makeup after {} phenotyping ({}):'.format(method, col_to_count))
    print(value_counts)
    print('Total objects: {}'.format(value_counts.sum()))

    # Return the final dataframe
    return df  #, phenotype_colnames
