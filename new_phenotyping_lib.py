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

def map_species_to_possibly_compound_phenotypes(df, phenotype_identification_file, full_marker_list, species_int_colname='species_int'):

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

    # Determine all the species integers in the dataframe
    species_int_not_in_id_file = df[species_int_colname].unique()

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
        curr_df_indexes = df[df[species_int_colname] == curr_species_int].index
        df.loc[curr_df_indexes, phenotype_colnames] = [int(x) for x in phenotype_str_list]  # new line
        df.loc[curr_df_indexes, 'phenotype_int'] = curr_phenotype_int
        num_updated_rows = num_updated_rows + len(curr_df_indexes)

        # Remove the current species integer from the list of all species integers in the dataframe
        species_int_not_in_id_file = species_int_not_in_id_file[species_int_not_in_id_file != curr_species_int]

    # Filter out species with integer IDs that do not appear to be present in the phenotype identification file
    if len(species_int_not_in_id_file) > 0:
        print(f'Filtering out species with integer IDs {species_int_not_in_id_file} because they do not appear to be present in the phenotype iD file {phenotype_identification_file}...')
        num_rows_before = len(df)
        df = df[~df[species_int_colname].isin(species_int_not_in_id_file)].copy()
        num_rows_after = len(df)
        print(f'Filtered out {num_rows_before - num_rows_after} rows')

    # Ensure the total number of rows modified equals the size of the dataframe itself
    assert num_updated_rows == len(df), 'ERROR: Not a one-to-one mapping of the rows'

    # Cast the new columns to integers since for some reason their default datatypes are floats
    cols_to_cast = phenotype_colnames + ['phenotype_int']
    for col_to_cast in cols_to_cast:
        df[col_to_cast] = df[col_to_cast].astype(int)

    # Return the final dataframe
    return df, phenotype_colnames

def map_to_common_species_id_and_get_pheno_name_mapping(df):

    # Extract the final phenotype columns (note the lowercase "p")
    phenotype_columns = [column for column in df.columns if column.startswith('phenotype ')]

    # Check that every row corresponds to one and only one phenotype
    assert df[phenotype_columns].sum(axis='columns').unique().tolist() == [1], 'ERROR: Non-one-to-one mapping from row to phenotype (#1)'

    # Get the unique combinations of species and phenotype IDs
    df_duplicates_dropped = df[['Species int', 'phenotype_int']].drop_duplicates(ignore_index=False)

    # Since a phenotype may correspond to more than one species, group by phenotype and get the corresponding list of species
    srs_agg = df_duplicates_dropped.groupby(by='phenotype_int')['Species int'].agg(list)

    # Initialize mapping dictionaries
    species_int_to_species_int = {}
    species_int_to_pheno_int = {}
    pheno_int_to_pheno_name = dict(zip(sorted(df['phenotype_int'].unique()), [column.removeprefix('phenotype ') for column in phenotype_columns[-1::-1]]))
    species_int_to_pheno_name = {}
    pheno_int_to_pheno_name = get_entity_id_to_list_mapper(df, entity_colname='phenotype_int', entity_column_prefix='phenotype ')
    assert set([len(x) for x in pheno_int_to_pheno_name.values()]) == {1}, 'ERROR: Non-one-to-one mapping from row to phenotype (#2)'
    pheno_int_to_pheno_name = dict(zip(pheno_int_to_pheno_name.keys(), [x[0] for x in pheno_int_to_pheno_name.values()]))

    # For each unique phenotype...
    for phenotype_int, species_list in srs_agg.items():

        # If the phenotype corresponds to more than one species...
        if len(species_list) > 1:

            # Choose as the "common" species identifier the one with the lowest identifier
            common_species_int = min(species_list)

            # For every species for the current phenotype, map to the common species
            for species_int in species_list:
                species_int_to_species_int[species_int] = common_species_int

        # If the phenotype corresponds to only a single species...
        else:

            # Choose the one species as the "common" species
            common_species_int = species_list[0]

            # Perform the corresponding trivial mapping definition
            species_int_to_species_int[common_species_int] = common_species_int

        # Map the common species corresponding to the current phenotype to that phenotype
        species_int_to_pheno_int[common_species_int] = phenotype_int

        # Map the common species corresponding to the current phenotype to the corresponding phenotype name
        species_int_to_pheno_name[common_species_int] = pheno_int_to_pheno_name[phenotype_int]

    # Map the species IDs to the common species ID for each phenotype
    df['Species int'] = df['Species int'].apply(lambda x: species_int_to_species_int[x])

    # Return the species-mapped dataframe and the phenotype name corresponding to each common species ID
    return df, species_int_to_pheno_name

def get_entity_id_to_list_mapper(df, entity_colname='Species int', entity_column_prefix='Phenotype '):
    # "Entity" can e.g. refer to phenotypes or markers
    # Another sample call: get_entity_id_to_list_mapper(df, entity_colname='phenotype_int', entity_column_prefix='phenotype ')

    # Get a list of unique entity IDs
    unique_entity_ids = df[entity_colname].unique()

    # Get the full list of possible entity units
    full_entity_list = [column.removeprefix(entity_column_prefix) for column in df.columns if column.startswith(entity_column_prefix)]

    # Get the number of possible entity units
    num_entities = len(full_entity_list)

    # Initialize the ID to list mapper
    entity_id_to_list_mapper = {}

    # For every unique entity ID...
    for entity_int in unique_entity_ids:

        # Get a binary string corresponding to the current entity ID of the full zero-padded length
        curr_id_str = format(entity_int, '0' + str(num_entities) + 'b')

        # Determine the locations within the string where the bits are "on"
        on_loc = [ix for ix, x in enumerate(curr_id_str) if x == '1']

        # Get the unit corresponding to each bit that is on
        positive_entity_list = [full_entity_list[curr_on_loc] for curr_on_loc in on_loc]

        # Add to the mapper the positive entity list corresponding to the current unique ID
        entity_id_to_list_mapper[entity_int] = positive_entity_list

    # Return the mapper
    return entity_id_to_list_mapper

def apply_phenotyping(csv_file_path_or_df, method, phenotype_identification_file, species_int_colname='species_int', remove_allneg_phenotypes=True):
    """Load a datafile and apply one of three phenotyping methods: Species, Marker, or Custom.
    """

    import utils

    # Apply some checks to the function parameters
    if phenotype_identification_file is not None:
        assert method == 'Custom', 'ERROR: The phenotype identification file is not None but it will not be used'
    if phenotype_identification_file is None:
        assert method in ['Species', 'Marker'], 'ERROR: The phenotype identification file is not specified but something other than the Species or Marker phenotyping method (i.e., Custom) has been specified'

    # Print what we're doing
    print('Applying "{}" phenotyping method'.format(method))

    if type(csv_file_path_or_df) == str:

        # Determine the field seperator
        sep = (',' if csv_file_path_or_df.split('.')[-1] == 'csv' else '\t')

        # Read in the datafile
        df = utils.downcast_dataframe_dtypes(pd.read_csv(csv_file_path_or_df, sep=sep))

        # From the detected datafile format, determine the coordinate columns and the marker information
        _, _, coord_cols, marker_prefix, _, markers_in_csv_file = dataset_formats.extract_datafile_metadata(csv_file_path_or_df)

        # Determine the marker columns from the datafile
        marker_cols = [marker_prefix + marker for marker in markers_in_csv_file]

        # Extract the relevant columns of the dataframe
        df = df[coord_cols + marker_cols]

    # Assume it's already a dataset_formats.py-transformed dataframe
    else:
        df = csv_file_path_or_df
        marker_cols = [column for column in df.columns if column.startswith('Phenotype ')]
        markers_in_csv_file = [column.removeprefix('Phenotype ') for column in marker_cols]

    # If the marker columns are not 1s and 0s, map their values to 1s and 0s
    marker_cols_first_row = df[marker_cols].iloc[0, :].to_list()  # get just the first row of marker values
    if (0 not in marker_cols_first_row) and (1 not in marker_cols_first_row):
        df[marker_cols] = df[marker_cols].map(lambda x: {'+': 1, '-': 0}[x[-1]])

    # Get the integer conversion of the markers lined up a base-two bits. This is much faster than the original code (<1 second instead of minutes)
    powers_of_two = 2**np.arange(df[marker_cols].columns.size)[::-1]
    df[species_int_colname] = df[marker_cols].dot(powers_of_two)

    # Drop objects that aren't positive for any markers of interest
    if remove_allneg_phenotypes:
        df = df[df[species_int_colname] != 0].copy()  # I probably make a copy to avoid a SettingWithCopyWarning
    else:
        any_allneg_present = (df[species_int_colname] == 0).any()

    # Print the initial species makeup of the dataframe prior to phenotyping
    value_counts = df[species_int_colname].value_counts()
    print('Dataframe makeup before {} phenotyping (species_int):'.format(method))
    print(value_counts)
    print('Total objects: {}'.format(value_counts.sum()))

    # If the user requests the Species phenotyping method...
    if method == 'Species':
        
        # Get a dictionary taking the species ID to the phenotype name, which in this case is the combination of all positive markers
        species_int_to_pheno_name0 = get_entity_id_to_list_mapper(df, entity_colname=species_int_colname, entity_column_prefix='Phenotype ')
        species_int_to_pheno_name = {}
        for species_int, positive_marker_list in zip(species_int_to_pheno_name0.keys(), species_int_to_pheno_name0.values()):
            positive_marker_list[-1] = positive_marker_list[-1] + '+'
            species_int_to_pheno_name[species_int] = '+ '.join(positive_marker_list)
    
    # If the user requests the Marker phenotyping method...
    elif method == 'Marker':

        # Decompound the species
        df = decompound_integer_field(df, species_int_colname, marker_cols)

        # Get a dictionary taking the species ID to the phenotype name, which in this case is the corresponding marker
        species_int_to_pheno_name = get_entity_id_to_list_mapper(df, entity_colname=species_int_colname, entity_column_prefix='Phenotype ')
        assert set([len(x) for x in species_int_to_pheno_name.values()]) == {1}, 'ERROR: Non-one-to-one mapping from row to marker'
        species_int_to_pheno_name = dict(zip(species_int_to_pheno_name.keys(), [x[0] for x in species_int_to_pheno_name.values()]))
    
    # If the user requests the Custom phenotyping method...
    elif method == 'Custom':
        
        # Use the TSV file to assign possibly more than one phenotype to each species
        df, phenotype_colnames = map_species_to_possibly_compound_phenotypes(df, phenotype_identification_file, markers_in_csv_file, species_int_colname=species_int_colname)

        # Decompound the possibly compound phenotypes so there is a single phenotype per row, likely duplicating coordinates in the process
        df = decompound_integer_field(df, 'phenotype_int', phenotype_colnames)

        # Currently the species IDs are those corresponding to a species phenotyping method; change them to ones corresponding to the custom phenotype assignments file
        df, species_int_to_pheno_name = map_to_common_species_id_and_get_pheno_name_mapping(df)

    # If any all-negative phenotypes are present and we didn't want to remove them, add the mapper for these species
    if (not remove_allneg_phenotypes) and any_allneg_present:
        species_int_to_pheno_name[0] = 'All negative'

    # Print the phenotype makeup of the dataframe after phenotyping
    if 'phenotype_int' in df.columns:
        col_to_count = 'phenotype_int'
    else:
        col_to_count = species_int_colname
    value_counts = df[col_to_count].value_counts()
    print('Dataframe makeup after {} phenotyping ({}):'.format(method, col_to_count))
    print(value_counts)
    print('Total objects: {}'.format(value_counts.sum()))

    # Return the final dataframe
    return df, np.array(markers_in_csv_file), species_int_to_pheno_name
