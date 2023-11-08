def extract_datafile_metadata(datafile_path):

    # Import relevant library
    import pandas as pd

    # Extract the field seperator using the filename
    sep = (',' if datafile_path.split('.')[-1] == 'csv' else '\t')

    # Efficiently get just the column names in the datafiles
    columns_list = pd.read_csv(datafile_path, nrows=0, sep=sep).columns.to_list()

    # If the file is in the newest, "HALO" format...
    if 'Image Location' in columns_list:
        file_format = 'HALO'
        relevant_column_str = 'Image Location'
        first_dapi_positive_col = [x for x in columns_list if 'DAPI Positive' in x][0]
        index_start = columns_list.index('YMax') + 1
        index_end = columns_list.index(first_dapi_positive_col)
        marker_prefix = ''
        marker_cols = columns_list[index_start:index_end]
        coord_cols = ['XMin', 'XMax', 'YMin', 'YMax']
        string_processing_func = lambda x: x.split('.tif')[0].split('_')[-1]

    # If the file is in the original GMB format...
    elif 'tag' in columns_list:
        file_format = 'Native'
        # relevant_column_str = 'tag'
        relevant_column_str = 'Slide ID'
        marker_prefix = 'Phenotype '
        marker_cols = [x for x in columns_list if x.startswith(marker_prefix)]
        coord_cols = ['Cell X Position', 'Cell Y Position']
        string_processing_func = lambda x: '{:>3s}'.format(x.split('-')[0])

    # If the file is in the newer GMB format...
    elif 'tag.x' in columns_list:
        file_format = 'GMBSecondGeneration'
        # relevant_column_str = 'tag.x'
        relevant_column_str = 'Slide ID'
        marker_prefix = 'Phenotype '
        marker_cols = [x for x in columns_list if x.startswith(marker_prefix)]
        coord_cols = ['Cell X Position', 'Cell Y Position']
        # string_processing_func = lambda x: '{:>3s}'.format(x.split('-')[0])
        string_processing_func = lambda x: x.removesuffix('-MACRO-POS')

    # If the file is of unknown format...
    else:
        print('WARNING: Unregistered datafile format')
        file_format, relevant_column_str, marker_prefix, marker_cols, coord_cols, string_processing_func = None, None, None, None, None, None

    # Obtain the actual markers from the marker columns and the marker prefix
    if marker_cols is not None:
        markers = [x.removeprefix(marker_prefix) for x in marker_cols]
    else:
        markers = None

    # Return the the datafile metadata
    return relevant_column_str, string_processing_func, coord_cols, marker_prefix, file_format, markers

class Native:
    """Class representing format of original Consolidated_data.txt file that Houssein sent us around 2017-2018

    The conversion methods (i.e., "adhere_to_XXXX_format()") are blank for this format as this is the format used natively in the time_cell_interaction_lib.py module.

    Sample instantiation:

        import dataset_format_conversion
        dataset_obj = dataset_format_conversion.Native(input_datafile='../data/Consolidated_data-original.txt', coord_units_in_microns=0.5, species_equivalents={36: 32, 20: 16, 10: 8, 33: 1, 18: 16, 17: 1, 9: 1})
        dataset_obj.process_dataset()

    I should dig up the email/conversation/notes in which Houssein told us the units of the original input datafile were half-microns.
    """

    def __init__(self, input_datafile, coord_units_in_microns, images_to_analyze, sep='\t', min_coord_spacing=None, species_equivalents={}, mapping_dict={}, roi_width=None, overlap=0, phenotype_identification_tsv_file=None):
        """Object initialization, which just consists of reading in the data from the datafile

        Args:
            input_datafile (str): Full pathname to the input datafile.
            coord_units_in_microns (float): Factor by which to multiply the coordinates in the input datafile in order to achieve units of microns.
            sep (str, optional): Separator in the input datafile. Defaults to '\t'.
            min_coord_spacing (float, optional): Minimum coordinate spacing after conversion to microns. If set to None, this value is automatically calculated when process_dataset() is called. Defaults to None.
            species_equivalents (dict, optional): Define equivalent species, which are identified by integer numbers, but sometimes differently-marked species should actually be the same species; specify that here. Unlike mapping_dict (a plotting input in the calling Jupyter notebook), this DOES affect the actual determination of the species and therefore the calculated data such as the P values. Species can be accessed via slices.plotting_map. Defaults to {}.
            mapping_dict (dict, optional): Dictionary for mapping markers to known species names in order to make the plots more clear. This does not affect the actual determination of the species; this is just labeling. Defaults to {}.
            roi_width (float, optional): Desired ROI width in microns for patching up the input slide. If None, no patching is performed. Defaults to None.
            overlap (float, optional): If patching is requested (roi_width!=None), amount in microns you want consecutive ROIs to overlap. Note that for the SIP code, in order to utilize all data most efficiently, this should equal twice the density radius (e.g., the "thickness" value in main.ipynb). Defaults to 0.
            phenotype_identification_tsv_file (str, optional): Name of the TSV file generated by pasting into a text editor columns from Excel in the order specified in the read_csv() method below. Defaults to None.

        Returns:
            Spatial analysis dataset object.
        """
        self.input_datafile = input_datafile
        self.sep = sep
        self.images_to_analyze = images_to_analyze
        self.data = self.read_datafile()
        self.coord_units_in_microns = coord_units_in_microns
        self.min_coord_spacing_ = min_coord_spacing
        self.species_equivalents = species_equivalents
        self.mapping_dict = mapping_dict
        self.roi_width = roi_width
        self.overlap = overlap
        self.phenotype_identification_tsv_file = phenotype_identification_tsv_file

    def read_datafile(self):
        """Read in the text-formatted datafile

        Returns:
            Pandas dataframe: Pandas dataframe of the imported text datafile.
        """

        # Variable definitions from attributes
        input_datafile = self.input_datafile
        sep = self.sep
        images_to_analyze = self.images_to_analyze

        # Import relevant libraries
        import os
        import pandas as pd

        # Import the text file using Pandas
        if os.path.exists(input_datafile):
            df = pd.read_csv(input_datafile, sep=sep)
            relevant_column_str, string_processing_func, _, _, _, _ = extract_datafile_metadata(input_datafile)
            df = df[df[relevant_column_str].apply(lambda x: string_processing_func(x).strip() in images_to_analyze)]
            print('Text file "{}" with separator "{}" has been successfully read and filtered to images {}'.format(input_datafile, sep, images_to_analyze))
        else:
            print('ERROR: Text file "{}" does not exist'.format(input_datafile))

        # Return the resulting Pandas dataframe
        return df

    def write_reformatted_datafile(self, pathname):
        """Write a tab-separated datafile

        Args:
            pathname (str): Full path to the datafile you want to generate.
        """

        # Variable definitions from attributes
        df = self.data

        # Import relevant libraries
        import os

        # Write the text file using Pandas if the desired pathname doesn't already exist
        if not os.path.exists(pathname):
            sep = '\t'
            curr_dirname = os.path.dirname(pathname)
            if not os.path.exists(curr_dirname):
                os.makedirs(curr_dirname)
            df.to_csv(path_or_buf=pathname, sep=sep)
            print('Text file "{}" with separator "{}" has been successfully written'.format(pathname, sep))
        else:
            print('WARNING: Text file "{}" has not been written because that pathname already exists'.format(pathname))

    def adhere_to_slide_id_format(self):
        """Ensure the "Slide ID" column of the data conforms to the required format
        """
        pass

    def adhere_to_tag_format(self):
        """Ensure the "tag" column of the data conforms to the required format
        """
        pass

    def adhere_to_cell_position_format(self):
        """Ensure the "Cell X Position" and "Cell Y Position" columns of the data conform to the required format
        """

        # Variable definitions from attributes
        srs_x_orig = self.data['Cell X Position']
        srs_y_orig = self.data['Cell Y Position']
        coord_units_in_microns = self.coord_units_in_microns

        # Convert coordinates to microns by multiplying by coord_units_in_microns
        srs_x_new = srs_x_orig * coord_units_in_microns
        srs_y_new = srs_y_orig * coord_units_in_microns

        # Attribute assignments from variables
        self.data['Cell X Position'] = srs_x_new
        self.data['Cell Y Position'] = srs_y_new
        self.coord_units_in_microns = 1  # implying 1 micron per coordinate unit

    def adhere_to_phenotype_format(self):
        """Ensure the "Phenotype XXXX" columns of the data conform to the required format
        """
        pass

    def trim_dataframe(self):
        """Replace the Pandas dataframe with a trimmed version containing just the necessary columns
        """

        # Variable definitions from attributes
        df = self.data

        # Extract just the columns to keep in the trimmed dataframe
        cols_to_keep = ['Slide ID', 'tag', 'Cell X Position', 'Cell Y Position'] + df.loc[0, :].filter(regex='^Phenotype ').index.tolist()

        # Attribute assignments from variables
        self.data = df.loc[:, cols_to_keep]

    def calculate_minimum_coordinate_spacing(self):
        """Calculate the minimum spacing (aside from zero) in the data coordinates after conversion to microns
        """

        # Variable definitions from attributes
        df = self.data
        min_coord_spacing = self.min_coord_spacing_

        # If min_coord_spacing hasn't been manually specified...
        if min_coord_spacing is None:

            # Import relevant libraries
            import numpy as np

            # Get the sorted x- and y-coordinates
            sorted_x = df['Cell X Position'].sort_values()
            sorted_y = df['Cell Y Position'].sort_values()

            # Get just the edges of the finite difference arrays, reset the indexes, take the differences, and extract the minimum value aside from zero
            spacing_x = (sorted_x[1:].reset_index(drop=True) - sorted_x[:-1].reset_index(drop=True)).sort_values().unique()[1]
            spacing_y = (sorted_y[1:].reset_index(drop=True) - sorted_y[:-1].reset_index(drop=True)).sort_values().unique()[1]

            # If the minimum coordinate spacings for the x- and y-coordinates don't agree, say so and set the global minimum to the smaller of the two
            if spacing_x != spacing_y:
                min_coord_spacing = np.min([spacing_x, spacing_y])
                print('WARNING: The x-coordinate spacing ({}) is different from the y-coordinate spacing ({}); using the smaller of the two: {}'.format(spacing_x, spacing_y, min_coord_spacing))

            # Otherwise, set the global minimum to the shared minimum value
            else:
                min_coord_spacing = spacing_x
                print('Calculated minimum coordinate spacing: {}'.format(min_coord_spacing))

        # Output to screen a snippet of the coordinates so the user can visually check the determined minimum coordinate spacing
        print(df.loc[:, ['Cell X Position', 'Cell Y Position']])
        print('The minimum coordinate spacing has been set to {}. If this does not look right from the data (e.g., from the sample of the data printed above), please manually set this value using the "min_coord_spacing" parameter in the call to time_cell_interaction_lib.preprocess_dataset()'.format(min_coord_spacing))

        # Attribute assignments from variables
        self.min_coord_spacing_ = min_coord_spacing

    def calculate_minimum_coordinate_spacing_per_roi(self):
        """Calculate the minimum spacing (aside from zero) in the data coordinates after conversion to microns

        Do this for each ROI (i.e., "tag" field), and choose the smallest over all ROIs.
        """

        # Variable definitions from attributes
        df = self.data
        min_coord_spacing = self.min_coord_spacing_

        # If min_coord_spacing hasn't been manually specified...
        if min_coord_spacing is None:

            # Set the global minimum coordinate spacing to an unreasonably high number
            min_coord_spacing = 4444

            # For each ROI in the overall dataframe...
            for roi in df['tag'].unique():

                # Get the current slice of the overall dataframe
                curr_df = df[df['tag'] == roi]

                # Calculate the minimum coordinate spacing in this slice of the overall dataframe
                min_coord_spacing_in_roi = get_min_coord_spacing_in_dataframe(curr_df)

                # Output the result
                print('Minimum spacing in ROI {}: {} microns'.format(roi, min_coord_spacing_in_roi))

                # Update the global minimum coordinate spacing if it's the smallest so far
                if min_coord_spacing_in_roi < min_coord_spacing:
                    min_coord_spacing = min_coord_spacing_in_roi

            # Round the calculated minimum coordinate spacing to 8 decimal places
            # round(min_spacing_holder.min(), 8)
            min_coord_spacing = round(min_coord_spacing, 8)

        # Output to screen a snippet of the coordinates so the user can visually check the determined minimum coordinate spacing
        print(df.loc[:, ['Cell X Position', 'Cell Y Position']])
        print('The minimum coordinate spacing has been set to {} by finding the minimum of the minimum coordinate spacings for each ROI (which can be significantly larger than the value when done for each slide). If this does not look right from the data (e.g., from the sample of the data printed above, but be careful; likely just use the automatically generated value here in a reasonable format), please manually set this value using the "min_coord_spacing" parameter in the configuration .yml file'.format(min_coord_spacing))

        # Attribute assignments from variables
        self.min_coord_spacing_ = min_coord_spacing

    def extra_processing(self):
        """Perform any extra processing
        """
        pass

    def process_dataset(self, write_new_datafile=True, new_datafile_suffix='-converted'):
        """Convert dataset to the format required for the SIP library

        Args:
            write_new_datafile (bool, optional): Whether to write the reformatted datafile to disk using new_datafile_suffix. Defaults to True.
            new_datafile_suffix (str, optional): Suffix to add to the original pathname as the name of the converted datafile. Defaults to '-converted'.
        """

        # Import relevant library
        import os

        # Variable definitions from attributes
        input_datafile = self.input_datafile

        # Apply all four sets of column reformattings
        self.adhere_to_slide_id_format()
        self.adhere_to_tag_format()
        self.adhere_to_cell_position_format()
        self.adhere_to_phenotype_format()

        # Retain just the necessary columns of data
        self.trim_dataframe()

        # Calculate and output the minimum coordinate spacing
        self.calculate_minimum_coordinate_spacing_per_roi()

        # Perform any additional processing for the dataset
        self.extra_processing()

        # Write the new dataframe to disk if desired
        if write_new_datafile:
            # output_datafile = add_suffix_to_pathname(input_datafile, new_datafile_suffix)
            extension = os.path.basename(input_datafile).split('.')[1]
            output_datafile = os.path.join('.', 'output', os.path.basename(input_datafile).split('.')[0] + new_datafile_suffix + '.' + extension)
            self.write_reformatted_datafile(output_datafile)

    def extract_original_parameters(self):
        """Extract attributes, some not actually needed, for compatibility with original way the SIP library worked before there was a dataset object

        Returns:
            tuple: coord_units_in_microns, min_coord_spacing, input_data_filename, species_equivalents, mapping_dict
        """

        # Variable definitions from attributes
        coord_units_in_microns = self.coord_units_in_microns
        min_coord_spacing = self.min_coord_spacing_
        input_data_filename = self.input_datafile
        species_equivalents = self.species_equivalents
        mapping_dict = self.mapping_dict
        phenotype_identification_tsv_file = self.phenotype_identification_tsv_file

        return coord_units_in_microns, min_coord_spacing, input_data_filename, species_equivalents, mapping_dict, phenotype_identification_tsv_file


class GMBSecondGeneration(Native):
    """Class representing format of the data the GMB group is using around Fall/Winter 2021

    Sample instantiation:

        import dataset_formats
        dataset_obj = dataset_formats.GMBSecondGeneration(input_datafile='../data/Consolidated_data-shanias_example_2.txt', coord_units_in_microns=0.250, min_coord_spacing=0.025)
        dataset_obj.process_dataset()

    Note that the units of coord_units_in_microns are um/pixel per the screenshare of Inform? Shania showed me on 12/2/21 and the units of the coordinates in input_datafile are pixels per an email Shania sent me on 11/30/21.
    """

    def adhere_to_slide_id_format(self):
        """Ensure the "Slide ID" column of the data conforms to the required format
        """

        # Variable definitions from attributes
        srs_slide_id_orig = self.data['Slide ID']

        # Import relevant libraries
        import string

        # Get a list of uppercase letters
        all_upper_letters = list(string.ascii_uppercase)

        # Determine the unique cases (i.e., people)
        unique_cases = srs_slide_id_orig.apply(lambda x: x.split('-')[0]).unique()
        print('Number of unique cases: {}'.format(len(unique_cases)))

        # Initialize some variables
        slide_id_keys = []
        slide_id_values = []
        case_id = 1

        # Loop over each unique case
        for unique_case in unique_cases:

            # Determine the unique samples (i.e., conditions, slides) for the current unique case
            unique_samples = srs_slide_id_orig.loc[srs_slide_id_orig.str.match(unique_case)].unique()
            num_unique_samples = len(unique_samples)
            print('Number of unique samples for case {} (case ID: {}): {}'.format(unique_case, case_id, num_unique_samples))

            # If there are more than 26 unique samples for the current case, we're out of letters, so throw an error and exit
            if len(unique_samples) > 26:
                print('ERROR: There are more than 26 samples for case {}; this is limited by the number of capital letters in the alphabet'.format(unique_case))
                exit()

            # Determine the array of sample IDs for the current case
            sample_ids = all_upper_letters[:num_unique_samples]

            # Loop over each unique sample
            for isample, unique_sample in enumerate(unique_samples):

                # Get the current sample ID
                sample_id = sample_ids[isample]

                # Determine the correctly formatted slide ID
                new_slide_id = str(case_id) + sample_id + '-' + unique_sample

                # Output what has been done
                print('  Sample {} has sample ID {}; new slide ID is {}'.format(unique_sample, sample_id, new_slide_id))

                # Append to lists that will define the dictionary that maps the old slide IDs to the new ones
                slide_id_keys.append(unique_sample)
                slide_id_values.append(new_slide_id)

            # Increment the case ID
            case_id = case_id + 1

        # Assmple the dictionary for mapping the old sample IDs to the new ones and assign the conversion to a variable
        slide_id_mapping = dict(zip(slide_id_keys, slide_id_values))
        srs_slide_id_new = srs_slide_id_orig.map(slide_id_mapping)

        # Attribute assignments from variables
        self.data['Slide ID'] = srs_slide_id_new

    def adhere_to_tag_format(self):
        """Ensure the "tag" column of the data conforms to the required format
        """

        # Variable definitions from attributes
        df = self.data

        # Combine two columns to create the ROI ID, i.e., "tag"
        srs_tag = df.apply(lambda x: x['Slide ID'] + '_' + x['tag.x'].split('_')[-1], axis='columns')

        # Attribute assignments from variables
        self.data['tag'] = srs_tag


class OMAL(Native):
    """Class representing format of the data the OMAL group is using around Winter/Spring 2022

    Sample instantiation:

        import dataset_formats
        dataset_obj = dataset_formats.OMAL(input_datafile='../../data/all_five_images_NOS2_and_COX2_and_CD8_Fused_ALW_05212021_entire_image.csv', coord_units_in_microns=0.325, sep=',')
        dataset_obj.process_dataset()

    Note that the units of coord_units_in_microns are um/pixel and the units of the coordinates in input_datafile are pixels.
    """

    def adhere_to_slide_id_format(self):
        """Ensure the "Slide ID" column of the data conforms to the required format
        """

        # Variable definitions from attributes
        image_location = self.data['Image Location']

        # Determine the image numbers
        # srs_imagenum = image_location.apply(lambda x: x.split('_')[-1].rstrip('.tif'))
        srs_imagenum = image_location.apply(lambda x: x.split('_')[-1].split('.')[0])  # make it not have to be a .tif extension but rather any extension

        # Get the unique image numbers
        unique_images = srs_imagenum.unique()

        # Calculate a dictionary mapper that maps the unique images to integers
        mapper = dict(zip(unique_images, [x + 1 for x in range(len(unique_images))]))

        # Attribute assignments
        self.data['Slide ID'] = srs_imagenum.apply(lambda x: '{}A-imagenum_{}'.format(mapper[x], x))

    def adhere_to_tag_format(self):
        """Ensure the "tag" column of the data conforms to the required format
        """

        # Variable definitions from attributes
        df = self.data
        coord_units_in_microns = self.coord_units_in_microns
        roi_width = self.roi_width
        overlap = self.overlap

        # Get the "Slide ID" column of the dataframe
        slide_id = df['Slide ID']

        # Determine the unique slides in the dataset
        unique_slides = slide_id.unique()

        # For each slide...
        for unique_slide in unique_slides:

            # Get the locations in the overall dataframe of the current slide
            curr_loc = slide_id == unique_slide

            # Get the current dataframe slice
            curr_df = df[curr_loc].copy()

            # Create the "tag" (i.e., ROI, i.e., patch) column
            # Note the coordinates used below are in pixels
            if roi_width is None:  # i.e., if we don't want to do any patching
                print('No patching will be done')
                df.loc[curr_loc, 'tag'] = curr_df['Slide ID'] + '_roi_' + '[{},{}]'.format(int((curr_df['XMin'].min() + curr_df['XMax'].max()) / 2 + 0.5), int((curr_df['YMin'].min() + curr_df['YMax'].max()) / 2 + 0.5))  # populate the "tag" column with the slide ID appended with the minimum x- and y-coordinates
            else:
                print('Patching will be performed...')
                df.loc[curr_loc, 'tag'] = break_up_slide_into_patches(curr_df, roi_width=(roi_width / coord_units_in_microns), overlap=(overlap / coord_units_in_microns))

        # Print out the current dataframe size
        print('Current dataframe length: {}'.format(len(df)))

        # If patching was performed, then the 'tag' column will be lists of strings, not strings
        if roi_width is not None:

            # Print out the expected length of the new dataframe
            print('Expected length of duplicated dataframe: {}'.format(df['tag'].apply(lambda x: len(x)).sum()))

            # Repeat objects that are in multiple ROIs
            df = duplicate_rows(df, col_to_repeat='tag')

            # Print the length of the new dataframe
            print('Length of duplicated dataframe: {}'.format(len(df)))

        # Sort the data by slide and then by ROI
        df = df.sort_values(by=['Slide ID', 'tag']).reset_index(drop=True)  # adding .reset_index(drop=True) to just *make sure* the SIP code in general doesn't assume the data index is already in sorted order

        # Overwrite the original dataframe with the one having the appended ROI column (I'd imagine this line is unnecessary)
        self.data = df

    def adhere_to_cell_position_format(self):
        """Ensure the "Cell X Position" and "Cell Y Position" columns of the data conform to the required format
        """

        # Variable definitions from attributes
        df = self.data
        coord_units_in_microns = self.coord_units_in_microns

        # Establish the x- and y-coordinate columns from their gates
        srs_x_orig = (df['XMin'] + df['XMax']) / 2
        srs_y_orig = (df['YMin'] + df['YMax']) / 2

        # Convert coordinates to microns by multiplying by coord_units_in_microns, creating the new columns for the x- and y-coordinates
        df['Cell X Position'] = srs_x_orig * coord_units_in_microns
        df['Cell Y Position'] = srs_y_orig * coord_units_in_microns

        # Attribute assignments from variables
        self.data = df
        self.coord_units_in_microns = 1  # implying 1 micron per coordinate unit (since now the coordinates are in microns)

    def adhere_to_phenotype_format(self):
        """Ensure the "Phenotype XXXX" columns of the data conform to the required format
        """

        # Variable definitions from attributes
        df = self.data
        input_datafile = self.input_datafile

        # Define the phenotypes of interest in the input dataset
        _, _, _, _, _, phenotype_columns = extract_datafile_metadata(input_datafile)

        # Rename the phenotype columns so that they are prepended with "Phenotype "
        df = df.rename(dict(zip(phenotype_columns, ['Phenotype {}'.format(x) for x in phenotype_columns])), axis='columns')

        # For each phenotype column, convert zeros and ones to -'s and +'s
        for col in df.filter(regex='^Phenotype\ '):
            df[col] = df[col].map({0: '-', 1: '+'})

        # Attribute assignments from variables
        self.data = df

    def extra_processing(self):
        """Delete any ROIs that contain only a single spatial coordinate, regardless of whether the objects are compound species or if there are simply different species located as the same coordinate

        In hindsight I should have also eliminated 1D ROIs instead of only 0D ROIs as I do here. I do this later though in time_cell_interaction_lib.py.
        """

        # Variable definitions from attributes
        df = self.data

        # Obtain just the phenotype columns and get the number of phenotypes
        df_phenotypes = df.filter(regex='^Phenotype\ ')
        num_phenotypes = df_phenotypes.shape[1]

        # Delete empty objects
        df_pared = df[df_phenotypes.apply(lambda x: ''.join(x), axis='columns') != '-' * num_phenotypes]

        # Get a Series containing the number of unique coordinates in every ROI, in which we first group by ROI and coordinates, and then group by just ROI, counting the number of unique coordinates within
        num_unique_coords_per_roi = df_pared.groupby(by=['tag', 'Cell X Position', 'Cell Y Position']).size().index.to_frame(index=False).groupby(by='tag').count().iloc[:, 0]

        # Get the set of ROIs containing a single unique coordinate
        rois_with_single_unique_coord = set(num_unique_coords_per_roi[num_unique_coords_per_roi == 1].index.to_list())

        # Drop these ROIs from the dataset
        print('Dropping {} ROIs with valid objects that have only a single unique spatial coordinate'.format(len(rois_with_single_unique_coord)))
        df = df.loc[df['tag'].apply(lambda x: x not in rois_with_single_unique_coord), :]

        # Add .reset_index(drop=True) to just *make sure* the SIP code in general doesn't assume the data index is already in sorted order
        df = df.reset_index(drop=True)

        # Attribute assignments from variables
        self.data = df


def add_suffix_to_pathname(input_datafile, new_datafile_suffix):
    """Add a suffix to the basename of a pathname, right before the file extension

    Note I could have just used os.path.splitext(input_datafile).

    Args:
        input_datafile (str): Pathname to which you want to add a suffix prior to the file extension
        new_datafile_suffix (str): Suffix you want to add to the inputted pathname

    Returns:
        str: Input pathname but with a suffix added prior to the file extension
    """

    # Get the full input pathname in reverse order character-by-character
    rev_pathname = input_datafile[::-1]

    # Determine the index location of the first "." in the reversed pathname
    dot_loc = rev_pathname.find('.')

    # From this location, determine the file extension and basename, in normal order
    extension = rev_pathname[:dot_loc][::-1]
    rest = rev_pathname[dot_loc + 1:][::-1]

    # Return the original pathname but with the suffix added prior to the file extension
    return rest + new_datafile_suffix + '.' + extension


def get_min_coord_spacing_in_dataframe(df, print_output=True):
    """Calculate the minimum coordinate spacing in a properly formatted Pandas dataframe

    Args:
        df (dataframe): Pandas dataframe containing "Cell X Position" and "Cell Y Position" fields
        print_output (bool, optional): Whether to print the calculated minimum coordinate spacing, including whether the x- and y-values are different. Defaults to True.

    Returns:
        float: Calculated minimum coordinate spacing
    """

    # Import relevant libraries
    import numpy as np

    def get_smallest_spacing_aside_from_zero(arr):
        # This should exit cleanly now when arr = [0], which should probably never happen though it might... yes I think it does because I don't think I've dropped the "trivial" ROIs yet
        tol = 1e-8
        large_num = 1e8
        arr = [[x for x in arr if x < tol][0]] + [x for x in arr if x > tol]  # account for their possibly being multiple essentially-zero spacings; just keep the first of these close-to-zeros
        nonzeros_in_smallest_spacings = [x for x in arr[:2] if x > tol]
        if len(nonzeros_in_smallest_spacings) >= 1:
            return nonzeros_in_smallest_spacings[0]
        else:
            # import sys
            print('WARNING: This ROI is likely "trivial." Here are its first ten smallest spacings: {}. Returning {}.'.format(arr[:10], large_num))
            # sys.exit()
            return large_num

    # Get the sorted x- and y-coordinates
    sorted_x = df['Cell X Position'].sort_values()
    sorted_y = df['Cell Y Position'].sort_values()

    # Get just the edges of the finite difference arrays, reset the indexes, take the differences, and extract the minimum value aside from zero
    # spacing_x = (sorted_x[1:].reset_index(drop=True) - sorted_x[:-1].reset_index(drop=True)).sort_values().unique()[1]
    # spacing_y = (sorted_y[1:].reset_index(drop=True) - sorted_y[:-1].reset_index(drop=True)).sort_values().unique()[1]
    spacing_x = get_smallest_spacing_aside_from_zero((sorted_x[1:].reset_index(drop=True) - sorted_x[:-1].reset_index(drop=True)).sort_values().unique())
    spacing_y = get_smallest_spacing_aside_from_zero((sorted_y[1:].reset_index(drop=True) - sorted_y[:-1].reset_index(drop=True)).sort_values().unique())

    # If the minimum coordinate spacings for the x- and y-coordinates don't agree, say so and set the global minimum to the smaller of the two
    if spacing_x != spacing_y:
        min_coord_spacing = np.min([spacing_x, spacing_y])
        if print_output:
            print('WARNING: The x-coordinate spacing ({}) is different from the y-coordinate spacing ({}); using the smaller of the two: {}'.format(spacing_x, spacing_y, min_coord_spacing))

    # Otherwise, set the global minimum to the shared minimum value
    else:
        min_coord_spacing = spacing_x
        if print_output:
            print('Calculated minimum coordinate spacing: {}'.format(min_coord_spacing))

    # Return the minimum coordinate spacing in the dataframe
    return min_coord_spacing


def calculate_roi_coords(min_coord, max_coord, roi_width, overlap=0):
    """Break up the coordinates of whole slide images into individual ROIs/patches.

    This function should be called independently for the x- and y-coordinates. Units should be consistent between all four arguments.

    Sample usage:
        x_roi_starts, x_roi_ends = calculate_roi_coords(min_coord=3, max_coord=16, roi_width=2, overlap=2*thickness)
        y_roi_starts, y_roi_ends = calculate_roi_coords(min_coord=-150, max_coord=300, roi_width=20, overlap=2*thickness)

    Args:
        min_coord (float): Smallest coordinate in the slide
        max_coord (float): Largest coordinate in the slide
        roi_width (float): Desired width of each ROI
        overlap (float, optional): Amount you want consecutive ROIs to overlap. Note that for the SIP code, in order to utilize all data most efficiently, this should equal twice the density radius (e.g., the "thickness" value in main.ipynb). Defaults to 0.

    Returns:
        1D numpy array: Beginning (smaller) coordinates of each patch
        1D numpy array: Ending (larger) coordinates of each patch
    """

    # Import relevant libraries
    import numpy as np

    # Calculate relevant scalars
    # See notebook notes on 4-15-22 for clarification
    num_rois = int(np.ceil((max_coord - min_coord - overlap) / (roi_width - overlap)))
    roi_span = num_rois * roi_width - (num_rois - 1) * overlap
    first_roi_start = min_coord - (roi_span - (max_coord - min_coord)) / 2
    first_roi_end = first_roi_start + roi_width
    roi_step_size = roi_width - overlap

    # Calculate relevant vectors
    roi_starts = np.arange(num_rois) * roi_step_size + first_roi_start
    roi_ends = np.arange(num_rois) * roi_step_size + first_roi_end

    # Return the start and end coordinates
    return roi_starts, roi_ends


def break_up_slide_into_patches(df, roi_width=1230.77, overlap=246.15):
    """Break up a slide's dataframe into patches/ROIs, returning a Pandas series with the appropriate ROI labels for each object

    Sample call:
        df.loc[curr_loc, 'tag'] = break_up_slide_into_patches(curr_df, roi_width=(roi_width / coord_units_in_microns), overlap=(overlap / coord_units_in_microns))

    Args:
        df (Pandas dataframe): Dataframe corresponding to a single slide containing coordinates (columns XMin, XMax, YMin, YMax) in pixels
        roi_width (float, optional): Desired ROI width in pixels. Defaults to 1230.77 ~ 400 / 0.325.
        overlap (float, optional): Desired patching overlap in pixels. Defaults to 246.15 ~ 80 / 0.325.

    Returns:
        Pandas series: Series corresponding to the tag/ROI column of the dataframe of the input slide
    """

    # Import relevant libraries
    import numpy as np

    # Get items for the ROI name
    max_num_digits = len(str(max([df['XMax'].max(), df['YMax'].max()])))  # get the maximum number of digits to format into the ROI names; assume the input coordinates are non-negative integers
    format_string = '[{:0' + str(max_num_digits) + 'd},{:0' + str(max_num_digits) + 'd}]'
    # slide_id = df.head(1)['Slide ID']  # get the current Slide ID; it should be same for every row in the input dataframe
    slide_id = df['Slide ID'].iloc[0]
    print('On slide {}...'.format(slide_id))

    # Get the patch coordinates for the entire slide
    x_roi_starts, x_roi_ends = calculate_roi_coords(min_coord=df['XMin'].min(), max_coord=df['XMax'].max(), roi_width=roi_width, overlap=overlap)
    y_roi_starts, y_roi_ends = calculate_roi_coords(min_coord=df['YMin'].min(), max_coord=df['YMax'].max(), roi_width=roi_width, overlap=overlap)

    # # Print out the minimum coordinates of every patch
    # print('Breaking up the slide into patches:')
    # print('  x minimum coordinates: {}'.format(x_roi_starts))
    # print('  y minimum coordinates: {}'.format(y_roi_starts))

    # Initialize four relevant columns in the dataframe
    # df.loc[:, 'tag'] = []
    # df.loc[:, 'roi_count'] = 0
    df['tag'] = [[]] * len(df)
    df['roi_count'] = 0
    xmid = df['XMin'] + (df['XMax'] - df['XMin']) / 2
    ymid = df['YMin'] + (df['YMax'] - df['YMin']) / 2

    # Loop over every ROI and assign an object to a ROI accordingly
    total_object_count = 0
    for xmin, xmax in zip(x_roi_starts, x_roi_ends):
        for ymin, ymax in zip(y_roi_starts, y_roi_ends):
            in_roi = (xmid >= xmin) & (xmid < xmax) & (ymid >= ymin) & (ymid < ymax)
            roi_name = slide_id + '_roi_' + format_string.format(int(np.round((xmin + xmax) / 2)), int(np.round((ymin + ymax) / 2)))  # the ROI name is the slide ID appended with the middle x- and y-coordinates of the current ROI
            # df.loc[in_roi, 'tag'] = roi_name
            # df.loc[in_roi, 'tag'] = df.loc[in_roi, 'tag'].apply(lambda x: x.append(roi_name), axis='rows')
            # df.loc[in_roi, 'tag'] = df.loc[in_roi, 'tag'].apply(lambda x: x.append(roi_name))
            df.loc[in_roi, 'tag'] = df.loc[in_roi, 'tag'].apply(lambda x: x + [roi_name])
            df.loc[in_roi, 'roi_count'] = df.loc[in_roi, 'roi_count'] + 1
            num_objects_in_roi = in_roi.sum()
            total_object_count = total_object_count + num_objects_in_roi
            print('Number of objects in ROI {}: {}'.format(roi_name, num_objects_in_roi))

    # Output some items to check
    print(len(df), total_object_count, df['roi_count'].sum())  # if there is overlap then the latter two numbers, which should be the same, should be at least as big as the first number. If there is no overlap, then all three numbers should be the same.
    print('Unique ROI counts: {}'.format(df['roi_count'].unique()))  # should be [1, 2, 4] if there is overlap; otherwise, just [1]

    # Return the ROI names
    return df['tag']


def duplicate_rows(df, col_to_repeat='tag'):
    """Duplicate rows of a dataframe based on a column containing iterables whose entries should be different for each duplicate

    Args:
        df (Pandas dataframe): Dataframe whose rows should be repeated with the differences in repeated rows specified by the entries in the col_to_repeat column
        col_to_repeat (str, optional): Name of the column containing iterables whose values should be the difference between the duplicated rows. Defaults to 'tag'.

    Returns:
        Pandas dataframe: New dataframe with repeated rows except for the col_to_repeat column
    """

    # Import relevant library
    import pandas as pd

    # Get the lengths of the entries in the column on which repetition should be performed
    col_lengths = df[col_to_repeat].apply(lambda x: len(x))

    # Get the unique lengths
    unique_col_lengths = col_lengths.unique()

    # Initialize a list from which the new dataframe will be built
    new_df_list = []

    # For each unique length...
    for unique_col_length in unique_col_lengths:

        # Get a view of the main dataframe with just the rows having the current column length
        curr_df = df[col_lengths == unique_col_length]

        # For each item in each row's iterable...
        for list_index in range(unique_col_length):

            # Make a copy of the view
            new_curr_df = curr_df.copy()

            # Replace the contents of the column of interest with the current item in the iterable
            new_curr_df[col_to_repeat] = curr_df[col_to_repeat].apply(lambda x: x[list_index])

            # Append the new dataframe to the master dataframe list
            new_df_list.append(new_curr_df)

    # Create the new Pandas dataframe from the master list
    new_df = pd.concat(new_df_list).reset_index(drop=True)

    # Return the new dataframe
    return new_df


def calculate_min_coord_spacing_per_slide(df):

    # Not actually used anywhere but is a good function.

    # Import relevant libraries
    # import pandas as pd
    import numpy as np

    # # Get the data in the appropriate format
    # datafile = '../data/all_five_images_NOS2_and_COX2_and_CD8_Fused_ALW_05212021_entire_image.csv'
    # df = pd.read_csv(datafile)
    # df['Cell X Position'] = (df['XMin'] + (df['XMax'] - df['XMin']) / 2) * 0.325
    # df['Cell Y Position'] = (df['YMin'] + (df['YMax'] - df['YMin']) / 2) * 0.325
    # df['Slide ID'] = df['Image Location']

    # Obtain the unique slides in the dataset
    unique_slides = df['Slide ID'].unique()

    # Initialize a minimum spacing holder
    min_spacing_holder = np.zeros((len(unique_slides), 2))

    # For each slide in the dataset...
    for islide, unique_slide in enumerate(unique_slides):

        # Store the dataset for the current slide
        df_curr_slide = df[df['Slide ID'] == unique_slide]

        # For each of the two coordinate columns...
        for icoord_col, coord_col in enumerate(['Cell X Position', 'Cell Y Position']):

            # Get the current coordinate column for the current slide
            srs_coords = df_curr_slide[coord_col]

            # Get the unique, sorted values in the coordinates (should never be zero because we're taking the unique values before we take the finite difference)
            coords_unique_sorted = srs_coords.sort_values().unique()

            # Store the minimum of the finite difference (which again is not zero)
            min_spacing_holder[islide, icoord_col] = (coords_unique_sorted[1:] - coords_unique_sorted[:-1]).min()

    # Calculate the minimum coordinate spacing
    detected_minimum_spacing = round(min_spacing_holder.min(), 8)

    # Print the output
    print('The detected minimum spacing is {}. Note that this can be 1/2 of the conversion from pixels to microns due to the calculation being based on the coordinate midpoints.'.format(detected_minimum_spacing))

    # # These are helpful outputs but not necessary
    # print(min_spacing_holder)
    # print(min_spacing_holder.min())
    # print(detected_minimum_spacing)

    # Return the calculated minimum coordinate spacing
    return detected_minimum_spacing