def full_df_to_just_coords_in_microns(df_full, coord_units_in_microns):
    """Create a function to take in the full set of columns, extract just the spatial bounds columns, calculate the midpoints, and convert to microns.

    Args:
        df_full (Pandas dataframe): Must contain the columns XMin, XMax, YMin, and YMax.
        coord_units_in_microns (float): Number of microns per coordinate unit in the input dataframe.

    Returns:
        Pandas dataframe: Contains the input dataframe, plus new xmid and ymid columns, all converted to microns.
    """

    # Get a copy of the input dataframe with just the coordinate bounds columns
    df_coords = df_full[['XMin', 'XMax', 'YMin', 'YMax']].copy()

    # Calculate the midpoints of the bounds columns, i.e., get the midpoint of each input object
    df_coords['xmid'] = (df_coords['XMax'] - df_coords['XMin']) / 2 + df_coords['XMin']
    df_coords['ymid'] = (df_coords['YMax'] - df_coords['YMin']) / 2 + df_coords['YMin']

    # Convert all columns to microns
    df_coords = df_coords.map(lambda x: x * coord_units_in_microns)

    # Return the new dataframe
    return df_coords

def load_annotation_data(annotations_dir, csv_files, phenotyping_method, phenotype_identification_file):
    # settings__phenotyping__method, settings__phenotyping__phenotype_identification_file
    """Load into dataframes and all the metadata for the files in the annotations directory.

    Args:
        annotations_dir (str, optional): Name of the directory in which the annotation data reside. Defaults to 'annotations'.
        region_type_prefix (str, optional): Prefix to the region ID in the annotation filenames. Defaults to '05212021_'.

    Returns:
        Pandas dataframe: Dataframe holding all the dataframes of the annotation data, plus their metadata (image ID, region type, CSV file path, and number of objects).
    """

    # Import relevant libraries
    import os
    import pandas as pd
    import dataset_formats
    import new_phenotyping_lib

    # Obtain the full list of CSV files
    # csv_files = sorted(os.listdir(annotations_dir))
    csv_files = sorted(csv_files)

    # Create a dataframe holding all metadata and the dataframes corresponding to all CSV files
    all_data = []
    for csv_file in csv_files:
        image_id = csv_file.split('__')[0]
        region_type = csv_file.split('.')[0].split('__')[-1]
        csv_file_path = os.path.join(annotations_dir, csv_file)
        _, _, _, marker_prefix, _, markers = dataset_formats.extract_datafile_metadata(csv_file_path)
        marker_column_names_list = [marker_prefix + x for x in markers]
        # df = pd.read_csv(csv_file_path)
        df = new_phenotyping_lib.apply_phenotyping(csv_file_path=csv_file_path, method=phenotyping_method, phenotype_identification_file=phenotype_identification_file)
        if marker_column_names_list != []:
            # marker_column_names_list = ast.literal_eval(marker_column_names_str)
            df = df[df[marker_column_names_list].apply(sum, axis='columns') != 0]
            df_length_marker = df[marker_column_names_list].sum().sum()
        else:
            df_length_marker = -1
        all_data.append({'image_id': image_id, 'region_type': region_type, 'csv_file_path': csv_file_path, 'dataframe': df, 'df_length': len(df), 'df_length_marker': df_length_marker})
    df = pd.DataFrame(all_data)

    # Return the resulting dataframe
    return df

def plot_annotation_data(df, coord_units_in_microns=0.325, alpha=0.4, buffer_frac=0.05, generate_single_figure=True, figsize=(10, 10)):
    """Either plot all annotation data for all images in a single figure (generate_single_figure=True, default), returning the overall figure handle, or send back the individual figures in a nested list.

    Args:
        df (Pandas dataframe): Dataframe holding all the annotation data, as read in e.g. by load_annotation_data().
        coord_units_in_microns (float, optional): Number of microns per coordinate unit in the input dataframe. Defaults to 0.325.
        alpha (float, optional): Transparency value for the plots of cells that can be on top of each other. Defaults to 0.4.
        buffer_frac (float, optional): Fraction of the "tight" axes configuration to add as a spatial buffer in the four edges of the plots. Defaults to 0.05.
        generate_single_figure (bool, optional): Whether to generate a single image holding all the figures together, or otherwise to send back a nested list of matplotlib figures. Defaults to True.

    Returns:
        Nested list of figure handles if generate_single_figure=False; single figure handle otherwise. If the former, also return the image IDs and the region types in the same order they appear in the nested list.
    """

    # Import relevant libraries
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np

    # Get the Seaborn color palette
    color_palette = sns.color_palette()

    # Get the image IDs and the region types in order of decreasing overall prevalence (for plotting purposes)
    image_ids = df['image_id'].unique()
    num_images = len(image_ids)
    # region_types_ordered = df.groupby(by='region_type')['df_length'].apply(sum).sort_values(ascending=False).index
    region_types_ordered = df.groupby(by='region_type')['df_length'].sum().sort_values(ascending=False).index
    num_region_types = len(region_types_ordered)

    # Create a new summary dataframe indexed by the image IDs and annotation types
    df_multiindex = df.set_index(keys=['image_id', 'region_type'])

    # Create a holder for all images, whether a list or a set of axes
    if not generate_single_figure:
        single_figure_holder = []
    else:
        fig, ax = plt.subplots(ncols=(num_region_types + 1), nrows=num_images, figsize=(10 * (num_region_types + 1), 10 * num_images))

    # For each image...
    for iimage, image_id in enumerate(image_ids):

        # Initialize unrealistic x and y limits
        x_range = np.array([100000, 0]) * coord_units_in_microns
        y_range = np.array([100000, 0]) * coord_units_in_microns

        # For the current image, get the axes for the composite plot
        if not generate_single_figure:
            fig_composite, ax_composite = plt.subplots(figsize=figsize)
        else:
            if num_images > 1:
                ax_composite = ax[iimage, 0]
            else:
                ax_composite = ax[0]

        # For each region type in the dataset...
        single_figure_holder_per_region = []
        for iregion_type, region_type in enumerate(region_types_ordered):

            # Generate the figure and axis in the multiple-figure case
            if not generate_single_figure:
                fig_single_region, ax_single_region = plt.subplots(figsize=figsize)
            else:
                if num_images > 1:
                    ax_single_region = ax[iimage, iregion_type + 1]
                else:
                    ax_single_region = ax[iregion_type + 1]

            # Grab the dataframe for the current image/region type
            try:

                # Get the current dataframe
                curr_df = df_multiindex.loc[image_id, region_type]['dataframe']

                # Get the number of datapoints in the current dataframe
                df_length = df_multiindex.loc[image_id, region_type]['df_length']
                df_length_marker = df_multiindex.loc[image_id, region_type]['df_length_marker']
                assert df_length == len(curr_df)

                # Get the midpoints for the current dataframe
                curr_df = full_df_to_just_coords_in_microns(curr_df, coord_units_in_microns)

                # Update the x and y limits
                x_range[0] = min(x_range[0], curr_df['xmid'].min())
                x_range[1] = max(x_range[1], curr_df['xmid'].max())
                y_range[0] = min(y_range[0], curr_df['ymid'].min())
                y_range[1] = max(y_range[1], curr_df['ymid'].max())

                # Plot the dataset for the current region type
                sns.scatterplot(data=curr_df, x='xmid', y='ymid', ax=ax_composite, alpha=alpha, label='{} ({}/{} cells)'.format(region_type, df_length, df_length_marker))
                sns.scatterplot(data=curr_df, x='xmid', y='ymid', ax=ax_single_region, alpha=alpha, label='{} ({}/{} cells)'.format(region_type, df_length, df_length_marker), color=color_palette[iregion_type])

            # If the current region type doesn't exist for the current image, skip that
            except KeyError:
                pass

            # Save the current region-specific figures to the figure holder
            if not generate_single_figure:
                single_figure_holder_per_region.append(fig_single_region)
                plt.close(fig_single_region)

        # Save the composite figure to the figure holder
        if not generate_single_figure:
            # single_figure_holder[iimage][0] = fig_composite
            single_figure_holder_per_region = [fig_composite] + single_figure_holder_per_region
            plt.close(fig_composite)

        # Update the x and y limits in the scatter plot to include a buffer, nominally 5% on either end of an axis
        x_buffer = (x_range[1] - x_range[0]) * buffer_frac
        x_range[0] = x_range[0] - x_buffer
        x_range[1] = x_range[1] + x_buffer
        y_buffer = (y_range[1] - y_range[0]) * buffer_frac
        y_range[0] = y_range[0] - y_buffer
        y_range[1] = y_range[1] + y_buffer

        # Set some axis properties for each region type in the dataset
        if not generate_single_figure:
            ax_composite = single_figure_holder_per_region[0].gca()
        else:
            if num_images > 1:
                ax_composite = ax[iimage, 0]
            else:
                ax_composite = ax[0]
        ax_composite.set_xlim(x_range)
        ax_composite.set_ylim(y_range)
        ax_composite.invert_yaxis()
        ax_composite.set_title('Image {} - combined'.format(image_id))
        ax_composite.set_xlabel('x-coordinates (microns)')
        ax_composite.set_ylabel('y-coordinates (microns)')
        ax_composite.set_aspect('equal')
        for iregion_type, region_type in enumerate(region_types_ordered):
            if not generate_single_figure:
                ax_single_region = single_figure_holder_per_region[iregion_type + 1].gca()
            else:
                if num_images > 1:
                    ax_single_region = ax[iimage, iregion_type + 1]
                else:
                    ax_single_region = ax[iregion_type + 1]
            ax_single_region.set_xlim(x_range)
            ax_single_region.set_ylim(y_range)
            ax_single_region.invert_yaxis()
            ax_single_region.set_title('Image {} - {}'.format(image_id, region_type))
            ax_single_region.set_xlabel('x-coordinates (microns)')
            ax_single_region.set_ylabel('y-coordinates (microns)')
            ax_single_region.set_aspect('equal')

        # Update the main figure holder
        if not generate_single_figure:
            single_figure_holder.append(single_figure_holder_per_region)

    # If we're not plotting all figures in a single image, return the nested list holding all the individual figures
    if not generate_single_figure:
        return single_figure_holder, image_ids, region_types_ordered
    else:
        return fig

def overlay_on_background_data(df_to_overlay, df_background, alpha=0.4, buffer_frac=0.05, figsize=(10, 10), annotation_color_index=0, analysis_color_index=4, df_background_length=None, df_background_length_marker=None):

    # Import relevant libraries
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np

    # Get the Seaborn color palette
    color_palette = sns.color_palette()

    # Initialize axes
    fig, ax = plt.subplots(figsize=figsize)

    # Run a quick check
    assert df_background_length == len(df_background)

    # Plot the dataset for the current region type
    sns.scatterplot(data=df_background, x='xmid', y='ymid', ax=ax, alpha=alpha, label='{} ({}/{} cells)'.format('annotation data', df_background_length, df_background_length_marker), color=color_palette[annotation_color_index])
    sns.scatterplot(data=df_to_overlay, x='Cell X Position', y='Cell Y Position', ax=ax, alpha=alpha, label='{} ({} cells)'.format('analysis data', len(df_to_overlay)), color=color_palette[analysis_color_index])

    # Update the x and y limits in the scatter plot to include a buffer, nominally 5% on either end of an axis
    x_range = np.array([
        min(df_background['xmid'].min(), df_to_overlay['Cell X Position'].min()),
        max(df_background['xmid'].max(), df_to_overlay['Cell X Position'].max())
    ])
    y_range = np.array([
        min(df_background['ymid'].min(), df_to_overlay['Cell Y Position'].min()),
        max(df_background['ymid'].max(), df_to_overlay['Cell Y Position'].max())
    ])
    x_buffer = (x_range[1] - x_range[0]) * buffer_frac
    x_range[0] = x_range[0] - x_buffer
    x_range[1] = x_range[1] + x_buffer
    y_buffer = (y_range[1] - y_range[0]) * buffer_frac
    y_range[0] = y_range[0] - y_buffer
    y_range[1] = y_range[1] + y_buffer

    # Set some axis properties
    ax.set_xlim(x_range)
    ax.set_ylim(y_range)
    ax.invert_yaxis()
    ax.set_title('Image overlaid on annotation')
    ax.set_xlabel('x-coordinates (microns)')
    ax.set_ylabel('y-coordinates (microns)')
    ax.set_aspect('equal')

    # Return the generated figure
    return fig

def recreate_data_property(df_data_by_roi):
    # Call like: recreate_data_property(slices.df_data_by_roi)

    # Import relevant library
    import pandas as pd

    # Create a list to hold one dictionary holding the data for each ROI
    roi_data_holder = []

    # For each ROI in df_data_by_roi...
    for roi_data in df_data_by_roi.iterrows():

        # Get a shortcut to the current Series
        curr_roi_data = roi_data[1]

        # Get the number of cells/objects in the current ROI
        num_cells = len(curr_roi_data['x_roi'])

        # Get or create the relevant arrays of data
        slide_name = [curr_roi_data['unique_slide']] * num_cells
        roi_name = [curr_roi_data['unique_roi']] * num_cells
        x_roi = curr_roi_data['x_roi']
        y_roi = curr_roi_data['y_roi']
        species_roi = curr_roi_data['species_roi']
        case_id = [curr_roi_data['unique_case']] * num_cells

        # Create a single dictionary holding all the data for the current ROI and append it to the dictionary holder
        roi_data_holder.append(pd.DataFrame({'Slide ID': slide_name, 'tag': roi_name, 'Cell X Position': x_roi, 'Cell Y Position': y_roi, 'Species int': species_roi, 'Case ID': case_id}))

    # Turn the list of dictionaries into a Dataframe
    df = pd.concat(roi_data_holder, ignore_index=True)

    # This is an old dataframe creator for the goal dataframe, originally used for function testing in main.ipynb
    # df1 = slices.data[df.columns].sort_values(df.columns.to_list()).reset_index(drop=True)

    # Sort the data by the columns in turn (and cast a column to int64), which is necessary only for my testing comparison to the goal dataframe (see commented out line above)
    df2 = df[df.columns].sort_values(df.columns.to_list()).reset_index(drop=True)
    df2['Species int'] = df2['Species int'].astype('int64')
    # df1.equals(df2)  # this was the check in main.ipynb to make sure we were faithfully recreating the "data" property of the TIMECellInteraction class

    # Return the result
    return df2

def transform_to_integer(data_to_transform, inverse_scale=1.0, neg_shift_after_scaling=0.0):
    """Transform input data to integers using specified scales and shifts.

    Args:
        data_to_transform (int, float, Series, list, ndarray): Possibly float data to transform to integers.
        inverse_scale (float, optional): 1 / scale factor to which first scale the input data. Defaults to 1.0.
        neg_shift_after_scaling (float, optional): Constant shift to subtract from the scaled data. Defaults to 0.0.

    Returns:
        int, Series, ndarray: Transformed-to-integer data.
    """

    # Import relevant libraries
    import numpy as np
    import pandas as pd

    # Convert the data to a numpy array, noting whether the input data is a scalar or vector and therefore whether to transform back to a scalar at the end of this function
    if type(data_to_transform) in [int, float]:
        data_to_transform = np.array([data_to_transform])
        convert_to_scalar_at_end = True
    elif type(data_to_transform) is pd.core.series.Series:
        convert_to_scalar_at_end = False
    else:
        data_to_transform = np.array(data_to_transform)
        convert_to_scalar_at_end = False

    # Perform the transform (scale, shift, cast)
    transformed_data = (data_to_transform / inverse_scale - neg_shift_after_scaling).astype(int)

    # If the original input data were scalar, convert back to scalar
    if convert_to_scalar_at_end:
        transformed_data = transformed_data[0]

    # Return the transformed dasta
    return transformed_data

def transform_annotation_data_to_integers(df_annotation, df_data_by_roi_for_selected_slide, annotation_microns_per_integer_unit=None):
    """Transform annotation min-max object coordinates to integers and create a corresponding 2D Boolean matrix that's True where an object exists.

    The transformation parameters--namely, the shifts--depend on the ROI mins-maxs in the analysis data so that the canvas for the image just fits both the analysis and annotation data.

    Args:
        df_annotation (dataframe): Annotation data dataframe for a given image and a given annotation region type (coordinates in microns, which are the same as in the analysis data below)
        df_data_by_roi_for_selected_slide (dataframe): Analysis data where one row corresponds to one ROI, for just the slide corresponding to that for the annotation data df_annotation
        TODO: add appropriate docs line for annotation_microns_per_integer_unit

    Returns:
        annot_int_arr, df_annotation, min_spacing_for_annotation, x_neg_shift_after_scaling, y_neg_shift_after_scaling
        ndarray: 2D array where indices are integer coordinates and values are True or False
        dataframe: Annotation dataframe with four new columns corresponding to the corresponding columns that were transformed
        float: Minimum coordinate spacing from the annotation extrema data, used to to generate integers from these data
        float: Value to negatively shift (after scaling) in the x direction as part of the integer transform
        float: Value to negatively shift (after scaling) in the y direction as part of the integer transform
    """
    
    # Import relevant library
    import numpy as np

    # If the user hasn't defined the annotation coordinate conversion from microns to integer units...
    if annotation_microns_per_integer_unit is None:
    
        # Columns holding the min and max coordinates for every object in the annotation file (which applies to a single slide, i.e., it corresponds to a slide)
        coord_cols = ['XMin', 'XMax', 'YMin', 'YMax']

        # Initialize a list holding the minimum spacing in each of these coordinate columns
        min_spacings = []

        # For each column of coordinate extrema...
        for coord_col in coord_cols:

            # Sort the values and convert to a numpy array
            tmp = df_annotation[coord_col].sort_values().to_numpy()

            # Add to the list the first non-trivial spacing between the values
            min_spacings.append(round(list(set(np.unique(tmp[1:] - tmp[:-1])) - set([0]))[0], ndigits=8))

        # Determine the minimum coordinate spacing for the current annotation by taking the minimum over all four columns. This is the number we'll use to scale all values in (from) microns (both in the annotations and in the real data)
        min_spacing_for_annotation = np.min(min_spacings)

    # ... if they have, then use what they defined
    else:

        min_spacing_for_annotation = annotation_microns_per_integer_unit

    # Get the minimum and maximum annotation x and y values in the current annotation, in particular, the minimum and maximum object gate coordinates
    x_range_annotation = np.array([df_annotation['XMin'].min(), df_annotation['XMax'].max()])
    y_range_annotation = np.array([df_annotation['YMin'].min(), df_annotation['YMax'].max()])

    # Get the minimum and maximum x and y values in the entire slide, in particular, the minimum and maximum ROI coordinates for the entire slide of data on which calculations were made
    x_range_slide = np.array([df_data_by_roi_for_selected_slide['x_range'].apply(min).min(), df_data_by_roi_for_selected_slide['x_range'].apply(max).max()])
    y_range_slide = np.array([df_data_by_roi_for_selected_slide['y_range'].apply(min).min(), df_data_by_roi_for_selected_slide['y_range'].apply(max).max()])

    # Obtain the initial shift to the scaled coordinates in order to get at least the annotation coordinates starting exactly at zero and ensuring they are all integers
    x_neg_shift_after_scaling = x_range_annotation[0] / min_spacing_for_annotation
    y_neg_shift_after_scaling = y_range_annotation[0] / min_spacing_for_annotation

    # In case the data coordinates are more negative than the annotation coordinates, then don't make the shift as quite as extreme so that the most negative of any data or annotation coordinate fall at zero exactly. This is the net, final shift to any scaled coordinates
    most_negative_x_coord = transform_to_integer(x_range_slide[0], inverse_scale=min_spacing_for_annotation, neg_shift_after_scaling=x_neg_shift_after_scaling)
    if most_negative_x_coord < 0:
        x_neg_shift_after_scaling = x_neg_shift_after_scaling + most_negative_x_coord
    most_negative_y_coord = transform_to_integer(y_range_slide[0], inverse_scale=min_spacing_for_annotation, neg_shift_after_scaling=y_neg_shift_after_scaling)
    if most_negative_y_coord < 0:
        y_neg_shift_after_scaling = y_neg_shift_after_scaling + most_negative_y_coord

    # Get ultimate width of the converted-to-integers array of coordinates just big enough to hold both the annotation and real data
    num_x_indices = np.max([transform_to_integer(x_range_annotation[1], inverse_scale=min_spacing_for_annotation, neg_shift_after_scaling=x_neg_shift_after_scaling),
                            transform_to_integer(x_range_slide[1], inverse_scale=min_spacing_for_annotation, neg_shift_after_scaling=x_neg_shift_after_scaling)
                            ]) + 1
    num_y_indices = np.max([transform_to_integer(y_range_annotation[1], inverse_scale=min_spacing_for_annotation, neg_shift_after_scaling=y_neg_shift_after_scaling),
                            transform_to_integer(y_range_slide[1], inverse_scale=min_spacing_for_annotation, neg_shift_after_scaling=y_neg_shift_after_scaling)
                            ]) + 1

    # Initialize the integer array holding all the annotation data
    annot_int_arr = np.full((num_x_indices, num_y_indices), False)

    # Add the relevant transformed columns to the annotation data
    df_annotation['XMin_int'] = transform_to_integer(df_annotation['XMin'], inverse_scale=min_spacing_for_annotation, neg_shift_after_scaling=x_neg_shift_after_scaling)
    df_annotation['XMax_int'] = transform_to_integer(df_annotation['XMax'], inverse_scale=min_spacing_for_annotation, neg_shift_after_scaling=x_neg_shift_after_scaling)
    df_annotation['YMin_int'] = transform_to_integer(df_annotation['YMin'], inverse_scale=min_spacing_for_annotation, neg_shift_after_scaling=y_neg_shift_after_scaling)
    df_annotation['YMax_int'] = transform_to_integer(df_annotation['YMax'], inverse_scale=min_spacing_for_annotation, neg_shift_after_scaling=y_neg_shift_after_scaling)

    # Set to true the indices corresponding to the span of each annotation object in the 2D annotation integer array
    # Note this allows us to plot the annotation data using, e.g., plt.imshow(annot_int_arr.T) (preferred) or plt.matshow(annot_int_arr.T)
    for curr_row in df_annotation.itertuples():
        annot_int_arr[curr_row.XMin_int:(curr_row.XMax_int+1), curr_row.YMin_int:(curr_row.YMax_int+1)] = True

    # Return the important variables
    return annot_int_arr, min_spacing_for_annotation, x_neg_shift_after_scaling, y_neg_shift_after_scaling

def add_raw_weighting_values_to_roi_data(df_data_by_roi_for_selected_slide, df_annotation, thickness, annot_int_arr, min_spacing_for_annotation, x_neg_shift_after_scaling, y_neg_shift_after_scaling, annotation_region_type):
    """Add, in a new dataframe, four columns to the by-ROI dataframe containing raw values that could aid in weighting the ROIs based on the selected annotation region type.

    TODO: Update args/returns below

    Args:
        df_data_by_roi_for_selected_slide (dataframe): Analysis data where one row corresponds to one ROI
        df_annotation (dataframe): Annotation data dataframe for a given image and a given annotation region type (coordinates in microns, which are the same as in the analysis data below)
        thickness (float): Analysis radius around each center in microns

    Returns:
        dataframe: Analysis data with four columns added that would help in weighting the ROIs per annotation region type
    """

    # Import relevant library
    import numpy as np
    import pandas as pd

    # Initialize a dataframe with the same index as that in df_data_by_roi_for_selected_slide; this is what we'll update in this function and simply return at the end
    valid_area_microns_sq_colname = 'valid_area_microns_sq_{}'.format(annotation_region_type)
    roi_integer_area_colname = 'roi_integer_area_{}'.format(annotation_region_type)
    num_ann_objects_within_roi_colname = 'num_ann_objects_within_roi_{}'.format(annotation_region_type)
    footprint_integer_area_colname = 'footprint_integer_area_{}'.format(annotation_region_type)
    df_weights_data = pd.concat([
        pd.Series(index=df_data_by_roi_for_selected_slide.index, dtype=np.float64, name=valid_area_microns_sq_colname),
        pd.Series(index=df_data_by_roi_for_selected_slide.index, dtype=np.int64, name=roi_integer_area_colname),
        pd.Series(index=df_data_by_roi_for_selected_slide.index, dtype=np.int64, name=num_ann_objects_within_roi_colname),
        pd.Series(index=df_data_by_roi_for_selected_slide.index, dtype=np.int64, name=footprint_integer_area_colname)
        ], axis='columns')

    # For every ROI in the dataset...
    for roi_name in df_data_by_roi_for_selected_slide['unique_roi']:

        # Store the current ROI data
        curr_data_by_roi = df_data_by_roi_for_selected_slide[df_data_by_roi_for_selected_slide['unique_roi'] == roi_name]  # get the one-row dataframe corresponding to the data in slices.df_data_by_roi for the current ROI
        curr_roi_index = curr_data_by_roi.index[0]

        # These are the same as the ranges calculated from the original dataset (slices.data) coordinates (doing the mins and maxs), shown below
        # According to the plots, which show consistent limits, these are in microns, which is further consistent with the fact that in dataset_formats.py I ensure the units of the datasets are always microns
        # The annotations dataframes in annotations.py are also in microns
        x_mid = curr_data_by_roi['x_roi'].iloc[0]
        y_mid = curr_data_by_roi['y_roi'].iloc[0]
        x_range_stored = curr_data_by_roi.iloc[0]['x_range']
        y_range_stored = curr_data_by_roi.iloc[0]['y_range']
        x_range_calc = np.array([x_mid.min(), x_mid.max()])
        y_range_calc = np.array([y_mid.min(), y_mid.max()])
        if (not (x_range_stored == x_range_calc).all()) or (not (y_range_stored == y_range_calc).all()):
            print('ERROR: Ranges are not the same between the calculated method and the stored method')

        # Run some other checks on the data in slices.df_data_by_roi
        if (not ((x_range_stored == np.array([curr_data_by_roi.iloc[0]['x_min_prior_to_decimation'], curr_data_by_roi.iloc[0]['x_max_prior_to_decimation']])).all())) or (not ((y_range_stored == np.array([curr_data_by_roi.iloc[0]['y_min_prior_to_decimation'], curr_data_by_roi.iloc[0]['y_max_prior_to_decimation']])).all())):
            print('ERROR')
        if (curr_data_by_roi.iloc[0]['width'] != (x_range_stored[1] - x_range_stored[0])) or (curr_data_by_roi.iloc[0]['height'] != (y_range_stored[1] - y_range_stored[0])):
            print('ERROR')

        # slices.thickness is in microns and is the radius
        # Thus, the VALID region within each ROI has these bounds
        radius = thickness
        x_range_valid = [x_range_stored[0] + radius, x_range_stored[1] - radius]
        y_range_valid = [y_range_stored[0] + radius, y_range_stored[1] - radius]

        # Ensure the valid regions of the ROI are physically existing regions (i.e., have non-zero area)
        if not ((x_range_valid[1] > x_range_valid[0]) and (y_range_valid[1] > y_range_valid[0])):
            # i.e., if there is no valid area
            # i.e., if (x_max <= x_min) or (y_max <= x_min)
            # print('WARNING: ROI {} has NO area'.format(tag))
            valid_area_microns = 0
            roi_integer_area = 0
            num_ann_objects_within_roi = 0
            footprint_integer_area = 0
        else:
            # i.e., if there is a real valid area
            # print('ROI {} has finite area'.format(tag))
            valid_area_microns = (x_range_valid[1] - x_range_valid[0]) * (y_range_valid[1] - y_range_valid[0])
            x_range_valid_int = transform_to_integer(x_range_valid, inverse_scale=min_spacing_for_annotation, neg_shift_after_scaling=x_neg_shift_after_scaling)
            y_range_valid_int = transform_to_integer(y_range_valid, inverse_scale=min_spacing_for_annotation, neg_shift_after_scaling=y_neg_shift_after_scaling)
            roi_integer_area = (x_range_valid_int[1] - x_range_valid_int[0] + 1) * (y_range_valid_int[1] - y_range_valid_int[0] + 1)  # integer area of the entire ROI
            num_ann_objects_within_roi = ((df_annotation['xmid'] >= x_range_valid[0]) & (df_annotation['xmid'] < x_range_valid[1]) & (df_annotation['ymid'] >= y_range_valid[0]) & (df_annotation['ymid'] < y_range_valid[1])).sum()
            footprint_integer_area = annot_int_arr[x_range_valid_int[0]:(x_range_valid_int[1] + 1), y_range_valid_int[0]:(y_range_valid_int[1] + 1)].sum()  # determine the integer footprint area of the annotation objects in the current ROI. If the scaling takes us back to pixels, then the units of these areas (this and roi_integer_area) are pixels

        # Add the calculated data to the ROI-by-ROI dataframe
        df_weights_data.loc[curr_roi_index, valid_area_microns_sq_colname] = valid_area_microns
        df_weights_data.loc[curr_roi_index, roi_integer_area_colname] = roi_integer_area
        df_weights_data.loc[curr_roi_index, num_ann_objects_within_roi_colname] = num_ann_objects_within_roi
        df_weights_data.loc[curr_roi_index, footprint_integer_area_colname] = footprint_integer_area

    # Return the modified dataframe
    return df_weights_data

def drop_duplicate_columns_in_weighting_data(df_data_by_roi, column_prefixes=['valid_area_microns_sq_', 'roi_integer_area_', 'num_ann_objects_within_roi_', 'footprint_integer_area_']):
    """Drop the duplicate columns containing raw weighting values likely recently added to the df_data_by_roi dataframe.

    Args:
        df_data_by_roi (Pandas dataframe): Original dataframe holding all the data by ROI.
        column_prefixes (list, optional): Holder of column prefixes to check one-by-one. Defaults to ['valid_area_microns_sq_', 'roi_integer_area_', 'num_ann_objects_within_roi_', 'footprint_integer_area_'].

    Returns:
        Pandas dataframe: Processed version of the input dataframe.
    """

    # Import relevant library
    import numpy as np

    # For each column prefix...
    for column_prefix in column_prefixes:

        # Determine the columns having that prefix in their name
        matching_columns = [column for column in df_data_by_roi.columns if column.startswith(column_prefix)]

        # Determine whether the entire columns (through all the rows, i.e., ROIs) are identical to each other
        columns_are_identical = (df_data_by_roi[matching_columns].apply(lambda x: np.abs(x.max() - x.min()) < 1e-8, axis='columns')).sum() == len(df_data_by_roi)

        # If the columns are identical...
        if columns_are_identical:

            # Drop all but the first column
            num_columns_before_drop = len(df_data_by_roi.columns)
            df_data_by_roi = df_data_by_roi.drop(matching_columns[1:], axis='columns')

            # Rename the first column to the prefix less the trailing underscore
            new_first_column_name = column_prefix.rstrip('_')
            df_data_by_roi = df_data_by_roi.rename({matching_columns[0]: new_first_column_name}, axis='columns')

            # Print what we've done
            print('All {} columns starting with {} are identical; {} of them have been dropped from the dataframe and the first one has been renamed to {}'.format(len(matching_columns), column_prefix, num_columns_before_drop - len(df_data_by_roi.columns), new_first_column_name))

    # Return the processed dataframe
    return df_data_by_roi

def check_raw_weights(df_data_by_roi, log_transform_data=False):
    """Perform checks on the two potential weights values num_ann_objects_within_roi_... and footprint_integer_area_....

    The fact that the data scale roughly linearly and that I checked one ROI by eye for num_ann_objects_within_roi tells us that the footprint_integer_area is likely also correct.

    Args:
        df_data_by_roi (Pandas datagrame): Original dataframe holding the main data one ROI per row.
        log_transform_data (bool, optional): Whether to plot (not propagate) the log of the weights data; otherwise, plot the weights data as is. Defaults to False.

    Returns:
        matplotlib.pyplot figure: Figure containing the final set of scatter plots of weights data to plot.
    """

    # Note I also visually checked the raw weighting values in the dataframe
    # I also visually checked that the plot sns.scatterplot(df_data_by_roi, x='valid_area_microns_sq', y='roi_integer_area') is linear as it should be

    # Import relevant libraries
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt

    # Check that the micron areas are zero if and only if the ROI integer areas are zero
    assert ((np.abs(df_data_by_roi['valid_area_microns_sq']) < 1e-8) == (np.abs(df_data_by_roi['roi_integer_area']) < 1e-8)).sum() == len(df_data_by_roi), 'The area in microns are not zero exactly when the area in integers are zero'

    # Check that other weighting data are all equal to zero when the ROI areas are zero
    assert abs(df_data_by_roi.loc[np.abs(df_data_by_roi['valid_area_microns_sq']) < 1e-8, 'valid_area_microns_sq':].map(np.abs).sum().sum()) < 1e-8, 'Some other annotation weighting data are not zero when the ROI areas are zero'

    # Check that there are no negative values anywhere
    assert df_data_by_roi.loc[:, 'valid_area_microns_sq':].map(lambda x: x < 0).sum().sum() == 0, 'There are some negative ROI weighting data'

    # Determine the lists of unique images and annotations
    slides = df_data_by_roi['unique_slide'].unique()
    annotations = [column.removeprefix('num_ann_objects_within_roi_') for column in df_data_by_roi.columns if column.startswith('num_ann_objects_within_roi_')]

    # Initialize a matplotlib figure and axes
    fig, ax = plt.subplots(nrows=len(slides), ncols=len(annotations), figsize=tuple(np.array([6.4, 4.8]) * np.array([len(slides), len(annotations)]).max()))

    # For every slide and every annotation...
    for islide, slide in enumerate(slides):
        for iannotation, annotation in enumerate(annotations):

            # Get the current axis
            if len(slides) > 1:
                curr_ax = ax[islide][iannotation]
            else:
                curr_ax = ax[iannotation]

            # Determine the names of the current annotation-specific columns
            num_ann_objects_within_roi_col = 'num_ann_objects_within_roi_{}'.format(annotation)
            footprint_print_integer_area_col = 'footprint_integer_area_{}'.format(annotation)

            # Make a copy of the weighting data dataframe for the current annotation in case we want to perform a transform
            curr_data = df_data_by_roi.loc[df_data_by_roi['unique_slide'] == slide, [num_ann_objects_within_roi_col, footprint_print_integer_area_col]].copy()

            # If we want to log transform the data
            if log_transform_data:
                curr_data = curr_data.map(np.log10)

            # Plot and label the weighting data for the current images and annotation
            sns.scatterplot(curr_data, x=num_ann_objects_within_roi_col, y=footprint_print_integer_area_col, ax=curr_ax)
            curr_ax.set_title('Image: {}, annotation: {} ({} transform)'.format(slide, annotation, ('log' if log_transform_data else 'no')))

    # Return the figure of scatter plots
    return fig

def plot_annotation_weights_on_annotation_data(df_data_by_roi_for_selected_slide, df_annotation, selected_full_slide_name='1A-imagenum_11431', selected_annotation='entire_image', weight_column_prefix='footprint_integer_area_', do_log_transform=False, scatterpoint_alpha=1.0, roi_weight_alpha=0.4, buffer_frac=0.05, scatterpoint_size=1, overplot_column_values=False, figsize=(20, 20), ax=None):

    # Import relevant libraries
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    import time_cell_interaction_lib as tci

    # Make a copy of the input dataframe so nothing is being defined in a view
    df_data_by_roi_for_selected_slide = df_data_by_roi_for_selected_slide.copy()

    # Determine the column name in df_data_by_roi to use as the weights, which were calculated from the annotation data
    column_for_roi_colors = '{}{}'.format(weight_column_prefix, selected_annotation)
    weights_colname = '{}_weights'.format(column_for_roi_colors)

    # Depending on whether the users selects to use log transformed data, set the transformation function
    if do_log_transform:
        transformation_func = np.log10
    else:
        transformation_func = lambda x: x

    # For the selected slide of data and for the selected annotation and column to use for the weights, apply the transform to the weights (either the identity or log)
    weights = df_data_by_roi_for_selected_slide[column_for_roi_colors].apply(transformation_func)

    # Force the weights into the range [0, 1]
    weights, vmin, vmax = standardize_weights_range(weights)

    # Add the weights to the dataframe
    df_data_by_roi_for_selected_slide[weights_colname] = weights

    # Add to the current dataframe a column containing a tuple of the minimum x- and y-coordinates, just as patches.Rectangle() requires
    df_data_by_roi_for_selected_slide['xy'] = df_data_by_roi_for_selected_slide[['x_roi', 'y_roi']].map(min).apply(tuple, axis='columns')

    # Initialize axes if they're not already sent in
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Plot the dataset for the current annotation region type
    sns.scatterplot(data=df_annotation, x='xmid', y='ymid', ax=ax, alpha=scatterpoint_alpha, s=scatterpoint_size)

    # Update the x and y limits in the scatter plot to include a buffer, nominally 5% on either end of an axis
    x_range = np.array([
        min(df_annotation['xmid'].min(), df_data_by_roi_for_selected_slide['x_range'].apply(np.min).min()),
        max(df_annotation['xmid'].max(), df_data_by_roi_for_selected_slide['x_range'].apply(np.max).max())
    ])
    y_range = np.array([
        min(df_annotation['ymid'].min(), df_data_by_roi_for_selected_slide['y_range'].apply(np.min).min()),
        max(df_annotation['ymid'].max(), df_data_by_roi_for_selected_slide['y_range'].apply(np.max).max())
    ])
    x_buffer = (x_range[1] - x_range[0]) * buffer_frac
    x_range[0] = x_range[0] - x_buffer
    x_range[1] = x_range[1] + x_buffer
    y_buffer = (y_range[1] - y_range[0]) * buffer_frac
    y_range[0] = y_range[0] - y_buffer
    y_range[1] = y_range[1] + y_buffer

    # Set some axis properties
    ax.set_xlim(x_range)
    ax.set_ylim(y_range)
    ax.invert_yaxis()
    ax.set_title('{}: {}{}'.format(selected_full_slide_name, ('log ' if do_log_transform else ''), column_for_roi_colors))
    ax.set_xlabel('x-coordinates (microns)')
    ax.set_ylabel('y-coordinates (microns)')
    ax.set_aspect('equal')

    # For each ROI...
    for _, row_data in df_data_by_roi_for_selected_slide.iterrows():

        # Get the scalar data value to plot
        roi_weight = row_data[weights_colname]

        # Overwrite the facecolor of the ROI to that corresponding to the datapoint, the minimum/maximum values of the colormap extremes, and the colormap name, appending the desired transparency to the final tuple
        facecolor = tci.get_properly_interpolated_rgba_from_data(val_to_map=roi_weight, vmin=vmin, vmax=vmax, colormap_name='rocket', print_byte_value=False)[:3] + (roi_weight_alpha,)

        # Plot the ROI outlines and the corresponding weight colors and optionally overplot the actual value being used
        ax.add_patch(patches.Rectangle(row_data['xy'], row_data['width'], row_data['height'], linewidth=1, edgecolor='k', facecolor=facecolor))
        if overplot_column_values:
            xmid = row_data['width'] / 2 + row_data['x_range'][0]
            ymid = row_data['height'] / 2 + row_data['y_range'][0]
            ax.text(x=xmid, y=ymid, s=int(row_data[column_for_roi_colors]), horizontalalignment='center', verticalalignment='center')

    # Return the final figure if an axis has not been sent in; otherwise, just return the axis
    if ax is None:
        return fig, df_data_by_roi_for_selected_slide

def plot_and_save_weight_permutations(df_data_by_roi_for_selected_slide, df_annotation, selected_annotation, dpi=150, single_plot=True, multiple_plot_figsize=(10, 10), multiple_plot_dpi=200):

    # Import relevant libraries
    import numpy as np
    import matplotlib.pyplot as plt

    # Constant
    four_k = (3840, 2160)  # px
    sharpness_factor = 2  # a value of 1 will make the final image four_k pixels in size. A larger value will make it proportionally larger. If we want the final image to be 4k pixels in size, then we can use e.g. Gimp to scale it down to the right number of pixels, increasing the ppi proportionally.

    # Nested lists holding the permutations of the possible types of weights
    weight_column_prefixes = ['num_ann_objects_within_roi_', 'footprint_integer_area_']
    do_log_transforms = [False, True]

    # Calculate the full figure size in inches
    full_figsize = np.array(four_k) / dpi  # this tuple basically determines the size of the text relative to the axes and should be set accordingly. The sharpness factor is what to increase to make it look better on the screen. 

    # Determine from the dataframe the current full slide name
    selected_full_slide_name = df_data_by_roi_for_selected_slide['unique_slide'].iloc[0]

    # Initialize the figure(s) and axes
    if single_plot:
        fig, ax = plt.subplots(nrows=2, ncols=2, figsize=full_figsize)
    else:
        fig_holder = []
        ax_holder = []
        weight_column_prefix_holder = []
        do_log_transform_holder = []
        for icol, weight_column_prefix in enumerate(weight_column_prefixes):
            for irow, do_log_transform in enumerate(do_log_transforms):
                curr_fig, curr_ax = plt.subplots(figsize=multiple_plot_figsize)
                fig_holder.append(curr_fig)
                ax_holder.append(curr_ax)
                weight_column_prefix_holder.append(weight_column_prefix)
                do_log_transform_holder.append(do_log_transform)

    # For each of two rows and columns...
    if not single_plot:
        ifigax = 0
    for icol, weight_column_prefix in enumerate(weight_column_prefixes):
        for irow, do_log_transform in enumerate(do_log_transforms):

            # Determine the current axis handle depending on whether a single plot with four subfigures is to be displayed or if four different figures are to be displayed
            if single_plot:
                curr_ax = ax[irow][icol]
            else:
                curr_ax = ax_holder[ifigax]

            # Plot on the appropriate axes the current permutation of the possible weight types
            plot_annotation_weights_on_annotation_data(df_data_by_roi_for_selected_slide, df_annotation, selected_full_slide_name=selected_full_slide_name, selected_annotation=selected_annotation, weight_column_prefix=weight_column_prefix, do_log_transform=do_log_transform, ax=curr_ax)

            # Increment the axes index if multiple figures are to be displayed
            if not single_plot:
                ifigax = ifigax + 1

    # Return the final figure(s) and scaled DPI value
    if single_plot:
        return fig, sharpness_factor * dpi
    else:
        return fig_holder, multiple_plot_dpi, weight_column_prefix_holder, do_log_transform_holder

def standardize_weights_range(weights):

    # Import relevant library
    import numpy as np

    # From the input weights determine the minima and maxima from the finite values
    are_weights_finite = weights.apply(np.isfinite)
    vmin = weights[are_weights_finite].min()
    vmax = weights[are_weights_finite].max()

    # Standardize the weights range accordingly
    if vmin == vmax:
        weights.loc[:] = 0.5
    else:
        weights = (weights - vmin) / (vmax - vmin)

    # Set the new min/max range
    vmin, vmax = 0, 1

    # Set to zero the weights that we don't want to use
    weights[~are_weights_finite] = 0

    # Return relevant values
    return weights, vmin, vmax

def average_pvals_over_annotation_new(df_data_by_roi, df_density_pvals_arrays, equal_weighting_for_all_rois=False, log_pval_range=(-50, 0)):

    # Import relevant libraries
    import numpy as np
    import pandas as pd
    # import sys

    # Constant
    tol = 1e-8

    # Combine the two main dataframes into one
    df_joined = df_data_by_roi.set_index('unique_roi').join(df_density_pvals_arrays.set_index('roi_name'))

    # Determine the number of possible species in the dataset
    # nspecies = df_joined.iloc[0]['log_dens_pvals_arr'].shape[0]
    srs_tmp = df_joined['log_dens_pvals_arr']
    nspecies = srs_tmp[srs_tmp.apply(lambda x: type(x) == np.ndarray)].iloc[0].shape[0]

    # Structures holding iterables of all possible parameters
    unique_slides = df_joined['unique_slide'].unique()
    annotation_types = [x.removeprefix('num_ann_objects_within_roi_') for x in df_joined.columns if x.startswith('num_ann_objects_within_roi_')]
    weight_column_prefixes = ['num_ann_objects_within_roi_', 'footprint_integer_area_']
    do_log_transforms = [False, True]

    # Initialize a list to hold all parameters and calculated quantities
    data_holder = []

    # Iterate over all possible parameters
    for unique_slide in unique_slides:
        for selected_annotation in annotation_types:
            for weight_column_prefix in weight_column_prefixes:
                for do_log_transform in do_log_transforms:

                    # Obtain a fresh copy of the original joined dataframe for the current slide
                    df_joined_slide = df_joined[df_joined['unique_slide'] == unique_slide].copy()

                    # Determine the column name in df_data_by_roi to use as the weights, which were calculated from the annotation data
                    column_for_roi_colors = '{}{}'.format(weight_column_prefix, selected_annotation)

                    # Depending on whether the users selects to use log transformed data, set the transformation function
                    if do_log_transform:
                        transformation_func = np.log10
                    else:
                        transformation_func = lambda x: x

                    # For the selected slide of data and for the selected annotation and column to use for the weights, apply the transform to the weights (either the identity or log)
                    weights = df_joined_slide[column_for_roi_colors].apply(transformation_func)

                    # If there's at least 1 ROI with a valid weight...
                    if weights.notna().sum() >= 1:

                        # Force the weights into the range [0, 1]
                        weights, vmin, vmax = standardize_weights_range(weights)

                        # Rename the weights column to "weights"
                        weights = weights.rename('weights')

                        # Run a check on the weights range
                        if not (np.array([weights.min(), weights.max(), vmin, vmax]) == [0, 1, 0, 1]).all():
                            print('WARNING: weights.min() and vmin are not both zero and/or weights.max() and vmax are not both 1; they are {} for parameter combination {}-{}-{}-{}'.format([weights.min(), weights.max(), vmin, vmax], unique_slide, selected_annotation, weight_column_prefix, do_log_transform))
                            print('Number of values in weights that are *not* NaN: {}'.format(df_joined_slide[column_for_roi_colors].apply(transformation_func).notna().sum()))
                            # sys.exit()

                        # Append the weights as a column to the dataframe for the current slide
                        df_joined_slide = pd.concat([df_joined_slide, weights], axis='columns')

                        # Initialize the intermediate data holder arrays
                        sum_holder = np.zeros((nspecies, nspecies, 2, 1))
                        num_not_nan_holder = np.zeros((nspecies, nspecies, 2, 1), dtype=int)
                        weight_holder = np.zeros((nspecies, nspecies, 2, 1))

                        # For every ROI in the slide...
                        for roi_index, row_data in df_joined_slide.iterrows():

                            # Get log_dens_pvals_arr for the current ROI
                            log_dens_pvals_arr = row_data['log_dens_pvals_arr']

                            # If there are valid density P value data for the current ROI...
                            if type(log_dens_pvals_arr) != float:  # if it's a float then I believe np.isnan(log_dens_pvals_arr). Otherwise, they are np.ndarrays of shape (3, 3, 2, 1) always

                                # ...well the implementation in Squidpy may return an array with all nans I believe, so the type()... check above isn't a general check for there being valid density P value data
                                if np.isnan(log_dens_pvals_arr).sum() == np.prod(log_dens_pvals_arr.shape):
                                    warning_str = 'WARNING: There is no valid density P value data for ROI {}. There may be other ROIs with no valid density P value data, but this ROI failed this particular check.'.format(roi_index)
                                    print(warning_str)
                                    import os
                                    logs_dir = os.path.join('.', 'output', 'logs')
                                    if not os.path.exists(logs_dir):
                                        os.mkdir(logs_dir)
                                    with open(os.path.join(logs_dir, 'squidpy.log'), 'w') as f:
                                        f.write('This is *probably* Squidpy-specific: {}\n'.format(warning_str))

                                # log_dens_pvals_arr does have at least one valid value...
                                else:

                                    # Get the density P values and weight for the current ROI
                                    if not equal_weighting_for_all_rois:
                                        weight = row_data['weights']
                                    else:
                                        weight = 1

                                    # Run some bounds checks on the density P values
                                    nanmin = np.nanmin(log_dens_pvals_arr)
                                    nanmax = np.nanmax(log_dens_pvals_arr)
                                    assert nanmin >= log_pval_range[0] - tol, 'ERROR: The log of the density P values can be < {} (e.g., {})'.format(log_pval_range[0], nanmin)
                                    assert nanmax <= log_pval_range[1] + tol, 'ERROR: The log of the density P values can be > {} (e.g., {})'.format(log_pval_range[1], nanmax)

                                    # For every center type, neighbor type, and P value type...
                                    for icenter in range(nspecies):
                                        for ineighbor in range(nspecies):
                                            for ileft_right in range(2):

                                                # Store the current P value
                                                element = log_dens_pvals_arr[icenter, ineighbor, ileft_right, 0]

                                                # If the current P value is not NaN (as I believe occurs e.g. when there are not enough valid centers)...
                                                if not np.isnan(element):

                                                    # Add to the running sums of the weighted P values, number of not-NaNs, and weights
                                                    sum_holder[icenter, ineighbor, ileft_right, 0] = sum_holder[icenter, ineighbor, ileft_right, 0] + weight * element
                                                    num_not_nan_holder[icenter, ineighbor, ileft_right, 0] = num_not_nan_holder[icenter, ineighbor, ileft_right, 0] + 1
                                                    weight_holder[icenter, ineighbor, ileft_right, 0] = weight_holder[icenter, ineighbor, ileft_right, 0] + weight

                        # Calculated the final weighted average
                        weighted_average = sum_holder / weight_holder

                        # Run some bounds checks on the weighted averaged density P values
                        nanmin = np.nanmin(weighted_average)
                        nanmax = np.nanmax(weighted_average)
                        assert nanmin >= log_pval_range[0] - tol, 'ERROR: The log of the averaged density P values can be < {} (e.g., {})'.format(log_pval_range[0], nanmin)
                        assert nanmax <= log_pval_range[1] + tol, 'ERROR: The log of the averaged density P values can be > {} (e.g., {})'.format(log_pval_range[1], nanmax)

                    # Note what happens when there are no valid weights data
                    else:
                        print('NOTE: Parameter combination {}-{}-{}-{} has no valid weights'.format(unique_slide, selected_annotation, weight_column_prefix, do_log_transform))
                        sum_holder, num_not_nan_holder, weight_holder, weighted_average = None, None, None, None

                    # Add the current set of parameters and derived values to the running list of data
                    data_holder.append({'slide': unique_slide, 'annotation': selected_annotation, 'weight_column_prefix': weight_column_prefix, 'do_log_transform': do_log_transform, 'sum_holder': sum_holder, 'num_not_nan_holder': num_not_nan_holder, 'weight_holder': weight_holder, 'weighted_average': weighted_average})

    # Return a dataframe version of the running list of data
    return pd.DataFrame(data_holder)

def average_over_rois_per_annotation_region_for_all_slides_and_annotations(df_data_by_roi, df_density_pvals_arrays, annotations_csv_files, phenotyping_method, phenotype_identification_file, marker_column_names_str='[\'NOS2\', \'COX2\', \'CD8\']', marker_column_names_list=[], annotation_coord_units_in_microns=0.325, alpha=0.4, axis_buffer_frac=0.05, figsize=(10.0, 10.0), annotation_microns_per_integer_unit=0.325, settings__analysis__thickness=40, save_figures=False, also_save_pixel_plot=False, equal_weighting_for_all_rois=False, webpage_dir='', pixel_plot_dpi=200, pixel_plot_figsize=(10, 10), downsample=8, log_pval_range=(-50, 0)):

    # Note other main data structure in time_cell_interaction_lib.py is slices.data

    # Import relevant libraries
    import time_cell_interaction_lib as tci
    import matplotlib.pyplot as plt
    import os

    # Load the full slides of analysis data from df_data_by_roi. If df.data is present from the original TCI library, that should be equivalent to the result (at least to within a few objects), df_analysis_depatched
    df_analysis_depatched, _ = tci.undo_patching_overlaps_and_decompounding(recreate_data_property(df_data_by_roi))

    # Load and plot the annotation data for all slides
    df_annotations = load_annotation_data(os.path.join('.', 'input', 'annotations'), annotations_csv_files, phenotyping_method, phenotype_identification_file)
    print('Plotting overall annotation data plot')
    fig_all_annotation_data = plot_annotation_data(df_annotations, generate_single_figure=True, coord_units_in_microns=annotation_coord_units_in_microns, alpha=alpha, buffer_frac=axis_buffer_frac, figsize=figsize)
    if save_figures:
        fig_all_annotation_data.savefig(os.path.join(webpage_dir, 'all_annotation_data.png'), bbox_inches='tight')
    plt.close(fig_all_annotation_data)

    # Determine the image IDs and region types present in the annotation data
    image_ids = df_annotations['image_id'].unique()
    region_types = df_annotations['region_type'].unique()

    # For every slide and annotation region type in the annotation data...
    for image_id in image_ids:
        for iregion_type, region_type in enumerate(region_types):

            # Get the full slide name in the analysis data according to the current annotation slide
            slide_id = [curr_slide_id for curr_slide_id in df_analysis_depatched['Slide ID'].unique() if image_id == curr_slide_id.split('_')[1]][0]

            # Obtain the analysis data for just the current slide (image_id)
            df_analysis = df_analysis_depatched[df_analysis_depatched['Slide ID'] == slide_id]

            # Obtain the annotation data for just the current slide and region type
            matching_annotations = (df_annotations['image_id'] == image_id) & (df_annotations['region_type'] == region_type)
            if matching_annotations.sum() == 1:
                df_annotation_slide_region = df_annotations.loc[matching_annotations, :]  # should be either a one-row dataframe or an empty dataframe
                srs_annotation_slide_region = df_annotation_slide_region.iloc[0]
                df_annotation = full_df_to_just_coords_in_microns(srs_annotation_slide_region['dataframe'], annotation_coord_units_in_microns)
                df_annotation_length = srs_annotation_slide_region['df_length']
                df_annotation_length_marker = srs_annotation_slide_region['df_length_marker']

                # Overlay the current analysis data on top of the corresponding annotation data. Note these plots show that the units are consistent between the main analysis data and the annotations data, all of which are in microns
                print('Plotting analysis data overlaid on the annotation data for image {} and region {}'.format(image_id, region_type))
                fig_analysis_overlaid_on_annotation = overlay_on_background_data(df_analysis, df_annotation, alpha=alpha, buffer_frac=axis_buffer_frac, figsize=figsize, annotation_color_index=iregion_type, analysis_color_index=len(region_types), df_background_length=df_annotation_length, df_background_length_marker=df_annotation_length_marker)
                if save_figures:
                    savedir = os.path.join(webpage_dir, 'analysis_overlaid_on_annotation')
                    if not os.path.exists(savedir):
                        os.makedirs(savedir)
                    fig_analysis_overlaid_on_annotation.savefig(os.path.join(savedir, 'analysis_overlaid_on_annotation-{}-{}.png'.format(image_id, region_type)), bbox_inches='tight')
                plt.close(fig_analysis_overlaid_on_annotation)

                # Transform the annotation data to integers (which could be pixels if that's what they are)
                annot_int_arr, min_spacing_for_annotation, x_neg_shift_after_scaling, y_neg_shift_after_scaling = transform_annotation_data_to_integers(df_annotation, df_data_by_roi[df_data_by_roi['unique_slide'] == slide_id], annotation_microns_per_integer_unit)

                # Potentially plot the "pixel plot"
                if save_figures and also_save_pixel_plot:  # making saving the pixel plot separate and optional since this takes 30-45 sec even on a laptop
                    print('Plotting a per-pixel "image" of the annotation data for image {} and region {}'.format(image_id, region_type))
                    fig_pixel_plot, ax_pixel_plot = plt.subplots(figsize=pixel_plot_figsize)
                    ax_pixel_plot.imshow(annot_int_arr.T[::downsample, ::downsample])
                    ax_pixel_plot.set_title('{} annotation for image {}'.format(region_type, image_id))
                    savedir = os.path.join(webpage_dir, 'pixel_plot')
                    if not os.path.exists(savedir):
                        os.makedirs(savedir)
                    fig_pixel_plot.savefig(os.path.join(savedir, 'pixel_plot-{}-{}.png'.format(image_id, region_type)), bbox_inches='tight', dpi=pixel_plot_dpi)
                    plt.close(fig_pixel_plot)

                # Add to df_data_by_roi the raw weights data for the current slide/annotation
                df_weights_data = add_raw_weighting_values_to_roi_data(df_data_by_roi[df_data_by_roi['unique_slide'] == slide_id], df_annotation, settings__analysis__thickness, annot_int_arr, min_spacing_for_annotation, x_neg_shift_after_scaling, y_neg_shift_after_scaling, region_type)
                df_data_by_roi.loc[df_weights_data.index, df_weights_data.columns] = df_weights_data

            else:
                print('WARNING: No matching annotations for image {} and annotation region type {}'.format(image_id, region_type))

    # Drop the just-added duplicate columns of raw weighting data
    df_data_by_roi = drop_duplicate_columns_in_weighting_data(df_data_by_roi, column_prefixes=['valid_area_microns_sq_', 'roi_integer_area_', 'num_ann_objects_within_roi_', 'footprint_integer_area_'])

    # Plot the weights data in scatterplots for both identity and log transforms
    for do_log_transform in [False, True]:
        print('Plotting scatter plots of the raw weights with {} transform'.format(('a log' if do_log_transform else 'an identity')))
        fig_raw_weights_check = check_raw_weights(df_data_by_roi, log_transform_data=do_log_transform)
        if save_figures:
            savedir = os.path.join(webpage_dir, 'raw_weights_check')
            if not os.path.exists(savedir):
                os.makedirs(savedir)
            fig_raw_weights_check.savefig(os.path.join(savedir, 'raw_weights_check-log_transform__{}.png'.format(do_log_transform)), bbox_inches='tight')
        plt.close(fig_raw_weights_check)

    # Overplot the annotation weights on the annotation data
    for image_id in image_ids:
        for iregion_type, region_type in enumerate(region_types):

            # Get the full slide name in the analysis data according to the current annotation slide
            slide_id = [curr_slide_id for curr_slide_id in df_analysis_depatched['Slide ID'].unique() if image_id == curr_slide_id.split('_')[1]][0]

            # Obtain the annotation data for just the current slide and region type
            matching_annotations = (df_annotations['image_id'] == image_id) & (df_annotations['region_type'] == region_type)
            if matching_annotations.sum() == 1:
                df_annotation_slide_region = df_annotations.loc[matching_annotations, :]  # should be either a one-row dataframe or an empty dataframe
                srs_annotation_slide_region = df_annotation_slide_region.iloc[0]
                df_annotation = full_df_to_just_coords_in_microns(srs_annotation_slide_region['dataframe'], annotation_coord_units_in_microns)

                # Overplot heatmaps of the weights by ROI on top of the slide data of annotations
                print('Plotting four total heatmaps of the annotation weights overlaid on the annotations for image {} and region {}'.format(image_id, region_type))
                fig_weight_heatmaps_on_annot, heatmap_dpi, weight_column_prefix_holder, do_log_transform_holder = plot_and_save_weight_permutations(df_data_by_roi[df_data_by_roi['unique_slide'] == slide_id], df_annotation, region_type, single_plot=False)
                if save_figures:
                    savedir = os.path.join(webpage_dir, 'weight_heatmaps_on_annot')
                    if not os.path.exists(savedir):
                        os.makedirs(savedir)
                    for curr_fig, curr_weight_column_prefix, curr_do_log_transform in zip(fig_weight_heatmaps_on_annot, weight_column_prefix_holder, do_log_transform_holder):
                        curr_fig.savefig(os.path.join(savedir, 'annotation_weights_on_annotation_data-{}-{}-{}-{}.png'.format(image_id, region_type, curr_weight_column_prefix, curr_do_log_transform)), dpi=heatmap_dpi, bbox_inches='tight')  # removing bbox_inches='tight' will make everything the expected number of pixels. Adding bbox_inches='tight' will make the number of pixels appear arbitrary and it would be very difficult to make sense of scales and such
                for curr_fig in fig_weight_heatmaps_on_annot:
                    plt.close(curr_fig)

    # Average the density P values over annotation region type
    df_heatmaps_averaged_over_region_types = average_pvals_over_annotation_new(df_data_by_roi, df_density_pvals_arrays, equal_weighting_for_all_rois=equal_weighting_for_all_rois, log_pval_range=log_pval_range)

    # Return the averaged data and the new df_data_by_roi dataframe containing the raw (not scaled) weights
    return df_heatmaps_averaged_over_region_types, df_data_by_roi

def get_annotation_plots_paths(top_plot_dir='../results/webpage/slices_1x40/real'):

    # Import relevant libraries
    import glob
    import os
    import pandas as pd

    # Initialize a list that will hold lists holding the paths and their metadata
    list_paths = []

    # Raw weights scatter plots
    glob_result = glob.glob(os.path.join(top_plot_dir, 'raw_weights_check', 'raw_weights_check-log_transform__*.png'))
    list_paths.append([{'plot_type': 'raw_weights_scatter_plots',
                        'image_id': None,
                        'region_type': None,
                        'weight_column_prefix': None,
                        'do_log_transform': path.split(os.path.sep)[-1].removesuffix('.png').split('__')[-1] == 'True',
                        'path': path} for path in glob_result])

    # All annotation data
    glob_result = glob.glob(os.path.join(top_plot_dir, 'all_annotation_data.png'))
    list_paths.append([{'plot_type': 'all_annotation_data',
                        'image_id': None,
                        'region_type': None,
                        'weight_column_prefix': None,
                        'do_log_transform': None,
                        'path': path} for path in glob_result])

    # Annotation weights heatmaps overlaid on annotation data
    glob_result = glob.glob(os.path.join(top_plot_dir, 'weight_heatmaps_on_annot', 'annotation_weights_on_annotation_data-*-*-*-*.png'))
    list_paths.append([{'plot_type': 'weights_heatmap_overlay',
                        'image_id': path.split(os.path.sep)[-1].removesuffix('.png').split('-')[1],
                        'region_type': path.split(os.path.sep)[-1].removesuffix('.png').split('-')[2],
                        'weight_column_prefix': path.split(os.path.sep)[-1].removesuffix('.png').split('-')[3],
                        'do_log_transform': path.split(os.path.sep)[-1].removesuffix('.png').split('-')[4] == 'True',
                        'path': path} for path in glob_result])

    # Per-pixel "images" of annotation data
    glob_result = glob.glob(os.path.join(top_plot_dir, 'pixel_plot', 'pixel_plot-*-*.png'))
    list_paths.append([{'plot_type': 'pixel_plot',
                        'image_id': path.split(os.path.sep)[-1].removesuffix('.png').split('-')[1],
                        'region_type': path.split(os.path.sep)[-1].removesuffix('.png').split('-')[2],
                        'weight_column_prefix': None,
                        'do_log_transform': None,
                        'path': path} for path in glob_result])

    # Analysis data overlaid on annotation data
    glob_result = glob.glob(os.path.join(top_plot_dir, 'analysis_overlaid_on_annotation', 'analysis_overlaid_on_annotation-*-*.png'))
    list_paths.append([{'plot_type': 'analysis_overlay',
                        'image_id': path.split(os.path.sep)[-1].removesuffix('.png').split('-')[1],
                        'region_type': path.split(os.path.sep)[-1].removesuffix('.png').split('-')[2],
                        'weight_column_prefix': None,
                        'do_log_transform': None,
                        'path': path} for path in glob_result])

    # Average density P value heatmaps
    glob_result = glob.glob(os.path.join(top_plot_dir, 'dens_pvals_per_annotation', 'density_pvals-real-*-*-*-*-slice_01_of_01-annotation_index_-1.png'))
    list_paths.append([{'plot_type': 'average_p_value_heatmap',
                        'image_id': path.split(os.path.sep)[-1].removesuffix('.png').split('-')[3].split('_')[-1],
                        'region_type': path.split(os.path.sep)[-1].removesuffix('.png').split('-')[4],
                        'weight_column_prefix': path.split(os.path.sep)[-1].removesuffix('.png').split('-')[5],
                        'do_log_transform': path.split(os.path.sep)[-1].removesuffix('.png').split('-')[6] == 'True',
                        'path': path} for path in glob_result])

    # Flatten the list of lists
    list_paths2 = []
    for curr_list in list_paths:
        for item in curr_list:
            list_paths2.append(item)

    # Convert to an appropriate Pandas dataframe
    df_paths = pd.DataFrame(list_paths2)
    df_paths = df_paths.sort_values(list(df_paths.columns[:-1])).reset_index(drop=True)

    # Return the final dataframe
    return df_paths
