# Import relevant libraries
import pandas as pd
import seaborn as sns
import streamlit as st
import plotly.graph_objects as go
import os
import dataset_formats
import plotly.express as px
import utils
import matplotlib.pyplot as plt
import numpy as np

# Import relevant libraries
import app_top_of_page as top
import streamlit_dataframe_editor as sde

# Function to load the data in a unified format
def load_data(input_datafile_path, coord_units_in_microns, dataset_format):
    transform = {'HALO': 'OMAL'}
    dataset_class = getattr(dataset_formats, (transform[dataset_format] if dataset_format in transform else dataset_format))  # done this way so that the format (e.g., “REEC”) can be select programmatically
    dataset_obj = dataset_class(input_datafile=input_datafile_path, coord_units_in_microns=coord_units_in_microns)
    dataset_obj.process_dataset(do_calculate_minimum_coordinate_spacing_per_roi=False, do_trimming=False, do_extra_processing=False)
    return dataset_obj.data

# Update the dependencies of the selectbox for the current analysis column
def update_dependencies_of_filtering_widgets():
    df = st.session_state['mg__df']
    column_for_filtering = st.session_state['mg__column_for_filtering']
    image_for_filtering = st.session_state['mg__selected_image']
    if image_for_filtering == 'All images':
        image_loc = df.index
    else:
        image_loc = df[df['Slide ID'] == image_for_filtering].index
    curr_series = df.loc[image_loc, column_for_filtering]
    if column_for_filtering in st.session_state['mg__all_numeric_columns']:
        st.session_state['mg__selected_column_type'] = 'numeric'
        st.session_state['mg__curr_column_range'] = (curr_series.min(), curr_series.max())
        st.session_state['mg__selected_value_range'] = st.session_state['mg__curr_column_range']  # initialize the selected range to the entire range
    else:
        st.session_state['mg__selected_column_type'] = 'categorical'
        st.session_state['mg__curr_column_unique_values'] = curr_series.unique()
        st.session_state['mg__selected_column_values'] = []  # initialize the selected values to nothing

# Set the column options to all columns of the desired type less columns that have already been used for filtering for the current phenotype
def update_column_options(key_for_column_list='mg__all_columns'):
    # old default: key_for_column_list='mg__all_numeric_columns'
    return [column for column in st.session_state[key_for_column_list] if column in (set(st.session_state[key_for_column_list]) - set(st.session_state['mg__de_current_phenotype'].reconstruct_edited_dataframe()['Column for filtering']))]

# Add to the current phenotype and update the previous parts of the app
def update_dependencies_of_button_for_adding_column_filter_to_current_phenotype(column_for_filtering, selected_min_val=None, selected_max_val=None, selected_column_values=None):

    # Get the current value of the current phenotype dataframe editor
    df_current_phenotype = st.session_state['mg__de_current_phenotype'].reconstruct_edited_dataframe()

    # Add the selected column filter to the current phenotype assignments dataframe and update the phenotype assignments dataframe with this new dataframe
    dataframe_to_add = pd.DataFrame(pd.Series({'Column for filtering': column_for_filtering, 'Minimum value': selected_min_val, 'Maximum value': selected_max_val, 'Selected values': selected_column_values}))
    new_df_contents = pd.concat([df_current_phenotype, dataframe_to_add.T]).reset_index(drop=True)
    st.session_state['mg__de_current_phenotype'].update_editor_contents(new_df_contents=new_df_contents)

    # Set the currently selected column as the first of the possible options
    st.session_state['mg__column_for_filtering'] = update_column_options(key_for_column_list='mg__all_numeric_columns')[0]

    # Since the column selection must have just changed, update its dependencies
    update_dependencies_of_filtering_widgets()

# Add to the phenotype assignments for the new dataset and update the previous parts of the app
def update_dependencies_of_button_for_adding_phenotype_to_new_dataset():

    if not st.session_state['mg__df_current_phenotype'].empty:
        # Get the current values of the two dataframe editors
        df_current_phenotype = st.session_state['mg__de_current_phenotype'].reconstruct_edited_dataframe()
        df_phenotype_assignments = st.session_state['mg__de_phenotype_assignments'].reconstruct_edited_dataframe()
        if 'Phenotype' in df_phenotype_assignments.columns:
            df_phenotype_assignments = df_phenotype_assignments.set_index('Phenotype')

        # Populate a dictionary of the column filters for the current phenotype
        curr_phenotype_dict = dict()
        for row in df_current_phenotype.itertuples(index=False):
            curr_col, curr_min, curr_max, curr_items = row
            if curr_items is None:  # numeric column filter
                curr_phenotype_dict[curr_col + ' [[min]]'] = curr_min
                curr_phenotype_dict[curr_col + ' [[max]]'] = curr_max
            else:  # categorical column filter
                curr_phenotype_dict[curr_col + ' [[items]]'] = curr_items

        # Update the contents of the phenotype assignments data editor
        dataframe_to_add = pd.DataFrame(pd.Series(curr_phenotype_dict, name=st.session_state['mg__current_phenotype_name']))
        new_df_contents = pd.concat([df_phenotype_assignments, dataframe_to_add.T]).rename_axis('Phenotype').reset_index(drop=False)
        st.session_state['mg__de_phenotype_assignments'].update_editor_contents(new_df_contents=new_df_contents)

        # Clear the current phenotype dataframe editor to its default value
        st.session_state['mg__de_current_phenotype'].reset_dataframe_content()

        # Set the currently selected column as the first of the possible options
        st.session_state['mg__column_for_filtering'] = update_column_options(key_for_column_list='mg__all_numeric_columns')[0]

        # Since the column selection must have just changed, update its dependencies
        update_dependencies_of_filtering_widgets()

# From the phenotype assignments, add one column per phenotype to the original dataframe containing pluses where all the phenotype criteria are met
def add_new_phenotypes_to_main_df(df, image_for_filtering):

    # Import relevant library
    from datetime import datetime

    if not st.session_state['mg__df_phenotype_assignments'].empty:

        # Reassign the *input* dataframe
        if image_for_filtering == 'All images':
            image_loc = df.index
            filtering_section_name = 'all_images'
        else:
            image_loc = df[df['Slide ID'] == image_for_filtering].index
            filtering_section_name = 'image_{}'.format(image_for_filtering)
        df = df.loc[image_loc, :]

        # Get the current values of the phenotype assignments data editor
        df_phenotype_assignments = st.session_state['mg__de_phenotype_assignments'].reconstruct_edited_dataframe().set_index('Phenotype')

        # Debugging output
        print('-- Phenotype criteria --')
        print()

        # For each set of phenotype assignments...
        for iphenotype, row in enumerate(df_phenotype_assignments.itertuples()):

            # Obtain the name of the current phenotype
            phenotype = row[0]

            # Debugging output
            print('Phenotype #{}: {}'.format(iphenotype + 1, phenotype))

            # Create a dataframe from the current row (i.e., phenotype) so that we can split and add columns
            curr_df = pd.Series(dict(zip(df_phenotype_assignments.columns, row[1:])), name=phenotype).dropna().to_frame().reset_index()

            # Add a "column" column containing the name of the filtering column
            curr_df['column'] = [x.split(' [[')[0] for x in curr_df['index']]

            # Drop the unnecessary "index" column
            curr_df = curr_df.drop(['index'], axis='columns')

            # Initialize a Series of booleans to True
            phenotype_bools = pd.Series([True] * len(df), index=df.index)

            # For each filtering column...
            object_count_holder = []
            for ifilter_col, (column, values) in enumerate(curr_df.groupby(by='column')[phenotype].apply(list).items()):  # note the groupby appears to result in the column's order of min then max

                # Set the booleans if the boolean column filter is numeric (values = [min, max])
                if len(values) == 2:
                    column_filter_bools = (df[column] >= values[0]) & (df[column] <= values[1])
                    print('  Column #{} (numeric): {}'.format(ifilter_col + 1, column))
                    print('    min: {}'.format(values[0]))
                    print('    max: {}'.format(values[1]))

                # Set the booleans if the boolean column filter is categorical (values = [items])
                else:
                    column_filter_bools = df[column].apply(lambda x: x in values[0])
                    print('  Column #{} (categorical): {}'.format(ifilter_col + 1, column))
                    print('    items: {}'.format(values[0]))

                # Debugging output
                curr_filter_column_count = column_filter_bools.sum()
                object_count_holder.append(curr_filter_column_count)
                print('    object count: {}'.format(curr_filter_column_count))

                # Determine where the current column values are within the specified range criterion
                phenotype_bools = phenotype_bools & column_filter_bools

            # Debugging output
            curr_phenotype_count = phenotype_bools.sum()
            assert curr_phenotype_count <= min(object_count_holder), 'ERROR: The object count for the total phenotype must be smaller than the smallest object count for its individual column filters'
            print('  Phenotype object count: {}'.format(curr_phenotype_count))

            # Add a column to the original dataframe with the new phenotype satisfying all of its filtering criteria
            st.session_state['mg__df'].loc[image_loc, 'Phenotype {}'.format(phenotype)] = phenotype_bools.apply(lambda x: ('+' if x else '-'))

        # Debugging output
        print('------------------------')

        # Save the gating table to disk
        gating_filename = 'gating_table_for_{}_for_datafile_{}-{}.csv'.format(filtering_section_name, st.session_state['mg__input_datafile_filename'], datetime.now().strftime("date%Y_%m_%d_time%H_%M_%S"))
        df_phenotype_assignments.to_csv(path_or_buf=os.path.join(os.path.join('.', 'output'), gating_filename), index=True)
        st.write('File {} written to disk'.format(gating_filename))

# Function to clear the session state as would be desired when loading a new dataset
def clear_session_state(keep_keys=[]):
    for key in (set(st.session_state.keys()) - set(keep_keys)):
        del st.session_state[key]

# Callback to run every time the intensity column multiselect is modified
def update_field_matching():

    # Get the current state of the multiselect
    selected_intensity_fields = st.session_state['mg__selected_intensity_fields']

    # Get the current state of the intensity-marker-matching dataframe
    curr_df = st.session_state['mg__de_field_matching'].reconstruct_edited_dataframe()

    # Create a new "blank" dataframe corresponding to the currently selected intensity columns in the multiselect
    new_df = pd.DataFrame({'Intensity field': selected_intensity_fields, 'Corresponding thresholded marker field': ['Select thresholded marker field 🔽' for _ in selected_intensity_fields]})

    # For every intensity field in the current matching dataframe...
    for intensity_field in curr_df['Intensity field']:

        # If the current intensity field is present in the new, multiselect-based dataframe...
        if intensity_field in new_df['Intensity field'].to_list():

            # Update the new dataframe with the corresponding existing value in the current dataframe
            curr_df_index = curr_df[curr_df['Intensity field'] == intensity_field].index
            new_df_index = new_df[new_df['Intensity field'] == intensity_field].index
            new_df.loc[new_df_index, 'Corresponding thresholded marker field'] = curr_df.loc[curr_df_index, 'Corresponding thresholded marker field']

    # Update the matching dataframe with the merged, "new" dataframe
    st.session_state['mg__de_field_matching'].update_editor_contents(new_df_contents=new_df)

# Callback particularly for the current phenotype definition in case that dataframe is edited, in particular one of the filter column names
# Note I can call this multiple places in this script because its contents are repeated many places
def basic_filter_column_updates():

    # Set the currently selected column as the first of the possible options
    st.session_state['mg__column_for_filtering'] = update_column_options(key_for_column_list='mg__all_numeric_columns')[0]

    # Since the column selection must have just changed, update its dependencies
    update_dependencies_of_filtering_widgets()

# Delete specified columns in the main dataframe
def delete_all_gated_phenotypes(new_phenotypes=[]):
    st.session_state['mg__df'] = st.session_state['mg__df'].drop(columns=new_phenotypes)

# Perform simple Z score normalization
def z_score_normalize(df, numeric_columns):

    # Copy the input dataframe as the output dataframe; this would be the whole function if no batch normalization were selected; this is essentially the identity transformation
    df_batch_normalized = df.copy()

    # Get the unique images in the dataframe
    unique_images = df['Slide ID'].unique()

    # For each image in the dataframe...
    for image_name in unique_images:

        # Get the locations of the data for the current image (this results in a boolean series)
        image_loc = df['Slide ID'] == image_name

        # Print what we're doing
        print('Z score normalizing image {} ({} rows)...'.format(image_name, image_loc.sum()))

        # Get just the numeric data for the current image
        curr_df_numeric = df.loc[image_loc, numeric_columns]

        # Z score normalize these data column by column, assigning the results to the output dataframe at the corresponding locations
        df_batch_normalized.loc[image_loc, numeric_columns] = (curr_df_numeric - curr_df_numeric.mean()) / curr_df_numeric.std()

    # Return the transformed dataset
    return df_batch_normalized

def main():
    '''
    Main function for running the page
    '''

    # Set page settings
    st.set_page_config(layout='wide', page_title='Multiaxial Gating')
    st.title('Multi-axial Gating')

    # Run streamlit-dataframe-editor library initialization tasks at the top of the page
    st.session_state = sde.initialize_session_state(st.session_state)

    # Run Top of Page (TOP) functions
    st.session_state = top.top_of_page_reqs(st.session_state)

    # Set the default dataframes to be edited
    default_df_current_phenotype     = pd.DataFrame(columns=['Column for filtering', 'Minimum value', 'Maximum value'])
    default_df_phenotype_assignments = pd.DataFrame()
    default_df_field_matching        = pd.DataFrame(columns=['Intensity field', 'Corresponding thresholded marker field'])

    # Constants
    input_directory = os.path.join('.', 'input')
    num_categorical_values_cutoff = 10
    tol = 1e-8

    # Set the options for input data filenames
    options_for_input_datafiles = [x for x in os.listdir(input_directory) if x.endswith(('.csv', '.tsv'))]

    # Initialize some things in the session state
    if 'mg__current_phenotype_name' not in st.session_state:
        st.session_state['mg__current_phenotype_name'] = ''
    if 'mg__de_current_phenotype' not in st.session_state:
        st.session_state['mg__de_current_phenotype'] = sde.DataframeEditor(df_name='mg__df_current_phenotype', default_df_contents=default_df_current_phenotype)
    if 'mg__de_phenotype_assignments' not in st.session_state:
        st.session_state['mg__de_phenotype_assignments'] = sde.DataframeEditor(df_name='mg__df_phenotype_assignments', default_df_contents=default_df_phenotype_assignments)
    if 'mg__de_field_matching' not in st.session_state:
        st.session_state['mg__de_field_matching'] = sde.DataframeEditor(df_name='mg__df_field_matching', default_df_contents=default_df_field_matching)
    if 'mg__input_datafile_filename' not in st.session_state:
        st.session_state['mg__input_datafile_filename'] = utils.get_first_element_or_none(options_for_input_datafiles)
    if 'mg__input_datafile_coordinate_units' not in st.session_state:
        st.session_state['mg__input_datafile_coordinate_units'] = 0.25
    if 'mg__do_batch_norm' not in st.session_state:
        st.session_state['mg__do_batch_norm'] = False
    if 'mg__selected_intensity_fields' not in st.session_state:
        st.session_state['mg__selected_intensity_fields'] = []
    if 'mg__selected_image' not in st.session_state:
        st.session_state['mg__selected_image'] = 'All images'

    # Create columns for the selections in the top part of the app
    topmost_input_columns = st.columns(2)

    # In the left half, create the data loading options and spinner
    with topmost_input_columns[0]:

        # Create columns for the data loading options (just one column is used)
        load_data_columns = st.columns([3, 1])

        # Set the load data options
        with load_data_columns[0]:
            st.selectbox('Filename:', options_for_input_datafiles, key='mg__input_datafile_filename', help='Input datafiles must be present in the "input" directory and have a .csv or .tsv extension.')
            input_datafilename = st.session_state['mg__input_datafile_filename']
            st.number_input('x-y coordinate units (microns):', min_value=0.0, key='mg__input_datafile_coordinate_units', help='E.g., if the coordinates in the input datafile were pixels, this number would be a conversion to microns in units of microns/pixel.', format='%.4f', step=0.0001)  # set the input datafile coordinate units in microns
            coord_units_in_microns = st.session_state['mg__input_datafile_coordinate_units']
            st.toggle(label='Perform batch normalization', key='mg__do_batch_norm')
            batch_normalization_func = (z_score_normalize if st.session_state['mg__do_batch_norm'] else lambda df, _: df)  # initially, create a very simple batch normalization option using the Z score, which appears to be justified in literature

        # Create columns for the data loading button and spinner
        data_butt_cols = st.columns([3, 1])

        # In the first column, create the data loading button
        with data_butt_cols[0]:
            MaG_load_hit = st.button('Load data', use_container_width=True, on_click=clear_session_state, kwargs={'keep_keys': ['mg__input_datafile_filename', 'mg__input_datafile_coordinate_units', 'mg__do_batch_norm']})

        # In the second column, create the data loading spinner
        with data_butt_cols[1]:
            if MaG_load_hit:
                with st.spinner('Loading Data'):
                    st.session_state['mg__df'] = load_data(os.path.join(input_directory, input_datafilename), coord_units_in_microns, dataset_formats.extract_datafile_metadata(os.path.join(input_directory, input_datafilename))[4])
                    unique_images = st.session_state['mg__df']['Slide ID'].unique()
                    st.session_state['mg__unique_images_short'] = [x.split('-imagenum_')[1] for x in unique_images]
                    st.session_state['mg__unique_image_dict'] = dict(zip(st.session_state['mg__unique_images_short'], unique_images))
                    phenotype_columns = [column for column in st.session_state['mg__df'].columns if column.startswith('Phenotype ')]
                    st.session_state['mg__df'] = st.session_state['mg__df'].rename(columns=dict(zip(phenotype_columns, [column.replace('Phenotype ', 'Phenotype_orig ') for column in phenotype_columns])))
                    st.session_state['mg__all_numeric_columns'] = st.session_state['mg__df'].select_dtypes(include='number').columns
                    st.session_state['mg__all_columns'] = st.session_state['mg__df'].columns
                    st.session_state['mg__column_config'] = {"Corresponding thresholded marker field": st.column_config.SelectboxColumn("Corresponding thresholded marker field", help="Tresholded marker field corresponding to the intensity at left", options=st.session_state['mg__all_columns'], required=True)}
                    st.session_state['mg__df_batch_normalized'] = batch_normalization_func(st.session_state['mg__df'], st.session_state['mg__all_numeric_columns'])
                    
        # Warn the user that they need to load the data at least once
        if 'mg__df' not in st.session_state:
            st.warning('Please load data from the selections above')

    # If the data have been loaded...
    if 'mg__df' in st.session_state:
    
        # Load the data and some resulting processed data
        df = st.session_state['mg__df']
        unique_images_short = st.session_state['mg__unique_images_short']
        unique_image_dict = st.session_state['mg__unique_image_dict']

        # In the right half, run the field matching
        with topmost_input_columns[1]:

            # Add expander, expanded by default just for the time being as sometimes otherwise it collapses unexpectedly
            with st.expander('Field matching (optional):', expanded=False):

                # Create two columns on the page
                cols_field_matching = st.columns(2)

                # In the first column, have a multiselect for the user to select raw intensity columns
                with cols_field_matching[0]:
                    st.multiselect('Select intensity fields:', st.session_state['mg__all_columns'], key='mg__selected_intensity_fields', on_change=update_field_matching)

                # In the second column, have an editable dataframe to allow the user to match thresholded marker columns to the selected raw intensity columns
                with cols_field_matching[1]:
                    st.session_state['mg__de_field_matching'].dataframe_editor(column_config=st.session_state['mg__column_config'], reset_data_editor_button=False)

        # Define the main columns
        main_columns = st.columns(3, gap='large')

        # Data column filter
        with main_columns[0]:

            # Column header
            st.header(':one: Column filter')

            # Initialize the filtering column selection
            if 'mg__column_for_filtering' not in st.session_state:
                st.session_state['mg__column_for_filtering'] = update_column_options(key_for_column_list='mg__all_numeric_columns')[0]

            # Allow the user to select either all images or just a single image
            st.selectbox('Image selection:', ['All images'] + df['Slide ID'].unique().tolist(), key='mg__selected_image', on_change=update_dependencies_of_filtering_widgets)

            # Have a dropdown for the column on which to perform a kernel density estimate
            st.selectbox(label='Column for filtering:', options=update_column_options(), key='mg__column_for_filtering', on_change=update_dependencies_of_filtering_widgets)
            column_for_filtering = st.session_state['mg__column_for_filtering']
            image_for_filtering = st.session_state['mg__selected_image']

            # Output information on the column
            if 'mg__selected_column_type' not in st.session_state:
                update_dependencies_of_filtering_widgets()
            if st.session_state['mg__selected_column_type'] == 'numeric':
                column_range = st.session_state['mg__curr_column_range']
                if np.abs(column_range[1] - column_range[0]) < tol:
                    trivial_column = True
                else:
                    trivial_column = False
            else:  # categorical column
                curr_column_unique_values = st.session_state['mg__curr_column_unique_values']
                st.write('Column\'s unique values: `{}`'.format(curr_column_unique_values))
                num_unique_values_in_curr_column = len(curr_column_unique_values)
                if num_unique_values_in_curr_column == 1:
                    trivial_column = True
                else:
                    trivial_column = False

            # If the column is trivial, don't allow it to be filtered on
            if trivial_column:
                st.info('Selected column is trivial; not available for filtering.', icon="ℹ️")

            # If the selected column contains more than one value...
            else:

                # Determine the x-y data to plot for the selected column, calculating the KDE for each column (and each image) only once ever
                if 'mg__kdes_or_hists_to_plot' not in st.session_state:
                    st.session_state['mg__kdes_or_hists_to_plot'] = pd.DataFrame(columns=st.session_state['mg__all_columns'], index=(['All images'] + df['Slide ID'].unique().tolist()), dtype='object')
                if st.session_state['mg__kdes_or_hists_to_plot'].loc[image_for_filtering, column_for_filtering] is np.nan:
                    if image_for_filtering == 'All images':
                        image_loc = df.index
                    else:
                        image_loc = df[df['Slide ID'] == image_for_filtering].index
                    if st.session_state['mg__selected_column_type'] == 'numeric':
                        line2d = sns.kdeplot(data=df.loc[image_loc, :], x=column_for_filtering).get_lines()[0]
                        curr_df = pd.DataFrame({'Value': line2d.get_xdata(), 'Density': line2d.get_ydata()})
                        st.session_state['mg__kdes_or_hists_to_plot'].loc[image_for_filtering, column_for_filtering] = [curr_df[(curr_df['Value'] >= column_range[0]) & (curr_df['Value'] <= column_range[1])]]  # needed because the KDE can extend outside the possible value range
                    else:
                        vc = df.loc[image_loc, column_for_filtering].value_counts(sort=False)
                        curr_df = pd.DataFrame({'Values': vc.index.to_list(), 'Counts': vc.to_list()})
                        st.session_state['mg__kdes_or_hists_to_plot'].loc[image_for_filtering, column_for_filtering] = [curr_df]
                kde_or_hist_to_plot_full = st.session_state['mg__kdes_or_hists_to_plot'].loc[image_for_filtering, column_for_filtering][0]

                # If the selected column is numeric...
                if st.session_state['mg__selected_column_type'] == 'numeric':

                    # Draw a range slider widget for selecting the range min and max
                    st.slider(label='Selected value range:', min_value=column_range[0], max_value=column_range[1], key='mg__selected_value_range')
                    selected_min_val, selected_max_val = st.session_state['mg__selected_value_range']

                    # Create a view of the full dataframe that is the selected subset
                    df_to_plot_selected = kde_or_hist_to_plot_full[(kde_or_hist_to_plot_full['Value'] >= selected_min_val) & (kde_or_hist_to_plot_full['Value'] <= selected_max_val)]

                    # Initialize the raw intensity cutoff number to None
                    intensity_cutoff = None

                    # Get the current field-matching dataframe (mapping intensity column to thresholded marker column)
                    srs_matched_marker_fields = st.session_state['mg__de_field_matching'].reconstruct_edited_dataframe().set_index('Intensity field').iloc[:, 0]

                    # If the selected column for performing filtering has a corresponding thresholded marker column identified...
                    if column_for_filtering in srs_matched_marker_fields.index.to_list():

                        # Read in the corresponding marker column name
                        marker_column = srs_matched_marker_fields.loc[column_for_filtering]  # string

                        # If the marker threshold column was actually set...
                        if marker_column != 'Select thresholded marker field 🔽':

                            # Get the current image indices
                            if image_for_filtering == 'All images':
                                image_loc = df.index
                            else:
                                image_loc = df[df['Slide ID'] == image_for_filtering].index

                            # Get the thresholded marker column values
                            srs_marker_column_values = df.loc[image_loc, marker_column]

                            # Set the indices of that series to the corresponding intensities
                            srs_marker_column_values.index = df.loc[image_loc, column_for_filtering]

                            # Sort the series by increasing intensity
                            srs_marker_column_values = srs_marker_column_values.sort_index()

                            # Determine whether the intensity values were deemed "positive" presumably by looking at the original image
                            if srs_marker_column_values.dtype == 'object':
                                positive_loc = srs_marker_column_values.apply(lambda x: x[-1] == '+')
                            else:
                                positive_loc = srs_marker_column_values == 1

                            # Get the lowest-intensity "positive" intensity/marker
                            intensity_cutoff = srs_marker_column_values[positive_loc].index[0]

                    # Plot the Plotly figure in Streamlit
                    fig = go.Figure()
                    fig.add_trace(go.Scatter(x=kde_or_hist_to_plot_full['Value'], y=kde_or_hist_to_plot_full['Density'], fill='tozeroy', mode='none', fillcolor='yellow', name='Full dataset', hovertemplate=' '))
                    fig.add_trace(go.Scatter(x=df_to_plot_selected['Value'], y=df_to_plot_selected['Density'], fill='tozeroy', mode='none', fillcolor='red', name='Selection', hoverinfo='skip'))
                    if intensity_cutoff is not None:
                        fig.add_vline(x=intensity_cutoff, line_color='green', line_width=3, line_dash="dash", annotation_text="Previous threshold: ~{}".format((intensity_cutoff)), annotation_font_size=18, annotation_font_color="green")
                    fig.update_layout(hovermode='x unified', xaxis_title='Column value', yaxis_title='Density')
                    fig.update_layout(legend=dict(yanchor="top", y=1.2, xanchor="left", x=0.01, orientation="h"))
                    
                    # Set Plotly chart in streamlit
                    st.plotly_chart(fig, use_container_width=True)

                    # Set the selection dictionary for the current filter to pass on to the current phenotype definition
                    selection_dict = {'column_for_filtering': column_for_filtering, 'selected_min_val': selected_min_val, 'selected_max_val': selected_max_val, 'selected_column_values': None}

                    # Whether to disable the add column button
                    add_column_button_disabled = False

                # If the selected column is categorical...
                else:

                    # If there's a tractable number of unique values in the selected column...
                    if num_unique_values_in_curr_column <= num_categorical_values_cutoff:

                        # Draw a multiselect widget for selecting the unique column values to use in the column filter
                        st.multiselect(label='Selected value items:', options=curr_column_unique_values, key='mg__selected_column_values')
                        selected_items = st.session_state['mg__selected_column_values']

                        # Plot in Streamlit the Matplotlib histogram of the possible values in the full dataset
                        # sns.histplot(data=kde_or_hist_to_plot_full, x=column_for_filtering, shrink=.8, ax=ax)
                        fig, ax = plt.subplots()
                        bar_colors = ['tab:red' if value in selected_items else 'y' for value in kde_or_hist_to_plot_full['Values']]
                        ax.bar(kde_or_hist_to_plot_full['Values'], kde_or_hist_to_plot_full['Counts'], color=bar_colors)
                        ax.set_xlabel('Value')
                        ax.set_ylabel('Count')
                        ax.set_title('Counts of values of column {}'.format(column_for_filtering))
                        ax.set_xticklabels(labels=ax.get_xticklabels(), rotation=45, ha='right')
                        st.pyplot(fig)

                        # Set the selection dictionary for the current filter to pass on to the current phenotype definition
                        selection_dict = {'column_for_filtering': column_for_filtering, 'selected_min_val': None, 'selected_max_val': None, 'selected_column_values': selected_items}

                        # Whether to disable the add column button
                        add_column_button_disabled = False

                    # If there are too many unique values in the column, do nothing
                    else:
                        st.info('There are too many unique values (more than {}) in the selected categorical column. Filtering not currently allowed.'.format(num_categorical_values_cutoff), icon="ℹ️")
                        selection_dict = dict()
                        add_column_button_disabled = True

                # Add the current column filter to the current phenotype assignment
                st.button(':star2: Add column filter to current phenotype :star2:', use_container_width=True, on_click=update_dependencies_of_button_for_adding_column_filter_to_current_phenotype, kwargs=selection_dict, disabled=add_column_button_disabled)

        # Current phenotype and phenotype assignments
        with main_columns[1]:

            # Column header
            st.header(':two: Current phenotype', help='Note you can refine non-list values in the following table by editing them directly or even deleting (or adding) whole rows.')

            # Output the dataframe holding the phenotype that's currently being built
            st.session_state['mg__de_current_phenotype'].dataframe_editor(reset_data_editor_button_text='Reset current phenotype definition', on_change=basic_filter_column_updates)

            # Choose a phenotype name
            st.text_input(label='Phenotype name:', key='mg__current_phenotype_name')

            # Add the current phenotype to the phenotype assignments table
            st.button(label=':star2: Add phenotype to assignments table :star2:',
                      use_container_width=True, 
                      on_click=update_dependencies_of_button_for_adding_phenotype_to_new_dataset)

            # Column header
            st.header(':three: Phenotype assignments', help='Note you can refine non-list values in the following table by editing them directly or even deleting whole rows.')

            # Output the dataframe holding the specifications for all phenotypes
            st.session_state['mg__de_phenotype_assignments'].dataframe_editor(reset_data_editor_button_text='Reset all phenotype definitions')

            # Generate the new dataset
            st.button(label=':star2: Append phenotype assignments to the dataset :star2:',
                      use_container_width=True,
                      on_click=add_new_phenotypes_to_main_df, args=(df, image_for_filtering))
            if image_for_filtering == 'All images':
                st.write('Clicking this button will apply the phenotype assignments to **all images in the dataset**')
            else:
                st.write('Clicking this button will apply the phenotype assignments to **just the image {}**'.format(image_for_filtering))

        # New dataset
        with main_columns[2]:

            # Column header
            st.header(':four: New dataset')

            # Print out a five-row sample of the main dataframe
            st.write('Augmented dataset sample:')
            st.dataframe(st.session_state['mg__df'].sample(5), hide_index=True)

            # Get a list of all new phenotypes
            new_phenotypes = [column for column in st.session_state['mg__df'].columns if column.startswith('Phenotype ')]

            # Print out the new phenotypes present
            st.write('There are {} gated phenotypes present in the augmented dataset: {}'.format(len(new_phenotypes), new_phenotypes))

            # If at least one phenotype has been assigned...
            if len(new_phenotypes) > 0:

                # Add an option to delete all generated phenotypes so far
                st.button('Delete all gated phenotypes', use_container_width=True, on_click=delete_all_gated_phenotypes, kwargs={'new_phenotypes': new_phenotypes})

                # Initialize the plot of an optional cell scatter plot to the first image in the dataset
                if 'mg__image_to_plot' not in st.session_state:
                    st.session_state['mg__image_to_plot'] = unique_images_short[0]
                if 'mg__phenotype_to_plot' not in st.session_state:
                    st.session_state['mg__phenotype_to_plot'] = new_phenotypes[0]

                # Generate widgets for the plotting parameters
                st.selectbox(label='Image to plot:', options=unique_images_short, key='mg__image_to_plot')
                st.selectbox(label='Phenotype to plot:', options=new_phenotypes, key='mg__phenotype_to_plot')

                # If the button is pressed
                if st.button('Plot the selected phenotype in the selected image'):
                    df_for_scatterplot = df.loc[df['Slide ID'] == unique_image_dict[st.session_state['mg__image_to_plot']], ['Cell X Position', 'Cell Y Position', st.session_state['mg__phenotype_to_plot']]]
                    fig = px.scatter(data_frame=df_for_scatterplot, x='Cell X Position', y='Cell Y Position', color=st.session_state['mg__phenotype_to_plot'])
                    fig.update_xaxes(scaleanchor='y')
                    st.plotly_chart(fig)

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)

# Call the main function
if __name__ == '__main__':
    main()
