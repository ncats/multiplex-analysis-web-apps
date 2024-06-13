# Import relevant libraries
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
import plotly.express as px
import os
import numpy as np
import utils
from scipy.stats import gaussian_kde
import app_top_of_page as top
import streamlit_dataframe_editor as sde
import random
import string
import image_filter


def generate_random_string(length=10):
    # Define the characters that will be used
    characters = string.ascii_letters + string.digits
    # Generate a random string of the specified length
    random_string = ''.join(random.choice(characters) for i in range(length))
    return random_string


def reset_x_axis_range(use_groups_for_plotting, kde_or_hist_to_plot_full):
    if not use_groups_for_plotting:
        st.session_state['mg__histogram_x_range'] = [kde_or_hist_to_plot_full['Value'].min(), kde_or_hist_to_plot_full['Value'].max()]
    else:
        st.session_state['mg__histogram_x_range'] = [min(kde_or_hist_to_plot_full[0]['Value'].min(), kde_or_hist_to_plot_full[1]['Value'].min()), max(kde_or_hist_to_plot_full[0]['Value'].max(), kde_or_hist_to_plot_full[1]['Value'].max())]
    st.session_state['mg__random_string'] = generate_random_string()


def plotly_chart_histogram_callback():
    plotly_chart_key = ('mg__plotly_chart_histogram_' + st.session_state['mg__random_string'] + '__do_not_persist')
    if plotly_chart_key in st.session_state:
        if st.session_state[plotly_chart_key]['selection']['box']:
            x_range = sorted(st.session_state[plotly_chart_key]['selection']['box'][0]['x'])
            if st.session_state['mg__histogram_box_selection_function'] == 'Positivity identification':
                st.session_state['mg__selected_value_range'] = tuple(x_range)
                st.session_state['mg__min_selection_value'] = x_range[0]
            else:
                st.session_state['mg__histogram_x_range'] = x_range
                st.session_state['mg__random_string'] = generate_random_string()


def plotly_chart_summary_callback():
    if 'mg__plotly_chart_summary__do_not_persist' in st.session_state:
        if st.session_state['mg__plotly_chart_summary__do_not_persist']['selection']['points']:
            selected_z_score = st.session_state['mg__plotly_chart_summary__do_not_persist']['selection']['points'][0]['x']
            df_summary_contents = st.session_state['mg__df_summary_contents'].set_index('Z score')
            df_selected_threshold = df_summary_contents.loc[selected_z_score, 'Threshold']
            selected_value_range = st.session_state['mg__selected_value_range']
            st.session_state['mg__selected_value_range'] = (df_selected_threshold, selected_value_range[1])
            st.session_state['mg__min_selection_value'] = df_selected_threshold


def df_summary_callback():
    if 'mg__df_summary__do_not_persist' in st.session_state:
        if st.session_state['mg__df_summary__do_not_persist']['selection']['rows']:
            df_selected_threshold = st.session_state['mg__df_summary_contents'].iloc[st.session_state['mg__df_summary__do_not_persist']['selection']['rows'][0]]['Threshold']
            selected_value_range = st.session_state['mg__selected_value_range']
            st.session_state['mg__selected_value_range'] = (df_selected_threshold, selected_value_range[1])
            st.session_state['mg__min_selection_value'] = df_selected_threshold


def generate_box_and_whisker(apply_another_filter, df, column_for_filtering, another_filter_column, values_on_which_to_filter, images_in_plotting_group_1, images_in_plotting_group_2, all_cells=True):

    # If we're ready to apply a filter, then create it
    if not apply_another_filter:
        filter_loc = pd.Series(True, index=df.index)
    else:
        filter_loc = df[another_filter_column].isin(values_on_which_to_filter)

    # Get the locations of the images in each group, including the optional filter just created
    image_loc_group_1 = df['Slide ID'].isin(images_in_plotting_group_1) & filter_loc
    image_loc_group_2 = df['Slide ID'].isin(images_in_plotting_group_2) & filter_loc

    # Get only the data for the column of interest for the first group of images with the filter applied
    ser_for_z_score = df.loc[image_loc_group_1, column_for_filtering]

    # From those data, get the values corresponding to -1 to 10 standard deviations above the mean
    z_scores = np.arange(-1, 11)
    thresholds = ser_for_z_score.mean() + z_scores * ser_for_z_score.std()

    # Initialize the positive percentages holders
    group_1_holder = []
    group_2_holder = []
    data_for_box_plot_holder = []

    # For each threshold...
    for threshold, z_score in zip(thresholds, z_scores):

        # Get the locations where the selected filtering column is at least the current threshold value
        positive_loc = df[column_for_filtering] >= threshold

        # If we want the positive percentage of all the cells in each group...
        if all_cells:

            # Get the positive percentage using the current threshold for each of the groups, for the entire dataset (not on a per-image basis)
            group_1_holder.append((image_loc_group_1 & positive_loc).sum() / image_loc_group_1.sum() * 100)
            group_2_holder.append((image_loc_group_2 & positive_loc).sum() / image_loc_group_2.sum() * 100)

        # If we want the positive percentage of the cells in each image, in each group...
        else:

            # Get the positive percentage using the current threshold for each the grooups, for each image
            # Each of these is a series with the image names as the index and the positive percentage as the values
            ser_group_1_pos_perc = df.loc[image_loc_group_1 & positive_loc, 'Slide ID'].value_counts() / df.loc[image_loc_group_1, 'Slide ID'].value_counts() * 100
            ser_group_2_pos_perc = df.loc[image_loc_group_2 & positive_loc, 'Slide ID'].value_counts() / df.loc[image_loc_group_2, 'Slide ID'].value_counts() * 100

            # Create two dataframes holding all the data in columns as will ultimately desired, adding them in turn to the box plot data holder
            # Doing it in this format for simple subsequent box plotting with plotly express
            df_group_1_pos_perc = ser_group_1_pos_perc.to_frame()
            df_group_1_pos_perc.columns = ['Positive %']
            df_group_2_pos_perc = ser_group_2_pos_perc.to_frame()
            df_group_2_pos_perc.columns = ['Positive %']
            df_group_1_pos_perc['Threshold'] = threshold
            df_group_2_pos_perc['Threshold'] = threshold
            df_group_1_pos_perc['Z score'] = z_score
            df_group_2_pos_perc['Z score'] = z_score
            df_group_1_pos_perc['Group'] = 'Baseline'
            df_group_2_pos_perc['Group'] = 'Signal'
            data_for_box_plot_holder.append(df_group_1_pos_perc)
            data_for_box_plot_holder.append(df_group_2_pos_perc)

            # Name each series the value of the current threshold
            ser_group_1_pos_perc.name = threshold
            ser_group_2_pos_perc.name = threshold

            # Append the series to the holders
            group_1_holder.append(ser_group_1_pos_perc)
            group_2_holder.append(ser_group_2_pos_perc)

    # If we want the positive percentage of all the cells in each group...
    if all_cells:
    
        # Create a plotly figure
        fig = go.Figure()

        # Plot positive percentage in whole dataset vs. threshold for group 1
        fig.add_trace(go.Scatter(x=z_scores, y=group_1_holder, mode='lines+markers', name='Baseline'))

        # Plot positive percentage in whole dataset vs. threshold for group 2
        fig.add_trace(go.Scatter(x=z_scores, y=group_2_holder, mode='lines+markers', name='Signal'))

        # Create the summary dataframe
        df_summary = pd.DataFrame({'Z score': z_scores, 'Threshold': thresholds, 'Positive percentage (all cells) for baseline group': group_1_holder, 'Positive percentage (all cells) for signal group': group_2_holder})

    # If we want the positive percentage of the cells in each image, in each group...
    else:

        # Create the dataframe holding all the data for the desired box plot
        df_box_plot = pd.concat(data_for_box_plot_holder, axis='rows')

        # Create a dataframe from the positive percentages for each image in each group
        df_group_1_pos_perc = pd.concat(group_1_holder, axis='columns')
        df_group_2_pos_perc = pd.concat(group_2_holder, axis='columns')

        # Calculate the average of each column for group 1 and group 2
        avg_group_1 = df_group_1_pos_perc.mean()
        avg_group_2 = df_group_2_pos_perc.mean()

        # Plot the positive percentage in each image vs. threshold for both groups
        # fig = go.Figure()
        # fig.add_trace(go.Scatter(x=avg_group_1.index, y=avg_group_1.values, mode='lines+markers', name='Baseline'))
        # fig.add_trace(go.Scatter(x=avg_group_2.index, y=avg_group_2.values, mode='lines+markers', name='Signal'))

        # Create the desired box plot
        fig = px.box(df_box_plot, x='Z score', y='Positive %', color='Group', points='all')

        # Create the summary dataframe
        df_summary = pd.DataFrame({'Z score': z_scores, 'Threshold': thresholds, 'Positive % (avg. over images) for baseline group': avg_group_1.values, 'Positive % (avg. over images) for signal group': avg_group_2.values})

    # Update the layout of the plot
    fig.update_layout(title='Positive percentage vs. baseline Z score',
                    xaxis_title='Baseline Z score',
                    yaxis_title='Positive percentage',
                    legend_title='Group')

    # Return the plot, the Z scores, and the thresholds
    return fig, df_summary


def reset_values_on_which_to_filter_another_column():
    st.session_state['mg__values_on_which_to_filter'] = []


def update_selected_value_range_from_minimum_selection_value():
    min_selection_value = st.session_state['mg__min_selection_value']
    selected_value_range = st.session_state['mg__selected_value_range']
    st.session_state['mg__selected_value_range'] = (min_selection_value, selected_value_range[1])


def update_minimum_selection_value_from_selected_value_range():
    selected_value_range = st.session_state['mg__selected_value_range']
    st.session_state['mg__min_selection_value'] = selected_value_range[0]


def calculate_histogram(ser):
    # calculate_histogram(df.loc[image_loc_group_1, column_for_filtering])
    if len(ser) != 0:
        vc = ser.value_counts(sort=False)
        df = pd.DataFrame({'Values': vc.index.to_list(), 'Counts': vc.to_list()})
    else:
        df = pd.DataFrame({'Values': [], 'Counts': []})
    return df


def calculate_kde(ser, kde_grid_size):
    # calculate_kde(df_batch_normalized.loc[image_loc_group_1, column_for_filtering], st.session_state['mg__kde_grid_size'])
    if len(ser) != 0:
        kde = gaussian_kde(ser)  # create a gaussian_kde object
        x = np.linspace(min(ser), max(ser), kde_grid_size)  # create a range of values over which to evaluate the KDE
        y = kde.evaluate(x)  # evaluate the KDE over the range of values
        df = pd.DataFrame({'Value': x, 'Density': y})
    else:
        df = pd.DataFrame({'Value': [], 'Density': []})
    return df


def calculate_kdes_on_same_grid(ser1, ser2, kde_grid_size):
    # calculate_kdes_on_same_grid(df_batch_normalized.loc[image_loc_group_1, column_for_filtering], df_batch_normalized.loc[image_loc_group_2, column_for_filtering], st.session_state['mg__kde_grid_size'])
    if (len(ser1) != 0) and (len(ser2) != 0):
        kde1 = gaussian_kde(ser1)  # create a gaussian_kde object
        kde2 = gaussian_kde(ser2)  # create a gaussian_kde object
        x = np.linspace(min(min(ser1), min(ser2)), max(max(ser1), max(ser2)), kde_grid_size)  # create a range of values over which to evaluate the KDE
        y1 = kde1.evaluate(x)  # evaluate the KDE over the range of values
        y2 = kde2.evaluate(x)  # evaluate the KDE over the range of values
        df1 = pd.DataFrame({'Value': x, 'Density': y1})
        df2 = pd.DataFrame({'Value': x, 'Density': y2})
    else:
        df1 = calculate_kde(ser1, kde_grid_size)
        df2 = calculate_kde(ser2, kde_grid_size)
    return df1, df2


# Add a function to reset the KDEs (numerical) and histograms (categorical) since we want to update this not just initially but also if the KDE resolution changes or they are calculated based on a different selection of images
def reset_kdes_and_hists(df):
    st.session_state['mg__kdes_or_hists_to_plot'] = pd.DataFrame(columns=st.session_state['mg__all_columns'], index=(['All images'] + df['Slide ID'].unique().tolist()), dtype='object')


# Update the dependencies of the selectbox for the current analysis column
def update_dependencies_of_filtering_widgets():
    df = st.session_state['mg__df_batch_normalized']
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
        st.session_state['mg__min_selection_value'] = st.session_state['mg__curr_column_range'][0]  # initialize the minimum selection value to the minimum of the range
        st.session_state['mg__histogram_x_range'] = list(st.session_state['mg__curr_column_range'])
        st.session_state['mg__random_string'] = generate_random_string()
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
        phenotype_name_changes = dict()
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
                    column_filter_bools = df[column].apply(lambda x: x in values[0]).astype('bool')
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
            old_phenotype_name = phenotype
            new_phenotype_name = phenotype.strip().replace('+', '(plus)').replace('-', '(dash)')
            if old_phenotype_name != new_phenotype_name:
                phenotype_name_changes[old_phenotype_name] = new_phenotype_name
            st.session_state['mg__df'].loc[image_loc, 'Phenotype {}'.format(new_phenotype_name)] = phenotype_bools.apply(lambda x: ('+' if x else '-'))  # since '+' and '-' are forbidden for the time being

        # Debugging output
        print('------------------------')

        # If any phenotype names were changed, output the changes
        if phenotype_name_changes:
            st.info(f'These phenotype transformations were made to avoid "+" and "-" characters: {phenotype_name_changes}')

        # Save the gating table to disk
        if st.session_state['input_metadata']['datafile_path'] is not None:
            datafile_name = os.path.splitext(os.path.basename(st.session_state['input_metadata']['datafile_path']))[0]
        else:
            datafile_name = 'from_memory'
        gating_filename = 'gating_table_for_{}_for_datafile_{}-{}.csv'.format(filtering_section_name, datafile_name, utils.get_timestamp())
        df_phenotype_assignments.to_csv(path_or_buf=os.path.join(os.path.join('.', 'output'), gating_filename), index=True)
        st.info('File {} written to disk'.format(gating_filename))


# Function to clear the session state as would be desired when loading a new dataset
def clear_session_state(keep_keys=[]):
    for key in (set([key for key in st.session_state.keys() if key.startswith('mg__')]) - set(keep_keys)):
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
    st.session_state['mg__df'].drop(columns=new_phenotypes, inplace=True)


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


# Allow a sample gating table to be loaded for quick testing
def load_sample_gating_table(csv_filename):
    st.session_state['mg__de_phenotype_assignments'].update_editor_contents(pd.read_csv(os.path.join('.', 'sample_phenotyping', csv_filename)), reset_key=False)


def main():
    '''
    Main function for running the page
    '''

    # Global variable
    st_key_prefix = 'mg__'

    # If 'input_dataset' isn't in the session state, print an error message and return
    if 'input_dataset' not in st.session_state:
        st.error('An input dataset has not yet been opened. Please do so using the "Open File" page in the sidebar.')
        return

    # Set the default dataframes to be edited
    default_df_current_phenotype     = pd.DataFrame(columns=['Column for filtering', 'Minimum value', 'Maximum value'])
    default_df_phenotype_assignments = pd.DataFrame()
    default_df_field_matching        = pd.DataFrame(columns=['Intensity field', 'Corresponding thresholded marker field'])

    # Constants
    num_categorical_values_cutoff = 10
    tol = 1e-8

    # Initialize some things in the session state
    if 'mg__current_phenotype_name' not in st.session_state:
        st.session_state['mg__current_phenotype_name'] = ''
    if 'mg__de_current_phenotype' not in st.session_state:
        st.session_state['mg__de_current_phenotype'] = sde.DataframeEditor(df_name='mg__df_current_phenotype', default_df_contents=default_df_current_phenotype)
    if 'mg__de_phenotype_assignments' not in st.session_state:
        st.session_state['mg__de_phenotype_assignments'] = sde.DataframeEditor(df_name='mg__df_phenotype_assignments', default_df_contents=default_df_phenotype_assignments)
    if 'mg__de_field_matching' not in st.session_state:
        st.session_state['mg__de_field_matching'] = sde.DataframeEditor(df_name='mg__df_field_matching', default_df_contents=default_df_field_matching)
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
            st.toggle(label='Perform batch normalization', key='mg__do_batch_norm')
            batch_normalization_func = (z_score_normalize if st.session_state['mg__do_batch_norm'] else lambda df, _: df)  # initially, create a very simple batch normalization option using the Z score, which appears to be justified in literature

        # Create columns for the data loading button and spinner
        data_butt_cols = st.columns([3, 1])

        # In the first column, create the data loading button
        with data_butt_cols[0]:
            MaG_load_hit = st.button('Load data', use_container_width=True, on_click=clear_session_state, kwargs={'keep_keys': ['mg__do_batch_norm', 'mg__df']})

        # Initialize the extra settings toggle to off
        if 'mg__extra_settings' not in st.session_state:
            st.session_state['mg__extra_settings'] = False
        extra_settings = st.checkbox('Show advanced settings', key='mg__extra_settings')

        # In the second column, create the data loading spinner
        with data_butt_cols[1]:
            if MaG_load_hit:
                with st.spinner('Loading Data'):
                    st.session_state['mg__df'] = st.session_state['input_dataset'].data  # not making a copy because we deliberately want to modify the in-memory dataframe. Note that all is ever done to the dataframe is renaming the original phenotype columns and then adding new ones.
                    unique_images = st.session_state['mg__df']['Slide ID'].unique()
                    st.session_state['mg__unique_images_short'] = ['-'.join(x.split('-')[1:]) for x in unique_images]
                    st.session_state['mg__unique_image_dict'] = dict(zip(st.session_state['mg__unique_images_short'], unique_images))
                    phenotype_columns = [column for column in st.session_state['mg__df'].columns if column.startswith('Phenotype ')]
                    st.session_state['mg__df'].rename(columns=dict(zip(phenotype_columns, [column.replace('Phenotype ', 'Phenotype_orig ') for column in phenotype_columns])), inplace=True)

                    srs_integer_columns = st.session_state['mg__df'].select_dtypes(include=['integer']).columns.to_series()
                    categorical_integer_columns = srs_integer_columns.loc[pd.Index([len(st.session_state['mg__df'][int_column].unique()) <= num_categorical_values_cutoff for int_column in srs_integer_columns])].to_list()
                    st.session_state['mg__all_numeric_columns'] = st.session_state['mg__df'].select_dtypes(include='number').columns.drop(categorical_integer_columns)
                    st.session_state['mg__all_columns'] = st.session_state['mg__df'].columns
                    st.session_state['mg__column_config'] = {"Corresponding thresholded marker field": st.column_config.SelectboxColumn("Corresponding thresholded marker field", help="Tresholded marker field corresponding to the intensity at left", options=st.session_state['mg__all_columns'], required=True)}
                    st.session_state['mg__df_batch_normalized'] = batch_normalization_func(st.session_state['mg__df'], st.session_state['mg__all_numeric_columns'])
                    
        # Warn the user that they need to load the data at least once
        if 'mg__df' not in st.session_state:
            st.warning('Please load data from the selections above.')
            return

    # Load the data and some resulting processed data
    df = st.session_state['mg__df']
    df_batch_normalized = st.session_state['mg__df_batch_normalized']
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

            # If extra settings are to be displayed, create widgets for selecting particular images to define two groups
            if st.session_state['mg__extra_settings']:

                # Instantiate the object
                image_selector = image_filter.ImageFilter(df, image_colname='Slide ID')

                # If the image filter is not ready (which means the filtering dataframe was not generated), return
                if not image_selector.ready:
                    return

                # Create two image filters
                st.session_state['mg__images_in_plotting_group_1'] = image_selector.select_images(key='baseline', color='blue')
                st.session_state['mg__images_in_plotting_group_2'] = image_selector.select_images(key='signal', color='red')

                # Define other keys which are no longer relevant but we are leaving in so the rest of the code does not break
                if 'mg__filter_on_another_column' not in st.session_state:
                    st.session_state['mg__filter_on_another_column'] = False
                if 'mg__values_on_which_to_filter' not in st.session_state:
                    reset_values_on_which_to_filter_another_column()

            else:
                st.session_state['mg__images_in_plotting_group_1'] = []
                st.session_state['mg__images_in_plotting_group_2'] = []
                st.session_state['mg__filter_on_another_column'] = False

            # If either list is not empty, then we know we want to plot 1 or 2 groups instead of all images or a single image
            if st.session_state['mg__images_in_plotting_group_1'] or st.session_state['mg__images_in_plotting_group_2']:
                use_groups_for_plotting = True
                groups_plotting_data = {'images_in_plotting_group_1': st.session_state['mg__images_in_plotting_group_1'], 'images_in_plotting_group_2': st.session_state['mg__images_in_plotting_group_2']}
            else:
                use_groups_for_plotting = False
                groups_plotting_data = None
            all_groups_plotting_data = {'use_groups_for_plotting': use_groups_for_plotting, 'groups_plotting_data': groups_plotting_data}

            # If we want to filter on another column and there exist values to which to filter, then we know we want to apply another filter
            if st.session_state['mg__filter_on_another_column'] and st.session_state['mg__values_on_which_to_filter']:
                apply_another_filter = True
                another_filter_data = {'another_filter_column': st.session_state['mg__another_filter_column'], 'values_on_which_to_filter': st.session_state['mg__values_on_which_to_filter']}
            else:
                apply_another_filter = False
                another_filter_data = None
            all_another_filter_data = {'apply_another_filter': apply_another_filter, 'another_filter_data': another_filter_data}

            # Reset the holder of the KDEs and histograms if some conditions are set. If the conditions are not met, since generating KDEs/histograms does not take trivial time, it's best to reuse their previously calculated values
            do_reset_kdes_and_hists = (
                ('mg__kdes_or_hists_to_plot' not in st.session_state) or
                ('mg__all_groups_plotting_data_previous' not in st.session_state) or
                (('mg__all_groups_plotting_data_previous' in st.session_state) and (all_groups_plotting_data != st.session_state['mg__all_groups_plotting_data_previous'])) or
                ('mg__all_another_filter_data_previous' not in st.session_state) or
                (('mg__all_another_filter_data_previous' in st.session_state) and (all_another_filter_data != st.session_state['mg__all_another_filter_data_previous']))
            )
            if do_reset_kdes_and_hists:
                reset_kdes_and_hists(df)

            # Store the current data for comparison on the next script top-to-bottom run
            st.session_state['mg__all_groups_plotting_data_previous'] = all_groups_plotting_data
            st.session_state['mg__all_another_filter_data_previous'] = all_another_filter_data

            # If no KDE or histogram has been calculated for the current image(s) and column for filtering, then calculate it. Otherwise, skip for efficiency since at this point there's no need to recalculate them (per the reset_kdes_and_hists(df) line above)
            if st.session_state['mg__kdes_or_hists_to_plot'].loc[image_for_filtering, column_for_filtering] is np.nan:

                # Initialize st.session_state['mg__kde_grid_size'] if the column of interest is numeric, since we're using it below
                if st.session_state['mg__selected_column_type'] == 'numeric':
                    if 'mg__kde_grid_size' not in st.session_state:
                        st.session_state['mg__kde_grid_size'] = 500

                # If we're ready to apply a filter, then create it
                if not apply_another_filter:
                    filter_loc = pd.Series(True, index=df.index)
                else:
                    filter_loc = df[st.session_state['mg__another_filter_column']].isin(st.session_state['mg__values_on_which_to_filter'])

                # If we want to plot all or one image...
                if not use_groups_for_plotting:
                    if image_for_filtering == 'All images':
                        image_loc = pd.Series(True, index=df.index) & filter_loc
                    else:
                        image_loc = (df['Slide ID'] == image_for_filtering) & filter_loc
                    if st.session_state['mg__selected_column_type'] == 'numeric':
                        st.session_state['mg__kdes_or_hists_to_plot'].loc[image_for_filtering, column_for_filtering] = [calculate_kde(df_batch_normalized.loc[image_loc, column_for_filtering], st.session_state['mg__kde_grid_size'])]
                    else:
                        st.session_state['mg__kdes_or_hists_to_plot'].loc[image_for_filtering, column_for_filtering] = [calculate_histogram(df.loc[image_loc, column_for_filtering])]

                # If we want to plot 1 or 2 groups...
                else:
                    image_loc_group_1 = df['Slide ID'].isin(st.session_state['mg__images_in_plotting_group_1']) & filter_loc
                    image_loc_group_2 = df['Slide ID'].isin(st.session_state['mg__images_in_plotting_group_2']) & filter_loc
                    if st.session_state['mg__selected_column_type'] == 'numeric':
                        curr_df_group_1, curr_df_group_2 = calculate_kdes_on_same_grid(df_batch_normalized.loc[image_loc_group_1, column_for_filtering], df_batch_normalized.loc[image_loc_group_2, column_for_filtering], st.session_state['mg__kde_grid_size'])
                    else:
                        curr_df_group_1 = calculate_histogram(df.loc[image_loc_group_1, column_for_filtering])
                        curr_df_group_2 = calculate_histogram(df.loc[image_loc_group_2, column_for_filtering])
                    st.session_state['mg__kdes_or_hists_to_plot'].loc[image_for_filtering, column_for_filtering] = [(curr_df_group_1, curr_df_group_2)]

            # Assign the KDE/histogram to a dataframe that will be saved
            kde_or_hist_to_plot_full = st.session_state['mg__kdes_or_hists_to_plot'].loc[image_for_filtering, column_for_filtering][0]

            # If the selected column is numeric...
            if st.session_state['mg__selected_column_type'] == 'numeric':

                # If extra settings are requested, add an option to update the KDE grid size from its default of 200
                if st.session_state['mg__extra_settings']:
                    st.number_input('KDE grid size:', min_value=1, max_value=1000, key='mg__kde_grid_size', on_change=reset_kdes_and_hists, args=(df,))
                    if 'mg__min_selection_value' not in st.session_state:
                        st.session_state['mg__min_selection_value'] = column_range[0]
                    st.number_input('Minimum selection value:', min_value=float(column_range[0]), max_value=float(column_range[1]), key='mg__min_selection_value', on_change=update_selected_value_range_from_minimum_selection_value)

                # Draw a range slider widget for selecting the range min and max
                st.slider(label='Selected value range:', min_value=float(column_range[0]), max_value=float(column_range[1]), key='mg__selected_value_range', on_change=update_minimum_selection_value_from_selected_value_range)
                selected_min_val, selected_max_val = st.session_state['mg__selected_value_range']

                # Create a view of the full dataframe that is the selected subset
                if not use_groups_for_plotting:
                    df_to_plot_selected = kde_or_hist_to_plot_full[(kde_or_hist_to_plot_full['Value'] >= selected_min_val) & (kde_or_hist_to_plot_full['Value'] <= selected_max_val)]
                else:  # kde_or_hist_to_plot_full is a tuple
                    df_to_plot_selected = (
                        kde_or_hist_to_plot_full[0][(kde_or_hist_to_plot_full[0]['Value'] >= selected_min_val) & (kde_or_hist_to_plot_full[0]['Value'] <= selected_max_val)],
                        kde_or_hist_to_plot_full[1][(kde_or_hist_to_plot_full[1]['Value'] >= selected_min_val) & (kde_or_hist_to_plot_full[1]['Value'] <= selected_max_val)]
                    )

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
                            image_loc = pd.Series(True, index=df.index)
                        else:
                            image_loc = df['Slide ID'] == image_for_filtering

                        # Get the thresholded marker column values
                        srs_marker_column_values = df_batch_normalized.loc[image_loc, marker_column]

                        # Set the indices of that series to the corresponding intensities
                        srs_marker_column_values.index = df_batch_normalized.loc[image_loc, column_for_filtering]

                        # Sort the series by increasing intensity
                        srs_marker_column_values = srs_marker_column_values.sort_index()

                        # Determine whether the intensity values were deemed "positive" presumably by looking at the original image
                        if srs_marker_column_values.dtype == 'object':
                            positive_loc = srs_marker_column_values.apply(lambda x: x[-1] == '+')
                        else:
                            positive_loc = srs_marker_column_values == 1

                        # Get the lowest-intensity "positive" intensity/marker
                        intensity_cutoff = srs_marker_column_values[positive_loc].index[0]

                if 'mg__histogram_box_selection_function' not in st.session_state:
                    st.session_state['mg__histogram_box_selection_function'] = 'Positivity identification'
                st.radio('Histogram box selection function:', ['Positivity identification', 'Zoom'], key='mg__histogram_box_selection_function')

                st.button('Reset x-axis zoom', on_click=reset_x_axis_range, args=(use_groups_for_plotting, kde_or_hist_to_plot_full))

                # Plot the Plotly figure in Streamlit
                fig = go.Figure()

                if not use_groups_for_plotting:
                    fig.add_trace(go.Scatter(x=kde_or_hist_to_plot_full['Value'], y=kde_or_hist_to_plot_full['Density'], fill='tozeroy', mode='markers', marker=dict(color='rgba(255, 0, 0, 0.25)', size=1), fillcolor='rgba(255, 0, 0, 0.25)', name='All selected images', hovertemplate=' '))
                    fig.add_trace(go.Scatter(x=df_to_plot_selected['Value'], y=df_to_plot_selected['Density'], fill='tozeroy', mode='none', fillcolor='rgba(255, 0, 0, 0.5)', name='Selection', hoverinfo='skip'))
                    if intensity_cutoff is not None:
                        fig.add_vline(x=intensity_cutoff, line_color='green', line_width=3, line_dash="dash", annotation_text="Previous threshold: ~{}".format((intensity_cutoff)), annotation_font_size=18, annotation_font_color="green")
                    fig.update_layout(hovermode='x unified', xaxis_title='Column value', yaxis_title='Density')
                    fig.update_layout(legend=dict(yanchor="top", y=1.2, xanchor="left", x=0.01, orientation="h"))
                else:
                    fig.add_trace(go.Scatter(x=kde_or_hist_to_plot_full[0]['Value'], y=kde_or_hist_to_plot_full[0]['Density'], fill='tozeroy', mode='markers', marker=dict(color='rgba(0, 255, 255, 0.25)', size=1), fillcolor='rgba(0, 255, 255, 0.25)', name='Baseline group', hovertemplate=' '))
                    fig.add_trace(go.Scatter(x=df_to_plot_selected[0]['Value'], y=df_to_plot_selected[0]['Density'], fill='tozeroy', mode='none', fillcolor='rgba(0, 255, 255, 0.5)', name='Baseline selection', hoverinfo='skip'))
                    fig.add_trace(go.Scatter(x=kde_or_hist_to_plot_full[1]['Value'], y=kde_or_hist_to_plot_full[1]['Density'], fill='tozeroy', mode='markers', marker=dict(color='rgba(255, 0, 0, 0.25)', size=1), fillcolor='rgba(255, 0, 0, 0.25)', name='Signal group', hovertemplate=' '))
                    fig.add_trace(go.Scatter(x=df_to_plot_selected[1]['Value'], y=df_to_plot_selected[1]['Density'], fill='tozeroy', mode='none', fillcolor='rgba(255, 0, 0, 0.5)', name='Signal selection', hoverinfo='skip'))
                    if intensity_cutoff is not None:
                        fig.add_vline(x=intensity_cutoff, line_color='green', line_width=3, line_dash="dash", annotation_text="Previous threshold: ~{}".format((intensity_cutoff)), annotation_font_size=18, annotation_font_color="green")
                    fig.update_layout(hovermode='x unified', xaxis_title='Column value', yaxis_title='Density')
                    fig.update_layout(legend=dict(yanchor="top", y=1.2, xanchor="left", x=0.01, orientation="h"))

                if 'mg__histogram_x_range' not in st.session_state:
                    reset_x_axis_range(use_groups_for_plotting, kde_or_hist_to_plot_full)

                fig.update_xaxes(range=st.session_state['mg__histogram_x_range'])

                # Set Plotly chart in streamlit
                st.plotly_chart(fig, on_select=plotly_chart_histogram_callback, key=('mg__plotly_chart_histogram_' + st.session_state['mg__random_string'] + '__do_not_persist'), selection_mode='box')

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

                    # Plot the plotly figure in Streamlit
                    fig = go.Figure()
                    if not use_groups_for_plotting:
                        bar_colors = ['rgba(255, 0, 0, 0.5)' if value in selected_items else 'rgba(255, 0, 0, 0.25)' for value in kde_or_hist_to_plot_full['Values']]  # define bar colors
                        fig.add_trace(go.Bar(x=kde_or_hist_to_plot_full['Values'], y=kde_or_hist_to_plot_full['Counts'], marker_color=bar_colors))  # create bar chart
                        fig.update_layout(xaxis_title='Value', yaxis_title='Count', title_text='Counts of values of column {}'.format(column_for_filtering))  # set labels and title
                        st.plotly_chart(fig)  # display plot
                    else:
                        bar_colors_group_1 = ['rgba(0, 255, 255, 0.5)' if value in selected_items else 'rgba(0, 255, 255, 0.25)' for value in kde_or_hist_to_plot_full[0]['Values']]  # define bar colors
                        bar_colors_group_2 = ['rgba(255, 0, 0, 0.5)' if value in selected_items else 'rgba(255, 0, 0, 0.25)' for value in kde_or_hist_to_plot_full[1]['Values']]
                        fig.add_trace(go.Bar(x=kde_or_hist_to_plot_full[0]['Values'], y=kde_or_hist_to_plot_full[0]['Counts'], marker_color=bar_colors_group_1, name='Baseline group'))  # create bar chart
                        fig.add_trace(go.Bar(x=kde_or_hist_to_plot_full[1]['Values'], y=kde_or_hist_to_plot_full[1]['Counts'], marker_color=bar_colors_group_2, name='Signal group'))
                        fig.update_layout(xaxis_title='Value', yaxis_title='Count', title_text='Counts of values of column {}'.format(column_for_filtering))  # set labels and title
                        st.plotly_chart(fig)  # display plot

                    # Set the selection dictionary for the current filter to pass on to the current phenotype definition
                    selection_dict = {'column_for_filtering': column_for_filtering, 'selected_min_val': None, 'selected_max_val': None, 'selected_column_values': selected_items}

                    # Whether to disable the add column button
                    add_column_button_disabled = False

                # If there are too many unique values in the column, do nothing
                else:
                    st.info('There are too many unique values (more than {}) in the selected categorical column. Filtering not currently allowed.'.format(num_categorical_values_cutoff), icon="ℹ️")
                    selection_dict = dict()
                    add_column_button_disabled = True

            # If applicable, allow the user to plot the desired box and whisker plot using the different thresholds
            if extra_settings and use_groups_for_plotting and (st.session_state['mg__selected_column_type'] == 'numeric'):
                if 'mg__plot_box_and_whisker' not in st.session_state:
                    st.session_state['mg__plot_box_and_whisker'] = False
                if st.toggle('Plot box and whisker', key='mg__plot_box_and_whisker'):
                    if 'mg__positive_percentage_per_image' not in st.session_state:
                        st.session_state['mg__positive_percentage_per_image'] = True
                    st.checkbox('Calculate positive percentages separately for each image', key='mg__positive_percentage_per_image')
                    fig, df_summary = generate_box_and_whisker(apply_another_filter, df_batch_normalized, column_for_filtering, st.session_state['mg__another_filter_column'], st.session_state['mg__values_on_which_to_filter'], st.session_state['mg__images_in_plotting_group_1'], st.session_state['mg__images_in_plotting_group_2'], all_cells=(not st.session_state['mg__positive_percentage_per_image']))
                    st.session_state['mg__df_summary_contents'] = df_summary
                    st.plotly_chart(fig, on_select=plotly_chart_summary_callback, key='mg__plotly_chart_summary__do_not_persist')
                    st.dataframe(df_summary, hide_index=True, key="mg__df_summary__do_not_persist", on_select=df_summary_callback, selection_mode=["single-row"])

            # Add the current column filter to the current phenotype assignment
            st.button(':star2: Add column filter to current phenotype :star2:', use_container_width=True, on_click=update_dependencies_of_button_for_adding_column_filter_to_current_phenotype, kwargs=selection_dict, disabled=add_column_button_disabled)

    # Current phenotype and phenotype assignments
    with main_columns[1]:

        # Column header
        st.header(':two: Current phenotype', help='Note you can refine non-list values in the following table by editing them directly or even deleting (or adding) whole rows.')

        # Output the dataframe holding the phenotype that's currently being built
        st.session_state['mg__de_current_phenotype'].dataframe_editor(reset_data_editor_button_text='Reset current phenotype definition', on_change=basic_filter_column_updates)

        # Choose a phenotype name
        st.write('It is fine to use "+" and "-" in the phenotype name below, but keep in mind they will be replaced with "(plus)" and "(dash)" in the downstream dataset.')
        st.text_input(label='Phenotype name:', key='mg__current_phenotype_name')

        # Add the current phenotype to the phenotype assignments table
        st.button(label=':star2: Add phenotype to assignments table :star2:',
                    use_container_width=True, 
                    on_click=update_dependencies_of_button_for_adding_phenotype_to_new_dataset)

        # Column header
        st.header(':three: Phenotype assignments', help='Note you can refine non-list values in the following table by editing them directly or even deleting whole rows.')

        # Output the dataframe holding the specifications for all phenotypes
        st.session_state['mg__de_phenotype_assignments'].dataframe_editor(reset_data_editor_button_text='Reset all phenotype definitions')

        # # Allow a sample gating table to be loaded
        # st.button('Load sample gating table', on_click=load_sample_gating_table, kwargs={'csv_filename': 'sample_gating_table.csv'})

        # Generate the new dataset
        st.button(label=':star2: Append phenotype assignments to the dataset :star2:',
                    use_container_width=True,
                    on_click=add_new_phenotypes_to_main_df, args=(df_batch_normalized, image_for_filtering))
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

        # if st.button('Save dataset to `output` folder'):
        #     st.session_state['mg__df'].to_csv('./output/saved_dataset.csv')

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
            if ('mg__phenotype_to_plot' not in st.session_state) or (st.session_state['mg__phenotype_to_plot'] not in new_phenotypes):
                st.session_state['mg__phenotype_to_plot'] = new_phenotypes[0] if new_phenotypes else None

            # Generate widgets for the plotting parameters
            st.selectbox(label='Image to plot:', options=unique_images_short, key='mg__image_to_plot')
            st.selectbox(label='Phenotype to plot:', options=new_phenotypes, key='mg__phenotype_to_plot')

            # If the button is pressed
            if st.toggle('Plot the selected phenotype in the selected image', help='This is relatively slow to keep updating with every page interaction, so feel free to toggle off if you\'d like faster page loading.'):
                df_for_scatterplot = df.loc[df['Slide ID'] == unique_image_dict[st.session_state['mg__image_to_plot']], ['Cell X Position', 'Cell Y Position', st.session_state['mg__phenotype_to_plot']]]
                fig = px.scatter(data_frame=df_for_scatterplot, x='Cell X Position', y='Cell Y Position', color=st.session_state['mg__phenotype_to_plot'], category_orders={st.session_state['mg__phenotype_to_plot']: ['-', '+']})
                fig.update_xaxes(scaleanchor='y')
                st.plotly_chart(fig)

# Call the main function
if __name__ == '__main__':

    # Set page settings
    st.set_page_config(layout='wide', page_title='Manual Phenotyping on Raw Intensities')
    st.title('Manual Phenotyping on Raw Intensities')

    # Run streamlit-dataframe-editor library initialization tasks at the top of the page
    st.session_state = sde.initialize_session_state(st.session_state)

    # Run Top of Page (TOP) functions
    st.session_state = top.top_of_page_reqs(st.session_state)

    main()

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)
