# Import relevant libraries
import streamlit as st
import app_top_of_page as top
import streamlit_dataframe_editor as sde
import radial_profiles
import time
import numpy as np
import plotly.express as px
from itertools import cycle, islice
import plotly.graph_objects as go
import pandas as pd
import scipy.spatial
from pages import multiaxial_gating
import utils

# Global variable
st_key_prefix = 'radial_profiles__'


def update_phenotyping_sub_options(phenotyping_method_options):
    if st.session_state[st_key_prefix + 'phenotyping_method'] == phenotyping_method_options[0]:
        st.session_state[st_key_prefix + 'apply_thresh_to_selected_group'] = False
        st.session_state[st_key_prefix + 'average_over_all_groups'] = False
    else:
        st.session_state[st_key_prefix + 'apply_thresh_to_selected_group'] = True
        st.session_state[st_key_prefix + 'average_over_all_groups'] = True


def plotly_mean_and_sem(dfs, df_names):

    # Create a Plotly figure
    fig = go.Figure()

    # For each dataframe and its name...
    for df, df_name in zip(dfs, df_names):

        # Calculate mean and SEM over the groups
        mean_values = df.mean(axis='columns')
        sem_values = df.sem(axis='columns')

        # Add a scatter plot
        fig.add_trace(go.Scatter(
            x=mean_values.index,
            y=mean_values,
            error_y=dict(
                type='data', # or 'percent' for percentage-based error bars
                array=sem_values,
                visible=True
            ),
            mode='markers+lines', # Use 'markers' or 'lines' or 'markers+lines'
            name=df_name
        ))

    # Customize layout
    fig.update_layout(
        title=f'Average Positive Percentage Mean with SEM as Error Bars',
        xaxis_title='Z score',
        yaxis_title='Average positive percentage mean',
        showlegend=True
    )

    # Return the plotly figure
    return fig


def get_box_and_whisker_data(df_grouped, df_thresholds, apply_thresh_to_selected_group, df, channel_for_phenotyping, column_identifying_baseline_signal, value_identifying_baseline, value_identifying_signal, row_selection, return_figure_and_summary=True):

    # If apply_thresh_to_selected_group (and not average_over_all_groups), the input into the function generate_box_and_whisker() corresponds to the selected group. The mean and std used are the defaults (None) for the function, i.e., those corresponding to the baseline group of the selected group.
    # If not apply_thresh_to_selected_group (and not average_over_all_groups), the input into the function generate_box_and_whisker() corresponds to the entire dataset. The mean and std (which, when corresponding to a selection from df_thresholds, correspond to the baseline group) used are those from the selected group. --> this is the original phenotyping method (like a single T=0 threshold)!
    # The thresholds (as in df_thresholds) are calculated from just baseline groups, whereas the groups of dataframes correspond to both the baseline and signal groups.
    # If average_over_all_groups, then it's as if we're selecting the groups in turn via a for loop. We don't want to transform the entire dataset every time; we only want to transform each part of the full dataset once. Thus, it generally doesn't make sense to average_over_all_groups and (not apply_thresh_to_selected_group), but it does make sense to average_over_all_groups and apply_thresh_to_selected_group, which essentially thresholds each group using its own baseline group, i.e., corresponds to the new phenotyping method.
    # Thus, here are the settings for the old and new phenotyping methods:
    #   Old: apply_thresh_to_selected_group=False, average_over_all_groups=False, DO have a particular group selected --> "Selected threshold applied to entire dataset"
    #   New: apply_thresh_to_selected_group=True, average_over_all_groups=True, DO NOT have a particular group selected --> "Group-specific threshold applied to each group"
    # Note that there are no other combinations of these settings that will transform each part of the full dataset once, so the above is complete!

    # Obtain the index and dataframe of the group identified by row_selection
    if isinstance(df_grouped, list):
        current_index = None
        df_selected = df_grouped[0][1]
    else:
        current_index = df_thresholds.iloc[row_selection].name
        df_selected = df_grouped.get_group(current_index)

    # Obtain the dataframe to actually transform as well as the mean and std to use for the transformation
    if apply_thresh_to_selected_group:
        df_transform = df_selected
        mean_for_zscore_calc = None
        std_for_zscore_calc = None
    else:
        df_transform = df
        ser_selected = df_thresholds.iloc[row_selection]
        mean_for_zscore_calc = ser_selected.loc['z score = 0']
        std_for_zscore_calc = ser_selected.loc['z score = 1'] - mean_for_zscore_calc

    # Obtain the data that one would get by generating a box and whisker plot
    return_values = multiaxial_gating.generate_box_and_whisker(
        df=df_transform,
        column_for_filtering=channel_for_phenotyping,
        apply_another_filter=False,
        another_filter_column=None,
        values_on_which_to_filter=None,
        images_in_plotting_group_1=df_transform.loc[df_transform[column_identifying_baseline_signal] == value_identifying_baseline, 'Slide ID'].unique(),
        images_in_plotting_group_2=df_transform.loc[df_transform[column_identifying_baseline_signal] == value_identifying_signal, 'Slide ID'].unique(),
        all_cells=False,
        mean_for_zscore_calc=mean_for_zscore_calc,
        std_for_zscore_calc=std_for_zscore_calc,
        return_figure_and_summary=return_figure_and_summary
        )

    # Return the desired values
    if return_figure_and_summary:
        return return_values
    else:
        return return_values, current_index


def delete_session_state_variable(key_suffix):
    del st.session_state[st_key_prefix + key_suffix]


def draw_single_image_scatter_plot(df, image_to_view, column_to_plot, values_to_plot, color_dict, xy_position_columns=['Cell X Position', 'Cell Y Position'], coordinate_scale_factor=1, annulus_spacing_um=250, use_coordinate_mins_and_maxs=False, xmin_col='Cell X Position', xmax_col='Cell X Position', ymin_col='Cell Y Position', ymax_col='Cell Y Position', units='microns', invert_y_axis=False, opacity=0.7):

    # Draw a header
    st.header('Single image scatter plot')

    # If the user wants to display the scatter plot, indicated by a toggle...
    if st_key_prefix + 'show_scatter_plot' not in st.session_state:
        st.session_state[st_key_prefix + 'show_scatter_plot'] = False
    if st.toggle('Show scatter plot', key=st_key_prefix + 'show_scatter_plot'):

        # Filter the DataFrame to include only the selected image
        df_selected_image_and_filter = df.loc[(df['Slide ID'] == image_to_view), xy_position_columns + [column_to_plot]]

        # Optionally scale the coordinates (probably not; should have been done in Datafile Unifier
        if coordinate_scale_factor != 1:
            df_selected_image_and_filter[xy_position_columns] = df_selected_image_and_filter[xy_position_columns] * coordinate_scale_factor

        # Calculate the radius edges of the annuli
        radius_edges, xy_mid = calculate_annuli_radius_edges(df_selected_image_and_filter, annulus_spacing_um=annulus_spacing_um, xy_position_columns=xy_position_columns)

        # Group the DataFrame for the selected image by unique value of the column to plot
        selected_image_grouped_by_value = df_selected_image_and_filter.groupby(column_to_plot)

        # Create the scatter plot
        fig = go.Figure()

        # Loop over the unique values in the column whose values to plot, in order of their frequency
        for value_to_plot in values_to_plot:

            # If the value exists in the selected image...
            if (value_to_plot in selected_image_grouped_by_value.groups) and (len(selected_image_grouped_by_value.groups[value_to_plot]) > 0):

                # Store the dataframe for the current value for the selected image
                df_group = selected_image_grouped_by_value.get_group(value_to_plot)

                # If value is a string, replace '(plus)' with '+' and '(dash)' with '-', since it could likely be a phenotype with those substitutions
                if isinstance(value_to_plot, str):
                    value_str_cleaned = value_to_plot.replace('(plus)', '+').replace('(dash)', '-')
                else:
                    value_str_cleaned = value_to_plot

                # Add the object index to the label
                df_group['hover_label'] = 'Index: ' + df_group.index.astype(str)

                # Works but doesn't scale the shapes
                if not use_coordinate_mins_and_maxs:
                    fig.add_trace(go.Scatter(x=df_group[xy_position_columns[0]], y=df_group[xy_position_columns[1]], mode='markers', name=value_str_cleaned, marker_color=color_dict[value_to_plot], hovertemplate=df_group['hover_label']))

                # Works really well
                else:
                    fig.add_trace(go.Bar(
                        x=((df_group[xmin_col] + df_group[xmax_col]) / 2),
                        y=df_group[ymax_col] - df_group[ymin_col],
                        width=df_group[xmax_col] - df_group[xmin_col],
                        base=df_group[ymin_col],
                        name=value_str_cleaned,
                        marker=dict(
                            color=color_dict[value_to_plot],
                            opacity=opacity,
                        ),
                        hovertemplate=df_group['hover_label']
                    ))

        # Plot circles of radii radius_edges[1:] centered at the midpoint of the coordinates
        for radius in radius_edges[1:]:
            fig.add_shape(
                type='circle',
                xref='x',
                yref='y',
                x0=xy_mid[xy_position_columns[0]] - radius,
                y0=xy_mid[xy_position_columns[1]] - radius,
                x1=xy_mid[xy_position_columns[0]] + radius,
                y1=xy_mid[xy_position_columns[1]] + radius,
                line=dict(
                    color='lime',
                    width=4,
                ),
                opacity=0.75,
            )

        # Update the layout
        fig.update_layout(
            xaxis=dict(
                scaleanchor="y",
                scaleratio=1,
            ),
            yaxis=dict(
                autorange=('reversed' if invert_y_axis else True),
            ),
            title=f'Scatter plot for {image_to_view}',
            xaxis_title=f'Cell X Position ({units})',
            yaxis_title=f'Cell Y Position ({units})',
            legend_title=column_to_plot,
            height=800,  # Set the height of the figure
            width=800,  # Set the width of the figure
        )

        # Plot the plotly chart in Streamlit
        st.plotly_chart(fig, use_container_width=True)

        # Write the analysis results for the selected image
        if st_key_prefix + 'df_analysis_results' in st.session_state:
            st.write(st.session_state[st_key_prefix + 'df_analysis_results'].loc[image_to_view])

    # We seem to need to render something on the page after rendering a plotly figure in order for the page to not automatically scroll back to the top when you go to the Previous or Next image... doesn't seem to always work at least in the scatter plotter
    st.write(' ')


def initialize_main_settings(df, unique_images):

    # Main settings section
    st.header('Main settings')

    # Define the main settings columns
    settings_columns_main = st.columns(3)

    # In the first column...
    with settings_columns_main[0]:

        # Store columns of certain types
        if st_key_prefix + 'categorical_columns' not in st.session_state:
            max_num_unique_values = 1000
            categorical_columns = []
            for col in df.columns:
                if df[col].nunique() <= max_num_unique_values:
                    categorical_columns.append(col)
            st.session_state[st_key_prefix + 'categorical_columns'] = categorical_columns
        if st_key_prefix + 'numeric_columns' not in st.session_state:
            st.session_state[st_key_prefix + 'numeric_columns'] = df.select_dtypes(include='number').columns
        categorical_columns = st.session_state[st_key_prefix + 'categorical_columns']
        numeric_columns = st.session_state[st_key_prefix + 'numeric_columns']

        # Choose a column to plot
        if st_key_prefix + 'column_to_plot' not in st.session_state:
            st.session_state[st_key_prefix + 'column_to_plot'] = categorical_columns[0]
        column_to_plot = st.selectbox('Select a column by which to color the points:', categorical_columns, key=st_key_prefix + 'column_to_plot')
        column_to_plot_has_changed = (st_key_prefix + 'column_to_plot_prev' not in st.session_state) or (st.session_state[st_key_prefix + 'column_to_plot_prev'] != column_to_plot)
        st.session_state[st_key_prefix + 'column_to_plot_prev'] = column_to_plot

        # Get some information about the images in the input dataset
        if st_key_prefix + 'ser_size_of_each_image' not in st.session_state:
            st.session_state[st_key_prefix + 'ser_size_of_each_image'] = df['Slide ID'].value_counts()  # calculate the number of objects in each image
        ser_size_of_each_image = st.session_state[st_key_prefix + 'ser_size_of_each_image']

        # Create an image selection selectbox
        if st_key_prefix + 'image_to_view' not in st.session_state:
            st.session_state[st_key_prefix + 'image_to_view'] = unique_images[0]
        image_to_view = st.selectbox('Select image to view:', unique_images, key=st_key_prefix + 'image_to_view')

        # Display the number of cells in the selected image
        st.write(f'Number of cells in image: {ser_size_of_each_image.loc[image_to_view]}')

        # Optionally navigate through the images using Previous and Next buttons
        cols = st.columns(2)
        with cols[0]:
            st.button('Previous image', on_click=go_to_previous_image, args=(unique_images,), disabled=(image_to_view == unique_images[0]), use_container_width=True)
        with cols[1]:
            st.button('Next image', on_click=go_to_next_image, args=(unique_images, ), disabled=(image_to_view == unique_images[-1]), use_container_width=True)

    # In the second column...
    with settings_columns_main[1]:

        # Optionally plot minimum and maximum coordinate fields
        if st_key_prefix + 'use_coordinate_mins_and_maxs' not in st.session_state:
            st.session_state[st_key_prefix + 'use_coordinate_mins_and_maxs'] = False
        use_coordinate_mins_and_maxs = st.checkbox('Use coordinate mins and maxs', key=st_key_prefix + 'use_coordinate_mins_and_maxs')
        settings_columns_refined = st.columns(2)
        if st_key_prefix + 'x_min_coordinate_column' not in st.session_state:
            st.session_state[st_key_prefix + 'x_min_coordinate_column'] = numeric_columns[0]
        if st_key_prefix + 'y_min_coordinate_column' not in st.session_state:
            st.session_state[st_key_prefix + 'y_min_coordinate_column'] = numeric_columns[0]
        if st_key_prefix + 'x_max_coordinate_column' not in st.session_state:
            st.session_state[st_key_prefix + 'x_max_coordinate_column'] = numeric_columns[0]
        if st_key_prefix + 'y_max_coordinate_column' not in st.session_state:
            st.session_state[st_key_prefix + 'y_max_coordinate_column'] = numeric_columns[0]
        with settings_columns_refined[0]:
            xmin_col = st.selectbox('Select a column for the minimum x-coordinate:', numeric_columns, key=st_key_prefix + 'x_min_coordinate_column', disabled=(not use_coordinate_mins_and_maxs))
        with settings_columns_refined[1]:
            xmax_col = st.selectbox('Select a column for the maximum x-coordinate:', numeric_columns, key=st_key_prefix + 'x_max_coordinate_column', disabled=(not use_coordinate_mins_and_maxs))
        with settings_columns_refined[0]:
            ymin_col = st.selectbox('Select a column for the minimum y-coordinate:', numeric_columns, key=st_key_prefix + 'y_min_coordinate_column', disabled=(not use_coordinate_mins_and_maxs))
        with settings_columns_refined[1]:
            ymax_col = st.selectbox('Select a column for the maximum y-coordinate:', numeric_columns, key=st_key_prefix + 'y_max_coordinate_column', disabled=(not use_coordinate_mins_and_maxs))
        units = ('coordinate units' if use_coordinate_mins_and_maxs else 'microns')

    # In the third column...
    with settings_columns_main[2]:

        # Add an option to invert the y-axis
        if st_key_prefix + 'invert_y_axis' not in st.session_state:
            st.session_state[st_key_prefix + 'invert_y_axis'] = False
        invert_y_axis = st.checkbox('Invert y-axis', key=st_key_prefix + 'invert_y_axis')

        # Choose the opacity of objects
        if st_key_prefix + 'opacity' not in st.session_state:
            st.session_state[st_key_prefix + 'opacity'] = 0.7
        opacity = st.number_input('Opacity:', min_value=0.0, max_value=1.0, step=0.1, key=st_key_prefix + 'opacity')

        # Define the colors for the values to plot
        if (st_key_prefix + 'color_dict' not in st.session_state) or column_to_plot_has_changed:
            reset_color_dict(df[column_to_plot])
        values_to_plot = st.session_state[st_key_prefix + 'values_to_plot']
        color_dict = st.session_state[st_key_prefix + 'color_dict']

        # Select a value whose color we want to modify
        if (st_key_prefix + 'value_to_change_color' not in st.session_state) or column_to_plot_has_changed:
            st.session_state[st_key_prefix + 'value_to_change_color'] = values_to_plot[0]
        value_to_change_color = st.selectbox('Value whose color to change:', values_to_plot, key=st_key_prefix + 'value_to_change_color')

        # Create a color picker widget for the selected value
        st.session_state[st_key_prefix + 'new_picked_color'] = color_dict[value_to_change_color]
        st.color_picker('Pick a new color:', key=st_key_prefix + 'new_picked_color', on_change=update_color_for_value, args=(value_to_change_color,))

        # Add a button to reset the colors to their default values
        st.button('Reset plotting colors to defaults', on_click=reset_color_dict, args=(df[column_to_plot],))
        color_dict = st.session_state[st_key_prefix + 'color_dict']

    # Return assigned variables
    return column_to_plot, image_to_view, use_coordinate_mins_and_maxs, xmin_col, xmax_col, ymin_col, ymax_col, units, invert_y_axis, opacity, color_dict, values_to_plot, categorical_columns


def initialize_radial_bin_calculation(df):

    # Radial bins calculation section
    st.header('Radial bins')

    # st.warning('Note, currently the image centroid is hardcoded!')

    # Get some information about the images in the input dataset
    if st_key_prefix + 'unique_images' not in st.session_state:
        st.session_state[st_key_prefix + 'unique_images'] = df['Slide ID'].unique()  # get the unique images in the dataset
    unique_images = st.session_state[st_key_prefix + 'unique_images']

    # Number input for the coordinate scale factor
    key = st_key_prefix + 'coordinate_scale_factor'
    if key not in st.session_state:
        st.session_state[key] = 1
    # coordinate_scale_factor = st.number_input('Coordinate scale factor:', min_value=0.0, key=key, help='This is only necessary if you have forgotten to scale the coordinates in the Datafile Unifier.')
    coordinate_scale_factor = st.session_state[key]

    # Number input for the annulus spacing
    key = st_key_prefix + 'annulus_spacing_um'
    if key not in st.session_state:
        st.session_state[key] = 250
    annulus_spacing_um = st.number_input('Annulus spacing (um):', min_value=0.0, key=key)

    # Multiselect for selection of coordinate columns
    xy_position_columns = ['Cell X Position', 'Cell Y Position']  # this is set in Open File so should always be the same, no need for a widget

    # Calculate radial bins
    if st.button('Calculate radial bins'):
        start_time = time.time()
        df = add_radial_bin_to_dataset(df, unique_images, coordinate_scale_factor=coordinate_scale_factor, annulus_spacing_um=annulus_spacing_um, xy_position_columns=xy_position_columns)
        st.session_state['input_dataset'].data = df
        st.write(f'Calculation of radial bins took {int(np.round(time.time() - start_time))} seconds')
        del st.session_state[st_key_prefix + 'categorical_columns']  # force the categorical columns to be recalculated since we just added one to the dataset

    # Return the necessary variables
    return df, unique_images, coordinate_scale_factor, annulus_spacing_um, xy_position_columns


def calculate_annuli_radius_edges(df, annulus_spacing_um=250, xy_position_columns=['Cell X Position', 'Cell Y Position']):

    # Get the x-y midpoint of the coordinates in df[xy_position_columns]
    xy_min = df[xy_position_columns].min()
    xy_max = df[xy_position_columns].max()
    xy_mid = (xy_min + xy_max) / 2


    # xy_mid.iloc[0] = 10130 / 2 * 0.32
    # xy_mid.iloc[1] = 10130 / 2 * 0.32


    # Get the radius edges that fit within the largest possible radius
    largest_possible_radius = np.linalg.norm(xy_max - xy_mid) + annulus_spacing_um
    num_intervals = largest_possible_radius // annulus_spacing_um  # calculate the number of intervals that fit within the largest_possible_radius
    end_value = (num_intervals + 1) * annulus_spacing_um  # calculate the end value for np.arange to ensure it does not exceed largest_possible_radius
    radius_edges = np.arange(0, end_value, annulus_spacing_um)  # generate steps

    # Return the necessary variables
    return radius_edges, xy_mid


# Function to determine the radial bin for every cell in the dataset
def add_radial_bin_to_dataset(df, unique_images, coordinate_scale_factor=1, annulus_spacing_um=250, xy_position_columns=['Cell X Position', 'Cell Y Position']):

    # For every image in the dataset...
    for current_image in unique_images:

        # Filter the DataFrame to the current image
        df_selected_image_and_filter = df.loc[(df['Slide ID'] == current_image), xy_position_columns]

        # Scale the coordinates. Generally unnecessary as this should have been done (converted to microns) in the unifier
        if coordinate_scale_factor != 1:
            df_selected_image_and_filter[xy_position_columns] = df_selected_image_and_filter[xy_position_columns] * coordinate_scale_factor

        # Calculate the radius edges of the annuli
        radius_edges, xy_mid = calculate_annuli_radius_edges(df_selected_image_and_filter, annulus_spacing_um=annulus_spacing_um, xy_position_columns=xy_position_columns)

        # Construct a KDTree for the current image
        kdtree = scipy.spatial.KDTree(df_selected_image_and_filter[xy_position_columns])

        # For every outer radius...
        prev_indices = None
        for radius in radius_edges[1:]:

            # Get the indices of the points in the image within the current radius
            curr_indices = kdtree.query_ball_point(xy_mid, radius)

            # Get the indices of the points in the current annulus (defined by the outer radius)
            if prev_indices is not None:
                annulus_indices = np.setdiff1d(curr_indices, prev_indices)
                # annulus_indices = curr_indices[~np.isin(curr_indices, prev_indices)]  # note copilot said this would be faster though may need to ensure curr_indices is a numpy array or else will get "TypeError: only integer scalar arrays can be converted to a scalar index"
            else:
                annulus_indices = curr_indices

            # Store the outer radius for all the cells in the current annulus in the current image
            if len(annulus_indices) > 0:
                df.loc[df_selected_image_and_filter.iloc[annulus_indices].index, 'Outer radius'] = radius

            # Store the current indices for the next iteration
            prev_indices = curr_indices

    # Make sure there are no NaNs in the "Outer radius" column using an assertion
    assert df['Outer radius'].isna().sum() == 0, 'There are NaNs in the "Outer radius" column but every cell should have been assigned a radial bin'

    # Return the dataframe with the new "Outer radius" column
    return df


# Function to calculate the percent positives in each annulus in each image
def calculate_percent_positives_for_entire_dataset(df, column_to_plot, unique_images, coordinate_scale_factor=1, annulus_spacing_um=250, xy_position_columns=['Cell X Position', 'Cell Y Position']):

    # For every image in the dataset...
    analysis_results_holder = []
    for current_image in unique_images:

        # Filter the DataFrame to the current image
        df_selected_image_and_filter = df.loc[(df['Slide ID'] == current_image), xy_position_columns + [column_to_plot]]

        # Scale the coordinates. Generally unnecessary as this should have been done (converted to microns) in the unifier
        if coordinate_scale_factor != 1:
            df_selected_image_and_filter[xy_position_columns] = df_selected_image_and_filter[xy_position_columns] * coordinate_scale_factor

        # Get the x-y midpoint of the coordinates in df_selected_image_and_filter[xy_position_columns]
        xy_min = df_selected_image_and_filter[xy_position_columns].min()
        xy_max = df_selected_image_and_filter[xy_position_columns].max()
        xy_mid = (xy_min + xy_max) / 2

        # Get the radius edges that fit within the largest possible radius
        largest_possible_radius = (xy_max - xy_mid).min()
        spacing_um = annulus_spacing_um
        num_intervals = largest_possible_radius // spacing_um  # calculate the number of intervals that fit within the largest_possible_radius
        end_value = (num_intervals + 1) * spacing_um  # calculate the end value for np.arange to ensure it does not exceed largest_possible_radius
        radius_edges = np.arange(0, end_value, spacing_um)  # generate steps

        # Construct a KDTree for the current image
        kdtree = scipy.spatial.KDTree(df_selected_image_and_filter[xy_position_columns])

        # For every outer radius...
        prev_indices = None
        percent_positives = []
        annulus_radius_strings = []
        for radius in radius_edges[1:]:

            # Get the indices of the points in the image within the current radius
            curr_indices = kdtree.query_ball_point(xy_mid, radius)

            # Get the indices of the points in the current annulus (defined by the outer radius)
            if prev_indices is not None:
                annulus_indices = np.setdiff1d(curr_indices, prev_indices)
                # annulus_indices = curr_indices[~np.isin(curr_indices, prev_indices)]  # note copilot said this would be faster though may need to ensure curr_indices is a numpy array or else will get "TypeError: only integer scalar arrays can be converted to a scalar index"
            else:
                annulus_indices = curr_indices

            # Get the series in the current image corresponding to the current annulus and positivity column
            ser_positivity_annulus = df_selected_image_and_filter.iloc[annulus_indices][column_to_plot]

            # Calculate the percent positive, denoted by "+"
            full_size = len(ser_positivity_annulus)
            if full_size != 0:
                percent_positives.append((ser_positivity_annulus == '+').sum() / full_size * 100)
            else:
                percent_positives.append(None)

            # Store the annulus radii as a string for the current annulus
            annulus_radius_strings.append(f'Annulus from {radius - annulus_spacing_um} to {radius} um')

            # Store the current indices for the next iteration
            prev_indices = curr_indices

        # Get a dictionary containing the calculation results for the current image
        percent_positives_dict = dict(zip(annulus_radius_strings, percent_positives))
        percent_positives_dict[f'Number of annuli of width {annulus_spacing_um} um'] = len(percent_positives)

        # Store the dictionary in the analysis results holder
        analysis_results_holder.append(percent_positives_dict)

    # Return a dataframe of the results
    return pd.DataFrame(analysis_results_holder, index=unique_images)


# Function to initialize the preprocessing section
def initialize_preprocessing(df):

    # Preprocessing section
    st.header('Preprocessing')

    # Checkbox for whether to run checks
    key = st_key_prefix + 'run_checks'
    if key not in st.session_state:
        st.session_state[key] = False
    run_checks = st.checkbox('Run checks', key=key)

    # Number input for the threshold for the RawIntNorm check
    key = st_key_prefix + 'perc_thresh_rawintnorm_column_check'
    if key not in st.session_state:
        st.session_state[key] = 0.01
    if run_checks:
        st.number_input('Threshold for the RawIntNorm column check (%):', min_value=0.0, max_value=100.0, key=key)
    perc_thresh_rawintnorm_column_check = st.session_state[key]

    # Number input to select the nuclear intensity channel
    key = st_key_prefix + 'nuclear_channel'
    if key not in st.session_state:
        st.session_state[key] = 1
    nuclear_channel = st.number_input('Nuclear channel:', min_value=1, key=key)

    # Checkbox for whether to apply the z-score filter
    key = st_key_prefix + 'do_z_score_filter'
    if key not in st.session_state:
        st.session_state[key] = True
    do_z_score_filter = st.checkbox('Do z-score filter', key=key)

    # Number input for the z-score filter threshold
    key = st_key_prefix + 'z_score_filter_threshold'
    if key not in st.session_state:
        st.session_state[key] = 3
    if do_z_score_filter:
        st.number_input('z-score filter threshold:', min_value=0.0, key=key)
    z_score_filter_threshold = st.session_state[key]

    # If dataset preprocessing is desired...
    if st.button('Preprocess dataset'):

        # Record the start time
        start_time = time.time()

        # Preprocess the dataset
        df = radial_profiles.preprocess_dataset(
            df,
            perc_thresh_rawintnorm_column_check=perc_thresh_rawintnorm_column_check,
            image_col='Slide ID',
            nuclear_channel=nuclear_channel,
            do_z_score_filter=do_z_score_filter,
            z_score_filter_threshold=z_score_filter_threshold,
            run_checks=run_checks
        )

        # Output the time taken
        st.write(f'Preprocessing took {int(np.round(time.time() - start_time))} seconds')

        # Calculate the memory usage of the transformed dataframe
        st.session_state['input_dataframe_memory_usage_bytes'] = df.memory_usage(deep=True).sum()

        # Update the preprocessing parameters
        st.session_state['input_metadata']['preprocessing'] = {
            'location': 'Radial Profiles app',
            'nuclear_channel': nuclear_channel,
            'do_z_score_filter': do_z_score_filter,
        }
        if do_z_score_filter:
            st.session_state['input_metadata']['preprocessing']['z_score_filter_threshold'] = z_score_filter_threshold

        # Display information about the new dataframe
        df.info()

        # In case df has been modified not-in-place in any way, reassign the input dataset as the modified df
        st.session_state['input_dataset'].data = df

    # Return the modified dataframe
    return df


# Function to update the color for a value
def update_color_for_value(value_to_change_color):
    st.session_state[st_key_prefix + 'color_dict'][value_to_change_color] = st.session_state[st_key_prefix + 'new_picked_color']


# Function to reset the color dictionary
def reset_color_dict(ser_to_plot):
    # Create a color sequence based on the frequency of the values to plot in the entire dataset
    values_to_plot = ser_to_plot.value_counts().index
    colors = list(islice(cycle(px.colors.qualitative.Plotly), len(values_to_plot)))
    st.session_state[st_key_prefix + 'values_to_plot'] = values_to_plot
    st.session_state[st_key_prefix + 'color_dict'] = dict(zip(values_to_plot, colors))  # map values to colors


def go_to_previous_image(unique_images):
    """
    Go to the previous image in the numpy array.

    Parameters:
    unique_images (numpy.ndarray): The unique images.

    Returns:
    None
    """

    # Get the current index in the unique images
    key = st_key_prefix + 'image_to_view'
    current_index = list(unique_images).index(st.session_state[key])

    # If we're not already at the first image, go to the previous image
    if current_index > 0:
        current_index -= 1
        st.session_state[key] = unique_images[current_index]


def go_to_next_image(unique_images):
    """
    Go to the next image in the numpy array.

    Parameters:
    unique_images (numpy.ndarray): The unique images.

    Returns:
    None
    """

    # Get the current index in the unique images
    key = st_key_prefix + 'image_to_view'
    current_index = list(unique_images).index(st.session_state[key])

    # If we're not already at the last image, go to the next image
    if current_index < len(unique_images) - 1:
        current_index += 1
        st.session_state[key] = unique_images[current_index]


def main():
    """
    Main function for the page.
    """

    # Ensure a dataset has been opened in the first place
    if 'input_dataset' not in st.session_state:
        st.warning('Please open a dataset from the Open File page at left.')
        return

    # Save a shortcut to the dataframe
    df = st.session_state['input_dataset'].data

    # Set up some columns
    columns = st.columns(3)

    # Set up preprocessing
    with columns[0]:
        df = initialize_preprocessing(df)

    # Set up calculation of radial bins
    with columns[1]:
        df, unique_images, coordinate_scale_factor, annulus_spacing_um, xy_position_columns = initialize_radial_bin_calculation(df)

    # Draw a divider
    st.divider()

    # Main settings section
    column_to_plot, image_to_view, use_coordinate_mins_and_maxs, xmin_col, xmax_col, ymin_col, ymax_col, units, invert_y_axis, opacity, color_dict, values_to_plot, categorical_columns = initialize_main_settings(df, unique_images)

    # Output/plutting section
    st.divider()

    # Draw a single image scatter plot
    draw_single_image_scatter_plot(df, image_to_view, column_to_plot, values_to_plot, color_dict, xy_position_columns=xy_position_columns, coordinate_scale_factor=coordinate_scale_factor, annulus_spacing_um=annulus_spacing_um, use_coordinate_mins_and_maxs=use_coordinate_mins_and_maxs, xmin_col=xmin_col, xmax_col=xmax_col, ymin_col=ymin_col, ymax_col=ymax_col, units=units, invert_y_axis=invert_y_axis, opacity=opacity)

    ####

    st.header('Adaptive phenotyping')

    columns = st.columns(3)

    with columns[0]:

        st.subheader(':one: Threshold calculation')

        # Select columns to use for grouping the threshold calculations
        key = st_key_prefix + 'columns_for_phenotype_grouping'
        if key not in st.session_state:
            st.session_state[key] = []
        columns_for_phenotype_grouping = st.multiselect('Columns for phenotype grouping:', categorical_columns, key=key)

        # Set the column name that describes the baseline field such as cell type
        key = st_key_prefix + 'column_identifying_baseline_signal'
        if key not in st.session_state:
            st.session_state[key] = categorical_columns[0]
        column_identifying_baseline_signal = st.selectbox('Column identifying baseline/signal:', categorical_columns, key=key, on_change=delete_session_state_variable, args=('value_identifying_baseline',))

        # Extract the baseline field value (such as a specific cell type) to be used for determining the thresholds
        available_baseline_signal_values = df[column_identifying_baseline_signal].unique()
        key = st_key_prefix + 'value_identifying_baseline'
        if key not in st.session_state:
            st.session_state[key] = available_baseline_signal_values[0]
        value_identifying_baseline = st.selectbox('Value identifying baseline:', available_baseline_signal_values, key=key)

        # Extract the available channels for performing phenotyping
        key = st_key_prefix + 'channel_for_phenotyping'
        if key not in st.session_state:
            st.session_state[key] = df.columns[0]
        channel_for_phenotyping = st.selectbox('Channel for phenotyping:', df.columns, key=key)

        # If adaptive thresholding is desired...
        if st.button('Calculate thresholds for phenotyping'):

            # Group the relevant subset of the dataframe by the selected variables
            if len(columns_for_phenotype_grouping) > 0:
                df_grouped = df[columns_for_phenotype_grouping + [column_identifying_baseline_signal, channel_for_phenotyping, 'Slide ID']].groupby(by=columns_for_phenotype_grouping)
            else:
                df_grouped = [(None, df)]

            # Define various z scores of interest to use for determining the thresholds
            z_scores = np.arange(-1, 11)

            # For every group of variables for which you want a different threshold...
            thresholds_outer = []
            for _, df_group in df_grouped:

                # Get the selected intensity data for just the baseline field value
                ser_baseline = df_group.loc[df_group[column_identifying_baseline_signal] == value_identifying_baseline, channel_for_phenotyping]

                # From that, calculate the mean and std
                mean_to_use = ser_baseline.mean()
                std_to_use = ser_baseline.std()

                # Determine the corresponding threshold for each z score
                thresholds_inner = []
                for z_score in z_scores:
                    thresholds_inner.append(mean_to_use + z_score * std_to_use)

                # Add to the main thresholds holder
                thresholds_outer.append(thresholds_inner)

            # Determine the multi-index for the dataframe
            if len(columns_for_phenotype_grouping) == 0:
                index = pd.Index([-1])
            elif len(columns_for_phenotype_grouping) == 1:
                index = pd.Index(list(df_grouped.groups.keys()), name=columns_for_phenotype_grouping[0])
            elif len(columns_for_phenotype_grouping) > 1:
                index = pd.MultiIndex.from_tuples(list(df_grouped.groups.keys()), names=columns_for_phenotype_grouping)
            index.name = 'Grouping'

            # Determine the columns index for the dataframe
            columns_index = pd.Index([f'z score = {z_score}' for z_score in z_scores])
            columns_index.name = 'Thresholds'

            # Set the dataframe of thresholds
            df_thresholds = pd.DataFrame(thresholds_outer, columns=columns_index, index=index).sort_index()
            st.session_state[st_key_prefix + 'df_thresholds'] = df_thresholds

            # Save the grouped data as well for plotting and phenotyping
            st.session_state[st_key_prefix + 'df_grouped'] = df_grouped

        # Make sure the phenotyping thresholds have been calculated
        key = st_key_prefix + 'df_thresholds'
        if key not in st.session_state:
            st.warning('Please calculate thresholds for phenotyping')
            return

        # Set a shortcut to the relevant data
        df_thresholds = st.session_state[key]
        df_grouped = st.session_state[st_key_prefix + 'df_grouped']

        # Display the thresholds, allowing the user to select a single row
        st.write('Calculated thresholds:')
        group_selection = st.dataframe(df_thresholds, on_select='rerun', selection_mode='single-row')

        # Create a new figure
        fig = go.Figure()

        # Loop through each column in df_thresholds to add a trace for each one
        for column in df_thresholds.columns:
            fig.add_trace(go.Scatter(
                # x=df_thresholds.index,  # Use the DataFrame index for the x-axis
                x=np.array(range(len(df_thresholds))),
                y=df_thresholds[column],  # Column values for the y-axis
                mode='markers+lines',  # Line plot
                name=column  # Use the column name as the trace name
            ))

        # Update the layout to add titles and adjust other aesthetics if needed
        fig.update_layout(
            title='Line Plot of Thresholds Dataframe Above',
            xaxis_title='Index',
            yaxis_title='Phenotyping channel',
            legend_title='Column'
        )

        # Plot the line plots
        st.plotly_chart(fig)

    with columns[1]:

        st.subheader(':two: Percent positives plotting')

        # Extract the value identifying the signal (such as a specific cell type) to be used for testing the thresholds
        key = st_key_prefix + 'value_identifying_signal'
        if key not in st.session_state:
            st.session_state[key] = available_baseline_signal_values[0]
        value_identifying_signal = st.selectbox('Value identifying signal:', available_baseline_signal_values, key=key)

        # Whether to apply the threshold to just the selected group
        key = st_key_prefix + 'apply_thresh_to_selected_group'
        if key not in st.session_state:
            st.session_state[key] = True
        apply_thresh_to_selected_group = st.checkbox('Apply threshold to each individual group (instead of to the entire dataset)', key=key)

        # Whether to average over all groups
        key = st_key_prefix + 'average_over_all_groups'
        if key not in st.session_state:
            st.session_state[key] = False
        average_over_all_groups = st.checkbox('Average over all groups', key=key)

        # If we want to generate a plot based on the currently selected group...
        if not average_over_all_groups:

            # Obtain the current selection
            row_selection_list = group_selection['selection']['rows']

            # If something is actually selected...
            if row_selection_list:

                # Get the box and whisker data
                fig, _ = get_box_and_whisker_data(df_grouped, df_thresholds, apply_thresh_to_selected_group, df, channel_for_phenotyping, column_identifying_baseline_signal, value_identifying_baseline, value_identifying_signal, row_selection_list[0], return_figure_and_summary=True)

                # Plot the box and whisker chart
                st.plotly_chart(fig)

        # If we want to generate a plot from averaging over all the groups...
        else:

            # Since this takes a non-trivial amount of time, hide the calculation behind a button
            if st.button('Calculate average positive percentages over all groups'):

                # Initialize the holders of the group indices and the average positive percentages
                index_holder = []
                ser_holder_baseline = []
                ser_holder_signal = []

                # For every group...
                for curr_row in range(len(df_thresholds)):

                    # Get the box and whisker data
                    df_summary, curr_index = get_box_and_whisker_data(df_grouped, df_thresholds, apply_thresh_to_selected_group, df, channel_for_phenotyping, column_identifying_baseline_signal, value_identifying_baseline, value_identifying_signal, curr_row, return_figure_and_summary=False)

                    # Append the desired data to the holders
                    df_summary = df_summary.drop('Threshold', axis='columns').set_index('Z score')  # the drop isn't necessary but it may make the operations marginally faster
                    index_holder.append(curr_index)
                    ser_holder_baseline.append(df_summary['Positive % (avg. over images) for baseline group'])
                    ser_holder_signal.append(df_summary['Positive % (avg. over images) for signal group'])

                # Calculate and save the average positive percentages for the baseline group
                df_baseline = pd.concat(ser_holder_baseline, axis='columns', keys=index_holder)
                df_baseline.columns.names = columns_for_phenotype_grouping if len(columns_for_phenotype_grouping) > 0 else ['No group']
                st.session_state[st_key_prefix + 'df_baseline'] = df_baseline

                # Calculate and save the average positive percentages for the signal group
                df_signal = pd.concat(ser_holder_signal, axis='columns', keys=index_holder)
                df_signal.columns.names = columns_for_phenotype_grouping if len(columns_for_phenotype_grouping) > 0 else ['No group']
                st.session_state[st_key_prefix + 'df_signal'] = df_signal

            # Ensure the average positive percentages have been calculated
            if st_key_prefix + 'df_baseline' not in st.session_state:
                st.warning('Please calculate average positive percentages over all groups')
                return

            # Get the difference in the average positive percentages between the baseline and signal groups
            df_baseline = st.session_state[st_key_prefix + 'df_baseline']
            df_signal = st.session_state[st_key_prefix + 'df_signal']
            df_diff = df_signal - df_baseline

            # Render some charts in Streamlit
            st.plotly_chart(plotly_mean_and_sem([df_baseline, df_signal], ['Baseline', 'Signal']))
            st.plotly_chart(plotly_mean_and_sem([df_diff], ['signal - baseline']))

    with columns[2]:

        st.subheader(':three: Phenotype generation')

        # Set the phenotype name
        key = st_key_prefix + 'phenotype_name'
        if key not in st.session_state:
            st.session_state[key] = ''
        phenotype_name = st.text_input('Phenotype name:', key=key)

        # Set the desired Z score
        key = st_key_prefix + 'desired_z_score'
        if key not in st.session_state:
            st.session_state[key] = 2.0
        desired_z_score = st.number_input('Desired Z score:', key=key)

        # Set the phenotyping method
        key = st_key_prefix + 'phenotyping_method'
        phenotyping_method_options = ["Selected threshold applied to entire dataset", "Group-specific threshold applied to each group"]
        if key not in st.session_state:
            st.session_state[key] = phenotyping_method_options[1]
        phenotyping_method = st.selectbox('Phenotyping method:', phenotyping_method_options, key=key, on_change=update_phenotyping_sub_options, args=(phenotyping_method_options,))

        # Obtain the current selection
        row_selection_list = group_selection['selection']['rows']

        # Determine whether the phenotyping button should be disabled
        phenotyping_button_disabled = False
        if (phenotyping_method == phenotyping_method_options[0]) and (not row_selection_list):
            st.warning(f'Phenotyping method "{phenotyping_method}" requires a group to be selected in the first column')
            phenotyping_button_disabled = True
        if (phenotyping_method == phenotyping_method_options[1]) and (row_selection_list):
            st.warning(f'Phenotyping method "{phenotyping_method}" does not actually require that no group be selected in the first column, but for clarity (since that method performs calculations for every group), we are enacting that requirement; please de-select the group selection in the first column')
            phenotyping_button_disabled = True

        # Render the button to perform the phenotyping
        if st.button('Perform phenotyping', disabled=phenotyping_button_disabled, on_click=update_phenotyping_sub_options, args=(phenotyping_method_options,)):

            # Perform the phenotyping method that applies the selected threshold to the entire dataset
            # Old: apply_thresh_to_selected_group=False, average_over_all_groups=False, DO have a particular group selected --> "Selected threshold applied to entire dataset"
            if phenotyping_method == phenotyping_method_options[0]:

                # Get the selected row of df_thresholds, which is a series
                ser_selected = df_thresholds.iloc[row_selection_list[0]]

                # Calculate the threshold from the mean/std from the selected group
                mean_for_zscore_calc = ser_selected.loc['z score = 0']
                std_for_zscore_calc = ser_selected.loc['z score = 1'] - mean_for_zscore_calc
                threshold = mean_for_zscore_calc + desired_z_score * std_for_zscore_calc

                # If the threshold for our desired Z score has already been calculated, enforce agreement
                if desired_z_score in np.arange(-1, 11, 1):
                    assert np.abs(threshold - ser_selected.loc[f'z score = {int(desired_z_score)}']) < 1e-8, 'The threshold calculated for the selected group does not match the threshold in the thresholds DataFrame'

                # Get the locations where the selected filtering column is at least the current threshold value and assign the integer version of this to a new series
                positive_loc = df[channel_for_phenotyping] >= threshold
                ser_phenotype = positive_loc.astype(int)

                # Output the acutally used threshold
                st.write(f'Threshold used for entire dataset: {threshold}')

            # Perform the phenotyping method that applies a group-specific threshold to each group
            # New: apply_thresh_to_selected_group=True, average_over_all_groups=True, DO NOT have a particular group selected --> "Group-specific threshold applied to each group"
            elif phenotyping_method == phenotyping_method_options[1]:

                # Initialize the phenotype column to all-negative
                ser_phenotype = pd.Series(-1, index=df.index)

                # For every group...
                thresholds = []
                for curr_row in range(len(df_thresholds)):

                    # Obtain the index and dataframe of the group identified by curr_row
                    if isinstance(df_grouped, list):
                        curr_index = None
                        curr_df = df_grouped[0][1]
                        curr_integer_indices_into_df = np.array(range(len(curr_df)))
                    else:
                        curr_index = df_thresholds.iloc[curr_row].name
                        curr_df = df_grouped.get_group(curr_index)
                        curr_integer_indices_into_df = df_grouped.indices[curr_index]  # this is just to ensure that the group is in the groupby object

                    # Get the locations of the images in the baseline group
                    images_in_plotting_group_1 = curr_df.loc[curr_df[column_identifying_baseline_signal] == value_identifying_baseline, 'Slide ID'].unique()
                    image_loc_group_1 = curr_df['Slide ID'].isin(images_in_plotting_group_1)

                    # Get only the data for the column of interest for the baseline images
                    ser_for_z_score = curr_df.loc[image_loc_group_1, channel_for_phenotyping]

                    # Use the mean and std from those data to calculate the desired threshold
                    mean_for_zscore_calc = ser_for_z_score.mean()
                    std_for_zscore_calc = ser_for_z_score.std()
                    threshold = mean_for_zscore_calc + desired_z_score * std_for_zscore_calc

                    # If the threshold for our desired Z score has already been calculated, enforce agreement
                    if desired_z_score in np.arange(-1, 11, 1):
                        assert np.abs(threshold - df_thresholds.iloc[curr_row].loc[f'z score = {int(desired_z_score)}']) < 1e-8, 'The threshold calculated for the current group does not match the threshold in the thresholds DataFrame'

                    # Get the locations where the selected filtering column is at least the current threshold value
                    positive_loc = curr_df[channel_for_phenotyping] >= threshold  # boolean series fitting the current group

                    # Update the phenotype assignments for the current group
                    ser_phenotype.iloc[curr_integer_indices_into_df] = positive_loc.astype(int)  # the LHS is the slice of ser_phenotype corresponding to the current group

                    # Store the threshold for the current group
                    thresholds.append(threshold)

                # Check that there are no values of -1 remaining in ser_phenotype
                assert -1 not in ser_phenotype, 'There are still cells that have not been assigned positivity, which shouldn\'t have happened'

                # Output the actually used thresholds
                thresholds = pd.Series(thresholds, index=df_thresholds.index)
                thresholds.name = 'Threshold'
                st.write('Thresholds used for each group in turn:')
                st.write(thresholds)

            # Add the phenotype column to the dataframe
            pheno_colname = f'Phenotype {phenotype_name}'
            df[pheno_colname] = utils.downcast_series_dtype(ser_phenotype)
            st.success(f'Phenotype column "{pheno_colname}" has been appended to the dataset')
            st.write(f'Number of cells in each phenotype group (0 = negative, 1 = positive):')
            st.write(df[pheno_colname].value_counts().reset_index(drop=True))

    ####


# Run the main function
if __name__ == '__main__':

    # Set page settings
    page_name = 'Radial Profiles - Calculations'
    st.set_page_config(layout='wide', page_title=page_name)
    st.title(page_name)

    # Run streamlit-dataframe-editor library initialization tasks at the top of the page
    st.session_state = sde.initialize_session_state(st.session_state)

    # Run Top of Page (TOP) functions
    st.session_state = top.top_of_page_reqs(st.session_state)

    # Call the main function
    main()

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)
