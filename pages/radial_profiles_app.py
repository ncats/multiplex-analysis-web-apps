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

# Global variable
st_key_prefix = 'radial_profiles__'


def initialize_preprocessing(df):

    st.header('Preprocessing')

    # On the first third vertical third of the page...
    with st.columns(3)[0]:
    
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

    return df


def update_color_for_value(value_to_change_color):
    st.session_state[st_key_prefix + 'color_dict'][value_to_change_color] = st.session_state[st_key_prefix + 'new_picked_color']


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

    df = initialize_preprocessing(df)

    st.divider()

    # Define the main settings columns
    settings_columns_main = st.columns(3)

    # In the first column...
    with settings_columns_main[0]:

        # Store columns of certain types
        if st_key_prefix + 'categorical_columns' not in st.session_state:
            max_num_unique_values = 1000
            categorical_columns = []
            for col in df.select_dtypes(include=('category', 'object')).columns:
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
        if st_key_prefix + 'unique_images' not in st.session_state:
            st.session_state[st_key_prefix + 'unique_images'] = df['Slide ID'].unique()  # get the unique images in the dataset
        if st_key_prefix + 'ser_size_of_each_image' not in st.session_state:
            st.session_state[st_key_prefix + 'ser_size_of_each_image'] = df['Slide ID'].value_counts()  # calculate the number of objects in each image
        unique_images = st.session_state[st_key_prefix + 'unique_images']
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

        # Optionally add another filter
        if st_key_prefix + 'add_another_filter' not in st.session_state:
            st.session_state[st_key_prefix + 'add_another_filter'] = False
        if st_key_prefix + 'column_to_filter_by' not in st.session_state:
            st.session_state[st_key_prefix + 'column_to_filter_by'] = categorical_columns[0]
        if st_key_prefix + 'values_to_filter_by' not in st.session_state:
            st.session_state[st_key_prefix + 'values_to_filter_by'] = []
        add_another_filter = st.checkbox('Add filter', key=st_key_prefix + 'add_another_filter')
        column_to_filter_by = st.selectbox('Select a column to filter by:', categorical_columns, key=st_key_prefix + 'column_to_filter_by', disabled=(not add_another_filter))
        values_to_filter_by = st.multiselect('Select values to filter by:', df[column_to_filter_by].unique(), key=st_key_prefix + 'values_to_filter_by', disabled=(not add_another_filter))

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

    # Draw a divider
    st.divider()

    # If the user wants to display the scatter plot, indicated by a toggle...
    if st_key_prefix + 'show_scatter_plot' not in st.session_state:
        st.session_state[st_key_prefix + 'show_scatter_plot'] = False
    if st.toggle('Show scatter plot', key=st_key_prefix + 'show_scatter_plot'):

        # Optionally set up another filter
        if add_another_filter:
            filter_loc = df[column_to_filter_by].isin(values_to_filter_by)
        else:
            filter_loc = pd.Series(True, index=df.index)

        # Filter the DataFrame to include only the selected image and filter
        df_selected_image_and_filter = df[(df['Slide ID'] == image_to_view) & filter_loc]



        #### This is only temporary and is due to not performing the conversion in the datafile unifier this time ####
        df_selected_image_and_filter[['Cell X Position', 'Cell Y Position']] = df_selected_image_and_filter[['Cell X Position', 'Cell Y Position']] * 0.32



        # Get the x-y midpoint of the coordinates in df_selected_image_and_filter[['Cell X Position', 'Cell Y Position']]
        xy_min = df_selected_image_and_filter[['Cell X Position', 'Cell Y Position']].min()
        xy_max = df_selected_image_and_filter[['Cell X Position', 'Cell Y Position']].max()
        xy_mid = (xy_min + xy_max) / 2

        # Get the radius edges that fit within the largest possible radius
        largest_possible_radius = (xy_max - xy_mid).min()
        spacing_um = 250
        num_intervals = largest_possible_radius // spacing_um  # calculate the number of intervals that fit within the largest_possible_radius
        end_value = (num_intervals + 1) * spacing_um  # calculate the end value for np.arange to ensure it does not exceed largest_possible_radius
        radius_edges = np.arange(0, end_value, spacing_um)  # generate steps


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
                    fig.add_trace(go.Scatter(x=df_group['Cell X Position'], y=df_group['Cell Y Position'], mode='markers', name=value_str_cleaned, marker_color=color_dict[value_to_plot], hovertemplate=df_group['hover_label']))

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
                x0=xy_mid['Cell X Position'] - radius,
                y0=xy_mid['Cell Y Position'] - radius,
                x1=xy_mid['Cell X Position'] + radius,
                y1=xy_mid['Cell Y Position'] + radius,
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

    # We seem to need to render something on the page after rendering a plotly figure in order for the page to not automatically scroll back to the top when you go to the Previous or Next image
    st.write(' ')


# Run the main function
if __name__ == '__main__':

    # Set page settings
    page_name = 'Radial Profiles'
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
