# Import relevant libraries
import streamlit as st
import app_top_of_page as top
import streamlit_dataframe_editor as sde
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
from itertools import cycle, islice


def update_color_for_value(value_to_change_color):
    st.session_state['rsp__color_dict'][value_to_change_color] = st.session_state['rsp__new_picked_color']


def reset_color_dict(ser_to_plot):
    # Create a color sequence based on the frequency of the values to plot in the entire dataset
    values_to_plot = ser_to_plot.value_counts().index
    colors = list(islice(cycle(px.colors.qualitative.Plotly), len(values_to_plot)))
    st.session_state['rsp__values_to_plot'] = values_to_plot
    st.session_state['rsp__color_dict'] = dict(zip(values_to_plot, colors))  # map values to colors


def go_to_previous_image(unique_images):
    """
    Go to the previous image in the numpy array.

    Parameters:
    unique_images (numpy.ndarray): The unique images.

    Returns:
    None
    """

    # Get the current index in the unique images
    current_index = list(unique_images).index(st.session_state['rsp__image_to_view'])

    # If we're not already at the first image, go to the previous image
    if current_index > 0:
        current_index -= 1
        st.session_state['rsp__image_to_view'] = unique_images[current_index]


def go_to_next_image(unique_images):
    """
    Go to the next image in the numpy array.

    Parameters:
    unique_images (numpy.ndarray): The unique images.

    Returns:
    None
    """

    # Get the current index in the unique images
    current_index = list(unique_images).index(st.session_state['rsp__image_to_view'])

    # If we're not already at the last image, go to the next image
    if current_index < len(unique_images) - 1:
        current_index += 1
        st.session_state['rsp__image_to_view'] = unique_images[current_index]


def draw_scatter_plot_with_options():

    # Define the main settings columns
    settings_columns_main = st.columns(3)

    # In the first column...
    with settings_columns_main[0]:

        # Allow user to select the dataframe containing the data to plot
        if 'rsp__data_to_plot' not in st.session_state:
            st.session_state['rsp__data_to_plot'] = 'Input data'
        data_to_plot = st.selectbox('Dataset containing plotting data:', ['Input data', 'Phenotyped data'], key='rsp__data_to_plot')
        input_dataset_has_changed = ('rsp__data_to_plot_prev' not in st.session_state) or (st.session_state['rsp__data_to_plot_prev'] != data_to_plot)
        st.session_state['rsp__data_to_plot_prev'] = data_to_plot

        # If they want to plot phenotyped data, ensure they've performed phenotyping
        if (data_to_plot == 'Input data') and ('input_dataset' not in st.session_state):
            st.warning('If you\'d like to plot the input data, please open a file first.')
            return

        # If they want to plot phenotyped data, ensure they've performed phenotyping
        if (data_to_plot == 'Phenotyped data') and (len(st.session_state['df']) == 1):
            st.warning('If you\'d like to plot the phenotyped data, please perform phenotyping first.')
            return
        
        # Set the shortcut to the dataframe of interest
        if data_to_plot == 'Input data':
            df = st.session_state['input_dataset'].data
        else:
            df = st.session_state['df']

        # Store columns of certain types
        if ('rsp__categorical_columns' not in st.session_state) or input_dataset_has_changed:
            max_num_unique_values = 1000
            categorical_columns = []
            for col in df.select_dtypes(include=('category', 'object')).columns:
                if df[col].nunique() <= max_num_unique_values:
                    categorical_columns.append(col)
            st.session_state['rsp__categorical_columns'] = categorical_columns
        if ('rsp__numeric_columns' not in st.session_state) or input_dataset_has_changed:
            st.session_state['rsp__numeric_columns'] = df.select_dtypes(include='number').columns
        categorical_columns = st.session_state['rsp__categorical_columns']
        numeric_columns = st.session_state['rsp__numeric_columns']

        # Choose a column to plot
        if ('rsp__column_to_plot' not in st.session_state) or input_dataset_has_changed:
            st.session_state['rsp__column_to_plot'] = categorical_columns[0]
        column_to_plot = st.selectbox('Select a column by which to color the points:', categorical_columns, key='rsp__column_to_plot')
        column_to_plot_has_changed = ('rsp__column_to_plot_prev' not in st.session_state) or (st.session_state['rsp__column_to_plot_prev'] != column_to_plot) or input_dataset_has_changed
        st.session_state['rsp__column_to_plot_prev'] = column_to_plot

        # Get some information about the images in the input dataset
        if input_dataset_has_changed:
            st.session_state['rsp__unique_images'] = df['Slide ID'].unique()  # get the unique images in the dataset
            st.session_state['rsp__ser_size_of_each_image'] = df['Slide ID'].value_counts()  # calculate the number of objects in each image
        unique_images = st.session_state['rsp__unique_images']
        ser_size_of_each_image = st.session_state['rsp__ser_size_of_each_image']

        # Create an image selection selectbox
        if 'rsp__image_to_view' not in st.session_state:
            st.session_state['rsp__image_to_view'] = unique_images[0]
        image_to_view = st.selectbox('Select image to view:', unique_images, key='rsp__image_to_view')

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
        if 'rsp__use_coordinate_mins_and_maxs' not in st.session_state:
            st.session_state['rsp__use_coordinate_mins_and_maxs'] = False
        use_coordinate_mins_and_maxs = st.checkbox('Use coordinate mins and maxs', key='rsp__use_coordinate_mins_and_maxs')
        settings_columns_refined = st.columns(2)
        if 'rsp__x_min_coordinate_column' not in st.session_state:
            st.session_state['rsp__x_min_coordinate_column'] = numeric_columns[0]
        if 'rsp__y_min_coordinate_column' not in st.session_state:
            st.session_state['rsp__y_min_coordinate_column'] = numeric_columns[0]
        if 'rsp__x_max_coordinate_column' not in st.session_state:
            st.session_state['rsp__x_max_coordinate_column'] = numeric_columns[0]
        if 'rsp__y_max_coordinate_column' not in st.session_state:
            st.session_state['rsp__y_max_coordinate_column'] = numeric_columns[0]
        with settings_columns_refined[0]:
            xmin_col = st.selectbox('Select a column for the minimum x-coordinate:', numeric_columns, key='rsp__x_min_coordinate_column', disabled=(not use_coordinate_mins_and_maxs))
        with settings_columns_refined[1]:
            xmax_col = st.selectbox('Select a column for the maximum x-coordinate:', numeric_columns, key='rsp__x_max_coordinate_column', disabled=(not use_coordinate_mins_and_maxs))
        with settings_columns_refined[0]:
            ymin_col = st.selectbox('Select a column for the minimum y-coordinate:', numeric_columns, key='rsp__y_min_coordinate_column', disabled=(not use_coordinate_mins_and_maxs))
        with settings_columns_refined[1]:
            ymax_col = st.selectbox('Select a column for the maximum y-coordinate:', numeric_columns, key='rsp__y_max_coordinate_column', disabled=(not use_coordinate_mins_and_maxs))
        units = ('coordinate units' if use_coordinate_mins_and_maxs else 'microns')

        # Optionally add another filter
        if ('rsp__add_another_filter' not in st.session_state) or input_dataset_has_changed:
            st.session_state['rsp__add_another_filter'] = False
        if ('rsp__column_to_filter_by' not in st.session_state) or input_dataset_has_changed:
            st.session_state['rsp__column_to_filter_by'] = categorical_columns[0]
        if ('rsp__values_to_filter_by' not in st.session_state) or input_dataset_has_changed:
            st.session_state['rsp__values_to_filter_by'] = []
        st.checkbox('Add filter', key='rsp__add_another_filter')
        st.selectbox('Select a column to filter by:', categorical_columns, key='rsp__column_to_filter_by', disabled=(not st.session_state['rsp__add_another_filter']))
        st.multiselect('Select values to filter by:', df[st.session_state['rsp__column_to_filter_by']].unique(), key='rsp__values_to_filter_by', disabled=(not st.session_state['rsp__add_another_filter']))
        add_another_filter = st.session_state['rsp__add_another_filter']
        column_to_filter_by = st.session_state['rsp__column_to_filter_by']
        values_to_filter_by = st.session_state['rsp__values_to_filter_by']

    # In the third column...
    with settings_columns_main[2]:

        # Add an option to invert the y-axis
        if 'rsp__invert_y_axis' not in st.session_state:
            st.session_state['rsp__invert_y_axis'] = False
        invert_y_axis = st.checkbox('Invert y-axis', key='rsp__invert_y_axis')

        # Choose the opacity of objects
        if 'rsp__opacity' not in st.session_state:
            st.session_state['rsp__opacity'] = 0.7
        opacity = st.number_input('Opacity:', min_value=0.0, max_value=1.0, step=0.1, key='rsp__opacity')

        # Define the colors for the values to plot
        if ('rsp__color_dict' not in st.session_state) or column_to_plot_has_changed:
            reset_color_dict(df[column_to_plot])
        values_to_plot = st.session_state['rsp__values_to_plot']
        color_dict = st.session_state['rsp__color_dict']

        # Select a value whose color we want to modify
        if ('rsp__value_to_change_color' not in st.session_state) or column_to_plot_has_changed:
            st.session_state['rsp__value_to_change_color'] = values_to_plot[0]
        value_to_change_color = st.selectbox('Value whose color to change:', values_to_plot, key='rsp__value_to_change_color')

        # Create a color picker widget for the selected value
        st.session_state['rsp__new_picked_color'] = color_dict[value_to_change_color]
        st.color_picker('Pick a new color:', key='rsp__new_picked_color', on_change=update_color_for_value, args=(value_to_change_color,))

        # Add a button to reset the colors to their default values
        st.button('Reset plotting colors to defaults', on_click=reset_color_dict, args=(df[column_to_plot],))
        color_dict = st.session_state['rsp__color_dict']

    # Draw a divider
    st.divider()

    # If the user wants to display the scatter plot, indicated by a toggle...
    if 'rsp__show_scatter_plot' not in st.session_state:
        st.session_state['rsp__show_scatter_plot'] = False
    if st.toggle('Show scatter plot', key='rsp__show_scatter_plot'):

        # Optionally set up another filter
        if add_another_filter:
            filter_loc = df[column_to_filter_by].isin(values_to_filter_by)
        else:
            filter_loc = pd.Series(True, index=df.index)

        # Filter the DataFrame to include only the selected image and filter
        df_selected_image_and_filter = df[(df['Slide ID'] == image_to_view) & filter_loc]

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

    # Return necessary variables
    return df, column_to_plot, values_to_plot, categorical_columns, unique_images


def main():
    """
    Main function for the page.
    """

    df, column_to_plot, values_to_plot, categorical_columns, unique_images = draw_scatter_plot_with_options()

    if 'rsp__get_percent_frequencies' not in st.session_state:
        st.session_state['rsp__get_percent_frequencies'] = False
    if st.toggle('Get percent frequencies of coloring column for entire dataset', key='rsp__get_percent_frequencies'):
        vc = df[column_to_plot].value_counts()
        st.dataframe((vc / vc.sum() * 100).astype(int).reset_index())

    if 'rsp__generate_box_and_whisker_plot' not in st.session_state:
        st.session_state['rsp__generate_box_and_whisker_plot'] = False
    if st.toggle('Generate box and whisker plot', key='rsp__generate_box_and_whisker_plot'):
        with st.columns(3)[0]:
            if 'rsp__box_and_whisker_plot_value' not in st.session_state:
                st.session_state['rsp__box_and_whisker_plot_value'] = values_to_plot[0]
            box_and_whisker_plot_value = st.selectbox('Select a value in the selected coloring column to analyze:', values_to_plot, key='rsp__box_and_whisker_plot_value')
            if 'rsp__column_identifying_trace' not in st.session_state:
                st.session_state['rsp__column_identifying_trace'] = categorical_columns[0]
            column_identifying_trace = st.selectbox('Select a column to identify different values/traces to plot:', categorical_columns, key='rsp__column_identifying_trace')
            match_loc = df[column_to_plot] == box_and_whisker_plot_value
            percent_match_holder = []
            trace_value_holder = []
            for image in unique_images:
                image_loc = df['Slide ID'] == image
                set_trace_values = set(df.loc[image_loc, column_identifying_trace])
                assert len(set_trace_values) == 1, 'There should only be one value for the column identifying the trace'
                trace_value_holder.append(set_trace_values.pop())
                percent_match_holder.append((image_loc & match_loc).sum() / image_loc.sum() * 100)
            df_boxplot = pd.DataFrame({'Image': unique_images, 'Percent': percent_match_holder, 'Trace': trace_value_holder})
            fig = px.box(df_boxplot, x='Trace', y='Percent', title=f'Box and whisker plot for {box_and_whisker_plot_value}', points='all')
            st.plotly_chart(fig, use_container_width=True)


# Run the main function
if __name__ == '__main__':

    # Set page settings
    page_name = 'Coordinate Scatter Plotter'
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
