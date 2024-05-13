# Import relevant libraries
import streamlit as st
import app_top_of_page as top
import streamlit_dataframe_editor as sde
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd


# TODO:
# Allow the user to customize the colors of the phenotypes


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


def main():
    """
    Main function for the page.
    """

    # Shortcut to the main, phenotyped dataframe
    df = st.session_state['df']

    # Get the unique images in the dataset
    unique_images = df['Slide ID'].unique()

    # Define the main settings columns
    settings_columns_main = st.columns(3)

    # In the first column...
    with settings_columns_main[0]:

        # Definitely calculate and optionally show the number of objects in each image
        ser_size_of_each_image = df['Slide ID'].value_counts()
        # with st.expander('Size of each image:', expanded=False):
        #     st.write(ser_size_of_each_image)

        # Create an image selection selectbox
        if 'rsp__image_to_view' not in st.session_state:
            st.session_state['rsp__image_to_view'] = unique_images[0]
        st.selectbox('Select image to view:', unique_images, key='rsp__image_to_view')

        # Display the number of cells in the selected image
        st.write(f'Number of cells in image: {ser_size_of_each_image.loc[st.session_state["rsp__image_to_view"]]}')

        # Optionally navigate through the images using Previous and Next buttons
        cols = st.columns(2)
        with cols[0]:
            st.button('Previous', on_click=go_to_previous_image, args=(unique_images,), disabled=(st.session_state['rsp__image_to_view'] == unique_images[0]), use_container_width=True)
        with cols[1]:
            st.button('Next', on_click=go_to_next_image, args=(unique_images, ), disabled=(st.session_state['rsp__image_to_view'] == unique_images[-1]), use_container_width=True)

        # Store columns of certain types
        categorical_columns = df.select_dtypes(include=('category', 'object')).columns
        numeric_columns = df.select_dtypes(include='number').columns

        # Choose a column to plot
        if 'rsp__column_to_plot' not in st.session_state:
            st.session_state['rsp__column_to_plot'] = categorical_columns[0]
        column_to_plot = st.selectbox('Select a column to plot:', categorical_columns, key='rsp__column_to_plot')

    # In the second column...
    with settings_columns_main[1]:

        # Optionally add another filter
        if 'rsp__add_another_filter' not in st.session_state:
            st.session_state['rsp__add_another_filter'] = False
        add_another_filter = st.checkbox('Add another filter', key='rsp__add_another_filter')
        if 'rsp__column_to_filter_by' not in st.session_state:
            st.session_state['rsp__column_to_filter_by'] = categorical_columns[0]
        column_to_filter_by = st.selectbox('Select a column to filter by:', categorical_columns, key='rsp__column_to_filter_by', disabled=(not st.session_state['rsp__add_another_filter']))
        if 'rsp__values_to_filter_by' not in st.session_state:
            st.session_state['rsp__values_to_filter_by'] = []
        values_to_filter_by = st.multiselect('Select values to filter by:', df[st.session_state['rsp__column_to_filter_by']].unique(), key='rsp__values_to_filter_by', disabled=(not st.session_state['rsp__add_another_filter']))

        # Add an option to invert the y-axis
        if 'rsp__invert_y_axis' not in st.session_state:
            st.session_state['rsp__invert_y_axis'] = False
        st.checkbox('Invert y-axis', key='rsp__invert_y_axis')

    # In the third column...
    with settings_columns_main[2]:
        
        # Choose the opacity of objects
        if 'rsp__opacity' not in st.session_state:
            st.session_state['rsp__opacity'] = 0.7
        st.number_input('Opacity:', min_value=0.0, max_value=1.0, step=0.1, key='rsp__opacity')

        # Optionally plot minimum and maximum coordinate fields
        if 'rsp__use_coordinate_mins_and_maxs' not in st.session_state:
            st.session_state['rsp__use_coordinate_mins_and_maxs'] = False
        st.checkbox('Use coordinate mins and maxs', key='rsp__use_coordinate_mins_and_maxs')
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
            st.selectbox('Select a column for the minimum x-coordinate:', numeric_columns, key='rsp__x_min_coordinate_column', disabled=(not st.session_state['rsp__use_coordinate_mins_and_maxs']))
        with settings_columns_refined[1]:
            st.selectbox('Select a column for the maximum x-coordinate:', numeric_columns, key='rsp__x_max_coordinate_column', disabled=(not st.session_state['rsp__use_coordinate_mins_and_maxs']))
        with settings_columns_refined[0]:
            st.selectbox('Select a column for the minimum y-coordinate:', numeric_columns, key='rsp__y_min_coordinate_column', disabled=(not st.session_state['rsp__use_coordinate_mins_and_maxs']))
        with settings_columns_refined[1]:
            st.selectbox('Select a column for the maximum y-coordinate:', numeric_columns, key='rsp__y_max_coordinate_column', disabled=(not st.session_state['rsp__use_coordinate_mins_and_maxs']))
        units = ('coordinate units' if st.session_state['rsp__use_coordinate_mins_and_maxs'] else 'microns')

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
        df_selected_image_and_filter = df[(df['Slide ID'] == st.session_state['rsp__image_to_view']) & filter_loc]

        # Create a color sequence based on the phenotype frequency in the entire dataset
        phenotypes = df[column_to_plot].value_counts().index
        colors = px.colors.qualitative.Plotly[:len(phenotypes)]  # get enough colors for all phenotypes
        color_dict = dict(zip(phenotypes, colors))  # map phenotypes to colors

        # Group the DataFrame for the selected image by phenotype
        selected_image_grouped_by_phenotype = df_selected_image_and_filter.groupby(column_to_plot)

        # Create the scatter plot
        fig = go.Figure()

        # Loop over the phenotypes in order of their frequency
        for phenotype in phenotypes:

            # If the phenotype exists in the selected image...
            if phenotype in selected_image_grouped_by_phenotype.groups:

                # Store the dataframe for the current phenotype for the selected image
                df_group = selected_image_grouped_by_phenotype.get_group(phenotype)

                # If phenotype is a string, replace '(plus)' with '+' and '(dash)' with '-'
                if isinstance(phenotype, str):
                    phenotype_str_cleaned = phenotype.replace('(plus)', '+').replace('(dash)', '-')

                # Add the object index to the label
                df_group['hover_label'] = 'Index: ' + df_group.index.astype(str)

                # Works but doesn't scale the shapes
                if not st.session_state['rsp__use_coordinate_mins_and_maxs']:
                    fig.add_trace(go.Scatter(x=df_group['Cell X Position'], y=df_group['Cell Y Position'], mode='markers', name=phenotype_str_cleaned, marker_color=color_dict[phenotype], hovertemplate=df_group['hover_label']))

                # Works really well
                else:
                    xmin_col = st.session_state['rsp__x_min_coordinate_column']
                    xmax_col = st.session_state['rsp__x_max_coordinate_column']
                    ymin_col = st.session_state['rsp__y_min_coordinate_column']
                    ymax_col = st.session_state['rsp__y_max_coordinate_column']
                    fig.add_trace(go.Bar(
                        x=((df_group[xmin_col] + df_group[xmax_col]) / 2),
                        y=df_group[ymax_col] - df_group[ymin_col],
                        width=df_group[xmax_col] - df_group[xmin_col],
                        base=df_group[ymin_col],
                        name=phenotype_str_cleaned,
                        marker=dict(
                            color=color_dict[phenotype],
                            opacity=st.session_state['rsp__opacity'],
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
                autorange=('reversed' if st.session_state['rsp__invert_y_axis'] else True),
            ),
            title=f'Scatter plot for {st.session_state["rsp__image_to_view"]}',
            xaxis_title=f'Cell X Position ({units})',
            yaxis_title=f'Cell Y Position ({units})',
            legend_title=column_to_plot,
            height=800,  # Set the height of the figure
            width=800,  # Set the width of the figure
        )

        # Plot the plotly chart in Streamlit
        st.plotly_chart(fig, use_container_width=True)

        
# Run the main function
if __name__ == '__main__':

    # Set page settings
    page_name = 'Scatter Plotter'
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
