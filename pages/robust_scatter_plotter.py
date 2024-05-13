# Import relevant libraries
import streamlit as st
import app_top_of_page as top
import streamlit_dataframe_editor as sde
import plotly.graph_objects as go
import plotly.express as px


# * Allow the user to customize the colors of the phenotypes
# * Improve tooltip over each cell
# * Add parameters for all remaining things that can be parametrized below
# * More?


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

    settings_columns = st.columns(3)

    # Shrink widgets to 1/3 of the page width
    with settings_columns[0]:

        # Definitely calculate and optionally show the number of objects in each image
        with st.expander('Size of each image:', expanded=False):
            ser_size_of_each_image = df['Slide ID'].value_counts()
            st.write(ser_size_of_each_image)

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

        # Add an option to invert the y-axis
        if 'rsp__invert_y_axis' not in st.session_state:
            st.session_state['rsp__invert_y_axis'] = False
        st.checkbox('Invert y-axis', key='rsp__invert_y_axis')

        if 'rsp__opacity' not in st.session_state:
            st.session_state['rsp__opacity'] = 0.7
        st.number_input('Opacity:', min_value=0.0, max_value=1.0, step=0.1, key='rsp__opacity')

    if 'rsp__use_coordinate_mins_and_maxs' not in st.session_state:
        st.session_state['rsp__use_coordinate_mins_and_maxs'] = False
    st.checkbox('Use coordinate mins and maxs', key='rsp__use_coordinate_mins_and_maxs')

    settings_columns = st.columns(3)

    numeric_columns = df.select_dtypes(include='number').columns
    if 'rsp__x_min_coordinate_column' not in st.session_state:
        st.session_state['rsp__x_min_coordinate_column'] = numeric_columns[0]
    if 'rsp__y_min_coordinate_column' not in st.session_state:
        st.session_state['rsp__y_min_coordinate_column'] = numeric_columns[0]
    if 'rsp__x_max_coordinate_column' not in st.session_state:
        st.session_state['rsp__x_max_coordinate_column'] = numeric_columns[0]
    if 'rsp__y_max_coordinate_column' not in st.session_state:
        st.session_state['rsp__y_max_coordinate_column'] = numeric_columns[0]
    with settings_columns[0]:
        st.selectbox('Select a column for the minimum x-coordinate:', numeric_columns, key='rsp__x_min_coordinate_column', disabled=(not st.session_state['rsp__use_coordinate_mins_and_maxs']))
    with settings_columns[1]:
        st.selectbox('Select a column for the maximum x-coordinate:', numeric_columns, key='rsp__x_max_coordinate_column', disabled=(not st.session_state['rsp__use_coordinate_mins_and_maxs']))
    with settings_columns[0]:
        st.selectbox('Select a column for the minimum y-coordinate:', numeric_columns, key='rsp__y_min_coordinate_column', disabled=(not st.session_state['rsp__use_coordinate_mins_and_maxs']))
    with settings_columns[1]:
        st.selectbox('Select a column for the maximum y-coordinate:', numeric_columns, key='rsp__y_max_coordinate_column', disabled=(not st.session_state['rsp__use_coordinate_mins_and_maxs']))

    # If the user wants to display the scatter plot, indicated by a toggle...
    if 'rsp__show_scatter_plot' not in st.session_state:
        st.session_state['rsp__show_scatter_plot'] = True
    if st.toggle('Show scatter plot', key='rsp__show_scatter_plot'):

        # Filter the DataFrame to include only the selected image
        df_selected_image = df[df['Slide ID'] == st.session_state['rsp__image_to_view']]

        # Create a color sequence based on the phenotype frequency in the entire dataset
        phenotypes = df['phenotype'].value_counts().index
        colors = px.colors.qualitative.Plotly[:len(phenotypes)]  # get enough colors for all phenotypes
        color_dict = dict(zip(phenotypes, colors))  # map phenotypes to colors

        # Group the DataFrame for the selected image by phenotype
        selected_image_grouped_by_phenotype = df_selected_image.groupby('phenotype')

        # Create the scatter plot
        fig = go.Figure()

        # Loop over the phenotypes in order of their frequency
        for phenotype in phenotypes:

            # If the phenotype exists in the selected image...
            if phenotype in selected_image_grouped_by_phenotype.groups:

                # Store the dataframe for the current phenotype for the selected image
                df_group = selected_image_grouped_by_phenotype.get_group(phenotype)

                phenotype_str_cleaned = phenotype.replace('(plus)', '+').replace('(dash)', '-')

                df_group['hover_label'] = 'Phenotype: ' + phenotype_str_cleaned + '<br>Index: ' + df_group.index.astype(str)

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

        units = ('coordinate units' if st.session_state['rsp__use_coordinate_mins_and_maxs'] else 'microns')
        
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
            legend_title='Phenotype',
            height=800,  # Set the height of the figure
            width=800,  # Set the width of the figure
        )

        # Plot the plotly chart in Streamlit
        st.plotly_chart(fig, use_container_width=True)

        
# Run the main function
if __name__ == '__main__':

    # Set page settings
    page_name = 'Robust Scatter Plotter'
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
