# Import relevant libraries
import streamlit as st
import app_top_of_page as top
import streamlit_dataframe_editor as sde
import plotly.graph_objects as go
import plotly.express as px


# * Plot ellipses
# * Make the main plot a larger height
# * Remove the (plus) and (dash) symbols
# * Allow the user to customize the colors of the phenotypes
# * Adjust tooltip over each cell
# * More?


def go_to_previous_image(unique_images):
    """
    Go to the previous image in the list.
    """

    current_index = list(unique_images).index(st.session_state['rsp__image_to_view'])

    if current_index > 0:
        current_index -= 1
        st.session_state['rsp__image_to_view'] = unique_images[current_index]


def go_to_next_image(unique_images):
    """
    Go to the next image in the list.
    """

    current_index = list(unique_images).index(st.session_state['rsp__image_to_view'])

    if current_index < len(unique_images) - 1:
        current_index += 1
        st.session_state['rsp__image_to_view'] = unique_images[current_index]


def main():
    """
    Main function for the page.
    """

    df = st.session_state['df']

    unique_images = df['Slide ID'].unique()

    with st.columns(3)[0]:

        if 'rsp__image_to_view' not in st.session_state:
            st.session_state['rsp__image_to_view'] = unique_images[0]
        st.selectbox('Select image to view:', unique_images, key='rsp__image_to_view')
        cols = st.columns(2)
        with cols[0]:
            st.button('Previous', on_click=go_to_previous_image, args=(unique_images,), disabled=(st.session_state['rsp__image_to_view'] == unique_images[0]), use_container_width=True)
        with cols[1]:
            st.button('Next', on_click=go_to_next_image, args=(unique_images, ), disabled=(st.session_state['rsp__image_to_view'] == unique_images[-1]), use_container_width=True)

    if 'rsp__show_scatter_plot' not in st.session_state:
        st.session_state['rsp__show_scatter_plot'] = False
    if st.toggle('Show scatter plot', key='rsp__show_scatter_plot'):

        # Filter the DataFrame to include only the selected image
        df_selected_image = df[df['Slide ID'] == st.session_state['rsp__image_to_view']]

        # Create a color sequence based on the phenotype frequency in the entire dataset
        phenotypes = df['phenotype'].value_counts().index
        colors = px.colors.qualitative.Plotly[:len(phenotypes)]  # Get enough colors for all phenotypes
        color_dict = dict(zip(phenotypes, colors))  # Map phenotypes to colors

        # Create the scatter plot
        fig = go.Figure()

        # Group the DataFrame for the selected image by phenotype
        grouped = df_selected_image.groupby('phenotype')

        # Loop over the phenotypes in order of their frequency
        for phenotype in phenotypes:
            if phenotype in grouped.groups:  # Check if the phenotype exists in the selected image
                group = grouped.get_group(phenotype)
                fig.add_trace(go.Scatter(x=group['Cell X Position'], y=group['Cell Y Position'], mode='markers', name=phenotype, marker_color=color_dict[phenotype]))

        # Update the layout once, after all traces have been added
        fig.update_layout(
            xaxis=dict(
                scaleanchor="y",
                scaleratio=1,
            ),
            title=f'Scatter plot for {st.session_state["rsp__image_to_view"]}',
            xaxis_title='Cell X Position (microns)',
            yaxis_title='Cell Y Position (microns)',
            legend_title='Phenotype'
        )

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
