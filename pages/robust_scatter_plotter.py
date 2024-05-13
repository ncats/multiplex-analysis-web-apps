# Import relevant libraries
import streamlit as st
import app_top_of_page as top
import streamlit_dataframe_editor as sde
import plotly.graph_objects as go
import plotly.express as px
import time
import pandas as pd


# * Plot ellipses
# * Remove the (plus) and (dash) symbols
# * Allow the user to customize the colors of the phenotypes
# * Adjust tooltip over each cell
# * Option to flip y-axis
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

    # st.write() a series containing the number of rows in each image
    if True:
        st.write(df['Slide ID'].value_counts())

    # df['Width'] = df['XMax'] - df['XMin']
    # df['Height'] = df['YMax'] - df['YMin']
    # df['XCenter'] = df['XMin'] + df['Width'] / 2
    # df['YCenter'] = df['YMin'] + df['Height'] / 2

    # start_time = time.time()

    # df['Path'] = 'M ' + (df['XCenter'] - df['Width'] / 2).astype(str) + ',' + df['YCenter'].astype(str) + ' a ' + (df['Width'] / 2).astype(str) + ',' + (df['Height'] / 2).astype(str) + ' 0 1,0 ' + df['Width'].astype(str) + ',0 a ' + (df['Width'] / 2).astype(str) + ',' + (df['Height'] / 2).astype(str) + ' 0 1,0 -' + df['Width'].astype(str) + ',0'

    # # # Create the components of the path as separate columns
    # # df['M'] = 'M'
    # # df['X1'] = df['XCenter'] - df['Width'] / 2
    # # df['Y1'] = df['YCenter']
    # # df['A1'] = 'a'
    # # df['Width1'] = df['Width'] / 2
    # # df['Height1'] = df['Height'] / 2
    # # df['Zero1'] = '0'
    # # df['One1'] = '1,0'
    # # df['Width2'] = df['Width']
    # # df['Zero2'] = ',0'
    # # df['A2'] = 'a'
    # # df['Width3'] = df['Width'] / 2
    # # df['Height2'] = df['Height'] / 2
    # # df['Zero3'] = '0'
    # # df['One2'] = '1,0'
    # # df['Width4'] = '-' + df['Width'].astype(str)
    # # df['Zero4'] = ',0'
    # # # Concatenate the columns to create the path
    # # df['Path'] = df[['M', 'X1', 'Y1', 'A1', 'Width1', 'Height1', 'Zero1', 'One1', 'Width2', 'Zero2', 'A2', 'Width3', 'Height2', 'Zero3', 'One2', 'Width4', 'Zero4']].astype(str).agg(' '.join, axis=1)

    # st.write(f'Creating the path took {time.time() - start_time:.2f} seconds.')

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
        # df_selected_image = df_selected_image.iloc[:100, :]

        # Create a color sequence based on the phenotype frequency in the entire dataset
        phenotypes = df['phenotype'].value_counts().index
        colors = px.colors.qualitative.Plotly[:len(phenotypes)]  # Get enough colors for all phenotypes
        color_dict = dict(zip(phenotypes, colors))  # Map phenotypes to colors

        # Create the scatter plot
        fig = go.Figure()

        # Group the DataFrame for the selected image by phenotype
        grouped = df_selected_image.groupby('phenotype')

        # Create a list of ellipses
        ellipses = []

        # Loop over the phenotypes in order of their frequency
        for phenotype in phenotypes:
            if phenotype in grouped.groups:  # Check if the phenotype exists in the selected image
                group = grouped.get_group(phenotype)

                # # Works but doesn't scale the shapes
                # fig.add_trace(go.Scatter(x=group['Cell X Position'], y=group['Cell Y Position'], mode='markers', name=phenotype, marker_color=color_dict[phenotype]))

                fig.add_trace(go.Bar(
                    x=group['XMin'],
                    y=group['YMax'] - group['YMin'],
                    width=group['XMax'] - group['XMin'],
                    base=group['YMin'],
                    name=phenotype,
                    marker=dict(
                        color=color_dict[phenotype],
                        opacity=0.5,
                        # line=dict(color='black', width=1)  # Add outlines
                    )
                    # marker_color=color_dict[phenotype],
                    # opacity=0.5,
                ))

                # # This never finished even with just 11K rows, maybe I should have tried again
                # for _, row in group.iterrows():
                #     fig.add_shape(
                #         type="path",
                #         path=row['Path'],
                #         line=dict(color=color_dict[phenotype]),
                #     )

                # # This works but is really slow, with type either "circle" or "rect"
                # # ----
                # # irow = 0
                # for _, row in group.iterrows():
                #     # print(irow)                    
                #     ellipse = dict(
                #         type="rect",
                #         xref="x",
                #         yref="y",
                #         x0=row['XMin'],
                #         y0=row['YMin'],
                #         x1=row['XMax'],
                #         y1=row['YMax'],
                #         # opacity=0.5,
                #         opacity=1,
                #         # fillcolor="blue",
                #         # line_color="blue",
                #         fillcolor=color_dict[phenotype],
                #         line_color=color_dict[phenotype],
                #     )
                #     ellipses.append(ellipse)
                #     # irow = irow + 1

        # Add the list of ellipses to the layout
        fig.update_layout(shapes=ellipses)
            # # ----

        # Update the layout once, after all traces have been added
        fig.update_layout(
            xaxis=dict(
                scaleanchor="y",
                scaleratio=1,
                # range=[df_selected_image['XMin'].min(), df_selected_image['XMax'].max()],  # Set the range of x-axis
            ),
            yaxis=dict(
                # range=[df_selected_image['YMin'].min(), df_selected_image['YMax'].max()],  # Set the range of y-axis
            ),
            title=f'Scatter plot for {st.session_state["rsp__image_to_view"]}',
            xaxis_title='Cell X Position (microns)',
            yaxis_title='Cell Y Position (microns)',
            legend_title='Phenotype',
            height=800,  # Set the height of the figure
            width=800,  # Set the width of the figure
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
