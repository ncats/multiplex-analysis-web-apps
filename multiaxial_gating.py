# Import relevant libraries
import pandas as pd
import seaborn as sns
import streamlit as st
import plotly.graph_objects as go
import os
import dataset_formats

# Function to load the data in a unified format
def load_data(input_datafile_path, coord_units_in_microns, dataset_format):
    dataset_class = getattr(dataset_formats, dataset_format)  # done this way so that the format (e.g., “REEC”) can be select programmatically
    dataset_obj = dataset_class(input_datafile=input_datafile_path, coord_units_in_microns=coord_units_in_microns)
    dataset_obj.process_dataset(do_calculate_minimum_coordinate_spacing_per_roi=False, do_trimming=False)
    return dataset_obj.data

# Update the dependencies of the selectbox for the current analysis column
def update_dependencies_of_column_to_plot():
    df = st.session_state['mg__df']
    column_to_plot = st.session_state['mg__column_to_plot']
    st.session_state['mg__curr_column_range'] = [df[column_to_plot].min(), df[column_to_plot].max()]
    st.session_state['mg__selected_min_val'] = st.session_state['mg__curr_column_range'][0]
    st.session_state['mg__selected_max_val'] = st.session_state['mg__curr_column_range'][1]

# Use a wide page layout
st.set_page_config(layout='wide')

# Constant
input_directory = os.path.join('.', 'input')

# Set datafile information --> turn into widgets soon
input_datafilename = 'measurementsEpCAMLy51MHCII-exported.csv'
coord_units_in_microns = 1

# Load the data only once, otherwise keep it in memory
if 'mg__df' not in st.session_state:
    st.session_state['mg__df'] = load_data(os.path.join(input_directory, input_datafilename), coord_units_in_microns, dataset_formats.extract_datafile_metadata(os.path.join(input_directory, input_datafilename))[4])
    unique_images = st.session_state['mg__df']['Slide ID'].unique()
    st.session_state['mg__unique_images_short'] = [x.split('-imagenum_')[1] for x in unique_images]
    st.session_state['mg__unique_image_dict'] = dict(zip(st.session_state['mg__unique_images_short'], unique_images))
df = st.session_state['mg__df']
unique_images_short = st.session_state['mg__unique_images_short']
unique_image_dict = st.session_state['mg__unique_image_dict']

# Define the main columns
main_columns = st.columns(3)

# Column filter
with main_columns[0]:

    # Column header
    st.header(':one: Column filter')

    # Have a dropdown for the column on which to perform a kernel density estimate
    st.selectbox(label='Column for filtering:', options=df.select_dtypes(include='number').columns, key='mg__column_to_plot', on_change=update_dependencies_of_column_to_plot)
    column_to_plot = st.session_state['mg__column_to_plot']

    # Initialize some things in the session state
    if 'mg__dfs_to_plot' not in st.session_state:
        st.session_state['mg__dfs_to_plot'] = dict()
    if 'mg__image_to_plot' not in st.session_state:
        st.session_state['mg__image_to_plot'] = unique_images_short[0]
    if 'mg__curr_column_range' not in st.session_state:
        update_dependencies_of_column_to_plot()
    if 'mg__df_current_phenotype' not in st.session_state:
        st.session_state['mg__df_current_phenotype'] = pd.DataFrame()
    if 'mg__df_phenotype_assignments' not in st.session_state:
        st.session_state['mg__df_phenotype_assignments'] = pd.DataFrame()
    if 'mg__current_phenotype_name' not in st.session_state:
        st.session_state['mg__current_phenotype_name'] = ''

    # Output information on the column range
    column_range = st.session_state['mg__curr_column_range']
    st.write('Column range: {}'.format(column_range))

    # Determine the x-y data to plot for the selected column, calculating the KDE for each column only once ever
    if column_to_plot not in list(st.session_state['mg__dfs_to_plot']):
        line2d = sns.kdeplot(df, x=column_to_plot).get_lines()[0]
        curr_df = pd.DataFrame({'Value': line2d.get_xdata(), 'Density': line2d.get_ydata()})
        st.session_state['mg__dfs_to_plot'][column_to_plot] = curr_df[(curr_df['Value'] >= column_range[0]) & (curr_df['Value'] <= column_range[1])]  # needed because the KDE can extend outside the possible value range
    df_to_plot_full = st.session_state['mg__dfs_to_plot'][column_to_plot]

    # Write some text boxes for the desired ranges for the current analysis column. Could make this stronger by enforcing the min_value and max_value parameters in each widget to correspond to the other so that the chosen max is never less than the chosen min, but initial attempts at this shows strange behavior and isn't worth debugging for now
    st.slider(label='Minimum value:', min_value=column_range[0], max_value=column_range[1], key='mg__selected_min_val')
    st.slider(label='Maximum value:', min_value=column_range[0], max_value=column_range[1], key='mg__selected_max_val')
    selected_min_val = st.session_state['mg__selected_min_val']
    selected_max_val = st.session_state['mg__selected_max_val']

    # Create a view of the full dataframe that is the selected subset
    df_to_plot_selected = df_to_plot_full[(df_to_plot_full['Value'] >= selected_min_val) & (df_to_plot_full['Value'] <= selected_max_val)]

    # Plot the Plotly figure in Streamlit
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=df_to_plot_full['Value'], y=df_to_plot_full['Density'], fill='tozeroy', mode='none', fillcolor='yellow', name='Full dataset', hovertemplate=' '))
    fig.add_trace(go.Scatter(x=df_to_plot_selected['Value'], y=df_to_plot_selected['Density'], fill='tozeroy', mode='none', fillcolor='red', name='Selection', hoverinfo='skip'))
    fig.update_layout(hovermode='x unified', xaxis_title='Column value', yaxis_title='Density')
    st.plotly_chart(fig)

    # If we want to add the current column filter to the current phenotype assignment...
    if st.button(':star2: Add column filter to current phenotype :star2:', use_container_width=True):
        st.session_state['mg__df_current_phenotype'] = pd.concat([st.session_state['mg__df_current_phenotype'], pd.DataFrame(pd.Series({'Column for filtering': column_to_plot, 'Minimum value': selected_min_val, 'Maximum value': selected_max_val})).T]).reset_index(drop=True)

    # # Optionally plot a cell scatter plot
    # st.selectbox(label='Image to plot:', options=unique_images_short, key='mg__image_to_plot')
    # image_to_plot = st.session_state['mg__image_to_plot']
    # if st.button('Update (or plot for the first time) the scatter plot of selected cells'):
    #     df_selected_image = df.loc[df['Slide ID'] == unique_image_dict[image_to_plot], ['Cell X Position', 'Cell Y Position', column_to_plot]].copy()
    #     df_selected_image['Label'] = 'Other'
    #     df_selected_image.loc[(df_selected_image[column_to_plot] >= selected_min_val) & (df_selected_image[column_to_plot] <= selected_max_val), 'Label'] = 'Selection'
    #     fig = px.scatter(data_frame=df_selected_image, x='Cell X Position', y='Cell Y Position', color='Label')
    #     fig.update_xaxes(scaleanchor='y')
    #     st.plotly_chart(fig)

# Current phenotype
with main_columns[1]:

    # Column header
    st.header(':two: Current phenotype')

    # Output the dataframe holding the phenotype that's currently being built
    st.dataframe(st.session_state['mg__df_current_phenotype'])

    # Choose a phenotype name
    st.text_input(label='Phenotype name:', key='mg__current_phenotype_name')

    # If we want to add the current phenotype to the new dataset...
    if st.button(label=':star2: Add phenotype to new dataset :star2:', use_container_width=True):
        pass

# New dataset
with main_columns[2]:

    # Column header
    st.header(':three: New dataset')

    st.write('bleh')
