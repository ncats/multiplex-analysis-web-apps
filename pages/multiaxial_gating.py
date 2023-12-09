# Import relevant libraries
import pandas as pd
import seaborn as sns
import streamlit as st
from st_pages import show_pages_from_config, add_indentation
import plotly.graph_objects as go
import os
import dataset_formats
import plotly.express as px
import streamlit_dataframe_editor as sde

# Function to load the data in a unified format
def load_data(input_datafile_path, coord_units_in_microns, dataset_format):
    dataset_class = getattr(dataset_formats, dataset_format)  # done this way so that the format (e.g., “REEC”) can be select programmatically
    dataset_obj = dataset_class(input_datafile=input_datafile_path, coord_units_in_microns=coord_units_in_microns)
    dataset_obj.process_dataset(do_calculate_minimum_coordinate_spacing_per_roi=False, do_trimming=False)
    return dataset_obj.data

# Update the dependencies of the selectbox for the current analysis column
def update_dependencies_of_column_for_filtering():
    df = st.session_state['mg__df']
    column_for_filtering = st.session_state['mg__column_for_filtering']
    st.session_state['mg__curr_column_range'] = [df[column_for_filtering].min(), df[column_for_filtering].max()]
    st.session_state['mg__selected_min_val'] = st.session_state['mg__curr_column_range'][0]
    st.session_state['mg__selected_max_val'] = st.session_state['mg__curr_column_range'][1]

# Set the column options to all numeric columns less columns that have already been used for filtering for the current phenotype
def update_column_options():
    return [column for column in st.session_state['mg__all_numeric_columns'] if column in (set(st.session_state['mg__all_numeric_columns']) - set(st.session_state['mg__de_current_phenotype'].reconstruct_edited_dataframe()['Column for filtering']))]

# Add to the current phenotype and update the previous parts of the app
def update_dependencies_of_button_for_adding_column_filter_to_current_phenotype(column_for_filtering, selected_min_val, selected_max_val):

    # Get the current value of the current phenotype dataframe editor
    df_current_phenotype = st.session_state['mg__de_current_phenotype'].reconstruct_edited_dataframe()

    # Add the selected column filter to the current phenotype assignments dataframe and update the phenotype assignments dataframe with this new dataframe
    new_df_contents = pd.concat([df_current_phenotype, pd.DataFrame(pd.Series({'Column for filtering': column_for_filtering, 'Minimum value': selected_min_val, 'Maximum value': selected_max_val})).T]).reset_index(drop=True)
    st.session_state['mg__de_current_phenotype'].update_editor_contents(new_df_contents=new_df_contents)

    # Set the currently selected column as the first of the possible options
    st.session_state['mg__column_for_filtering'] = update_column_options()[0]

    # Since the column selection must have just changed, update its dependencies
    update_dependencies_of_column_for_filtering()

# Add to the phenotype assignments for the new dataset and update the previous parts of the app
def update_dependencies_of_button_for_adding_phenotype_to_new_dataset():

    # Get the current values of the two dataframe editors
    df_current_phenotype = st.session_state['mg__de_current_phenotype'].reconstruct_edited_dataframe()
    df_phenotype_assignments = st.session_state['mg__de_phenotype_assignments'].reconstruct_edited_dataframe()
    if 'Phenotype' in df_phenotype_assignments.columns:
        df_phenotype_assignments = df_phenotype_assignments.set_index('Phenotype')

    # Populate a dictionary of the column filteres for the current phenotype
    curr_phenotype_dict = dict()
    for row in df_current_phenotype.itertuples(index=False):
        curr_col, curr_min, curr_max = row
        curr_phenotype_dict[curr_col + ' [[min]]'] = curr_min
        curr_phenotype_dict[curr_col + ' [[max]]'] = curr_max

    # Update the contents of the phenotype assignments data editor
    new_df_contents = pd.concat([df_phenotype_assignments, pd.DataFrame(pd.Series(curr_phenotype_dict, name=st.session_state['mg__current_phenotype_name'])).T]).rename_axis('Phenotype').reset_index(drop=False)
    st.session_state['mg__de_phenotype_assignments'].update_editor_contents(new_df_contents=new_df_contents)

    # Clear the current phenotype dataframe editor to its default value
    st.session_state['mg__de_current_phenotype'].reset_dataframe_content()

    # Set the currently selected column as the first of the possible options
    st.session_state['mg__column_for_filtering'] = update_column_options()[0]

    # Since the column selection must have just changed, update its dependencies
    update_dependencies_of_column_for_filtering()

# From the phenotype assignments, add one column per phenotype to the original dataframe containing pluses where all the phenotype criteria are met
def add_new_phenotypes_to_main_df(df):

    # Get the current values of the phenotype assignments data editor
    df_phenotype_assignments = st.session_state['mg__de_phenotype_assignments'].reconstruct_edited_dataframe().set_index('Phenotype')

    # For each set of phenotype assignments...
    for row in df_phenotype_assignments.itertuples():

        # Obtain the name of the current phenotype
        phenotype = row[0]

        # Create a dataframe from the current row so that we can split and add columns
        curr_df = pd.Series(dict(zip(df_phenotype_assignments.columns, row[1:])), name=phenotype).dropna().to_frame().reset_index()

        # Add a "column" column containing the name of the filtering column
        curr_df['column'] = [x.split(' [[')[0] for x in curr_df['index']]

        # Drop the unnecessary "index" column
        curr_df = curr_df.drop(['index'], axis='columns')

        # Initialize a Series of booleans to True
        phenotype_bools = pd.Series([True] * len(df))

        # For each filtering column...
        for column, value_range in curr_df.groupby(by='column')[phenotype].apply(lambda x: list(x)).items():  # note the groupby appears to result in the column's order of min then max

            # Determine where the current column values are within the specified range criterion
            phenotype_bools = phenotype_bools & (df[column] >= value_range[0]) & (df[column] <= value_range[1])

        # Add a column to the original dataframe with the new phenotype satisfying all of its filtering criteria
        st.session_state['mg__df']['Phenotype {}'.format(phenotype)] = phenotype_bools.apply(lambda x: ('+' if x else '-'))

def main():

    # Set page settings
    st.set_page_config(layout='wide', page_title='Multiaxial Gating')
    st.title('Multiaxial Gating')

    # Apply pages order and indentation
    add_indentation()
    show_pages_from_config()

    # Run streamlit-dataframe-editor library initialization tasks at the top of the page
    st.session_state = sde.initialize_session_state(st.session_state)

    # Set the default dataframes to be edited
    default_df_current_phenotype = pd.DataFrame(columns=['Column for filtering', 'Minimum value', 'Maximum value'])
    default_df_phenotype_assignments = pd.DataFrame()

    # Initialize some things in the session state
    if 'mg__dfs_to_plot' not in st.session_state:
        st.session_state['mg__dfs_to_plot'] = dict()
    if 'mg__current_phenotype_name' not in st.session_state:
        st.session_state['mg__current_phenotype_name'] = ''
    if 'mg__de_current_phenotype' not in st.session_state:
        st.session_state['mg__de_current_phenotype'] = sde.DataframeEditor(df_name='mg__df_current_phenotype', default_df_contents=default_df_current_phenotype)
    if 'mg__de_phenotype_assignments' not in st.session_state:
        st.session_state['mg__de_phenotype_assignments'] = sde.DataframeEditor(df_name='mg__df_phenotype_assignments', default_df_contents=default_df_phenotype_assignments)

    # Constant
    input_directory = os.path.join('.', 'input')

    # Set datafile information --> turn into widgets soon
    # input_datafilename = 'measurementsEpCAMLy51MHCII-exported.csv'
    input_datafilename = 'measurementsthymus-exported.csv'
    coord_units_in_microns = 1

    # Load the data
    if st.button('Load data'):
        st.session_state['mg__df'] = load_data(os.path.join(input_directory, input_datafilename), coord_units_in_microns, dataset_formats.extract_datafile_metadata(os.path.join(input_directory, input_datafilename))[4])
        unique_images = st.session_state['mg__df']['Slide ID'].unique()
        st.session_state['mg__unique_images_short'] = [x.split('-imagenum_')[1] for x in unique_images]
        st.session_state['mg__unique_image_dict'] = dict(zip(st.session_state['mg__unique_images_short'], unique_images))
        st.session_state['mg__all_numeric_columns'] = st.session_state['mg__df'].select_dtypes(include='number').columns
        phenotype_columns = [column for column in st.session_state['mg__df'].columns if column.startswith('Phenotype ')]
        st.session_state['mg__df'] = st.session_state['mg__df'].rename(columns=dict(zip(phenotype_columns, [column.replace('Phenotype ', 'Phenotype_orig ') for column in phenotype_columns])))
    
    # Warn the user that they need to load the data at least once
    if 'mg__df' not in st.session_state:
        st.warning('You must click the "Load data" button above before you can use this app!')

    # If the data have been loaded...
    else:
    
        # Load the data and some resulting processed data
        df = st.session_state['mg__df']
        unique_images_short = st.session_state['mg__unique_images_short']
        unique_image_dict = st.session_state['mg__unique_image_dict']

        # Define the main columns
        main_columns = st.columns(3, gap='medium')

        # Data column filter
        with main_columns[0]:

            # Column header
            st.header(':one: Column filter')

            # Have a dropdown for the column on which to perform a kernel density estimate
            st.selectbox(label='Column for filtering:', options=update_column_options(), key='mg__column_for_filtering', on_change=update_dependencies_of_column_for_filtering)
            column_for_filtering = st.session_state['mg__column_for_filtering']

            # Output information on the column range
            if 'mg__curr_column_range' not in st.session_state:
                update_dependencies_of_column_for_filtering()
            column_range = st.session_state['mg__curr_column_range']
            st.write('Column range: {}'.format(column_range))

            # Determine the x-y data to plot for the selected column, calculating the KDE for each column only once ever
            if column_for_filtering not in list(st.session_state['mg__dfs_to_plot']):
                line2d = sns.kdeplot(df, x=column_for_filtering).get_lines()[0]
                curr_df = pd.DataFrame({'Value': line2d.get_xdata(), 'Density': line2d.get_ydata()})
                st.session_state['mg__dfs_to_plot'][column_for_filtering] = curr_df[(curr_df['Value'] >= column_range[0]) & (curr_df['Value'] <= column_range[1])]  # needed because the KDE can extend outside the possible value range
            df_to_plot_full = st.session_state['mg__dfs_to_plot'][column_for_filtering]

            # Write some text boxes for the desired ranges for the current analysis column. Could make this stronger by enforcing the min_value and max_value parameters in each widget to correspond to the other so that the chosen max is never less than the chosen min, but initial attempts at this shows strange behavior and isn't worth debugging for now
            # Implement range slider here!
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

            # Add the current column filter to the current phenotype assignment
            st.button(':star2: Add column filter to current phenotype :star2:', use_container_width=True, on_click=update_dependencies_of_button_for_adding_column_filter_to_current_phenotype, args=(column_for_filtering, selected_min_val, selected_max_val))

        # Current phenotype and phenotype assignments
        with main_columns[1]:

            # Column header
            st.header(':two: Current phenotype', help='Note you can refine values in the following table by editing them directly or even deleting (or adding) whole rows.')

            # Output the dataframe holding the phenotype that's currently being built
            st.session_state['mg__de_current_phenotype'].dataframe_editor()

            # Choose a phenotype name
            st.text_input(label='Phenotype name:', key='mg__current_phenotype_name')

            # Add the current phenotype to the phenotype assignments table
            st.button(label=':star2: Add phenotype to assignments table :star2:', use_container_width=True, on_click=update_dependencies_of_button_for_adding_phenotype_to_new_dataset)

            # Column header
            st.header(':three: Phenotype assignments', help='Note you can refine values in the following table by editing them directly or even deleting whole rows.')

            # Output the dataframe holding the specifications for all phenotypes
            st.session_state['mg__de_phenotype_assignments'].dataframe_editor()

            # Generate the new dataset
            st.button(label=':star2: Generate the new dataset from the phenotype assignments :star2:', use_container_width=True, on_click=add_new_phenotypes_to_main_df, args=(df,))

        # New dataset
        with main_columns[2]:

            # Column header
            st.header(':four: New dataset')

            # Print out a sample of the main dataframe
            st.write('Augmented dataset sample:')
            st.dataframe(st.session_state['mg__df'].head(5))

            # Get a list of all new phenotypes
            new_phenotypes = [column for column in st.session_state['mg__df'].columns if column.startswith('Phenotype ')]

            # If at least one phenotype has been assigned...
            if len(new_phenotypes) > 0:

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

                # Optionally run some checks
                if st.toggle(label='Do checks'):

                    # Write out the detected cutoffs to use for basic validation
                    phenotypes_orig = ['Phenotype_orig MHCII', 'Phenotype_orig Ly51', 'Phenotype_orig EpCAM']
                    intensities = ['MHC II (CH3 Fluor Cy3): Membrane: Mean', 'Ly51 (CH4 Fluor Cy5): Membrane: Mean', 'EpCAM (CH2 Fluor GFP): Membrane: Mean']
                    for ipheno in range(len(phenotypes_orig)):
                        tmp = st.session_state['mg__df'][phenotypes_orig[ipheno]]
                        tmp.index = st.session_state['mg__df'][intensities[ipheno]]
                        tmp = tmp.sort_index()
                        st.write('Intensity cutoff for {}: {}'.format(phenotypes_orig[ipheno], tmp[tmp == '+'].index[0]))

                    # Write the numbers of positive markers to compare the new phenotypes with the original phenotypes
                    for phenotype_orig in phenotypes_orig:
                        st.write(st.session_state['mg__df'][phenotype_orig].value_counts())
                        st.write(st.session_state['mg__df'][phenotype_orig.replace('_orig', '').replace('MHCII', 'MHC II') + ' new'].value_counts())

        # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
        st.session_state = sde.finalize_session_state(st.session_state)

if __name__ == '__main__':
    main()
