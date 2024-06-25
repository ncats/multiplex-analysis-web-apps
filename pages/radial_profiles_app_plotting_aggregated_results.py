# Import relevant libraries
import streamlit as st
import app_top_of_page as top
import streamlit_dataframe_editor as sde
import image_filter
import pandas as pd
import scipy.stats
import numpy as np
import plotly.graph_objects as go


# Global variable
st_key_prefix = 'radial_profiles_plotting__'


# Calculate the 95% confidence interval of a series
def get_confidence_interval(ser, type='bootstrap'):
    if type == 'normal':
        # This is a common approach but works well primarily when the sample size is large (usually n > 30) and the data distribution is not heavily skewed.
        mean = ser.mean()
        margin_of_error = ser.sem() * 1.96
        return mean - margin_of_error, mean + margin_of_error
    elif type == 'bootstrap':
        # Largely distribution-independent but sample sizes less than 10 should be interpreted with caution
        res = scipy.stats.bootstrap((ser,), np.mean, confidence_level=0.95)  # calculate the bootstrap confidence interval
        confidence_interval = res.confidence_interval  # extract the confidence interval
        return confidence_interval.low, confidence_interval.high
    else:
        return ser.mean()


def main():
    """
    Main function for the page.
    """

    # Ensure analysis has been run
    key = 'radial_profiles__' + 'df_analysis_results'
    if key not in st.session_state:
        st.warning('You must run the analysis on the previous page before plotting the results here.')
        return

    # Save a shortcut to the analysis results dataframe
    df_analysis_results = st.session_state[key]

    # Save a shortcut to the main dataframe
    df = st.session_state['input_dataset'].data

    # Get the columns with a single unique value in each image
    key = st_key_prefix + 'columns_with_single_unique_value'
    if key not in st.session_state:
        grouped_nunique = df.groupby('Slide ID').agg(lambda x: x.nunique())  # group by 'Slide ID' and calculate nunique for each column within each group
        max_unique_values = grouped_nunique.max()  # find the maximum number of unique values per column across all groups
        columns_with_single_unique_value = max_unique_values[max_unique_values == 1].index  # create a mask for columns where the maximum number of unique values is 1
        st.session_state[key] = columns_with_single_unique_value
    columns_with_single_unique_value = st.session_state[key]

    # Instantiate the image selector
    image_selector = image_filter.ImageFilter(df, image_colname='Slide ID', st_key_prefix=st_key_prefix, possible_filtering_columns=columns_with_single_unique_value)

    # If the image filter is not ready (which means the filtering dataframe was not generated), return
    if not image_selector.ready:
        st.warning('Please prepare the filtering data first')
        return

    # Create the image filter, returning the selected images and the dataframe on which the selection was performed, saving shortcuts to these values
    st.session_state[st_key_prefix + 'images_in_plotting_group_1'], st.session_state[st_key_prefix + 'df_masked_group_1'] = image_selector.select_images(key='group 1', color='blue', return_df_masked=True)
    images_in_plotting_group_1 = st.session_state[st_key_prefix + 'images_in_plotting_group_1']
    df_masked_group_1 = st.session_state[st_key_prefix + 'df_masked_group_1']

    # If no images were selected, return
    if len(images_in_plotting_group_1) == 0:
        st.warning('No images selected in group 1')
        return
    
    # Sort the selected images by their timepoint (probably unnecessary), saving the resulting images and times
    df_masked_group_1_selected_sorted = df_masked_group_1.loc[images_in_plotting_group_1].sort_values('T')
    image_names_in_time_order = df_masked_group_1_selected_sorted.index
    times = df_masked_group_1_selected_sorted['T']

    # Extract the analysis results (percentages) for the selected images in order
    df_analysis_results_selected = df_analysis_results.loc[image_names_in_time_order]

    # Get the smallest number of annuli across all selected images
    minimum_number_of_annuli = df_analysis_results_selected.iloc[:, -1].min()

    # Combine the analysis dataframe with the corresponding times
    df_combined = pd.concat([df_analysis_results_selected.iloc[:, :minimum_number_of_annuli], times], axis='columns')

    # Group by the time column
    grouped_by_T = df_combined.groupby('T')

    # Calculate the mean of each group
    df_means = grouped_by_T.mean()

    # Extract just the outer radii from the column names and set those to the new column names
    outer_radii = [float(column.split(' ')[4]) for column in df_means.columns]
    df_means.columns = outer_radii

    # Sort by index and columns
    df_means = df_means.sort_index()
    df_means = df_means[sorted(df_means.columns)]

    # Name the index and column accordingly
    df_means.index.name = 'Time T'
    df_means.columns.name = 'Outer radius (um)'

    # Write the means
    st.write(df_means)

    # Calculate the 95% confidence interval for each group
    df_confidence_interval = grouped_by_T.agg(lambda ser: get_confidence_interval(ser, type='bootstrap'))

    # Set the columns to just the outer radii
    df_confidence_interval.columns = outer_radii

    # Sort by index and columns
    df_confidence_interval = df_confidence_interval.sort_index()
    df_confidence_interval = df_confidence_interval[sorted(df_confidence_interval.columns)]

    # Name the index and column accordingly
    df_confidence_interval.index.name = 'Time T'
    df_confidence_interval.columns.name = 'Outer radius (um)'

    # Write the confidence intervals
    st.write(df_confidence_interval)

    # Create a heatmap of df_means
    heatmap = go.Heatmap(
        x=df_means.columns,  # Column names as x-axis labels
        y=df_means.index,    # Index as y-axis labels
        z=df_means.values,    # DataFrame values as heatmap values
        colorscale='Jet'
    )

    # Create a figure and add the heatmap
    fig = go.Figure(data=[heatmap])

    # Customize layout
    fig.update_layout(
        title='Heatmap of df_means',
        xaxis_title='X Axis Title',
        yaxis_title='Y Axis Title'
    )

    st.plotly_chart(fig)


# Run the main function
if __name__ == '__main__':

    # Set page settings
    page_name = 'Radial Profiles Plots (Aggregated Results)'
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
