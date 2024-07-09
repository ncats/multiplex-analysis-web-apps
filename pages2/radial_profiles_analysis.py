# Import relevant libraries
import streamlit as st
import app_top_of_page as top
import streamlit_dataframe_editor as sde
import image_filter
import pandas as pd
import scipy.stats
import numpy as np
import plotly.graph_objects as go
import plotly.express as px

# Global variable
st_key_prefix = 'radial_profiles_plotting__'


def get_line_plots(percent_positives, well_id_indices, plot_confidence_intervals, unique_vals_for_series, position_in_percent_positives, series_name, unique_vals_for_x, xaxis_title, alpha=0.1, ci_type='bootstrap'):

    # Potentially permute the axes of percent_positives
    if position_in_percent_positives == 1:
        percent_positives_transposed = percent_positives
    else:
        percent_positives_transposed = np.transpose(percent_positives, axes=(0, 2, 1))

    # Get the default colors to use for both the main lines and the shaded confidence intervals
    colors_to_use = get_default_colors(len(unique_vals_for_series))
    colors_to_use_with_alpha = get_default_colors(len(unique_vals_for_series), alpha=alpha)

    # Initialize the plotly figure
    fig = go.Figure()

    # For each series...
    method_used_holder = []
    bootstrap_flag_holder = []
    for series_val_index, series_val in enumerate(unique_vals_for_series):

        # Get the percent positives for the current series
        curr_percent_positives = percent_positives_transposed[well_id_indices, series_val_index, :]

        # Calculate the confidence intervals for the current series
        confidence_intervals, bootstrap_flag, method_used = get_confidence_intervals(curr_percent_positives, ci_type=ci_type)
        method_used_holder.append(method_used)
        bootstrap_flag_holder.append(bootstrap_flag)

        # Optionally plot the confidence intervals
        if plot_confidence_intervals:

            # Plot the lower bound of the confidence interval
            fig.add_trace(go.Scatter(x=unique_vals_for_x, y=confidence_intervals[0, :], mode='lines', line=dict(width=0), showlegend=False))

            # Plot the upper bound of the confidence interval with filling to the next Y (the lower bound)
            fig.add_trace(go.Scatter(x=unique_vals_for_x, y=confidence_intervals[1, :], mode='lines', fill='tonexty', fillcolor=colors_to_use_with_alpha[series_val_index], line=dict(width=0), showlegend=True, name=f'{series_name} {series_val} CI'))  # color format example: 'rgba(0,100,80,0.2)'

    # For each series...
    for series_val_index, series_val in enumerate(unique_vals_for_series):

        # Get the percent positives for the current series
        curr_percent_positives = percent_positives_transposed[well_id_indices, series_val_index, :]

        # Calculate the means for the current series
        y = np.nanmean(curr_percent_positives, axis=0)

        # Plot the means
        fig.add_trace(go.Scatter(x=unique_vals_for_x, y=y, mode='lines+markers', name=f'{series_name} {series_val}', line=dict(color=colors_to_use[series_val_index]), marker=dict(color=colors_to_use[series_val_index])))

    # Update the layout of the figure
    fig.update_layout(title=f'Percent positive averaged over all selected wells', xaxis_title=xaxis_title, yaxis_title='Percent positive (%)')

    # Display a warning if we couldn't calculate potentially desired bootstrap confidence intervals for some of the data
    if np.any(bootstrap_flag_holder):
        st.warning('Normal confidence intervals were calculated for some of the data')
        with st.expander('Expand to see which data used normal confidence due to there being fewer than two non-NaN well_ids:', expanded=False):
            st.dataframe(pd.DataFrame(method_used_holder, index=unique_vals_for_series, columns=unique_vals_for_x))

    # Return the plotly figure
    return fig


def hex_to_rgb(hex_color):
    # Remove the '#' character and convert the remaining string to an integer using base 16
    # Then extract each color component
    hex_color = hex_color.lstrip('#')
    r, g, b = int(hex_color[0:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)
    return (r, g, b)


def get_default_colors(num_colors, alpha=None):

    # Get the color sequence
    color_sequence = px.colors.qualitative.Plotly

    # Return a list of the first 15 colors, cycling through color_sequence if necessary
    hex_colors = [color_sequence[i % len(color_sequence)] for i in range(num_colors)]

    if alpha is not None:
        return [f'rgba{hex_to_rgb(hex_color) + (alpha,)}' for hex_color in hex_colors]
    else:
        return hex_colors


def get_confidence_intervals(array2d, ci_type='bootstrap'):
    size_of_second_dim = array2d.shape[1]
    ci = np.ones((2, size_of_second_dim)) * np.nan
    method_used = []
    bootstrap_flag = False
    for i in range(size_of_second_dim):
        curr_set_of_data = array2d[:, i]
        curr_set_of_data = curr_set_of_data[~np.isnan(curr_set_of_data)]  # select out just the non-nan values in curr_set_of_data
        if (ci_type == 'bootstrap') and (len(curr_set_of_data) < 2):
            bootstrap_flag = True
            ci_type_to_use = 'normal'
        else:
            ci_type_to_use = ci_type
        method_used.append(ci_type_to_use)
        curr_ci = get_confidence_interval(pd.Series(curr_set_of_data), ci_type=ci_type_to_use)
        ci[:, i] = curr_ci
    return ci, bootstrap_flag, method_used


# Calculate the 95% confidence interval of a series
def get_confidence_interval(ser, ci_type='bootstrap'):
    assert ci_type in ['normal', 'bootstrap'], 'ci_type must be either "normal" or "bootstrap"'
    if ci_type == 'normal':
        # This is a common approach but works well primarily when the sample size is large (usually n > 30) and the data distribution is not heavily skewed.
        mean = ser.mean()
        margin_of_error = ser.sem() * 1.96
        return mean - margin_of_error, mean + margin_of_error
    elif ci_type == 'bootstrap':
        # Largely distribution-independent but sample sizes less than 10 should be interpreted with caution
        res = scipy.stats.bootstrap((ser,), np.mean, confidence_level=0.95)  # calculate the bootstrap confidence interval
        confidence_interval = res.confidence_interval  # extract the confidence interval
        return confidence_interval.low, confidence_interval.high


def main():
    """
    Main function for the page.
    """

    ####

    phenotype_columns = [column for column in st.session_state['input_dataset'].data.columns if column.startswith('Phenotype ')]

    df = st.session_state['input_dataset'].data[['Slide ID', 'T', 'REEC', 'well_id', 'Outer radius'] + phenotype_columns]

    st.write(df.sample(100))

    unique_well_ids = sorted(df['well_id'].unique())
    unique_time_vals = sorted(df['T'].unique())
    unique_outer_radii = sorted(df['Outer radius'].unique())

    st.write('well_id', df['well_id'].dtype, unique_well_ids)
    st.write('T', df['T'].dtype, unique_time_vals)
    st.write('Outer radius', df['Outer radius'].dtype, unique_outer_radii)

    percent_positives = np.ones((len(unique_well_ids), len(unique_time_vals), len(unique_outer_radii))) * np.nan

    tot_len = 0
    # ser_analysis_holder = []
    for well_id, df_group in df.groupby('well_id'):
        well_id_loc = unique_well_ids.index(well_id)
        st.write(well_id)
        for (time_val, outer_radius), df_group2 in df_group.groupby(by=['T', 'Outer radius']):
            time_val_loc = unique_time_vals.index(time_val)
            outer_radius_loc = unique_outer_radii.index(outer_radius)
            ser = df_group2['Phenotype SORE6']
            percent_positive = (ser == 1).sum() / len(ser) * 100
            # ser_analysis_holder.append({'well_id': well_id, 'T': time_val, 'Outer radius': outer_radius, 'Percent positive': percent_positive})
            percent_positives[well_id_loc, time_val_loc, outer_radius_loc] = percent_positive
            tot_len += len(df_group2)
    st.write(tot_len)
    # st.write(pd.DataFrame(ser_analysis_holder))

    # Get the number of nans in percent_positives
    st.write(np.isnan(percent_positives).sum())

    # Print out for what indices in percent_positives there are nans
    for i in range(len(unique_well_ids)):
        for j in range(len(unique_time_vals)):
            for k in range(len(unique_outer_radii)):
                if np.isnan(percent_positives[i, j, k]):
                    st.write(unique_well_ids[i], unique_time_vals[j], unique_outer_radii[k])

    well_id_indices = range(len(unique_well_ids))

    plot_confidence_intervals = True

    # For a given list of well_id indices, create a plotly lineplot on the percent_positives array using the radii on the x-axis and the percent positives on the y-axis
    st.plotly_chart(get_line_plots(percent_positives, well_id_indices, plot_confidence_intervals, unique_time_vals, 1, 'Time', unique_outer_radii, 'Outer radius (um)', alpha=0.1, ci_type='bootstrap'))

    # For a given list of well_id indices, create a plotly lineplot on the percent_positives array using the time values on the x-axis and the percent positives on the y-axis
    st.plotly_chart(get_line_plots(percent_positives, well_id_indices, plot_confidence_intervals, unique_outer_radii, 2, 'Outer radius', unique_time_vals, 'Time', alpha=0.1, ci_type='bootstrap'))

    ####

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
    df_confidence_interval = grouped_by_T.agg(lambda ser: get_confidence_interval(ser, ci_type='bootstrap'))

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
    page_name = 'Radial Profiles - Plotting Aggregated Results'
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
