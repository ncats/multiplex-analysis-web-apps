# Import relevant libraries
import streamlit as st
import pandas as pd
import scipy.stats
import numpy as np
import plotly.graph_objects as go
import plotly.express as px

# Global variable
st_key_prefix = 'radial_profiles_plotting__'


def calculate_significant_differences_between_groups(percent_positives, well_id_indices1, well_id_indices2, unique_time_vals, unique_outer_radii, confidence_level=0.95):

    # Initialize the flags array
    flags = np.ones((len(unique_time_vals), len(unique_outer_radii))) * np.nan

    # For every unique time and outer radius...
    for itime in range(len(unique_time_vals)):
        for iradius in range(len(unique_outer_radii)):

            # Get the percent positives for the two groups
            percent_positives_group1 = percent_positives[well_id_indices1, itime, iradius]
            percent_positives_group2 = percent_positives[well_id_indices2, itime, iradius]

            # Remove any NaNs
            percent_positives_group1 = percent_positives_group1[~np.isnan(percent_positives_group1)]
            percent_positives_group2 = percent_positives_group2[~np.isnan(percent_positives_group2)]

            # If there are fewer than two non-NaN values in either group, skip this iteration
            if (len(percent_positives_group1) < 2) or (len(percent_positives_group2) < 2):
                continue

            # Calculate the confidence intervals for the two groups
            res = scipy.stats.bootstrap((percent_positives_group1, percent_positives_group2), lambda arr1, arr2: arr2.mean(axis=0) - arr1.mean(axis=0), confidence_level=confidence_level)

            # Determine the flag
            if res.confidence_interval.low > 0:
                flags[itime, iradius] = 1
            elif res.confidence_interval.high < 0:
                flags[itime, iradius] = -1
            else:
                flags[itime, iradius] = 0

    # Create the heatmap
    heatmap = go.Heatmap(
        x=unique_outer_radii,
        y=unique_time_vals,
        z=flags,
        colorscale='Jet',
        zmin=-1,
        zmax=1
    )

    # Create a figure and add the heatmap
    fig = go.Figure(data=[heatmap])

    # Customize layout
    fig.update_layout(
        title='Flags Indicating Significant Differences (group2 - group1)',
        xaxis_title='Outer radius (um)',
        yaxis_title='Timepoint'
    )

    # Return the figure
    return fig


def get_heatmap(percent_positives, well_id_indices, unique_time_vals, unique_outer_radii):

    # Create the heatmap
    heatmap = go.Heatmap(
        x=unique_outer_radii,
        y=unique_time_vals,
        z=np.nanmean(percent_positives[well_id_indices, :, :], axis=0),
        colorscale='Jet',
        zmin=0,
        zmax=100
    )

    # Create a figure and add the heatmap
    fig = go.Figure(data=[heatmap])

    # Customize layout
    fig.update_layout(
        title='Average Percent Positive Over Selected Wells',
        xaxis_title='Outer radius (um)',
        yaxis_title='Timepoint'
    )

    # Return the figure
    return fig


def calculate_percent_positives(df, phenotype_column_for_analysis, unique_well_ids, unique_time_vals, unique_outer_radii):

    # Create a 3D array to hold the percent positives
    percent_positives = np.ones((len(unique_well_ids), len(unique_time_vals), len(unique_outer_radii))) * np.nan

    # For each well...
    for well_id, df_group in df.groupby('well_id'):

        # Get the location of the well_id in the unique_well_ids list
        well_id_loc = unique_well_ids.index(well_id)

        # For each unique combination of time and outer radius...
        for (time_val, outer_radius), df_group2 in df_group.groupby(by=['T', 'Outer radius']):

            # Get the location of the time_val and outer_radius in their respective lists
            time_val_loc = unique_time_vals.index(time_val)
            outer_radius_loc = unique_outer_radii.index(outer_radius)

            # Calculate the percent positive for the current group
            ser = df_group2[phenotype_column_for_analysis]
            percent_positives[well_id_loc, time_val_loc, outer_radius_loc] = (ser == 1).sum() / len(ser) * 100

    # Return the percent positives
    return percent_positives


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
        st.write('⚠️ Normal confidence intervals were calculated for some of the data instead of the selected bootstrap method.')
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

    # Optionally add an alpha value to each color and if so return the rbga values; otherwise return the hex values
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

    # Ensure the user has loaded a dataset
    if 'input_dataset' not in st.session_state:
        st.warning('Please open a dataset from the Open File page at left.')
        return

    # Save a shortcut to the dataframe
    df = st.session_state['input_dataset'].data

    # Obtain the phenotype columns
    phenotype_columns = [column for column in df.columns if column.startswith('Phenotype ')]

    # Ensure phenotype columns exist
    if len(phenotype_columns) == 0:
        st.warning('No phenotype columns found in the dataset. Please run the Adaptive Phenotyping page at left.')
        return

    # Make sure all of the columns in ['T', 'REEC', 'well_id', 'Outer radius'] are present in df
    if not all([col in df.columns for col in ['T', 'REEC', 'well_id', 'Outer radius']]):
        st.warning('The columns "T", "REEC", "well_id", and "Outer radius" must be present in the dataset.')
        return

    # Keep only the necessary columns
    df = df[['Slide ID', 'T', 'REEC', 'well_id', 'Outer radius'] + phenotype_columns]

    # Get the unique values of particular columns of interest
    unique_well_ids = sorted(df['well_id'].unique())
    unique_time_vals = sorted(df['T'].unique())
    unique_outer_radii = sorted(df['Outer radius'].unique())

    with st.columns(3)[0]:

        # Select a phenotype column on which to perform the analysis
        key = st_key_prefix + 'phenotype_column_for_analysis'
        if key not in st.session_state:
            st.session_state[key] = phenotype_columns[0]
        phenotype_column_for_analysis = st.selectbox('Select a phenotype column on which to perform the analysis:', phenotype_columns, key=key)

        # Calculate the percent positives
        if st.button('Calculate the percent positives'):

            # Save the percent_positives array to the session state
            st.session_state[st_key_prefix + 'percent_positives'] = calculate_percent_positives(df, phenotype_column_for_analysis, unique_well_ids, unique_time_vals, unique_outer_radii)

        # Ensure the percent positives have been calculated
        key = st_key_prefix + 'percent_positives'
        if key not in st.session_state:
            st.warning('Please calculate the percent positives.')
            return

        # Save a shortcut to the percent positives
        percent_positives = st.session_state[key]

        # Get the number of nans in percent_positives and if there are any, print out where they are
        if np.isnan(percent_positives).sum() > 0:
            st.write('⚠️ There are NaNs in the percent_positives array.')
            with st.expander('Expand to see where the NaNs are located (well ID, time, outer radius):', expanded=False):
                for i in range(len(unique_well_ids)):
                    for j in range(len(unique_time_vals)):
                        for k in range(len(unique_outer_radii)):
                            if np.isnan(percent_positives[i, j, k]):
                                st.write(unique_well_ids[i], unique_time_vals[j], unique_outer_radii[k])

        # Obtain the wells and their properties (just the REEC for now). This takes 0.10 to 0.15 seconds
        df_to_select = df[['well_id', 'REEC']].drop_duplicates().sort_values(['well_id', 'REEC'])

        # Allow the user to select more than one set of wells
        key = st_key_prefix + 'select_two_groups_of_wells'
        if key not in st.session_state:
            st.session_state[key] = False
        select_two_groups_of_wells = st.checkbox('Select two groups of wells', key=key)

        # Initialize columns appropriately
        if select_two_groups_of_wells:
            group1_column, group2_column = st.columns(2)
        else:
            group1_column = st.columns(1)[0]

        # Get the user selection from this dataframe
        with group1_column:
            if not select_two_groups_of_wells:
                st.write('Select the well(s) whose average percent positives we will plot:')
            else:
                st.write('Select the first group of wells:')
            selected_rows = st.dataframe(df_to_select, on_select='rerun', hide_index=True, key='group1_well_selection__do_not_persist')['selection']['rows']
            ser_selected_well_ids = df_to_select.iloc[selected_rows]['well_id']

        # If the user wants to select two groups of wells, allow them to select the second group
        if select_two_groups_of_wells:
            with group2_column:
                st.write('Select the second group of wells:')
                selected_rows = st.dataframe(df_to_select, on_select='rerun', hide_index=True, key='group2_well_selection__do_not_persist')['selection']['rows']
                ser_selected_well_ids2 = df_to_select.iloc[selected_rows]['well_id']

        # Ensure at least one well is selected
        if len(ser_selected_well_ids) == 0:
            if not select_two_groups_of_wells:
                st.warning('No wells selected. Please select them from the left side of the table above.')
            else:
                st.warning('No wells selected for the first group. Please select them from the left side of the left table above.')
            return

        # If the user wants to select two groups of wells, ensure at least one well is selected for the second group
        if select_two_groups_of_wells and len(ser_selected_well_ids2) == 0:
            st.warning('No wells selected for the second group. Please select them from the left side of the right table above.')
            return

        # Obtain the indices of the selected wells in unique_well_ids
        well_id_indices = [unique_well_ids.index(selected_well_id) for selected_well_id in ser_selected_well_ids]

        # If the user wants to select two groups of wells, obtain the indices of the selected wells in unique_well_ids
        if select_two_groups_of_wells:
            well_id_indices2 = [unique_well_ids.index(selected_well_id) for selected_well_id in ser_selected_well_ids2]

        # Checkbox for whether to plot confidence intervals
        key = st_key_prefix + 'plot_confidence_intervals'
        if key not in st.session_state:
            st.session_state[key] = True
        plot_confidence_intervals = st.checkbox('Plot confidence intervals (at least two wells must be selected)', key=key)

        # Set some widget defaults if they don't exist
        key = st_key_prefix + 'ci_type'
        if key not in st.session_state:
            st.session_state[key] = 'bootstrap'
        ci_type = st.session_state[key]
        key = st_key_prefix + 'alpha'
        if key not in st.session_state:
            st.session_state[key] = 0.2
        alpha = st.session_state[key]

        # Allow the user to customize their values if they matter
        if plot_confidence_intervals:

            # Radio button to select the type of confidence interval
            ci_type = st.radio('Select the type of confidence interval:', ['bootstrap', 'normal'], key=st_key_prefix + 'ci_type')

            # Number input for the alpha value
            alpha = st.number_input('Alpha value:', min_value=0.0, max_value=1.0, key=st_key_prefix + 'alpha')

    # Button to generate line plots
    if st.button('Generate line plots'):

        # For a given list of well_id indices, create a plotly lineplot on the percent_positives array using the radii on the x-axis and the percent positives on the y-axis
        st.plotly_chart(get_line_plots(percent_positives, well_id_indices, plot_confidence_intervals, unique_time_vals, 1, 'Time', unique_outer_radii, 'Outer radius (um)', alpha=alpha, ci_type=ci_type))

        # For a given list of well_id indices, create a plotly lineplot on the percent_positives array using the time values on the x-axis and the percent positives on the y-axis
        st.plotly_chart(get_line_plots(percent_positives, well_id_indices, plot_confidence_intervals, unique_outer_radii, 2, 'Bin', unique_time_vals, 'Time', alpha=alpha, ci_type=ci_type))

    # Button to generate heatmaps
    if st.button('Generate heatmaps'):

        # Display the heatmap
        st.plotly_chart(get_heatmap(percent_positives, well_id_indices, unique_time_vals, unique_outer_radii))

    # If the user wants to compare two groups of wells, allow them to set the confidence level
    if select_two_groups_of_wells:
        with st.columns(3)[0]:
            key = st_key_prefix + 'confidence_level'
            if key not in st.session_state:
                st.session_state[key] = 0.95
            confidence_level = st.number_input('Confidence level:', min_value=0.0, max_value=1.0, key=key, format='%.2f')

    # Button to compare differences between the two groups of wells
    if select_two_groups_of_wells and st.button('Assess whether the two group means are significantly different'):

        # Calculate the flags and display the heatmap
        st.plotly_chart(calculate_significant_differences_between_groups(percent_positives, well_id_indices, well_id_indices2, unique_time_vals, unique_outer_radii, confidence_level=confidence_level))


# Run the main function
if __name__ == '__main__':
    main()
