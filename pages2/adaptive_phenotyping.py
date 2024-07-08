# Import relevant libraries
import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from pages2 import multiaxial_gating
import utils

# Global variable
st_key_prefix = 'adaptive_phenotyping__'


def update_phenotyping_sub_options(phenotyping_method_options):
    if st.session_state[st_key_prefix + 'phenotyping_method'] == phenotyping_method_options[0]:
        st.session_state[st_key_prefix + 'apply_thresh_to_selected_group'] = False
        st.session_state[st_key_prefix + 'average_over_all_groups'] = False
    else:
        st.session_state[st_key_prefix + 'apply_thresh_to_selected_group'] = True
        st.session_state[st_key_prefix + 'average_over_all_groups'] = True


def plotly_mean_and_sem(dfs, df_names):

    # Create a Plotly figure
    fig = go.Figure()

    # For each dataframe and its name...
    for df, df_name in zip(dfs, df_names):

        # Calculate mean and SEM over the groups
        mean_values = df.mean(axis='columns')
        sem_values = df.sem(axis='columns')

        # Add a scatter plot
        fig.add_trace(go.Scatter(
            x=mean_values.index,
            y=mean_values,
            error_y=dict(
                type='data', # or 'percent' for percentage-based error bars
                array=sem_values,
                visible=True
            ),
            mode='markers+lines', # Use 'markers' or 'lines' or 'markers+lines'
            name=df_name
        ))

    # Customize layout
    fig.update_layout(
        title=f'Average Positive Percentage Mean with SEM as Error Bars',
        xaxis_title='Z score',
        yaxis_title='Average positive percentage mean',
        showlegend=True
    )

    # Return the plotly figure
    return fig


def get_box_and_whisker_data(df_grouped, df_thresholds, apply_thresh_to_selected_group, df, channel_for_phenotyping, column_identifying_baseline_signal, value_identifying_baseline, value_identifying_signal, row_selection, return_figure_and_summary=True):

    # If apply_thresh_to_selected_group (and not average_over_all_groups), the input into the function generate_box_and_whisker() corresponds to the selected group. The mean and std used are the defaults (None) for the function, i.e., those corresponding to the baseline group of the selected group.
    # If not apply_thresh_to_selected_group (and not average_over_all_groups), the input into the function generate_box_and_whisker() corresponds to the entire dataset. The mean and std (which, when corresponding to a selection from df_thresholds, correspond to the baseline group) used are those from the selected group. --> this is the original phenotyping method (like a single T=0 threshold)!
    # The thresholds (as in df_thresholds) are calculated from just baseline groups, whereas the groups of dataframes correspond to both the baseline and signal groups.
    # If average_over_all_groups, then it's as if we're selecting the groups in turn via a for loop. We don't want to transform the entire dataset every time; we only want to transform each part of the full dataset once. Thus, it generally doesn't make sense to average_over_all_groups and (not apply_thresh_to_selected_group), but it does make sense to average_over_all_groups and apply_thresh_to_selected_group, which essentially thresholds each group using its own baseline group, i.e., corresponds to the new phenotyping method.
    # Thus, here are the settings for the old and new phenotyping methods:
    #   Old: apply_thresh_to_selected_group=False, average_over_all_groups=False, DO have a particular group selected --> "Selected threshold applied to entire dataset"
    #   New: apply_thresh_to_selected_group=True, average_over_all_groups=True, DO NOT have a particular group selected --> "Group-specific threshold applied to each group"
    # Note that there are no other combinations of these settings that will transform each part of the full dataset once, so the above is complete!

    # Obtain the index and dataframe of the group identified by row_selection
    if isinstance(df_grouped, list):
        current_index = None
        df_selected = df_grouped[0][1]
    else:
        current_index = df_thresholds.iloc[row_selection].name
        df_selected = df_grouped.get_group(current_index)

    # Obtain the dataframe to actually transform as well as the mean and std to use for the transformation
    if apply_thresh_to_selected_group:
        df_transform = df_selected
        mean_for_zscore_calc = None
        std_for_zscore_calc = None
    else:
        df_transform = df
        ser_selected = df_thresholds.iloc[row_selection]
        mean_for_zscore_calc = ser_selected.loc['z score = 0']
        std_for_zscore_calc = ser_selected.loc['z score = 1'] - mean_for_zscore_calc

    # Obtain the data that one would get by generating a box and whisker plot
    return_values = multiaxial_gating.generate_box_and_whisker(
        df=df_transform,
        column_for_filtering=channel_for_phenotyping,
        apply_another_filter=False,
        another_filter_column=None,
        values_on_which_to_filter=None,
        images_in_plotting_group_1=df_transform.loc[df_transform[column_identifying_baseline_signal] == value_identifying_baseline, 'Slide ID'].unique(),
        images_in_plotting_group_2=df_transform.loc[df_transform[column_identifying_baseline_signal] == value_identifying_signal, 'Slide ID'].unique(),
        all_cells=False,
        mean_for_zscore_calc=mean_for_zscore_calc,
        std_for_zscore_calc=std_for_zscore_calc,
        return_figure_and_summary=return_figure_and_summary
        )

    # Return the desired values
    if return_figure_and_summary:
        return return_values
    else:
        return return_values, current_index


def main():
    """
    Main function for the page.
    """

    # Ensure a dataset has been opened in the first place
    if 'input_dataset' not in st.session_state:
        st.warning('Please open a dataset from the Open File page at left.')
        return
    
    # Get some necessary variables from the session state
    df = st.session_state['input_dataset'].data

    # Store columns of certain types
    if st_key_prefix + 'categorical_columns' not in st.session_state:
        st.session_state[st_key_prefix + 'categorical_columns'] = utils.get_categorical_columns_including_numeric(df, max_num_unique_values=1000)
    categorical_columns = st.session_state[st_key_prefix + 'categorical_columns']

    # Initialize three columns
    columns = st.columns(3)

    # In the first column...
    with columns[0]:

        st.subheader(':one: Threshold calculation')

        # Select columns to use for grouping the threshold calculations
        key = st_key_prefix + 'columns_for_phenotype_grouping'
        if key not in st.session_state:
            st.session_state[key] = []
        columns_for_phenotype_grouping = st.multiselect('Columns for phenotype grouping:', categorical_columns, key=key)

        # Optionally force-update the list of categorical columns
        st.button('Update phenotype grouping columns 💡', help='If you don\'t see the column you want to group, click this button to update the list of potential phenotype grouping columns.', on_click=lambda: st.session_state.pop(st_key_prefix + 'categorical_columns', None))

        # Set the column name that describes the baseline field such as cell type
        key = st_key_prefix + 'column_identifying_baseline_signal'
        if key not in st.session_state:
            st.session_state[key] = categorical_columns[0]
        column_identifying_baseline_signal = st.selectbox('Column identifying baseline/signal:', categorical_columns, key=key, on_change=lambda: st.session_state.pop(st_key_prefix + 'value_identifying_baseline', None))

        # Extract the baseline field value (such as a specific cell type) to be used for determining the thresholds
        available_baseline_signal_values = df[column_identifying_baseline_signal].unique()
        key = st_key_prefix + 'value_identifying_baseline'
        if key not in st.session_state:
            st.session_state[key] = available_baseline_signal_values[0]
        value_identifying_baseline = st.selectbox('Value identifying baseline:', available_baseline_signal_values, key=key)

        # Extract the available channels for performing phenotyping
        key = st_key_prefix + 'channel_for_phenotyping'
        if key not in st.session_state:
            st.session_state[key] = df.columns[0]
        channel_for_phenotyping = st.selectbox('Channel for phenotyping:', df.columns, key=key)

        # If adaptive thresholding is desired...
        if st.button('Calculate thresholds for phenotyping'):

            # Group the relevant subset of the dataframe by the selected variables
            if len(columns_for_phenotype_grouping) > 0:
                df_grouped = df[columns_for_phenotype_grouping + [column_identifying_baseline_signal, channel_for_phenotyping, 'Slide ID']].groupby(by=columns_for_phenotype_grouping)
            else:
                df_grouped = [(None, df)]

            # Define various z scores of interest to use for determining the thresholds
            z_scores = np.arange(-1, 11)

            # For every group of variables for which you want a different threshold...
            thresholds_outer = []
            for _, df_group in df_grouped:

                # Get the selected intensity data for just the baseline field value
                ser_baseline = df_group.loc[df_group[column_identifying_baseline_signal] == value_identifying_baseline, channel_for_phenotyping]

                # From that, calculate the mean and std
                mean_to_use = ser_baseline.mean()
                std_to_use = ser_baseline.std()

                # Determine the corresponding threshold for each z score
                thresholds_inner = []
                for z_score in z_scores:
                    thresholds_inner.append(mean_to_use + z_score * std_to_use)

                # Add to the main thresholds holder
                thresholds_outer.append(thresholds_inner)

            # Determine the multi-index for the dataframe
            if len(columns_for_phenotype_grouping) == 0:
                index = pd.Index([-1])
            elif len(columns_for_phenotype_grouping) == 1:
                index = pd.Index(list(df_grouped.groups.keys()), name=columns_for_phenotype_grouping[0])
            elif len(columns_for_phenotype_grouping) > 1:
                index = pd.MultiIndex.from_tuples(list(df_grouped.groups.keys()), names=columns_for_phenotype_grouping)
            index.name = 'Grouping'

            # Determine the columns index for the dataframe
            columns_index = pd.Index([f'z score = {z_score}' for z_score in z_scores])
            columns_index.name = 'Thresholds'

            # Set the dataframe of thresholds
            df_thresholds = pd.DataFrame(thresholds_outer, columns=columns_index, index=index).sort_index()
            st.session_state[st_key_prefix + 'df_thresholds'] = df_thresholds

            # Save the grouped data as well for plotting and phenotyping
            st.session_state[st_key_prefix + 'df_grouped'] = df_grouped

        # Make sure the phenotyping thresholds have been calculated
        key = st_key_prefix + 'df_thresholds'
        if key not in st.session_state:
            st.warning('Please calculate thresholds for phenotyping')
            return

        # Set a shortcut to the relevant data
        df_thresholds = st.session_state[key]
        df_grouped = st.session_state[st_key_prefix + 'df_grouped']

        # Display the thresholds, allowing the user to select a single row
        st.write('Calculated thresholds:')
        group_selection = st.dataframe(df_thresholds, on_select='rerun', selection_mode='single-row')

        # Create a new figure
        fig = go.Figure()

        # Loop through each column in df_thresholds to add a trace for each one
        for column in df_thresholds.columns:
            fig.add_trace(go.Scatter(
                # x=df_thresholds.index,  # Use the DataFrame index for the x-axis
                x=np.array(range(len(df_thresholds))),
                y=df_thresholds[column],  # Column values for the y-axis
                mode='markers+lines',  # Line plot
                name=column  # Use the column name as the trace name
            ))

        # Update the layout to add titles and adjust other aesthetics if needed
        fig.update_layout(
            title='Line Plot of Thresholds Dataframe Above',
            xaxis_title='Index',
            yaxis_title='Phenotyping channel',
            legend_title='Column'
        )

        # Plot the line plots
        st.plotly_chart(fig)

    # In the second column...
    with columns[1]:

        st.subheader(':two: Percent positives plotting')

        # Extract the value identifying the signal (such as a specific cell type) to be used for testing the thresholds
        key = st_key_prefix + 'value_identifying_signal'
        if key not in st.session_state:
            st.session_state[key] = available_baseline_signal_values[0]
        value_identifying_signal = st.selectbox('Value identifying signal:', available_baseline_signal_values, key=key)

        # Whether to apply the threshold to just the selected group
        # In the future, this should be made a radio button probably, individual group vs. entire dataset
        key = st_key_prefix + 'apply_thresh_to_selected_group'
        if key not in st.session_state:
            st.session_state[key] = True
        apply_thresh_to_selected_group = st.checkbox('Apply threshold to each individual group (instead of to the entire dataset)', key=key)

        # Whether to average over all groups
        key = st_key_prefix + 'average_over_all_groups'
        if key not in st.session_state:
            st.session_state[key] = False
        average_over_all_groups = st.checkbox('Average over all groups', key=key)

        # If we want to generate a plot based on the currently selected group...
        if not average_over_all_groups:

            # Obtain the current selection
            row_selection_list = group_selection['selection']['rows']

            # If something is actually selected...
            if row_selection_list:

                # Get the box and whisker data
                fig, _ = get_box_and_whisker_data(df_grouped, df_thresholds, apply_thresh_to_selected_group, df, channel_for_phenotyping, column_identifying_baseline_signal, value_identifying_baseline, value_identifying_signal, row_selection_list[0], return_figure_and_summary=True)

                # Plot the box and whisker chart
                st.plotly_chart(fig)

        # If we want to generate a plot from averaging over all the groups...
        else:

            # Since this takes a non-trivial amount of time, hide the calculation behind a button
            if st.button('Calculate average positive percentages over all groups'):

                # Initialize the holders of the group indices and the average positive percentages
                index_holder = []
                ser_holder_baseline = []
                ser_holder_signal = []

                # For every group...
                for curr_row in range(len(df_thresholds)):

                    # Get the box and whisker data
                    df_summary, curr_index = get_box_and_whisker_data(df_grouped, df_thresholds, apply_thresh_to_selected_group, df, channel_for_phenotyping, column_identifying_baseline_signal, value_identifying_baseline, value_identifying_signal, curr_row, return_figure_and_summary=False)

                    # Append the desired data to the holders
                    df_summary = df_summary.drop('Threshold', axis='columns').set_index('Z score')  # the drop isn't necessary but it may make the operations marginally faster
                    index_holder.append(curr_index)
                    ser_holder_baseline.append(df_summary['Positive % (avg. over images) for baseline group'])
                    ser_holder_signal.append(df_summary['Positive % (avg. over images) for signal group'])

                # Calculate and save the average positive percentages for the baseline group
                df_baseline = pd.concat(ser_holder_baseline, axis='columns', keys=index_holder)
                df_baseline.columns.names = columns_for_phenotype_grouping if len(columns_for_phenotype_grouping) > 0 else ['No group']
                st.session_state[st_key_prefix + 'df_baseline'] = df_baseline

                # Calculate and save the average positive percentages for the signal group
                df_signal = pd.concat(ser_holder_signal, axis='columns', keys=index_holder)
                df_signal.columns.names = columns_for_phenotype_grouping if len(columns_for_phenotype_grouping) > 0 else ['No group']
                st.session_state[st_key_prefix + 'df_signal'] = df_signal

            # Ensure the average positive percentages have been calculated
            if st_key_prefix + 'df_baseline' not in st.session_state:
                st.warning('Please calculate average positive percentages over all groups')
                return

            # Get the difference in the average positive percentages between the baseline and signal groups
            df_baseline = st.session_state[st_key_prefix + 'df_baseline']
            df_signal = st.session_state[st_key_prefix + 'df_signal']
            df_diff = df_signal - df_baseline

            # Render some charts in Streamlit
            st.plotly_chart(plotly_mean_and_sem([df_baseline, df_signal], ['Baseline', 'Signal']))
            st.plotly_chart(plotly_mean_and_sem([df_diff], ['signal - baseline']))

    # In the third column...
    with columns[2]:

        st.subheader(':three: Phenotype generation')

        # Set the phenotype name
        key = st_key_prefix + 'phenotype_name'
        if key not in st.session_state:
            st.session_state[key] = ''
        phenotype_name = st.text_input('Phenotype name:', key=key)

        # Set the desired Z score
        key = st_key_prefix + 'desired_z_score'
        if key not in st.session_state:
            st.session_state[key] = 2.0
        desired_z_score = st.number_input('Desired Z score:', key=key)

        # Set the phenotyping method
        key = st_key_prefix + 'phenotyping_method'
        phenotyping_method_options = ["Selected threshold applied to entire dataset", "Group-specific threshold applied to each group"]
        if key not in st.session_state:
            st.session_state[key] = phenotyping_method_options[1]
        phenotyping_method = st.selectbox('Phenotyping method:', phenotyping_method_options, key=key, on_change=update_phenotyping_sub_options, args=(phenotyping_method_options,))

        # Obtain the current selection
        row_selection_list = group_selection['selection']['rows']

        # Determine whether the phenotyping button should be disabled
        phenotyping_button_disabled = False
        if (phenotyping_method == phenotyping_method_options[0]) and (not row_selection_list):
            st.warning(f'Phenotyping method "{phenotyping_method}" requires a group to be selected in the first column')
            phenotyping_button_disabled = True
        if (phenotyping_method == phenotyping_method_options[1]) and (row_selection_list):
            st.warning(f'Phenotyping method "{phenotyping_method}" does not actually require that no group be selected in the first column, but for clarity (since that method performs calculations for every group), we are enacting that requirement; please de-select the group selection in the first column')
            phenotyping_button_disabled = True

        # Render the button to perform the phenotyping
        if st.button('Perform phenotyping', disabled=phenotyping_button_disabled, on_click=update_phenotyping_sub_options, args=(phenotyping_method_options,)):

            # Perform the phenotyping method that applies the selected threshold to the entire dataset
            # Here: apply_thresh_to_selected_group=False, average_over_all_groups=False, DO have a particular group selected --> "Selected threshold applied to entire dataset"
            if phenotyping_method == phenotyping_method_options[0]:

                # Get the selected row of df_thresholds, which is a series
                ser_selected = df_thresholds.iloc[row_selection_list[0]]

                # Calculate the threshold from the mean/std from the selected group
                mean_for_zscore_calc = ser_selected.loc['z score = 0']
                std_for_zscore_calc = ser_selected.loc['z score = 1'] - mean_for_zscore_calc
                threshold = mean_for_zscore_calc + desired_z_score * std_for_zscore_calc

                # If the threshold for our desired Z score has already been calculated, enforce agreement
                if desired_z_score in np.arange(-1, 11, 1):
                    assert np.abs(threshold - ser_selected.loc[f'z score = {int(desired_z_score)}']) < 1e-8, 'The threshold calculated for the selected group does not match the threshold in the thresholds DataFrame'

                # Get the locations where the selected filtering column is at least the current threshold value and assign the integer version of this to a new series
                positive_loc = df[channel_for_phenotyping] >= threshold
                ser_phenotype = positive_loc.astype(int)

                # Output the acutally used threshold
                st.write(f'Threshold used for entire dataset: {threshold}')

            # Perform the phenotyping method that applies a group-specific threshold to each group
            # Here: apply_thresh_to_selected_group=True, average_over_all_groups=True, DO NOT have a particular group selected --> "Group-specific threshold applied to each group"
            elif phenotyping_method == phenotyping_method_options[1]:

                # Initialize the phenotype column to all-negative
                ser_phenotype = pd.Series(-1, index=df.index)

                # For every group...
                thresholds = []
                for curr_row in range(len(df_thresholds)):

                    # Obtain the index and dataframe of the group identified by curr_row
                    if isinstance(df_grouped, list):
                        curr_df = df_grouped[0][1]
                        curr_integer_indices_into_df = np.array(range(len(curr_df)))
                    else:
                        curr_index = df_thresholds.iloc[curr_row].name
                        curr_df = df_grouped.get_group(curr_index)
                        curr_integer_indices_into_df = df_grouped.indices[curr_index]

                    # Get the locations of the images in the baseline group
                    images_in_plotting_group_1 = curr_df.loc[curr_df[column_identifying_baseline_signal] == value_identifying_baseline, 'Slide ID'].unique()
                    image_loc_group_1 = curr_df['Slide ID'].isin(images_in_plotting_group_1)

                    # Get only the data for the column of interest for the baseline images
                    ser_for_z_score = curr_df.loc[image_loc_group_1, channel_for_phenotyping]

                    # Use the mean and std from those data to calculate the desired threshold
                    mean_for_zscore_calc = ser_for_z_score.mean()
                    std_for_zscore_calc = ser_for_z_score.std()
                    threshold = mean_for_zscore_calc + desired_z_score * std_for_zscore_calc

                    # If the threshold for our desired Z score has already been calculated, enforce agreement
                    if desired_z_score in np.arange(-1, 11, 1):
                        assert np.abs(threshold - df_thresholds.iloc[curr_row].loc[f'z score = {int(desired_z_score)}']) < 1e-8, 'The threshold calculated for the current group does not match the threshold in the thresholds DataFrame'

                    # Get the locations where the selected filtering column is at least the current threshold value
                    positive_loc = curr_df[channel_for_phenotyping] >= threshold  # boolean series fitting the current group

                    # Update the phenotype assignments for the current group
                    ser_phenotype.iloc[curr_integer_indices_into_df] = positive_loc.astype(int)  # the LHS is the slice of ser_phenotype corresponding to the current group

                    # Store the threshold for the current group
                    thresholds.append(threshold)

                # Check that there are no values of -1 remaining in ser_phenotype
                assert -1 not in ser_phenotype, 'There are still cells that have not been assigned positivity, which shouldn\'t have happened'

                # Output the actually used thresholds
                thresholds = pd.Series(thresholds, index=df_thresholds.index)
                thresholds.name = 'Threshold'
                st.write('Thresholds used for each group in turn:')
                st.write(thresholds)

            # Add the phenotype column to the dataframe
            pheno_colname = f'Phenotype {phenotype_name}'
            df[pheno_colname] = utils.downcast_series_dtype(ser_phenotype)
            st.success(f'Phenotype column "{pheno_colname}" has been appended to (or modified in) the dataset')
            st.write(f'Number of cells in each phenotype group (0 = negative, 1 = positive):')
            st.write(df[pheno_colname].value_counts().reset_index(drop=True))


# Run the main function
if __name__ == '__main__':

    # Call the main function
    main()
