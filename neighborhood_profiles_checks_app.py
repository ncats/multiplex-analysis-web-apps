# Import relevant libraries
import neighborhood_profiles_checks
import pandas as pd
import os
import numpy as np
import streamlit as st
import plotly.express as px
import plotly.graph_objects as go


def plotly_scatter_plot(df_to_plot, x_colname='x coord', y_colname='y coord', label_colname='my_label', unique_labels=[0, 1], plot_title='My Plot', opacity_colname=None):

    # Define a color scale
    color_scale = px.colors.qualitative.Plotly

    # Map the 'cluster_label' values to colors
    color_mapping = {label: color_scale[i % len(color_scale)] for i, label in enumerate(unique_labels)}

    # Create a new column in the DataFrame that contains the colors of the markers
    df_to_plot.loc[:, 'color'] = df_to_plot[label_colname].map(color_mapping)

    # Create the scatter plot
    fig = go.Figure()

    # For each label, add a trace to the figure
    for label in unique_labels:
        df_subset = df_to_plot[df_to_plot[label_colname] == label]
        if opacity_colname is not None:
            fig.add_trace(go.Scatter(x=df_subset[x_colname], y=df_subset[y_colname], mode='markers', marker=dict(color=df_subset['color'], opacity=df_subset[opacity_colname]), name=label))
        else:
            fig.add_trace(go.Scatter(x=df_subset[x_colname], y=df_subset[y_colname], mode='markers', marker=dict(color=df_subset['color']), name=label))

    # Set the aspect ratio and plot title
    fig.update_layout(xaxis={"scaleanchor": "y", "scaleratio": 1}, title=plot_title)

    # Return the figure
    return fig


def draw_plots(df_image_labels, umap_x_colname='UMAP_1_20230327_152849', umap_y_colname='UMAP_2_20230327_152849', property_colnames=['XMin', 'XMax', 'YMin', 'YMax'], image_colname='ShortName', binary_colname='Survival_5yr', spatial_x_colname='CentroidX', spatial_y_colname='CentroidY'):
    
    # Get shortcuts to variables in the session state
    df_by_bin = st.session_state['df_by_bin']
    df_by_cell = st.session_state['df_by_cell']
    unique_images = st.session_state['unique_images']

    # Write a header
    st.header('Plots by dataset')

    # For plotting purposes below (to get a discrete legend instead of a colorbar), we need to convert the 'cluster_label' column to categorical, and also save their unique values so we can preserve plotting order
    df_by_bin['cluster_label'] = df_by_bin['cluster_label'].astype('category')
    df_by_cell['cluster_label'] = df_by_cell['cluster_label'].astype('category')
    unique_cluster_labels_by_bin = set([cluster_label for cluster_label in df_by_bin['cluster_label'].unique() if not np.isnan(cluster_label)])
    unique_cluster_labels_by_cell = set([cluster_label for cluster_label in df_by_cell['cluster_label'].unique() if not np.isnan(cluster_label)])
    assert unique_cluster_labels_by_bin == unique_cluster_labels_by_cell, f'The cluster labels are not the same for the bins ({unique_cluster_labels_by_bin}) and cells ({unique_cluster_labels_by_cell}).'
    unique_cluster_labels = list(unique_cluster_labels_by_bin)

    # Create two columns for the per-dataset plots
    col1, col2 = st.columns(2)

    # Plots by bin
    with col1:
        col1.plotly_chart(plotly_scatter_plot(df_by_bin, x_colname=umap_x_colname, y_colname=umap_y_colname, label_colname='cluster_label', unique_labels=unique_cluster_labels, plot_title='UMAP by Bin for Whole Dataset'))
        col1.plotly_chart(px.line(df_by_bin.groupby('cluster_label', observed=True)[property_colnames].mean().reset_index().melt(id_vars='cluster_label', var_name='column', value_name='value'), x='column', y='value', color='cluster_label', markers=True, title='Property Means by Bin for Whole Dataset', category_orders={'cluster_label': unique_cluster_labels}))  # get the neighbor vectors for each cluster averaged over the histogram bins falling in that cluster

    # Plots by cell
    with col2:
        col2.plotly_chart(plotly_scatter_plot(df_by_cell, x_colname=umap_x_colname, y_colname=umap_y_colname, label_colname='cluster_label', unique_labels=unique_cluster_labels, plot_title='UMAP by Cell for Whole Dataset'))
        col2.plotly_chart(px.line(df_by_cell.groupby('cluster_label', observed=True)[property_colnames].mean().reset_index().melt(id_vars='cluster_label', var_name='column', value_name='value'), x='column', y='value', color='cluster_label', markers=True, title='Property Means by Cell for Whole Dataset', category_orders={'cluster_label': unique_cluster_labels}))  # get the neighbor vectors for each cluster averaged over the cells falling in that cluster

    # Write a header
    st.header('Plots by image')

    # Use one column for the image selection
    col1, _ = st.columns(2)
    with col1:

        # Allow the user to select an image by dropdown
        if 'image_to_plot' not in st.session_state:
            if unique_images:
                st.session_state['image_to_plot'] = unique_images[0]
            else:
                st.session_state['image_to_plot'] = None
        st.selectbox('Image to plot:', unique_images, key='image_to_plot')

        # Add previous and next buttons to select the image by button
        current_index = unique_images.index(st.session_state['image_to_plot'])
        def update_image_to_plot(new_index):
            st.session_state['image_to_plot'] = unique_images[new_index]
        button_col1, button_col2 = st.columns(2)
        with button_col1:
            st.button('Previous', on_click=update_image_to_plot, kwargs={'new_index': current_index - 1}, disabled=(current_index == 0), use_container_width=True)
        with button_col2:
            st.button('Next', on_click=update_image_to_plot, kwargs={'new_index': current_index + 1}, disabled=(current_index == len(unique_images) - 1), use_container_width=True)

        # Get a shortcut to the selected image
        selected_image = st.session_state['image_to_plot']

        # Get the current image binary label
        current_image_label = df_image_labels.loc[selected_image, binary_colname]

    # Create two columns for the per-image plots
    col1, col2 = st.columns(2)

    # Plots by bin
    with col1:

        # Filter the binned data to the selected image
        image_in_set = pd.Series([selected_image in images for images in df_by_bin['unique_images']], index=df_by_bin.index)
        num_bins_with_cluster_labels = df_by_bin['cluster_label'].loc[image_in_set].notnull().sum()
        df_by_bin_filtered = df_by_bin[image_in_set]
        df_by_bin_filtered.loc[:, 'cluster_label'] = df_by_bin_filtered['cluster_label'].cat.remove_unused_categories()

        # Write some image information
        st.write(f'In the selected image, there are {image_in_set.sum()} bins present, {num_bins_with_cluster_labels} of which have been assigned a cluster label.')
        st.write(f'Label for image {selected_image}: {current_image_label}.')

        # Draw the plots
        st.plotly_chart(plotly_scatter_plot(df_by_bin_filtered, x_colname=spatial_x_colname, y_colname=spatial_y_colname, label_colname='cluster_label', unique_labels=unique_cluster_labels, plot_title=f'Spatial by Bin for {selected_image} (this is meaningless; don\'t read into it)'))
        df_by_bin['opacity'] = 0.05
        df_by_bin.loc[image_in_set, 'opacity'] = 1
        st.plotly_chart(plotly_scatter_plot(df_by_bin, x_colname=umap_x_colname, y_colname=umap_y_colname, label_colname='cluster_label', unique_labels=unique_cluster_labels, plot_title=f'UMAP by Bin for {selected_image}', opacity_colname='opacity'))
        st.plotly_chart(px.line(df_by_bin_filtered.groupby('cluster_label', observed=True)[property_colnames].mean().reset_index().melt(id_vars='cluster_label', var_name='column', value_name='value'), x='column', y='value', color='cluster_label', markers=True, title=f'Property Means by Bin for {selected_image}', category_orders={'cluster_label': unique_cluster_labels}))  # get the neighbor vectors for each cluster averaged over the histogram bins falling in that cluster

    # Plots by cell
    with col2:

        # Filter the cell data to the selected image
        cell_in_image = df_by_cell[image_colname] == selected_image
        num_cells_with_cluster_labels = df_by_cell['cluster_label'].loc[cell_in_image].notnull().sum()
        df_by_cell_filtered = df_by_cell[cell_in_image]
        df_by_cell_filtered.loc[:, 'cluster_label'] = df_by_cell_filtered['cluster_label'].cat.remove_unused_categories()

        # Write some image information
        st.write(f'In the selected image, there are {cell_in_image.sum()} cells present, {num_cells_with_cluster_labels} of which have been assigned a cluster label.')
        st.write(f'Label for image {selected_image}: {current_image_label}.')

        # Draw the plots
        st.plotly_chart(plotly_scatter_plot(df_by_cell_filtered, x_colname=spatial_x_colname, y_colname=spatial_y_colname, label_colname='cluster_label', unique_labels=unique_cluster_labels, plot_title=f'Spatial by Cell for {selected_image}'))
        df_by_cell['opacity'] = 0.05
        df_by_cell.loc[cell_in_image, 'opacity'] = 1
        st.plotly_chart(plotly_scatter_plot(df_by_cell, x_colname=umap_x_colname, y_colname=umap_y_colname, label_colname='cluster_label', unique_labels=unique_cluster_labels, plot_title=f'UMAP by Cell for {selected_image}', opacity_colname='opacity'))  # note that not all cells in the current image (e.g., cell_in_image.sum()=1000) actually have valid cluster labels (e.g., df_by_cell.loc[cell_in_image, 'cluster_label'].isin(unique_cluster_labels).sum())
        st.plotly_chart(px.line(df_by_cell_filtered.groupby('cluster_label', observed=True)[property_colnames].mean().reset_index().melt(id_vars='cluster_label', var_name='column', value_name='value'), x='column', y='value', color='cluster_label', markers=True, title=f'Property Means by Cell for {selected_image}', category_orders={'cluster_label': unique_cluster_labels}))  # get the neighbor vectors for each cluster averaged over the cells falling in that cluster


def get_score_and_prediction(value_counts, actual_label=None):

    # Convert the value counts to fractions
    vc_frac = value_counts / value_counts.sum()

    # If we want an accuracy score, which we can get when we provide the true value...
    if actual_label is not None:

        # If the standard two-value, binary case...
        if len(vc_frac) == 2:
            predicted_label = vc_frac.index[0]  # value counts is in decreasing order by default so the higher one is first, i.e., index=0
            score = vc_frac.loc[actual_label]

        # If we only have a single value...
        elif len(vc_frac) == 1:
            missing_label = (set([0, 1]) - set(vc_frac.index)).pop()
            vc_frac = pd.concat([vc_frac, pd.Series([0], index=[missing_label])])  # the existing counts frac must be 1, so the one for the missing label must be 0 since they must add to 1
            predicted_label = vc_frac.index[0]  # value counts is in decreasing order by default so the higher one is first, i.e., index=0
            score = vc_frac.loc[actual_label]

        # If there are no values...
        elif len(vc_frac) == 0:
            predicted_label = None
            score = 0  # technically this is undefined but we'll just say it's 0 because we really want to penalize not being able to make a prediction

        # If something unexpected happened...
        else:
            raise ValueError(f'The value counts series has too many values: {vc_frac}')
        
    # If we really just want the prediction and don't know the true value...
    else:  # here, score is somewhat meaningless (maybe means confidence?) since we don't have the true value to evaluate against

        # If value counts has any number of values at all...
        if len(vc_frac) in {1, 2}:
            predicted_label = vc_frac.index[0]  # value counts is in decreasing order by default so the higher one is first, i.e., index=0
            score = vc_frac.iloc[0]

        # If there are no values...
        elif len(vc_frac) == 0:
            predicted_label = None
            score = None

        # If something unexpected happened...
        else:
            raise ValueError(f'The value counts series has too many values: {vc_frac}')

    # Return the predicted label and score
    return predicted_label, score


def get_predictions(umap_test_colname, diff_cutoff_frac, image_name, df_by_bin, df_by_cell, df_image_labels, binary_colname, image_colname, debug=False):

    # Get the actual, true label for the current image and store it and other external information in a dictionary
    actual_label = df_image_labels.loc[image_name, binary_colname]
    predictions_dict = {'umap_test_colname': umap_test_colname, 'diff_cutoff_frac': diff_cutoff_frac, 'image_name': image_name, 'actual_label': actual_label}

    # Filter the binned data to the selected image
    image_in_set = pd.Series([image_name in images for images in df_by_bin['unique_images']], index=df_by_bin.index)
    num_bins_with_cluster_labels = df_by_bin['cluster_label'].loc[image_in_set].notnull().sum()
    df_by_bin_filtered = df_by_bin[image_in_set]

    # Get the estimate using the by-bin analysis
    vc = df_by_bin_filtered['cluster_label'].value_counts()
    assert vc.sum() == num_bins_with_cluster_labels
    estimated_label, score = get_score_and_prediction(vc, actual_label=actual_label)
    if debug:
        print(f'For image {image_name}, by bin: actual_label={actual_label}, estimated_label={estimated_label}, score={score}')
    predictions_dict['estimated_label_by_bin'] = estimated_label
    predictions_dict['score_by_bin'] = score

    # Filter the cell data to the selected image
    cell_in_image = df_by_cell[image_colname] == image_name
    num_cells_with_cluster_labels = df_by_cell['cluster_label'].loc[cell_in_image].notnull().sum()
    df_by_cell_filtered = df_by_cell[cell_in_image]

    # Get the estimate using the by-cell analysis
    vc = df_by_cell_filtered['cluster_label'].value_counts()
    assert vc.sum() == num_cells_with_cluster_labels
    estimated_label, score = get_score_and_prediction(vc, actual_label=actual_label)
    if debug:
        print(f'For image {image_name}, by cell: actual_label={actual_label}, estimated_label={estimated_label}, score={score}')
    predictions_dict['estimated_label_by_cell'] = estimated_label
    predictions_dict['score_by_cell'] = score

    # Return the dictionary containing the current prediction
    return predictions_dict


# Main function definition
def main():

    # Define some constants
    data_filename = 'Combo_CSVfiles_20230327_152849.csv'
    image_colname = 'ShortName'
    spatial_x_colname = 'CentroidX'
    spatial_y_colname = 'CentroidY'
    property_colnames = ['XMin', 'XMax', 'YMin', 'YMax']
    binary_colnames = ['Outcome', 'Survival_5yr']
    num_umap_bins = 200
    diff_cutoff_frac_default = 0.375  # Outcome: 0.15
    min_cells_per_bin = 1
    plot_manual_histogram_diff = False
    workflow_options = ['Figure visualization', 'Prediction']

    # Generate three columns for the settings and umap differences
    col1, col2, col3 = st.columns(3)

    # In the first column...
    with col1:

        # Click a button to load the data
        if st.button('Load data'):
            st.session_state['df'] = pd.read_csv(os.path.join('.', 'input', data_filename))
            st.session_state['unique_images'] = list(st.session_state['df'][image_colname].unique())

        # Ensure the data has been loaded
        if 'df' not in st.session_state:
            st.warning('Please load the data.')
            return
        else:
            df = st.session_state['df']
            unique_images = st.session_state['unique_images']

        # Click a button to load the densities from disk
        densities_path = os.path.join('.', 'output', 'npc_densities.pkl')
        densities_file_exists = os.path.exists(densities_path)
        if not densities_file_exists:
            load_densities_help = f'Densities file {densities_path} not found'
        else:
            load_densities_help = None
        if st.button('Load densities', disabled=(not densities_file_exists), help=load_densities_help):
            densities = pd.read_pickle(densities_path)
            df.drop(columns=[column for column in densities.columns if column in df.columns], inplace=True)
            len_df_orig = len(df)
            df = pd.concat([df, densities], axis='columns')
            assert len(df) == len_df_orig, f'Length of df after concatenation is {len(df)} but should be {len_df_orig}.'
            st.session_state['df'] = df
            
        # Click a button to calculate the densities
        if densities_file_exists:
            calculate_densities_help = f'WARNING: Clicking this will overwrite the current densities file {densities_path}'
            button_text = '⚠️ Calculate densities'
        else:
            calculate_densities_help = None
            button_text = 'Calculate densities'
        if st.button(button_text, help=calculate_densities_help):
            densities = neighborhood_profiles_checks.run_density_calculation()
            densities.to_pickle(densities_path)
            df.drop(columns=[column for column in densities.columns if column in df.columns], inplace=True)
            len_df_orig = len(df)
            df = pd.concat([df, densities], axis='columns')
            assert len(df) == len_df_orig, f'Length of df after concatenation is {len(df)} but should be {len_df_orig}.'
            st.session_state['df'] = df

        # Get the density columns
        density_columns = [col for col in st.session_state['df'].columns if col.startswith('spatial_umap_density of ')]
        if len(density_columns) == 0:
            st.warning('Please load or calculate densities.')
            return
        else:
            df = st.session_state['df']

        # Create a dropdown for the user to choose the binary label of interest
        if 'npc__binary_colname' not in st.session_state:
            st.session_state['npc__binary_colname'] = binary_colnames[0]
        binary_colname = st.selectbox('Binary label:', binary_colnames, key='npc__binary_colname')

        # Get the actual binary value/label for each image
        df_image_labels = df[[image_colname, binary_colname]].groupby(image_colname).agg(set)
        assert df_image_labels[binary_colname].apply(len).max() == 1, f'There are images with multiple binary labels: {df_image_labels[df_image_labels[binary_colname].apply(len) > 1]}'
        df_image_labels[binary_colname] = df_image_labels[binary_colname].apply(lambda x: list(x)[0])

        # Allow user to select the fractions of the size of the smallest image to use for UMAP training and "testing"
        if 'npc__frac_train' not in st.session_state:
            st.session_state['npc__frac_train'] = 0.1
        frac_train = st.number_input('Fraction of the number of cells in the smallest image to use for training the UMAP:', key='npc__frac_train', format='%.2f', min_value=0.0, max_value=0.5)
        if 'npc__frac_test' not in st.session_state:
            st.session_state['npc__frac_test'] = 0.1
        frac_test = st.number_input('Fraction of the number of cells in the smallest image to use for generating the 2D histogram:', key='npc__frac_test', format='%.2f', min_value=0.0, max_value=0.5)

        # Create a dropdown for the user to choose the workflow
        if 'workflow' not in st.session_state:
            st.session_state['workflow'] = workflow_options[0]
        st.selectbox('Workflow:', workflow_options, key='workflow')

        # Write a divider
        st.divider()

    # If we want to run the figure checks visualization workflow...
    if st.session_state['workflow'] == workflow_options[0]:

        # In the first column...
        with col1:

            # Add columns to partition the dataframe into train and test sets for the UMAP, roughly similar to Giraldo et. al. 2021
            if st.button('Partition dataset into train and test sets for the UMAP'):
                st.session_state['df'] = neighborhood_profiles_checks.get_umap_train_and_test_sets(df, frac_train=frac_train, frac_test=frac_test, image_colname=image_colname, num_umap_test_sets=1)

            # Ensure partitioning has been performed
            if 'umap_train' not in st.session_state['df']:
                st.warning('Please partition the dataset into train and test sets for the UMAP.')
                return
            else:
                df = st.session_state['df']
            
            # Click a button to load the UMAP from disk
            umap_path = os.path.join('.', 'output', 'npc_umap.pkl')
            umap_file_exists = os.path.exists(umap_path)
            if not umap_file_exists:
                load_umap_help = f'UMAP file {umap_path} not found'
            else:
                load_umap_help = None
            if st.button('Load UMAP', disabled=(not umap_file_exists), help=load_umap_help):
                umap = pd.read_pickle(umap_path)
                df.drop(columns=[column for column in umap.columns if column in df.columns], inplace=True)
                len_df_orig = len(df)
                df = pd.concat([df, umap], axis='columns')
                assert len(df) == len_df_orig, f'Length of df after concatenation is {len(df)} but should be {len_df_orig}.'
                st.session_state['df'] = df
                
            # Click a button to calculate the UMAP
            if umap_file_exists:
                calculate_umap_help = f'WARNING: Clicking this will overwrite the current UMAP file {umap_path}'
                button_text = '⚠️ Calculate UMAP'
            else:
                calculate_umap_help = None
                button_text = 'Calculate UMAP'
            if st.button(button_text, help=calculate_umap_help):
                umap = neighborhood_profiles_checks.run_umap_calculation()
                umap.to_pickle(umap_path)
                df.drop(columns=[column for column in umap.columns if column in df.columns], inplace=True)
                len_df_orig = len(df)
                df = pd.concat([df, umap], axis='columns')
                assert len(df) == len_df_orig, f'Length of df after concatenation is {len(df)} but should be {len_df_orig}.'
                st.session_state['df'] = df

            # Get the UMAP columns
            umap_column_options = [col for col in st.session_state['df'].columns if col.startswith('UMAP_')]
            if len(umap_column_options) == 0:
                st.warning('Please load UMAP, calculate UMAP, or load a dataset that contains UMAP columns.')
                return
            else:
                df = st.session_state['df']
            
            # Create a dropdown for the user to choose the UMAP columns
            if 'npc__umap_x_colname' not in st.session_state:
                st.session_state['npc__umap_x_colname'] = umap_column_options[0]
            umap_x_colname = st.selectbox('UMAP x column:', umap_column_options, key='npc__umap_x_colname')
            if 'npc__umap_y_colname' not in st.session_state:
                st.session_state['npc__umap_y_colname'] = umap_column_options[1]
            umap_y_colname = st.selectbox('UMAP y column:', umap_column_options, key='npc__umap_y_colname')

            # Write a number_input widget for diff_cutoff_frac
            if 'diff_cutoff_frac' not in st.session_state:
                st.session_state['diff_cutoff_frac'] = diff_cutoff_frac_default
            diff_cutoff_frac = st.number_input('Difference cutoff fraction:', key='diff_cutoff_frac', format='%.3f')

            # Click a button to calculate a dictionary of clusters as keys and list of bin tuples as values using the normalized histogram differences between two conditions for a single set of UMAP "test" data
            if st.button('Calculate the difference clusters'):
                st.session_state['npc__clusters'], st.session_state['npc__edges_x'], st.session_state['npc__edges_y'], st.session_state['npc__fig_d_diff'], st.session_state['npc__fig_d_diff_manual'] = neighborhood_profiles_checks.calculate_difference_clusters(df, diff_cutoff_frac, umap_x_colname=umap_x_colname, umap_y_colname=umap_y_colname, binary_colname=binary_colname, num_umap_bins=num_umap_bins, plot_manual_histogram_diff=plot_manual_histogram_diff, plot_diff_matrix=True, umap_test_colname='umap_test_0')

            # Ensure the data has been processed
            if 'npc__clusters' not in st.session_state:
                st.warning('Please calculate the difference clusters.')
                return
            else:
                clusters = st.session_state['npc__clusters']
                edges_x = st.session_state['npc__edges_x']
                edges_y = st.session_state['npc__edges_y']
                fig_d_diff = st.session_state['npc__fig_d_diff']
                fig_d_diff_manual = st.session_state['npc__fig_d_diff_manual']
            
            # Click a button to perform UMAP binning
            if st.button('Perform UMAP binning'):

                # Determine the bins of the UMAP x and y coordinates of an input dataframe
                st.session_state['df_test_with_bins'] = neighborhood_profiles_checks.perform_umap_binning(df[df['umap_test_0']], edges_x, edges_y, image_colname=image_colname, spatial_x_colname=spatial_x_colname, spatial_y_colname=spatial_y_colname, umap_x_colname=umap_x_colname, umap_y_colname=umap_y_colname, property_colnames=property_colnames)

            # Ensure UMAP binning has been performed
            if 'df_test_with_bins' not in st.session_state:
                st.warning('Please perform UMAP binning.')
                return
            else:
                df_test_with_bins = st.session_state['df_test_with_bins']

            # Click a button to assign the cluster labels
            if st.button('Assign cluster labels'):

                # Given the UMAP bins of each cell and the labels assigned to the bins, assign the labels to a per-cell dataframe and to a transformed dataframe grouped by bin, i.e., assign the label to each bin in that transformed dataframe
                st.session_state['df_by_bin'], st.session_state['df_by_cell'] = neighborhood_profiles_checks.assign_cluster_labels(df_test_with_bins, clusters, image_colname=image_colname, min_cells_per_bin=min_cells_per_bin)

        # In the second column, plot Giraldo's difference histogram
        with col2:
            if fig_d_diff is not None:
                st.write('Difference matrix:')
                st.pyplot(fig_d_diff)

        # In the third column, optionally plot the manual difference histogram
        with col3:
            if plot_manual_histogram_diff:
                if fig_d_diff_manual is not None:
                    st.write('Difference matrix - manual')
                    st.pyplot(fig_d_diff_manual)

        # Ensure we're ready to display the figure checks
        with col1:
            if 'df_by_bin' not in st.session_state:
                st.warning('Please assign cluster labels.')
                return
        
        # Draw the plots
        draw_plots(df_image_labels=df_image_labels, umap_x_colname=umap_x_colname, umap_y_colname=umap_y_colname, property_colnames=property_colnames, image_colname=image_colname, binary_colname=binary_colname, spatial_x_colname=spatial_x_colname, spatial_y_colname=spatial_y_colname)

    # If we want to run the prediction workflow...
    else:

        # In the first column, draw widgets for the prediction workflow
        with col1:
            if 'npc__num_umap_test_sets' not in st.session_state:
                st.session_state['npc__num_umap_test_sets'] = 10
            num_umap_test_sets = st.number_input('Number of UMAP "test" sets to generate:', min_value=1, key='npc__num_umap_test_sets')
            if 'npc__cutoff_frac_start' not in st.session_state:
                st.session_state['npc__cutoff_frac_start'] = 0.05  # Outcome: 0.075
            cutoff_frac_start = st.number_input('Cutoff fraction start:', min_value=0.0, max_value=1.0, format='%.3f', key='npc__cutoff_frac_start')
            if 'npc__cutoff_frac_end' not in st.session_state:
                st.session_state['npc__cutoff_frac_end'] = 0.5  # Outcome: 0.3
            cutoff_frac_end = st.number_input('Cutoff fraction end:', min_value=0.0, max_value=1.0, format='%.3f', key='npc__cutoff_frac_end')
            if 'npc__cutoff_frac_step' not in st.session_state:
                st.session_state['npc__cutoff_frac_step'] = 0.025  # Outcome: 0.025
            cutoff_frac_step = st.number_input('Cutoff fraction step:', min_value=0.0, max_value=1.0, format='%.3f', key='npc__cutoff_frac_step')

        # If we're ready to generate predictions...
        if st.button('Generate predictions'):
            with st.spinner('Generating predictions...'):

                # Add columns to partition the dataframe into train and test sets for the UMAP, roughly similar to Giraldo et. al. 2021
                st.session_state['df'] = neighborhood_profiles_checks.get_umap_train_and_test_sets(df, frac_train=frac_train, frac_test=frac_test, image_colname=image_colname, num_umap_test_sets=num_umap_test_sets)
                df = st.session_state['df']

                # TODO: Run the UMAP on the train set based on the following code block
                # umap = neighborhood_profiles_checks.run_umap_calculation()
                # umap.to_pickle(umap_path)
                # df.drop(columns=[column for column in umap.columns if column in df.columns], inplace=True)
                # len_df_orig = len(df)
                # df = pd.concat([df, umap], axis='columns')
                # assert len(df) == len_df_orig, f'Length of df after concatenation is {len(df)} but should be {len_df_orig}.'
                # st.session_state['df'] = df

                # For now we're just going to use the existing UMAP columns
                umap_column_options = [col for col in st.session_state['df'].columns if col.startswith('UMAP_')]
                df = st.session_state['df']
                umap_x_colname = umap_column_options[0]
                umap_y_colname = umap_column_options[1]

                # Iterate over possible difference cutoff fractions
                predictions_holder = []
                if cutoff_frac_step != 0:
                    cutoff_frac_range = np.arange(cutoff_frac_start, cutoff_frac_end + cutoff_frac_step, cutoff_frac_step)
                else:
                    cutoff_frac_range = [cutoff_frac_start]
                for diff_cutoff_frac in cutoff_frac_range:

                    # For every UMAP "test" i.e. calculation set...
                    for umap_test_colname in [column for column in df.columns if column.startswith('umap_test_')]:

                        # Calculate a dictionary of clusters as keys and list of bin tuples as values using the normalized histogram differences between two conditions for the current set of UMAP "test" data
                        clusters, edges_x, edges_y, fig_d_diff, fig_d_diff_manual = neighborhood_profiles_checks.calculate_difference_clusters(df, diff_cutoff_frac, umap_x_colname=umap_x_colname, umap_y_colname=umap_y_colname, binary_colname=binary_colname, num_umap_bins=num_umap_bins, plot_manual_histogram_diff=False, plot_diff_matrix=False, umap_test_colname=umap_test_colname)

                        # Determine the bins of the UMAP x and y coordinates of an input dataframe
                        df_test_with_bins = neighborhood_profiles_checks.perform_umap_binning(df[df[umap_test_colname]], edges_x, edges_y, image_colname=image_colname, spatial_x_colname=spatial_x_colname, spatial_y_colname=spatial_y_colname, umap_x_colname=umap_x_colname, umap_y_colname=umap_y_colname, property_colnames=property_colnames)

                        # Given the UMAP bins of each cell and the labels assigned to the bins, assign the labels to a per-cell dataframe and to a transformed dataframe grouped by bin, i.e., assign the label to each bin in that transformed dataframe
                        df_by_bin, df_by_cell = neighborhood_profiles_checks.assign_cluster_labels(df_test_with_bins, clusters, image_colname=image_colname, min_cells_per_bin=min_cells_per_bin)

                        # For each image in the dataset, predict the binary label
                        for image_name in unique_images:
                            predictions_dict = get_predictions(umap_test_colname=umap_test_colname, diff_cutoff_frac=diff_cutoff_frac, image_name=image_name, df_by_bin=df_by_bin, df_by_cell=df_by_cell, df_image_labels=df_image_labels, binary_colname=binary_colname, image_colname=image_colname, debug=False)
                            predictions_holder.append(predictions_dict)

                # Create a dataframe holding all the predictions and save it to the session state
                df_predictions = pd.DataFrame(predictions_holder)
                st.session_state['df_predictions'] = df_predictions

        # If predictions have been made...
        if 'df_predictions' in st.session_state:

            # Get a shortcut to the predictions dataframe
            df_predictions = st.session_state['df_predictions']

            # Add to it whether the predictions were correct
            df_predictions['bin_label_correct'] = df_predictions['actual_label'] == df_predictions['estimated_label_by_bin']
            df_predictions['cell_label_correct'] = df_predictions['actual_label'] == df_predictions['estimated_label_by_cell']

            # Write the predictions to the screen
            st.write(df_predictions)

            # Get the average scores by diff_cutoff_frac (as the index), write it to screen, and plot it
            df_scores_vs_cutoff = df_predictions.groupby(['umap_test_colname', 'diff_cutoff_frac'])[['score_by_bin', 'score_by_cell']].mean().reset_index().groupby('diff_cutoff_frac')[['score_by_bin', 'score_by_cell']].mean()
            st.write(df_scores_vs_cutoff)
            st.plotly_chart(px.line(df_scores_vs_cutoff.reset_index().melt(id_vars='diff_cutoff_frac', var_name='Column', value_name='Value'), x='diff_cutoff_frac', y='Value', color='Column', markers=True))

            # Write to screen the number of times each image was predicted correctly by cell
            st.write(df_predictions.groupby('image_name')['cell_label_correct'].sum())


# Main script block
if __name__ == '__main__':

    # Set page settings
    page_name = 'Neighborhood Profiles Checker'
    st.set_page_config(layout='wide', page_title=page_name)
    st.title(page_name)

    # Call the main function
    main()
