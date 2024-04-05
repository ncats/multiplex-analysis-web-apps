# Import relevant libraries
import neighborhood_profiles_checks
import pandas as pd
import os
import numpy as np
import PlottingTools_orig as PlottingTools
import streamlit as st
import matplotlib.pyplot as plt
import plotly.express as px


def load_data(data_filename='Combo_CSVfiles_20230327_152849.csv', image_colname='ShortName'):
    with st.spinner('Loading data...'):
        st.session_state['df'] = pd.read_csv(os.path.join('.', 'input', data_filename))
        st.session_state['unique_images'] = list(st.session_state['df'][image_colname].unique())


def process_data(image_colname='ShortName', umap_x_colname='UMAP_1_20230327_152849', umap_y_colname='UMAP_2_20230327_152849', binary_colname='Survival_5yr', number_of_samples_frac=0.1, num_umap_bins=200, plot_manual_histogram_diff=False):

    with st.spinner('Processing data...'):

        # Get shortcuts to session state variables
        df = st.session_state['df']
        diff_cutoff_frac = st.session_state['diff_cutoff_frac']

        # Get an equal number of cells from each image
        num_cells_from_each_sample = int(round(df.groupby(image_colname).size().min() * number_of_samples_frac, -3))
        indices = df.groupby(image_colname).apply(lambda x: x.sample(n=num_cells_from_each_sample, replace=False).index, include_groups=False).explode().values

        # Assign the sample to 'umap_test'
        df['umap_test'] = False
        df.loc[indices, 'umap_test'] = True

        # Get a universal set of edges for the UMAPs
        edges_x = np.linspace(df[umap_x_colname].min(), df[umap_x_colname].max(), num_umap_bins + 1)
        edges_y = np.linspace(df[umap_y_colname].min(), df[umap_y_colname].max(), num_umap_bins + 1)

        # Get subsets of the full data that are the entire test set and both binary subsets of the test set
        df_test = df[df['umap_test']]
        binary_values = sorted(list(df[binary_colname].unique()))
        df_binary_subset_0 = df_test[df_test[binary_colname] == binary_values[0]]
        df_binary_subset_1 = df_test[df_test[binary_colname] == binary_values[1]]

        # Get the 2D histograms for each condition
        d_binary_subset_0 = PlottingTools.plot_2d_density(df_binary_subset_0[umap_x_colname], df_binary_subset_0[umap_y_colname], bins=[edges_x, edges_y], return_matrix=True)
        d_binary_subset_1 = PlottingTools.plot_2d_density(df_binary_subset_1[umap_x_colname], df_binary_subset_1[umap_y_colname], bins=[edges_x, edges_y], return_matrix=True)

        # Get the difference between the histograms for the binary subsets
        d_diff = d_binary_subset_1 - d_binary_subset_0

        # "Mask" the difference matrix based on a cutoff
        cutoff = np.abs(d_diff).max() * diff_cutoff_frac
        d_diff[d_diff > cutoff] = 1
        d_diff[d_diff < -cutoff] = -1
        d_diff[(d_diff >= -cutoff) & (d_diff <= cutoff)] = 0

        # Get the clusters based only on the cutoff, deliberately not performing any clustering
        clusters = {0: [tuple([x[1], x[0]]) for x in np.transpose(np.where(np.isclose(d_diff, -1)))], 1: [tuple([x[1], x[0]]) for x in np.transpose(np.where(np.isclose(d_diff, 1)))]}

        # Plot the difference matrix
        fig_d_diff, ax_d_diff = plt.subplots()
        PlottingTools.plot_2d_density(d_diff, bins=[edges_x, edges_y], n_pad=30, circle_type='arch', cmap=plt.get_cmap('bwr'), ax=ax_d_diff)

        # Optionally plot the difference matrix manually
        if plot_manual_histogram_diff:
            fig, ax = plt.subplots()
            c = ax.pcolormesh(edges_x, edges_y, d_diff, cmap='bwr')
            fig.colorbar(c, ax=ax)
            st.session_state['fig_d_diff_manual'] = fig

        # Save the processed data to the session state
        st.session_state['processed_data'] = {'edges_x': edges_x, 'edges_y': edges_y, 'clusters': clusters}

        # Save the difference plot to the session state
        st.session_state['fig_d_diff'] = fig_d_diff


def run_checks(image_colname='ShortName', spatial_x_colname='CentroidX', spatial_y_colname='CentroidY', umap_x_colname='UMAP_1_20230327_152849', umap_y_colname='UMAP_2_20230327_152849', property_colnames=['XMin', 'XMax', 'YMin', 'YMax'], min_cells_per_bin=1):

    with st.spinner('Running checks...'):

        # Get shortcuts to the original and processed data
        df = st.session_state['df']
        edges_x = st.session_state['processed_data']['edges_x']
        edges_y = st.session_state['processed_data']['edges_y']
        clusters = st.session_state['processed_data']['clusters']
        
        # Perform binning
        df_test = neighborhood_profiles_checks.perform_binning(df, edges_x, edges_y, image_colname=image_colname, spatial_x_colname=spatial_x_colname, spatial_y_colname=spatial_y_colname, umap_x_colname=umap_x_colname, umap_y_colname=umap_y_colname, property_colnames=property_colnames)

        # Generate the checks as plotly figures
        df_by_bin, df_by_cell = neighborhood_profiles_checks.assign_cluster_labels(df_test, clusters, image_colname=image_colname, min_cells_per_bin=min_cells_per_bin)

        # Save the figures to the session state
        st.session_state['df_by_bin'] = df_by_bin
        st.session_state['df_by_cell'] = df_by_cell


def draw_plots(df_image_labels, umap_x_colname='UMAP_1_20230327_152849', umap_y_colname='UMAP_2_20230327_152849', property_colnames=['XMin', 'XMax', 'YMin', 'YMax'], image_colname='ShortName', binary_colname='Survival_5yr', spatial_x_colname='CentroidX', spatial_y_colname='CentroidY'):
    
    # Get shortcuts to variables in the session state
    df_by_bin = st.session_state['df_by_bin']
    df_by_cell = st.session_state['df_by_cell']
    unique_images = st.session_state['unique_images']

    # Write a header
    st.header('Plots by dataset')

    # For plotting purposes below, we need to convert the 'cluster_label' column to categorical, and also save their unique values so we can preserve plotting order
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
        col1.plotly_chart(px.scatter(df_by_bin, x=umap_x_colname, y=umap_y_colname, color='cluster_label', title='UMAP by Bin for Whole Dataset', category_orders={'cluster_label': unique_cluster_labels}).update_layout(xaxis={"scaleanchor": "y", "scaleratio": 1}))
        col1.plotly_chart(px.line(df_by_bin.groupby('cluster_label')[property_colnames].mean().reset_index().melt(id_vars='cluster_label', var_name='column', value_name='value'), x='column', y='value', color='cluster_label', markers=True, title='Property Means by Bin for Whole Dataset', category_orders={'cluster_label': unique_cluster_labels}))  # get the neighbor vectors for each cluster averaged over the histogram bins falling in that cluster

    # Plots by cell
    with col2:
        col2.plotly_chart(px.scatter(df_by_cell, x=umap_x_colname, y=umap_y_colname, color='cluster_label', title='UMAP by Cell for Whole Dataset', category_orders={'cluster_label': unique_cluster_labels}).update_layout(xaxis={"scaleanchor": "y", "scaleratio": 1}))
        col2.plotly_chart(px.line(df_by_cell.groupby('cluster_label')[property_colnames].mean().reset_index().melt(id_vars='cluster_label', var_name='column', value_name='value'), x='column', y='value', color='cluster_label', markers=True, title='Property Means by Cell for Whole Dataset', category_orders={'cluster_label': unique_cluster_labels}))  # get the neighbor vectors for each cluster averaged over the cells falling in that cluster

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
        df_by_bin_filtered['cluster_label'] = df_by_bin_filtered['cluster_label'].cat.remove_unused_categories()

        # Write some image information
        st.write(f'In the selected image, there are {image_in_set.sum()} bins present, {num_bins_with_cluster_labels} of which have been assigned a cluster label.')
        st.write(f'Label for image {selected_image}: {current_image_label}.')

        # Draw the plots
        st.plotly_chart(px.scatter(df_by_bin_filtered, x=spatial_x_colname, y=spatial_y_colname, color='cluster_label', title=f'Spatial by Bin for {selected_image} (this is meaningless; don\'t read into it)', category_orders={'cluster_label': unique_cluster_labels}).update_layout(xaxis={"scaleanchor": "y", "scaleratio": 1}))
        st.plotly_chart(px.scatter(df_by_bin_filtered, x=umap_x_colname, y=umap_y_colname, color='cluster_label', title=f'UMAP by Bin for {selected_image}', category_orders={'cluster_label': unique_cluster_labels}).update_layout(xaxis={"scaleanchor": "y", "scaleratio": 1}))
        st.plotly_chart(px.line(df_by_bin_filtered.groupby('cluster_label')[property_colnames].mean().reset_index().melt(id_vars='cluster_label', var_name='column', value_name='value'), x='column', y='value', color='cluster_label', markers=True, title=f'Property Means by Bin for {selected_image}', category_orders={'cluster_label': unique_cluster_labels}))  # get the neighbor vectors for each cluster averaged over the histogram bins falling in that cluster

    # Plots by cell
    with col2:

        # Filter the cell data to the selected image
        cell_in_image = df_by_cell[image_colname] == selected_image
        num_cells_with_cluster_labels = df_by_cell['cluster_label'].loc[cell_in_image].notnull().sum()
        df_by_cell_filtered = df_by_cell[cell_in_image]
        df_by_cell_filtered['cluster_label'] = df_by_cell_filtered['cluster_label'].cat.remove_unused_categories()

        # Write some image information
        st.write(f'In the selected image, there are {cell_in_image.sum()} cells present, {num_cells_with_cluster_labels} of which have been assigned a cluster label.')
        st.write(f'Label for image {selected_image}: {current_image_label}.')

        # Draw the plots
        st.plotly_chart(px.scatter(df_by_cell_filtered, x=spatial_x_colname, y=spatial_y_colname, color='cluster_label', title=f'Spatial by Cell for {selected_image}', category_orders={'cluster_label': unique_cluster_labels}).update_layout(xaxis={"scaleanchor": "y", "scaleratio": 1}))
        st.plotly_chart(px.scatter(df_by_cell_filtered, x=umap_x_colname, y=umap_y_colname, color='cluster_label', title=f'UMAP by Cell for {selected_image}', category_orders={'cluster_label': unique_cluster_labels}).update_layout(xaxis={"scaleanchor": "y", "scaleratio": 1}))
        st.plotly_chart(px.line(df_by_cell_filtered.groupby('cluster_label')[property_colnames].mean().reset_index().melt(id_vars='cluster_label', var_name='column', value_name='value'), x='column', y='value', color='cluster_label', markers=True, title=f'Property Means by Cell for {selected_image}', category_orders={'cluster_label': unique_cluster_labels}))  # get the neighbor vectors for each cluster averaged over the cells falling in that cluster


# Main function definition
def main():

    # Define some constants
    data_filename = 'Combo_CSVfiles_20230327_152849.csv'
    image_colname = 'ShortName'
    spatial_x_colname = 'CentroidX'
    spatial_y_colname = 'CentroidY'
    umap_x_colname = 'UMAP_1_20230327_152849'
    umap_y_colname = 'UMAP_2_20230327_152849'
    property_colnames = ['XMin', 'XMax', 'YMin', 'YMax']
    binary_colname = 'Survival_5yr'
    number_of_samples_frac = 0.1
    num_umap_bins = 200
    diff_cutoff_frac_default = 0.2
    min_cells_per_bin = 1
    plot_manual_histogram_diff = False

    # Generate three columns for the settings and umap differences
    col1, col2, col3 = st.columns(3)

    # In the first column...
    with col1:

        # Click a button to load the data
        st.button('Load data', on_click=load_data, kwargs={'data_filename': data_filename, 'image_colname': image_colname})

        # Ensure the data has been loaded
        if 'df' not in st.session_state:
            st.warning('Please load the data.')
            return

        # Write a number_input widget for diff_cutoff_frac
        if 'diff_cutoff_frac' not in st.session_state:
            st.session_state['diff_cutoff_frac'] = diff_cutoff_frac_default
        st.number_input('Difference cutoff fraction:', key='diff_cutoff_frac')

        # Click a button to process the data
        st.button('Process the data', on_click=process_data, kwargs={'image_colname': image_colname, 'umap_x_colname': umap_x_colname, 'umap_y_colname': umap_y_colname, 'binary_colname': binary_colname, 'number_of_samples_frac': number_of_samples_frac, 'num_umap_bins': num_umap_bins, 'plot_manual_histogram_diff': plot_manual_histogram_diff})

        # Ensure the data has been processed
        if 'processed_data' not in st.session_state:
            st.warning('Please process the data.')
            return
        
        # Click a button to run the checks
        st.button('Run checks', on_click=run_checks, kwargs={'image_colname': image_colname, 'spatial_x_colname': spatial_x_colname, 'spatial_y_colname': spatial_y_colname, 'umap_x_colname': umap_x_colname, 'umap_y_colname': umap_y_colname, 'property_colnames': property_colnames, 'min_cells_per_bin': min_cells_per_bin})

    # In the second column, plot Giraldo's difference histogram
    with col2:
        if 'fig_d_diff' in st.session_state:
            st.write('Difference matrix:')
            st.pyplot(st.session_state['fig_d_diff'])

    # In the third column, optionally plot the manual difference histogram
    with col3:
        if plot_manual_histogram_diff:
            if 'fig_d_diff_manual' in st.session_state:
                st.write('Difference matrix - manual')
                st.pyplot(st.session_state['fig_d_diff_manual'])

    # Ensure we're ready to display the figure checks
    with col1:
        if 'df_by_bin' not in st.session_state:
            st.warning('Please run the checks.')
            return
    
    # Get the actual binary value/label for each image
    df_image_labels = st.session_state['df'][[image_colname, binary_colname]].groupby(image_colname).agg(set)
    assert df_image_labels[binary_colname].apply(len).max() == 1, f'There are images with multiple binary labels: {df_image_labels[df_image_labels[binary_colname].apply(len) > 1]}'
    df_image_labels[binary_colname] = df_image_labels[binary_colname].apply(lambda x: list(x)[0])

    # Draw the plots
    draw_plots(df_image_labels=df_image_labels, umap_x_colname=umap_x_colname, umap_y_colname=umap_y_colname, property_colnames=property_colnames, image_colname=image_colname, binary_colname=binary_colname, spatial_x_colname=spatial_x_colname, spatial_y_colname=spatial_y_colname)

    # Get shortcuts to variables in the session state
    df_by_bin = st.session_state['df_by_bin']
    df_by_cell = st.session_state['df_by_cell']

    # difference_cutoff_fractions = np.linspace(0.1, 0.9, 9)

    # Loop through all the images in the dataset
    for image_name in st.session_state['unique_images']:

        # Get the actual, true label for the current image
        actual_label = df_image_labels.loc[image_name, binary_colname]

        # Filter the binned data to the selected image
        image_in_set = pd.Series([image_name in images for images in df_by_bin['unique_images']], index=df_by_bin.index)
        num_bins_with_cluster_labels = df_by_bin['cluster_label'].loc[image_in_set].notnull().sum()
        df_by_bin_filtered = df_by_bin[image_in_set]
        df_by_bin_filtered['cluster_label'] = df_by_bin_filtered['cluster_label'].cat.remove_unused_categories()

        # Get the estimate using the by-bin analysis
        vc = df_by_bin_filtered['cluster_label'].value_counts()
        assert vc.sum() == num_bins_with_cluster_labels
        scores_by_bin = vc / vc.sum()
        if vc.sum() == 0:
            estimated_label = None
            score = None
        else:
            estimated_label = scores_by_bin.index[0]
            score = scores_by_bin.iloc[0]
        st.write(f'For image {image_name}, by bin: actual_label={actual_label}, estimated_label={estimated_label}, score={score}')

        # Filter the cell data to the selected image
        cell_in_image = df_by_cell[image_colname] == image_name
        num_cells_with_cluster_labels = df_by_cell['cluster_label'].loc[cell_in_image].notnull().sum()
        df_by_cell_filtered = df_by_cell[cell_in_image]
        df_by_cell_filtered['cluster_label'] = df_by_cell_filtered['cluster_label'].cat.remove_unused_categories()

        # Get the estimate using the by-cell analysis
        vc = df_by_cell_filtered['cluster_label'].value_counts()
        assert vc.sum() == num_cells_with_cluster_labels
        scores_by_cell = vc / vc.sum()
        if vc.sum() == 0:
            estimated_label = None
            score = None
        else:
            estimated_label = scores_by_cell.index[0]
            score = scores_by_cell.iloc[0]
        st.write(f'For image {image_name}, by cell: actual_label={actual_label}, estimated_label={estimated_label}, score={score}')


# Main script block
if __name__ == '__main__':

    # Set page settings
    page_name = 'Neighborhood Profiles Checker'
    st.set_page_config(layout='wide', page_title=page_name)
    st.title(page_name)

    # Call the main function
    main()
