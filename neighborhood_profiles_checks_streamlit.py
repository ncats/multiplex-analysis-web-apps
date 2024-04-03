# Import relevant libraries
import neighborhood_profiles_checks
import pandas as pd
import os
import numpy as np
import PlottingTools_orig as PlottingTools
import streamlit as st
import matplotlib.pyplot as plt


# Main function definition
def main():

    # Define some constants
    data_filename = 'Combo_CSVfiles_20230327_152849.csv'
    spatial_x_colname = 'CentroidX'
    spatial_y_colname = 'CentroidY'
    umap_x_colname = 'UMAP_1_20230327_152849'
    umap_y_colname = 'UMAP_2_20230327_152849'
    property_colnames = ['XMin', 'XMax', 'YMin', 'YMax']
    number_of_samples_frac = 0.1
    num_umap_bins = 200
    diff_cutoff_frac = 0.2
    min_cells_per_bin = 1

    # Click a button to load the data
    if st.button('Load data'):
        st.session_state['df'] = pd.read_csv(os.path.join('.', 'input', data_filename))

    # Click a button to process the data
    if st.button('Process the data'):

        # Ensure the data has been loaded
        if 'df' not in st.session_state:
            st.warning('Please load the data first.')
            return

        # Get a shortcut to the original dataframe
        df = st.session_state['df']

        # Get an equal number of cells from each image
        num_cells_from_each_sample = int(round(df.groupby('ShortName').size().min() * number_of_samples_frac, -3))
        indices = df.groupby('ShortName').apply(lambda x: x.sample(n=num_cells_from_each_sample, replace=False).index, include_groups=False).explode().values

        # Assign the sample to 'umap_test'
        df['umap_test'] = False
        df.loc[indices, 'umap_test'] = True

        # Get a universal set of edges for the UMAPs
        edges_x = np.linspace(df[umap_x_colname].min(), df[umap_x_colname].max(), num_umap_bins + 1)
        edges_y = np.linspace(df[umap_y_colname].min(), df[umap_y_colname].max(), num_umap_bins + 1)

        # Get subsets of the full data that are (1) the entire test set, (2) the alive cells in the test set, and (3) the dead cells in the test set
        df_test = df[df['umap_test']]
        df_alive = df_test[df_test['Survival_5yr'] == 1]  # assume "alive" means Survival_5yr == 1
        df_dead = df_test[df_test['Survival_5yr'] == 0]

        # Get the 2D histograms for each condition
        d_alive = PlottingTools.plot_2d_density(df_alive[umap_x_colname], df_alive[umap_y_colname], bins=[edges_x, edges_y], return_matrix=True)
        d_dead = PlottingTools.plot_2d_density(df_dead[umap_x_colname], df_dead[umap_y_colname], bins=[edges_x, edges_y], return_matrix=True)

        # Get the difference between the alive and dead histograms
        d_diff = d_alive - d_dead

        # "Mask" the difference matrix based on a cutoff
        cutoff = np.abs(d_diff).max() * diff_cutoff_frac
        d_diff[d_diff > cutoff] = 1
        d_diff[d_diff < -cutoff] = -1
        d_diff[(d_diff >= -cutoff) & (d_diff <= cutoff)] = 0

        # Get the clusters based only on the cutoff
        clusters = {0: [tuple(x) for x in np.transpose(np.where(d_diff == -1))], 1: [tuple(x) for x in np.transpose(np.where(d_diff == 1))]}

        # Plot the difference matrix
        fig_d_diff, ax_d_diff = plt.subplots()
        PlottingTools.plot_2d_density(d_diff, bins=[edges_x, edges_y], n_pad=30, circle_type='arch', cmap=plt.get_cmap('bwr').copy(), ax=ax_d_diff)

        # Save the processed data to the session state
        st.session_state['processed_data'] = {'df': df, 'edges_x': edges_x, 'edges_y': edges_y, 'clusters': clusters}

        # Save the difference plot to the session state
        st.session_state['fig_d_diff'] = fig_d_diff

    # Click a button to run the checks
    if st.button('Run checks'):

        # Ensure the data has been processed
        if 'processed_data' not in st.session_state:
            st.warning('Please process the data first.')
            return
        
        # Get shortcuts to the processed data
        df = st.session_state['processed_data']['df']
        edges_x = st.session_state['processed_data']['edges_x']
        edges_y = st.session_state['processed_data']['edges_y']
        clusters = st.session_state['processed_data']['clusters']
        
        # Perform binning
        df_test = neighborhood_profiles_checks.perform_binning(df, edges_x, edges_y, pd.Series(True, index=df.index), spatial_x_colname=spatial_x_colname, spatial_y_colname=spatial_y_colname, umap_x_colname=umap_x_colname, umap_y_colname=umap_y_colname, property_colnames=property_colnames)

        # Generate the checks as plotly figures
        fig_property_means_by_bin, fig_property_means_by_cell, fig_umap_by_bin, fig_spatial_by_bin, fig_umap_by_cell, fig_spatial_by_cell = neighborhood_profiles_checks.generate_figures(df_test, clusters, spatial_x_colname=spatial_x_colname, spatial_y_colname=spatial_y_colname, umap_x_colname=umap_x_colname, umap_y_colname=umap_y_colname, property_colnames=property_colnames, min_cells_per_bin=min_cells_per_bin)

        # Save the figures to the session state
        st.session_state['figs_checks'] = [fig_property_means_by_bin, fig_property_means_by_cell, fig_umap_by_bin, fig_spatial_by_bin, fig_umap_by_cell, fig_spatial_by_cell]

    # If the figure checks have been generated, display them
    if 'fig_d_diff' in st.session_state:
        st.write('Difference matrix')
        st.pyplot(st.session_state['fig_d_diff'])
    if 'figs_checks' in st.session_state:
        col1, col2 = st.columns(2)
        col1.write('Property means by bin')
        col1.plotly_chart(st.session_state['figs_checks'][0])
        col2.write('Property means by cell')
        col2.plotly_chart(st.session_state['figs_checks'][1])
        col1.write('UMAP by bin')
        col1.plotly_chart(st.session_state['figs_checks'][2])
        col2.write('UMAP by cell')
        col2.plotly_chart(st.session_state['figs_checks'][4])
        col1.write('Spatial by bin')
        col1.plotly_chart(st.session_state['figs_checks'][3])
        col2.write('Spatial by cell')
        col2.plotly_chart(st.session_state['figs_checks'][5])


# Main script block
if __name__ == '__main__':

    # Set page settings
    page_name = 'Neighborhood Profiles Checker'
    st.set_page_config(layout='wide', page_title=page_name)
    st.title(page_name)

    # Call the main functino
    main()
