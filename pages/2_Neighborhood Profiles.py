import streamlit as st
from streamlit_javascript import st_javascript
import time

# Import relevant libraries
import nidap_dashboard_lib as ndl   # Useful functions for dashboards connected to NIDAP
import basic_phenotyper_lib as bpl  # Useful functions for phenotyping collections of cells

def main():
    '''
    Main function for running the page
    '''
    # Use the whole page width
    st.set_page_config(
        page_title="Neighborhood Profiles",
        layout="wide"
    )

    for key, val in st.session_state.items():
        if (not key.endswith('__do_not_persist')) and (not key.startswith('FormSubmitter:')):
            st.session_state[key] = val

    st.header('Neighborhood Profiles\nNCATS-NCI-DMAP')

    ### SIDE BAR ORGANIZATION ###
    with st.sidebar:
        url = st_javascript("await fetch('').then(r => window.parent.location.href)")
        st.write(f'''[Open app in new Tab]({url})\n (MS Edge/ Google Chrome)''')

    clust_minmax = [1, 40]
    cellCountsCol, UMAPCol, ClustCol = st.columns(3)
    with cellCountsCol:
        with st.form('Cell Counts/Areas'):
            submitCountsCell = st.form_submit_button('Perform Cell Counts/Areas')
            if submitCountsCell:
                countTSt = time.time()
                st.session_state.spatial_umap = bpl.setup_Spatial_UMAP(st.session_state.df,
                                                                        st.session_state.marker_pre,
                                                                        st.session_state.phenoOrder,
                                                                        st.session_state.cpu_pool_size)
                st.write('Done Calculating Cell Counts/Areas')

                # Record time elapsed
                countElapsed = time.time() - countTSt
                st.session_state.bc.set_value_df('time_to_run_counts', countElapsed)
    with UMAPCol:
        with st.form('UMAP Settings'):
            UMAPStyle = 'Densities' #st.radio('Choose SpatialUMAP measure', ['Densities', 'Proportions'], horizontal=True)
            submitUMAP = st.form_submit_button('Perform UMAP')
            
            # UMAP Submit Button
            if submitUMAP:
                UMAPTSt = time.time()
                st.session_state.spatial_umap = bpl.perform_spatialUMAP(st.session_state.spatial_umap, UMAPStyle)
                st.write('Done Calculating Spatial UMAP')
                
                # Record time elapsed
                UMAPElapsed = time.time() - UMAPTSt
                st.session_state.bc.set_value_df('time_to_run_UMAP', UMAPElapsed)

                # List of possible UMAP Lineages as defined by the completed UMAP
                st.session_state.umapPheno = [st.session_state.defLineageOpt]
                st.session_state.umapPheno.extend(st.session_state.spec_summ['phenotype'])
                st.session_state.umapMarks = [st.session_state.defLineageOpt]
                st.session_state.umapMarks.extend(st.session_state.spatial_umap.markers)
                st.session_state.umapMarks.extend(['Other'])

                # List of possible outcome variables as defined by the config yaml files
                st.session_state.umapOutcomes = [st.session_state.defumapOutcomes]
                st.session_state.umapOutcomes.extend(st.session_state.outcomes)
                st.session_state.inciOutcomes = [st.session_state.definciOutcomes]
                st.session_state.inciOutcomes.extend(st.session_state.outcomes)

                # Perform possible cluster variations with the completed UMAP
                st.session_state.clust_range, st.session_state.wcss = bpl.measure_possible_clust(st.session_state.spatial_umap, clust_minmax)
                st.session_state.umapCompleted = True
    with ClustCol:
        with st.form('Clustering Settings'):
            st.slider('Number of K-means clusters', min_value=clust_minmax[0], max_value=clust_minmax[1], key = 'slider_clus_val')
            submitCluster = st.form_submit_button('Perform Clustering')
            if submitCluster:
                ClustTSt = time.time()
                st.session_state.spatial_umap = bpl.perform_clusteringUMAP(st.session_state.spatial_umap, st.session_state.slider_clus_val)
                st.session_state.selected_nClus = st.session_state.slider_clus_val
                st.write('Done Calculating Clusters')

                # Record time elapsed
                ClustElapsed = time.time() - ClustTSt
                st.session_state.bc.set_value_df('time_to_run_cluster', ClustElapsed)

    ### Clustering Meta Analysis and Description ###
    with st.expander('Cluster Meta-Analysis'):
        wcssCol1, wcssCol2 = st.columns(2)
        with wcssCol1:
            st.markdown('''The within-cluster sum of squares (WCSS) is a measure of the 
                        variability of the observations within each cluster. In general, 
                        a cluster that has a small sum of squares is more compact than a 
                        cluster that has a large sum of squares. Clusters that have higher 
                        values exhibit greater variability of the observations within the 
                        cluster.''')
        with wcssCol2:
            if st.session_state.umapCompleted:
                elbowFig = bpl.draw_wcss_elbow_plot(st.session_state.clust_range, st.session_state.wcss, st.session_state.selected_nClus)
                st.pyplot(elbowFig)

    ### Visualizations ###
    uScatCol, uNeighPCol = st.columns(2)

    # Scatterplot Figure Column
    with uScatCol:
        # Print a column header
        st.header('Clusters Plot')

        if st.session_state.umapCompleted:
            st.session_state = ndl.setFigureObjs_UMAP(st.session_state)

            visOpCol1 , visOpCol2 = st.columns(2)
            with visOpCol1:
                clustOPheno = st.radio(
                        'Plot Colors by: ',
                            ('Clusters', 'Phenotype'),
                            horizontal = True, index = 0, key = 'clustOPheno')

            if clustOPheno == 'Clusters':
                st.pyplot(st.session_state.seabornFig_clust)
            else:
                st.pyplot(st.session_state.phenoFig)

    # Neighborhood Profiles Figure Column
    with uNeighPCol:
        st.header('Neighborhood Profiles')
        if 'spatial_umap' in st.session_state:
            if hasattr(st.session_state.spatial_umap, "densMeansDict"):
                selNeighFig = st.selectbox('Select a cluster to view', list(range(st.session_state.selected_nClus)))

                NeiProFig = bpl.neighProfileDraw(st.session_state.spatial_umap, selNeighFig)
                st.pyplot(fig=NeiProFig)

                with st.form('Exporting Neighborhood Profiles'):
                    selectProjOutPNG_U = st.selectbox(
                                'Select your Unstructured Dataset',
                                (st.session_state.OutputPNGPaths),
                                key = 'selectProjOutPNG_U')

                    st.session_state.NeighborProPng = st.text_input('Image File name', st.session_state.NeighborProPng, key = 'NeighborProPng_FileName')

                    submitSaveNeig = st.form_submit_button('Save Neighborhood Profiles')
                    if submitSaveNeig:
                        ndl.save_png_dataset(st.session_state.fiol, selectProjOutPNG_U, st.session_state.NeighborProPng, NeiProFig)
                        st.write('Export Successful')


if __name__ == '__main__':
    main()
