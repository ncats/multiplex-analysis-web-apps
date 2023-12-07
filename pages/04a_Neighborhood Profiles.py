'''
This is the python script which produces the NEIGHBORHOOD PROFILES PAGE
'''
import time
import pandas as pd
from datetime import datetime
import streamlit as st
from st_pages import show_pages_from_config, add_indentation
from streamlit_extras.add_vertical_space import add_vertical_space 
from streamlit_extras.app_logo import add_logo

# Import relevant libraries
import nidap_dashboard_lib as ndl   # Useful functions for dashboards connected to NIDAP
import basic_phenotyper_lib as bpl  # Useful functions for phenotyping collections of cells
import app_top_of_page as top

def init_spatial_umap():
    st.session_state.bc.startTimer()
    with st.spinner('Calculating Cell Counts and Areas'):
        st.session_state.spatial_umap = bpl.setup_Spatial_UMAP(st.session_state.df,
                                                               st.session_state.marker_multi_sel,
                                                               st.session_state.phenoOrder,
                                                               st.session_state.cpu_pool_size)
    st.write('Done Calculating Cell Counts and Areas')

    # Record time elapsed
    st.session_state.bc.set_value_df('time_to_run_counts', st.session_state.bc.elapsedTime())

    st.session_state.cell_counts_completed = True

def apply_umap(UMAPStyle):
    clust_minmax = [1, 40]
    st.session_state.bc.startTimer()
    with st.spinner('Calculating UMAP'):
        st.session_state.spatial_umap = bpl.perform_spatialUMAP(st.session_state.spatial_umap, UMAPStyle)
    st.write('Done Calculating Spatial UMAP')
    
    # Record time elapsed
    st.session_state.bc.printElapsedTime(msg = f'Performing UMAP')
    st.session_state.bc.set_value_df('time_to_run_UMAP', st.session_state.bc.elapsedTime())

    # List of possible UMAP Lineages as defined by the completed UMAP
    st.session_state.umapPheno = [st.session_state.defLineageOpt]
    st.session_state.umapPheno.extend(st.session_state.pheno_summ['phenotype'])
    st.session_state.umapMarks = [st.session_state.defLineageOpt]
    st.session_state.umapMarks.extend(st.session_state.spatial_umap.markers)
    st.session_state.umapMarks.extend(['Other'])

    # List of possible outcome variables as defined by the config yaml files
    st.session_state.umapOutcomes = [st.session_state.defumapOutcomes]
    st.session_state.umapOutcomes.extend(st.session_state.outcomes)
    st.session_state.inciOutcomes = [st.session_state.definciOutcomes]
    st.session_state.inciOutcomes.extend(st.session_state.outcomes)

    st.session_state.df_umap = st.session_state.spatial_umap.cells.loc[st.session_state.spatial_umap.cells['umap_test'], :]

    # Perform possible cluster variations with the completed UMAP
    st.session_state.bc.startTimer()
    with st.spinner('Calculating Possible Clusters'):
        st.session_state.clust_range, st.session_state.wcss = bpl.measure_possible_clust(st.session_state.spatial_umap, clust_minmax)
    st.session_state.bc.printElapsedTime(msg = f'Calculating posstible clusters')

    st.session_state.wcss_calc_completed = True
    st.session_state.umapCompleted = True

def set_clusters():
    st.session_state.bc.startTimer()
    st.session_state.spatial_umap = bpl.perform_clusteringUMAP(st.session_state.spatial_umap, st.session_state.slider_clus_val)
    st.session_state.selected_nClus = st.session_state.slider_clus_val
    st.write('Done Calculating Clusters')

    # Record time elapsed
    st.session_state.bc.printElapsedTime(msg = f'Setting Clusters')
    st.session_state.bc.set_value_df('time_to_run_cluster', st.session_state.bc.elapsedTime())

    st.session_state.clustering_completed = True

def slide_id_callback():
    # st.session_state['idxSlide ID'] = st.session_state['uniSlide ID_short'].index(st.session_state['selSlide ID_short'])
    idx =  st.session_state['idxSlide ID'] = st.session_state['uniSlide ID_short'].index(st.session_state['selSlide ID_short'])
    st.session_state['selSlide ID'] = st.session_state['uniSlide ID'][idx]

def main():
    '''
    Main function for running the page
    '''
    # Use the whole page width
    st.set_page_config(
        page_title="Neighborhood Profiles",
        layout="wide"
    )

    # Remove key values from session_state that should not persist
    for key, val in st.session_state.items():
        if (not key.endswith('__do_not_persist')) and (not key.startswith('FormSubmitter:')):
            st.session_state[key] = val

    # Apply pages order and indentation
    add_indentation()
    show_pages_from_config()

    # Sidebar organization
    with st.sidebar:
        st.write('**:book: [Documentation](https://ncats.github.io/multiplex-analysis-web-apps)**')

    # Add logo to page
    add_logo('app_images/mawa_logo-width315.png', height=150)

    # Run Top of Page (TOP) functions
    st.session_state = top.check_for_platform(st.session_state)

    if 'init' not in st.session_state:
        settings_yaml_file = 'config_files/OMAL_REEC.yml'
        # Initialize session_state values for streamlit processing
        st.session_state = ndl.init_session_state(st.session_state, settings_yaml_file)

    st.header('Neighborhood Profiles\nNCATS-NCI-DMAP')

    clust_minmax = [1, 40]
    neiProCols = st.columns([1, 1, 2])
    with neiProCols[0]:
        cellCountsButt = st.button('Perform Cell Counts/Areas Analysis')
        umapButt       = st.button('Perform UMAP Analysis')
        st.slider('Number of K-means clusters', min_value=clust_minmax[0], max_value=clust_minmax[1], key = 'slider_clus_val')
        clustButt      = st.button('Perform Clustering Analysis')

    with neiProCols[1]:
        if cellCountsButt:
            if st.session_state.phenotyping_completed:
                init_spatial_umap()
        if umapButt:
            if st.session_state.cell_counts_completed:
                apply_umap(UMAPStyle = 'Densities')
        if clustButt:
            if st.session_state.umapCompleted: 
                set_clusters()

    with neiProCols[2]:
        if st.session_state.umapCompleted:
            ### Clustering Meta Analysis and Description ###
            with st.expander('Cluster Meta-Analysis'):
                wcss_cols = st.columns(2)
                with wcss_cols[0]:
                    st.markdown('''The within-cluster sum of squares (WCSS) is a measure of the 
                                variability of the observations within each cluster. In general, 
                                a cluster that has a small sum of squares is more compact than a 
                                cluster that has a large sum of squares. Clusters that have higher 
                                values exhibit greater variability of the observations within the 
                                cluster.''')
                with wcss_cols[1]:
                    if st.session_state.umapCompleted:
                        elbowFig = bpl.draw_wcss_elbow_plot(st.session_state.clust_range, st.session_state.wcss, st.session_state.selected_nClus)
                        st.pyplot(elbowFig)

    if not st.session_state.phenotyping_completed:
        st.warning('Step 0: Please complete phentoyping analysis (See Phenotyping Page)', icon="⚠️")
    elif not st.session_state.cell_counts_completed:
        st.warning('Step 1: Please complete Cell Counts and Areas analysis', icon="⚠️")
    elif not st.session_state.umapCompleted:
        st.warning('Step 2: Please run UMAP analysis', icon="⚠️")
    elif not st.session_state.clustering_completed:
        st.warning('Step 3: Please run clustering analysis', icon="⚠️")
    else:
        add_vertical_space(2)

    ### Visualizations ###
    uScatCol, uNeighPCol = st.columns(2)

    # Scatterplot Figure Column
    with uScatCol:
        # Print a column header
        st.header('Clusters Plot')

        clustOPheno = st.radio('Plot Colors by: ',
                               ('Clusters', 'Phenotype'),
                               horizontal = True, index = 0, key = 'clustOPheno')
        
        imageProgCol = st.columns([3, 1, 1, 2])
        with imageProgCol[0]:
            st.selectbox('Slide ID',
                         (st.session_state['uniSlide ID_short']),
                         key = 'selSlide ID_short',
                         on_change=slide_id_callback)
        with imageProgCol[1]:
            add_vertical_space(2)
            # st.button('←', on_click=slide_id_prog_left_callback, disabled=st.session_state.prog_left_disabeled)
        with imageProgCol[2]:
            add_vertical_space(2)
            # st.button('→', on_click=slide_id_prog_right_callback, disabled=st.session_state.prog_right_disabeled)
        with imageProgCol[3]:
            add_vertical_space(2)
            st.write(f'Image {st.session_state["idxSlide ID"]+1} of {st.session_state["numSlide ID"]}')

        if st.session_state.umapCompleted:
            st.session_state.df_umap_filt = st.session_state.df_umap.loc[st.session_state.df_umap['Slide ID'] == st.session_state['selSlide ID'], :]
            st.session_state = ndl.setFigureObjs_UMAP(st.session_state)

            if clustOPheno == 'Clusters':
                st.pyplot(st.session_state.seabornFig_clust)
            else:
                st.pyplot(st.session_state.phenoFig)

            clust_fig_col = st.columns([2, 1])
            with clust_fig_col[0]:
                st.text_input('.png file suffix (Optional)', key = 'cluster_scatter_suffix')
            with clust_fig_col[1]:
                add_vertical_space(2)
                if st.button('Append Export List', key = 'appendexportbutton_clusterscatter__do_not_persist'):
                    ndl.save_png(st.session_state.seabornFig_clust, 'Cluster Scatterplot', st.session_state.cluster_scatter_suffix)
                    st.toast(f'Added {st.session_state.cluster_scatter_suffix} to export list')

    # Neighborhood Profiles Figure Column
    with uNeighPCol:
        st.header('Neighborhood Profiles')
        if 'spatial_umap' in st.session_state:
            selNeighFig = st.selectbox('Select a cluster to view', 
                                       list(range(st.session_state.selected_nClus)))
            if hasattr(st.session_state.spatial_umap, 'dens_df'):
                
                NeiProFig = bpl.neighProfileDraw(st.session_state.spatial_umap, selNeighFig)
                st.pyplot(fig=NeiProFig)

                neigh_prof_col = st.columns([2, 1])
                with neigh_prof_col[0]:
                    st.text_input('.png file suffix (Optional)', key = 'neigh_prof_line_suffix')
                with neigh_prof_col[1]:
                    add_vertical_space(2)
                    if st.button('Append Export List', key = 'appendexportbutton_neighproline__do_not_persist'):
                        
                        ndl.save_png(NeiProFig, 'Neighborhood Profiles', st.session_state.neigh_prof_line_suffix)
                        st.toast(f'Added {st.session_state.neigh_prof_line_suffix} to export list')

if __name__ == '__main__':
    main()
