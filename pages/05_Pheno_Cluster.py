# Import relevant libraries
import streamlit as st
import app_top_of_page as top
import streamlit_dataframe_editor as sde
import streamlit as st 
import pandas as pd
import anndata as ad
import scanpy as sc
import seaborn as sns
import os
import matplotlib.pyplot as plt
import phenograph
import parc
from utag import utag
import numpy as np
import scanpy.external as sce
import plotly.express as px

# make adata from unifier table
# hard coded column names, nned to think about a better way
def phenocluster__make_adata(df):
    colNames = list(df)
    object_index = colNames.index("Object")
    area_index = colNames.index("area")
    extracted_elements = colNames[object_index + 1 : area_index]
    #start_index = colNames.index("area") - 1 
    #end_index = len(colNames)
    end_index = colNames.index("Image")
    metaCols =  colNames[0 : end_index]
    mat = df[extracted_elements]
    meta = df[metaCols]
    #st.write(list(meta))
    #st.write(list(mat))
    adata = ad.AnnData(mat)
    adata.obs = meta
    adata.layers["counts"] = adata.X.copy()
    return adata

# scanpy clustering
def RunNeighbClust(adata, n_neighbors):
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=None)
    sc.tl.leiden(adata,resolution=1)
    adata.obs['Cluster'] = adata.obs['leiden']
    sc.tl.umap(adata)
    adata.obsm['spatial'] = np.array(adata.obs[["Centroid X (µm)_(standardized)", "Centroid Y (µm)_(standardized)"]])
    return adata

# phenograph clustering
def RunPhenographClust(adata, n_neighbors, k):
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=None)
    communities, graph, Q = phenograph.cluster(adata.X, clustering_algo = 'leiden', k=k)
    adata.obs['Cluster'] = communities
    adata.obs['Cluster'] = adata.obs['Cluster'].astype(str)
    sc.tl.umap(adata)
    adata.obsm['spatial'] = np.array(adata.obs[["Centroid X (µm)_(standardized)", "Centroid Y (µm)_(standardized)"]])
    return adata

# parc clustering
def run_parc_clust(adata, n_neighbors):
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=None)
    parc_results = parc.PARC(adata.X)
    parc_results.run_PARC()
    adata.obs['Cluster'] = parc_results.labels
    adata.obs['Cluster'] = adata.obs['Cluster'].astype(str)
    sc.tl.umap(adata)
    adata.obsm['spatial'] = np.array(adata.obs[["Centroid X (µm)_(standardized)", "Centroid Y (µm)_(standardized)"]])
    return adata

# utag clustering
# need to make image selection based on the variable
def run_utag_clust(adata, n_neighbors, resolutions):
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=None)
    sc.tl.umap(adata)
    adata.obsm['spatial'] = np.array(adata.obs[["Centroid X (µm)_(standardized)", "Centroid Y (µm)_(standardized)"]])
    utag_results = utag(adata,
        slide_key="Image ID_(standardized)",
        max_dist=15,
        normalization_mode='l1_norm',
        apply_clustering=True,
        clustering_method = 'leiden', 
        resolutions = resolutions
    )
    curClusterCol = 'UTAG Label_leiden_'  + str(resolutions[0])
    utag_results.obs['Cluster'] = utag_results.obs[curClusterCol]
    adata.obs['Cluster'] = utag_results.obs[curClusterCol]
    return utag_results

# plot umaps
def phenocluster__plotly_umaps(adata, umap_cur_col, umap_cur_groups, umap_color_col):
    with phenocluster__col2:
        subcol1, subcol2 = st.columns(2)
        for i, umap_cur_group in enumerate(umap_cur_groups):
            if umap_cur_group == "All":
                subDat = adata
            else:
                subDat = adata[adata.obs[umap_cur_col] == umap_cur_group]
            umap_coords = subDat.obsm['X_umap']
            df = pd.DataFrame(umap_coords, columns=['UMAP_1', 'UMAP_2'])
            clustersList = list(subDat.obs[umap_color_col] )
            df[umap_color_col] = clustersList
            df[umap_color_col] = df[umap_color_col].astype(str)
            # Create the seaborn plot
            fig = px.scatter(df, 
             x="UMAP_1", 
             y="UMAP_2", 
             color=umap_color_col, 
             title="UMAP " + umap_cur_group
             #color_discrete_sequence=px.colors.sequential.Plasma
             )
            fig.update_traces(marker=dict(size=3)) # Adjust the size of the dots
            fig.update_layout(
                title=dict(
                    text="UMAP " + umap_cur_group,
                    x=0.5, # Center the title
                    xanchor='center',
                    yanchor='top'
                ),
                legend=dict(
                    orientation="h",
                    yanchor="top",
                    y=-0.2,
                    xanchor="right",
                    x=1
                )
            )
            if i % 2 == 0:
                subcol1.plotly_chart(fig, use_container_width=True)
            else:
                subcol2.plotly_chart(fig, use_container_width=True)

# plot spatial
def spatial_plots_cust_2(adata, umap_cur_col, umap_cur_groups, umap_color_col):
    with phenocluster__col3:
        subcol3, subcol4 = st.columns(2)
        for i, umap_cur_group in enumerate(umap_cur_groups):
            if umap_cur_group == "All":
                subDat = adata
            else:
                subDat = adata[adata.obs[umap_cur_col] == umap_cur_group]
            umap_coords = subDat.obs[['Centroid X (µm)_(standardized)', 'Centroid Y (µm)_(standardized)']]
            df = pd.DataFrame(umap_coords).reset_index().drop('index', axis = 1)
            clustersList = list(subDat.obs[umap_color_col] )
            df[umap_color_col] = clustersList
            df[umap_color_col] = df[umap_color_col].astype(str)
            fig = px.scatter(df, 
             x="Centroid X (µm)_(standardized)", 
             y="Centroid Y (µm)_(standardized)", 
             color=umap_color_col, 
             title="Spatial " + umap_cur_group
             #color_discrete_sequence=px.colors.sequential.Plasma
             )
            fig.update_traces(marker=dict(size=3)) # Adjust the size of the dots
            fig.update_layout(
                title=dict(
                    text="Spatial " + umap_cur_group,
                    x=0.5, # Center the title
                    xanchor='center',
                    yanchor='top'
                ),
                legend=dict(
                    orientation="h",
                    yanchor="top",
                    y=-0.2,
                    xanchor="right",
                    x=1
                )
            )
            if i % 2 == 0:
                subcol3.plotly_chart(fig, use_container_width=True)
            else:
                subcol4.plotly_chart(fig, use_container_width=True)

# make Unaps and Spatial Plots
def make_all_plots():
    # make umaps
        phenocluster__plotly_umaps(st.session_state['phenocluster__clustering_adata'], 
            st.session_state['phenocluster__umap_cur_col'], 
            st.session_state['phenocluster__umap_cur_groups'],
            st.session_state['phenocluster__umap_color_col'])
    # make spatial plots
        spatial_plots_cust_2(st.session_state['phenocluster__clustering_adata'], 
            st.session_state['phenocluster__umap_cur_col'], 
            st.session_state['phenocluster__umap_cur_groups'],
            st.session_state['phenocluster__umap_color_col'])
    

    # default session state values
def phenocluster__default_session_state():
    if 'phenocluster__cluster_method' not in st.session_state:
        st.session_state['phenocluster__cluster_method'] = "phenograph"

    if 'phenocluster__n_neighbors_state' not in st.session_state:
        st.session_state['phenocluster__n_neighbors_state'] = 10

    if 'phenocluster__phenograph_k' not in st.session_state:
        st.session_state['phenocluster__phenograph_k'] = 30

    if 'phenocluster__umap_cur_col' not in st.session_state:
        st.session_state['phenocluster__umap_cur_col'] = "Image"

    if 'phenocluster__umap_cur_groups' not in st.session_state:
        st.session_state['phenocluster__umap_cur_groups'] = ["All"]
    
# main
def main():
    """
    Main function for the page.
    """
    #st.write(st.session_state['unifier__df'].head())
    
    # set default values
    phenocluster__default_session_state()
    
    # make layout with columns
    adata = phenocluster__make_adata(st.session_state['unifier__df'])
    
    # options
    with phenocluster__col1:
        clusteringMethods = ['phenograph', 'neighb', 'parc', 'utag']
        selected_clusteringMethod = st.selectbox('Select Clustering method:', clusteringMethods, 
                                                key='clusteringMethods_dropdown') 

        # Update session state on every change
        st.session_state['phenocluster__cluster_method'] = selected_clusteringMethod

        # default widgets
        st.session_state['phenocluster__n_neighbors_state']  = st.number_input(label = "Input number of neighbors", 
                                value=st.session_state['phenocluster__n_neighbors_state'])
        
        if st.session_state['phenocluster__cluster_method'] == "phenograph":
            st.session_state['phenocluster__phenograph_k'] = st.number_input(label = "Phenograph k", 
                                value=st.session_state['phenocluster__phenograph_k'])
        
        # add options if clustering has been run
        if st.button('Run Clustering'):
            if st.session_state['phenocluster__cluster_method'] == "phenograph":
                with st.spinner('Wait for it...'):
                    st.session_state['phenocluster__clustering_adata'] = RunPhenographClust(adata=adata, n_neighbors=st.session_state['phenocluster__n_neighbors_state'], 
                                                                                            k=st.session_state['phenocluster__phenograph_k'])
            elif st.session_state['phenocluster__cluster_method'] == "neighb":
                with st.spinner('Wait for it...'):
                    st.session_state['phenocluster__clustering_adata'] = RunNeighbClust(adata=adata, 
                                                                                        n_neighbors=st.session_state['phenocluster__n_neighbors_state'])
            #st.session_state['phenocluster__clustering_adata'] = adata
            elif st.session_state['phenocluster__cluster_method'] == "parc":
                with st.spinner('Wait for it...'):
                    st.session_state['phenocluster__clustering_adata'] = run_parc_clust(adata=adata, 
                                                                                        n_neighbors=st.session_state['phenocluster__n_neighbors_state'])
            elif st.session_state['phenocluster__cluster_method'] == "utag":
                with st.spinner('Wait for it...'):
                    st.session_state['phenocluster__clustering_adata'] = run_utag_clust(adata=adata, 
                                                                                        n_neighbors=st.session_state['phenocluster__n_neighbors_state'], resolutions=[1])
            # save clustering result
            st.session_state['phenocluster__clustering_adata'].write("input/clust_dat.h5ad")
                    
            # set default values for umap color column
            if 'phenocluster__umap_color_col' not in st.session_state:
                st.session_state['phenocluster__umap_color_col'] = "Cluster"
                
            
            # umap
        if 'phenocluster__clustering_adata' in st.session_state:
            #st.write(pd.unique(st.session_state['phenocluster__clustering_adata'].obs["Cluster"]))
            
            st.session_state['phenocluster__umeta_columns'] = list(st.session_state['phenocluster__clustering_adata'].obs.columns)
            st.session_state['phenocluster__umap_color_col_index'] = st.session_state['phenocluster__umeta_columns'].index(st.session_state['phenocluster__umap_color_col'])
            #st.write(st.session_state['phenocluster__umap_color_col_index'])
            
            # select column for umap coloring
            st.session_state['phenocluster__umap_color_col'] = st.selectbox('Select column for UMAP coloring:', 
                                                                    st.session_state['phenocluster__umeta_columns'],
                                                                    index=st.session_state['phenocluster__umap_color_col_index']
                                                                    )
            
            # select column for umap subsetting
            st.session_state['phenocluster__umap_cur_col'] = st.selectbox('Select column for UMAP subsetting:', 
                                                                    st.session_state['phenocluster__umeta_columns'], key='phenocluster__umap_col_dropdown_subset'
                                                                    )
            
            # list of available subsetting options
            umap_cur_groups=  ["All"] + list(pd.unique(st.session_state['phenocluster__clustering_adata'].obs[st.session_state['phenocluster__umap_cur_col']]))
            umap_sel_groups = st.multiselect('Select groups to be plotted on UMAP',
                                                                            options = umap_cur_groups)
            st.session_state['phenocluster__umap_cur_groups'] = umap_sel_groups
            
            st.button('Make Plots' , on_click=make_all_plots)
            
            
    

# Run the main function
if __name__ == '__main__':

    # Set page settings
    page_name = 'Phenotype Clustering Andrei'
    st.set_page_config(layout='wide', page_title=page_name)
    st.title(page_name)

    phenocluster__col1, phenocluster__col2, phenocluster__col3 = st.columns([1, 3, 3])
    
    # Run streamlit-dataframe-editor library initialization tasks at the top of the page
    st.session_state = sde.initialize_session_state(st.session_state)

    # Run Top of Page (TOP) functions
    st.session_state = top.top_of_page_reqs(st.session_state)

    # Call the main function
    main()

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)