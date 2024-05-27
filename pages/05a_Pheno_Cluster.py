# Import relevant libraries
import streamlit as st
import subprocess

try:
    import parc
    print("Successfully imported parc!")
except ImportError:
    print("Installing parc...")
    subprocess.run("pip install parc", shell=True)
    try:
        import parc
    except ImportError:
       print("Failed to import parc.")
        
from ast import arg
from pyparsing import col
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
#from utag import utag
import numpy as np
import scanpy.external as sce
import plotly.express as px
import time
from utag.segmentation import utag

def phenocluster__make_adata(df, x_cols, meta_cols):
    mat = df[x_cols]
    meta = df[meta_cols]
    adata = ad.AnnData(mat)
    adata.obs = meta
    adata.layers["counts"] = adata.X.copy()
    #adata.write("input/clust_dat.h5ad")
    return adata

# scanpy clustering
def RunNeighbClust(adata, n_neighbors, metric, resolution, random_state, n_principal_components):
    if n_principal_components > 0:
        sc.pp.pca(adata, n_comps=n_principal_components)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, metric=metric, n_pcs=n_principal_components)
    sc.tl.leiden(adata,resolution=resolution, random_state=random_state, n_iterations=5, flavor="igraph")
    adata.obs['Cluster'] = adata.obs['leiden']
    #sc.tl.umap(adata)
    adata.obsm['spatial'] = np.array(adata.obs[["Centroid X (µm)_(standardized)", "Centroid Y (µm)_(standardized)"]])
    return adata

# phenograph clustering
def RunPhenographClust(adata, n_neighbors, clustering_algo, min_cluster_size, 
                       primary_metric, resolution_parameter, nn_method, random_seed, n_principal_components):
    #sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=0)
    if n_principal_components == 0:
        communities, graph, Q = phenograph.cluster(adata.X, clustering_algo=clustering_algo, k=n_neighbors, 
                                                min_cluster_size=min_cluster_size, primary_metric=primary_metric, 
                                                resolution_parameter=resolution_parameter, nn_method=nn_method,
                                                seed=random_seed, n_iterations=5)
    else:
        sc.pp.pca(adata, n_comps=n_principal_components)
        communities, graph, Q = phenograph.cluster(adata.obsm['X_pca'], clustering_algo=clustering_algo, k=n_neighbors, 
                                                min_cluster_size=min_cluster_size, primary_metric=primary_metric, 
                                                resolution_parameter=resolution_parameter, nn_method=nn_method,
                                                seed=random_seed, n_iterations=5)
    adata.obs['Cluster'] = communities
    adata.obs['Cluster'] = adata.obs['Cluster'].astype(str)
    #sc.tl.umap(adata)
    adata.obsm['spatial'] = np.array(adata.obs[["Centroid X (µm)_(standardized)", "Centroid Y (µm)_(standardized)"]])
    return adata

# parc clustering
def run_parc_clust(adata, n_neighbors, dist_std_local, jac_std_global, small_pop, 
                   random_seed, resolution_parameter, hnsw_param_ef_construction, n_principal_components):
    #sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=0)
    if n_principal_components == 0:
        parc_results = parc.PARC(adata.X, dist_std_local=dist_std_local, jac_std_global=jac_std_global, 
                                small_pop=small_pop, random_seed=random_seed, knn=n_neighbors,
                                resolution_parameter=resolution_parameter, 
                                hnsw_param_ef_construction=hnsw_param_ef_construction,
                                partition_type="RBConfigurationVP",
                                n_iter_leiden=5)
    else:
        sc.pp.pca(adata, n_comps=n_principal_components)
        parc_results = parc.PARC(adata.obsm['X_pca'], dist_std_local=dist_std_local, jac_std_global=jac_std_global, 
                                small_pop=small_pop, random_seed=random_seed, knn=n_neighbors,
                                resolution_parameter=resolution_parameter, 
                                hnsw_param_ef_construction=hnsw_param_ef_construction,
                                partition_type="RBConfigurationVP",
                                n_iter_leiden=5)
    parc_results.run_PARC()
    adata.obs['Cluster'] = parc_results.labels
    adata.obs['Cluster'] = adata.obs['Cluster'].astype(str)
    #sc.tl.umap(adata)
    adata.obsm['spatial'] = np.array(adata.obs[["Centroid X (µm)_(standardized)", "Centroid Y (µm)_(standardized)"]])
    return adata

# utag clustering
# need to make image selection based on the variable
def run_utag_clust(adata, n_neighbors, resolutions, clustering_method, max_dist, n_principal_components):
    #sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=0)
    #sc.tl.umap(adata)
    adata.obsm['spatial'] = np.array(adata.obs[["Centroid X (µm)_(standardized)", "Centroid Y (µm)_(standardized)"]])
    utag_results = utag(adata,
        slide_key="Image ID_(standardized)",
        max_dist=max_dist,
        normalization_mode='l1_norm',
        apply_clustering=True,
        clustering_method = clustering_method, 
        resolutions = resolutions,
        leiden_kwargs={"n_iterations": 5, "random_state": 42},
        pca_kwargs = {"n_comps": n_principal_components}
    )

                                
    curClusterCol = 'UTAG Label_leiden_'  + str(resolutions[0])
    utag_results.obs['Cluster'] = utag_results.obs[curClusterCol]
    adata.obs['Cluster'] = utag_results.obs[curClusterCol]
        
    curClusterCol = 'UTAG Label_leiden_'  + str(resolutions[0])
    utag_results.obs['Cluster'] = utag_results.obs[curClusterCol]
    adata.obs['Cluster'] = utag_results.obs[curClusterCol]
    return utag_results

def phenocluster__scanpy_umap(adata, n_neighbors, metric, n_principal_components):
    if n_principal_components > 0:
        sc.pp.pca(adata, n_comps=n_principal_components)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, metric=metric, n_pcs=n_principal_components)
    sc.tl.umap(adata)
    st.session_state['phenocluster__clustering_adata'] = adata

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
    with phenocluster__col2:
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

# make Umaps and Spatial Plots
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
    
    if 'phenocluster__subset_data' not in st.session_state:
        st.session_state['phenocluster__subset_data'] = False
    
    if 'phenocluster__cluster_method' not in st.session_state:
        st.session_state['phenocluster__cluster_method'] = "phenograph"
        
    if 'phenocluster__resolution' not in st.session_state:
        st.session_state['phenocluster__resolution'] = 1.0

    # phenograph options
    if 'phenocluster__n_neighbors_state' not in st.session_state:
        st.session_state['phenocluster__n_neighbors_state'] = 10

    if 'phenocluster__phenograph_k' not in st.session_state:
        st.session_state['phenocluster__phenograph_k'] = 30
    
    if 'phenocluster__phenograph_clustering_algo' not in st.session_state:
        st.session_state['phenocluster__phenograph_clustering_algo'] = 'louvain'
    
    if 'phenocluster__phenograph_min_cluster_size' not in st.session_state:
        st.session_state['phenocluster__phenograph_min_cluster_size'] = 10
    
    if 'phenocluster__metric' not in st.session_state:
        st.session_state['phenocluster__metric'] = 'euclidean'
        
    if 'phenocluster__phenograph_nn_method' not in st.session_state:
        st.session_state['phenocluster__phenograph_nn_method'] = 'kdtree'
        
    if 'phenocluster__n_principal_components' not in st.session_state:
        st.session_state['phenocluster__n_principal_components'] = 10
        
    # parc options
    # dist_std_local, jac_std_global, small_pop, random_seed, resolution_parameter, hnsw_param_ef_construction
    if 'phenocluster__parc_dist_std_local' not in st.session_state:
        st.session_state['phenocluster__parc_dist_std_local'] = 3
    
    if 'phenocluster__parc_jac_std_global' not in st.session_state:
        st.session_state['phenocluster__parc_jac_std_global'] = 0.15
        
    if 'phenocluster__parc_small_pop' not in st.session_state:
        st.session_state['phenocluster__parc_small_pop'] = 50
        
    if 'phenocluster__random_seed' not in st.session_state:
        st.session_state['phenocluster__random_seed'] = 42
        
    if 'phenocluster__hnsw_param_ef_construction' not in st.session_state:
        st.session_state['phenocluster__hnsw_param_ef_construction'] = 150
        
    # utag options
    #clustering_method ["leiden", "parc"]; resolutions; max_dist = 20
    if 'phenocluster__utag_clustering_method' not in st.session_state:
        st.session_state['phenocluster__utag_clustering_method'] = 'leiden'
        
    if 'phenocluster__utag_max_dist' not in st.session_state:
        st.session_state['phenocluster__utag_max_dist'] = 20

    # umap options
    #if 'phenocluster__umap_cur_col' not in st.session_state:
        #st.session_state['phenocluster__umap_cur_col'] = "Image"
        
    if 'phenocluster__umap_color_col' not in st.session_state:
        st.session_state['phenocluster__umap_color_col'] = "Cluster"

    if 'phenocluster__umap_cur_groups' not in st.session_state:
        st.session_state['phenocluster__umap_cur_groups'] = ["All"]
    
    # differential intensity options    
    if 'phenocluster__de_col' not in st.session_state:
        st.session_state['phenocluster__de_col'] = "Cluster"
    
    if 'phenocluster__de_sel_group' not in st.session_state:
        st.session_state['phenocluster__de_sel_groups'] = ["All"]
        
    if 'phenocluster__plot_diff_intensity_method' not in st.session_state:
        st.session_state['phenocluster__plot_diff_intensity_method'] = "Rank Plot"
        
    if 'phenocluster__plot_diff_intensity_n_genes' not in st.session_state:
        st.session_state['phenocluster__plot_diff_intensity_n_genes'] = 10
    
  
# subset data set
def phenocluster__subset_data(adata, subset_col, subset_vals):
    adata_subset = adata[adata.obs[subset_col].isin(subset_vals)]
    st.session_state['phenocluster__clustering_adata'] = adata_subset

# clusters differential expression
def phenocluster__diff_expr(adata, phenocluster__de_col, phenocluster__de_sel_groups):
    sc.tl.rank_genes_groups(adata, groupby = phenocluster__de_col, method="wilcoxon", layer="counts")
    
    if "All" in phenocluster__de_sel_groups:
        phenocluster__de_results = sc.get.rank_genes_groups_df(adata, group=None)
    else:
        phenocluster__de_results = sc.get.rank_genes_groups_df(adata, group=phenocluster__de_sel_groups)
    with phenocluster__col2:
        st.dataframe(phenocluster__de_results, use_container_width=True)
    

# main
def main():
    """
    Main function for the page.
    """
    #st.write(st.session_state['unifier__df'].head())
    
    # set default values
    phenocluster__default_session_state()
    
    
    
    # make layout with columns    
    # options
    
    with phenocluster__col_0[0]:
        
        numeric_cols = st.multiselect('Select numeric columns for clustering:', options = st.session_state['input_dataset'].data.columns, 
                    key='phenocluster__X_cols')
        meta_columns = st.multiselect('Select columns for metadata:', options = st.session_state['input_dataset'].data.columns, 
                key='phenocluster__meta_cols')
        
        items_to_add = ['Centroid X (µm)_(standardized)', 'Centroid Y (µm)_(standardized)']
        
        #meta_columns = st.session_state['phenocluster__meta_cols']
        # Add the new items if they don't already exist in the list
        for item in items_to_add:
            if item not in meta_columns:
                meta_columns.append(item)
        
        if st.button('Submuit columns'):
            st.session_state['phenocluster__clustering_adata'] = phenocluster__make_adata(st.session_state['input_dataset'].data, 
                                            numeric_cols,
                                            meta_columns)
    
    if 'phenocluster__clustering_adata' in st.session_state:
    
        with phenocluster__col1:
            
            # subset data
            st.checkbox('Subset Data:', key='phenocluster__subset_data')
            if st.session_state['phenocluster__subset_data'] == True:
                st.session_state['phenocluster__subset_options'] = list(st.session_state['phenocluster__clustering_adata'].obs.columns)
                phenocluster__subset_col = st.selectbox('Select column for subsetting:', st.session_state['phenocluster__subset_options'])
                st.session_state["phenocluster__subset_col"] = phenocluster__subset_col 
                st.session_state['phenocluster__subset_values_options'] = list(pd.unique(st.session_state['phenocluster__clustering_adata'].obs[st.session_state["phenocluster__subset_col"]]))
                phenocluster__subset_vals = st.multiselect('Select value for subsetting:', options = st.session_state['phenocluster__subset_values_options'], key='phenocluster__subset_vals_1')
                st.session_state["phenocluster__subset_vals"] = phenocluster__subset_vals 
                if st.button('Subset Data'):
                    phenocluster__subset_data(st.session_state['phenocluster__clustering_adata'],
                                            st.session_state["phenocluster__subset_col"],
                                            st.session_state["phenocluster__subset_vals"])
                
            # st.button('Subset Data' , on_click=phenocluster__subset_data, args = [st.session_state['phenocluster__clustering_adata'],
            #                                                                         st.session_state["phenocluster__subset_col"],
            #                                                                         st.session_state["phenocluster__subset_vals"]
                                                                                    
            # ]
            # )
                            
            clusteringMethods = ['phenograph', 'scanpy', 'parc', 'utag']
            selected_clusteringMethod = st.selectbox('Select Clustering method:', clusteringMethods, 
                                                    key='clusteringMethods_dropdown') 

            # Update session state on every change
            st.session_state['phenocluster__cluster_method'] = selected_clusteringMethod

            # default widgets
            
            st.number_input(label = "Number of Principal Components", key='phenocluster__n_principal_components', step = 1)
            
            st.session_state['phenocluster__n_neighbors_state']  = st.number_input(label = "K Nearest Neighbors", 
                                    value=st.session_state['phenocluster__n_neighbors_state'])
            
            st.number_input(label = "Clustering resolution", key='phenocluster__resolution')
            
            if st.session_state['phenocluster__cluster_method'] == "phenograph":
                # st.session_state['phenocluster__phenograph_k'] = st.number_input(label = "Phenograph k", 
                #                     value=st.session_state['phenocluster__phenograph_k'])
                st.selectbox('Phenograph clustering algorithm:', ['louvain', 'leiden'], key='phenocluster__phenograph_clustering_algo')
                st.number_input(label = "Phenograph min cluster size", key='phenocluster__phenograph_min_cluster_size', step = 1)
                st.selectbox('Distance metric:', ['euclidean', 'manhattan', 'correlation', 'cosine'], key='phenocluster__metric')
                st.selectbox('Phenograph nn method:', ['kdtree', 'brute'], key='phenocluster__phenograph_nn_method')
            
            elif st.session_state['phenocluster__cluster_method'] == "scanpy":
                st.selectbox('Distance metric:', ['euclidean', 'manhattan', 'correlation', 'cosine'], key='phenocluster__metric')
            
            elif st.session_state['phenocluster__cluster_method'] == "parc":
                # make parc specific widgets
                st.number_input(label = "Parc dist std local", key='phenocluster__parc_dist_std_local', step = 1)
                st.number_input(label = "Parc jac std global", key='phenocluster__parc_jac_std_global', step = 0.01)
                st.number_input(label = "Minimum cluster size to be considered a separate population",
                                key='phenocluster__parc_small_pop', step = 1)
                st.number_input(label = "Random seed", key='phenocluster__random_seed', step = 1)
                st.number_input(label = "HNSW exploration factor for construction", 
                                key='phenocluster__hnsw_param_ef_construction', step = 1)
            elif st.session_state['phenocluster__cluster_method'] == "utag":
                # make utag specific widgets
                #st.selectbox('UTAG clustering method:', ['leiden', 'parc'], key='phenocluster__utag_clustering_method')
                st.number_input(label = "UTAG max dist", key='phenocluster__utag_max_dist', step = 1)
            
            # add options if clustering has been run
            if st.button('Run Clustering'):
                start_time = time.time()
                if st.session_state['phenocluster__cluster_method'] == "phenograph":
                    with st.spinner('Wait for it...'):
                        st.session_state['phenocluster__clustering_adata'] = RunPhenographClust(adata=st.session_state['phenocluster__clustering_adata'], n_neighbors=st.session_state['phenocluster__n_neighbors_state'],
                                                                                                clustering_algo=st.session_state['phenocluster__phenograph_clustering_algo'],
                                                                                                min_cluster_size=st.session_state['phenocluster__phenograph_min_cluster_size'],
                                                                                                primary_metric=st.session_state['phenocluster__metric'],
                                                                                                resolution_parameter=st.session_state['phenocluster__resolution'],
                                                                                                nn_method=st.session_state['phenocluster__phenograph_nn_method'],
                                                                                                random_seed=st.session_state['phenocluster__random_seed'],
                                                                                                n_principal_components=st.session_state['phenocluster__n_principal_components']
                                                                                                )
                elif st.session_state['phenocluster__cluster_method'] == "scanpy":
                    with st.spinner('Wait for it...'):
                        st.session_state['phenocluster__clustering_adata'] = RunNeighbClust(adata=st.session_state['phenocluster__clustering_adata'], 
                                                                                            n_neighbors=st.session_state['phenocluster__n_neighbors_state'],
                                                                                            metric=st.session_state['phenocluster__metric'],
                                                                                            resolution=st.session_state['phenocluster__resolution'],
                                                                                            random_state=st.session_state['phenocluster__random_seed'],
                                                                                            n_principal_components=st.session_state['phenocluster__n_principal_components']
                                                                                            )
                #st.session_state['phenocluster__clustering_adata'] = adata
                elif st.session_state['phenocluster__cluster_method'] == "parc":
                    with st.spinner('Wait for it...'):                  
                        st.session_state['phenocluster__clustering_adata'] = run_parc_clust(adata=st.session_state['phenocluster__clustering_adata'], 
                                                                                            n_neighbors=st.session_state['phenocluster__n_neighbors_state'],
                                                                                            dist_std_local=st.session_state['phenocluster__parc_dist_std_local'],
                                                                                            jac_std_global= st.session_state['phenocluster__parc_jac_std_global'],
                                                                                            small_pop=st.session_state['phenocluster__parc_small_pop'],
                                                                                            random_seed=st.session_state['phenocluster__random_seed'],
                                                                                            resolution_parameter=st.session_state['phenocluster__resolution'],
                                                                                            hnsw_param_ef_construction=st.session_state['phenocluster__hnsw_param_ef_construction'],
                                                                                            n_principal_components=st.session_state['phenocluster__n_principal_components']
                                                                                            )
                elif st.session_state['phenocluster__cluster_method'] == "utag":
                    phenocluster__utag_resolutions = [st.session_state['phenocluster__resolution']]
                    with st.spinner('Wait for it...'):
                        st.session_state['phenocluster__clustering_adata'] = run_utag_clust(adata=st.session_state['phenocluster__clustering_adata'], 
                                                                                            n_neighbors=st.session_state['phenocluster__n_neighbors_state'], 
                                                                                            resolutions=phenocluster__utag_resolutions,
                                                                                            clustering_method=st.session_state['phenocluster__utag_clustering_method'],
                                                                                            max_dist=st.session_state['phenocluster__utag_max_dist'],
                                                                                            n_principal_components=st.session_state['phenocluster__n_principal_components']
                                                                                            )
                # save clustering result
                #st.session_state['phenocluster__clustering_adata'].write("input/clust_dat.h5ad")
                end_time = time.time()
                execution_time = end_time - start_time
                rounded_time = round(execution_time, 2)
                st.write('Execution time: ', rounded_time, 'seconds')       
                    
                
                # umap
            if 'Cluster' in st.session_state['phenocluster__clustering_adata'].obs.columns:
                st.session_state['input_dataset'].data["Phenotype_Cluster"] = st.session_state['phenocluster__clustering_adata'].obs["Cluster"]
                #st.write(pd.unique(st.session_state['phenocluster__clustering_adata'].obs["Cluster"]))

                st.session_state['phenocluster__umeta_columns'] = list(st.session_state['phenocluster__clustering_adata'].obs.columns)
                st.session_state['phenocluster__umap_color_col_index'] = st.session_state['phenocluster__umeta_columns'].index(st.session_state['phenocluster__umap_color_col'])
                #st.write(st.session_state['phenocluster__umap_color_col_index'])
                
                # select column for umap coloring
                st.session_state['phenocluster__umap_color_col'] = st.selectbox('Select column for groups coloring:', 
                                                                        st.session_state['phenocluster__umeta_columns'],
                                                                        index=st.session_state['phenocluster__umap_color_col_index']
                                                                        )
                
                # select column for umap subsetting
                st.session_state['phenocluster__umap_cur_col'] = st.selectbox('Select column to subset plots:', 
                                                                        st.session_state['phenocluster__umeta_columns'], key='phenocluster__umap_col_dropdown_subset'
                                                                        )
                
                # list of available subsetting options
                umap_cur_groups=  ["All"] + list(pd.unique(st.session_state['phenocluster__clustering_adata'].obs[st.session_state['phenocluster__umap_cur_col']]))
                umap_sel_groups = st.multiselect('Select groups to be plotted',
                                                                                options = umap_cur_groups)
                st.session_state['phenocluster__umap_cur_groups'] = umap_sel_groups
                
                st.button('Make Spatial Plots' , on_click=spatial_plots_cust_2, args = [st.session_state['phenocluster__clustering_adata'], 
                st.session_state['phenocluster__umap_cur_col'], 
                st.session_state['phenocluster__umap_cur_groups'],
                st.session_state['phenocluster__umap_color_col']
                ]
                        )
                
                st.button("Compute UMAPs", on_click=phenocluster__scanpy_umap, args = [st.session_state['phenocluster__clustering_adata'],
                                                                                    st.session_state['phenocluster__n_neighbors_state'],
                                                                                    st.session_state['phenocluster__metric'],
                                                                                    st.session_state['phenocluster__n_principal_components']
                                                                                    ]
                        )
                
                st.button('Plot UMAPs' , on_click=phenocluster__plotly_umaps, args = [st.session_state['phenocluster__clustering_adata'], 
                st.session_state['phenocluster__umap_cur_col'], 
                st.session_state['phenocluster__umap_cur_groups'],
                st.session_state['phenocluster__umap_color_col'],
                ]
                        )
            
                
            
    

# Run the main function
if __name__ == '__main__':

    # Set page settings
    page_name = 'Unsupervised Phenotype Clustering'
    st.set_page_config(layout='wide', page_title=page_name)
    st.title(page_name)
    phenocluster__col_0 = st.columns(1)
    phenocluster__col1, phenocluster__col2 = st.columns([1, 6])
    
    # Run streamlit-dataframe-editor library initialization tasks at the top of the page
    st.session_state = sde.initialize_session_state(st.session_state)

    # Run Top of Page (TOP) functions
    st.session_state = top.top_of_page_reqs(st.session_state)

    # Call the main function
    main()

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)
    
# need to make differential expression on another page 
