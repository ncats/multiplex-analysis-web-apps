# Import relevant libraries
import streamlit as st
import hnswlib
import parc
#from parc import PARC
import annoy
import sklearn_ann
from ast import arg
from pyparsing import col
import pandas as pd
import anndata as ad
import scanpy as sc
import seaborn as sns
import os
import matplotlib.pyplot as plt
import phenograph
import numpy as np
import scanpy.external as sce
import plotly.express as px
import time
#from pynndescent import PyNNDescentTransformer
from scipy.sparse import csr_matrix
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn_ann.utils import TransformerChecksMixin
import typing as tp
import anndata
from tqdm import tqdm
import parmap
import typing as tp
import scipy
import squidpy as sq
import leidenalg
import igraph as ig
from scipy import stats
from igraph.community import _community_leiden
community_leiden = _community_leiden
import framework.utils as framework_utils
import framework.analysis_framework as analysis_framework

def phenocluster__make_adata(df, x_cols, meta_cols, 
                             z_normalize, normalize_total, 
                             log_normalize, select_high_var_features, n_features):
    
    print(select_high_var_features)
    print(n_features)
    
    mat = df[x_cols]
    meta = df[meta_cols]
    adata = ad.AnnData(mat)
    adata.obs = meta
    adata.layers["counts"] = adata.X.copy()
    #adata.write("input/clust_dat.h5ad")
    if normalize_total:
        sc.pp.normalize_total(adata)
    if log_normalize:
        sc.pp.log1p(adata)
    if z_normalize:
        sc.pp.scale(adata)
        
    if select_high_var_features:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_features, flavor='cell_ranger')
        adata = adata[:, adata.var.highly_variable].copy()
    print(adata.shape)
    return adata

def phenocluster__scanpy_umap(adata, n_neighbors, metric, n_principal_components):
    if "X_pca" not in adata.obsm:
            if n_principal_components > 0:
                sc.pp.pca(adata, n_comps=n_principal_components)
    if "neighbors" not in adata.uns:
        print("Finding nearest neigbours")         
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, metric=metric, n_pcs=n_principal_components)
    sc.tl.umap(adata)
    st.session_state['phenocluster__clustering_adata'] = adata

# plot umaps
def phenocluster__plotly_umaps(adata, umap_cur_col, umap_cur_groups, umap_color_col, plot_col):
    with plot_col:
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
                ),
                xaxis=dict(
                    scaleanchor="y",
                    scaleratio=1)
            )
            if i % 2 == 0:
                subcol1.plotly_chart(fig, use_container_width=True)
            else:
                subcol2.plotly_chart(fig, use_container_width=True)

# plot spatial
def spatial_plots_cust_2(adata, umap_cur_col, umap_cur_groups, umap_color_col, plot_col):
    with plot_col:
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
                    yanchor='top',
                ),
                legend=dict(
                    orientation="h",
                    yanchor="top",
                    y=-0.2,
                    xanchor="right",
                    x=1
                ),
                xaxis=dict(
                    scaleanchor="y",
                    scaleratio=1)
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
        
    if 'phenocluster__n_jobs' not in st.session_state:
        st.session_state['phenocluster__n_jobs'] = 7
        
    if 'phenocluster__n_iterations' not in st.session_state:
        st.session_state['phenocluster__n_iterations'] = 5
        
    if 'phenocluster__n_features' not in st.session_state:
        st.session_state['phenocluster__n_features'] = 0

    # phenograph options
    if 'phenocluster__n_neighbors_state' not in st.session_state:
        st.session_state['phenocluster__n_neighbors_state'] = 30
    
    if 'phenocluster__phenograph_clustering_algo' not in st.session_state:
        st.session_state['phenocluster__phenograph_clustering_algo'] = 'leiden'
    
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
    
def phenocluster_select_high_var_features(adata, n_features):
    sc.pp.highly_variable_features(adata, n_top_features=n_features)
    adata = adata[:, adata.var.highly_variable]
    print(adata)
    return adata
   
def phenocluster__add_clusters_to_input_df():
    if "phenocluster__phenotype_cluster_cols" in st.session_state:
        cur_df = st.session_state['input_dataset'].data
        cur_df = cur_df.drop(columns=st.session_state["phenocluster__phenotype_cluster_cols"])
        st.session_state['input_dataset'].data = cur_df
    print(pd.unique(st.session_state['phenocluster__clustering_adata'].obs['Cluster']))
    st.session_state['input_dataset'].data["Phenotype_Cluster"] = 'Phenotype ' + str(st.session_state['phenocluster__clustering_adata'].obs["Cluster"])
    print(st.session_state['input_dataset'].data["Phenotype_Cluster"])
    dummies = pd.get_dummies(st.session_state['phenocluster__clustering_adata'].obs["Cluster"], prefix='Phenotype Cluster').astype(int)
    #dummies = dummies.replace({1: '+', 0: '-'})
    cur_df = pd.concat([st.session_state['input_dataset'].data, dummies], axis=1)
    st.session_state['input_dataset'].data = cur_df
    new_cluster_cols = list(dummies.columns)
    st.session_state["phenocluster__phenotype_cluster_cols"] = new_cluster_cols
    print(st.session_state["phenocluster__phenotype_cluster_cols"])

# check that only numeric columns are included in the adata.X
def phenocluster__check_input_dat(input_dat, numeric_cols):
    for cur_col in numeric_cols:
        if pd.api.types.is_numeric_dtype(input_dat[cur_col]):
            pass
        else:
            st.error("Column " + cur_col + " is not numeric. Only numeric columns can be included in the matrix",
                     icon="🚨")          


# main
def main():
    """
    Main function for the page.
    """
    #st.write(st.session_state['unifier__df'].head())
    phenocluster__col_0, phenocluster__col_0a = st.columns([10,1])

    phenocluster__col1, phenocluster__col2 = st.columns([2, 6])
    # set default values
    phenocluster__default_session_state()
    
    
    # make layout with columns    
    # options
    
    with phenocluster__col_0:
        st.multiselect('Select numeric columns for clustering:', options = st.session_state['input_dataset'].data.columns, 
                    key='phenocluster__X_cols')
        
        numeric_cols = st.session_state['phenocluster__X_cols']
        phenocluster__check_input_dat(input_dat=st.session_state['input_dataset'].data, numeric_cols=numeric_cols)
        
        st.multiselect('Select columns for metadata:', options = st.session_state['input_dataset'].data.columns, 
            key='phenocluster__meta_cols')
                
                
        meta_columns = st.session_state['phenocluster__meta_cols']
        #Add the new items if they don't already exist in the list
        items_to_add = ['Centroid X (µm)_(standardized)', 'Centroid Y (µm)_(standardized)']
        for item in items_to_add:
            if item not in st.session_state['phenocluster__meta_cols']:
                st.session_state['phenocluster__meta_cols'].append(item)
                
        st.toggle("Z-score normalize columns", key='phenocluster__zscore_normalize')
        st.toggle("Normalize total intensity", key='phenocluster__normalize_total_intensity')
        st.toggle("Log normalize", key='phenocluster__log_normalize')
        st.toggle("Select high variance features", key='phenocluster__select_high_var_features')
        if st.session_state['phenocluster__select_high_var_features'] == True:
            st.number_input(label = "Number of features", key='phenocluster__n_features', step = 1)
            
        
        if st.button('Submit columns'):
            st.session_state['phenocluster__clustering_adata'] = phenocluster__make_adata(st.session_state['input_dataset'].data, 
                                            numeric_cols,
                                            meta_columns,
                                            z_normalize = st.session_state['phenocluster__zscore_normalize'],
                                            normalize_total = st.session_state['phenocluster__normalize_total_intensity'],
                                            log_normalize = st.session_state['phenocluster__log_normalize'],
                                            select_high_var_features = st.session_state['phenocluster__select_high_var_features'],
                                            n_features = st.session_state['phenocluster__n_features']
                                            )
        
    if 'phenocluster__clustering_adata' in st.session_state:
    
        with phenocluster__col1:
            
            # subset data
            st.checkbox('Subset Data', key='phenocluster__subset_data', help = '''Subset data based on a variable''')
            if st.session_state['phenocluster__subset_data'] == True:
                st.session_state['phenocluster__subset_options'] = list(st.session_state['phenocluster__clustering_adata'].obs.columns)
                phenocluster__subset_col = st.selectbox('Select column for subsetting:', st.session_state['phenocluster__subset_options'])
                st.session_state["phenocluster__subset_col"] = phenocluster__subset_col 
                st.session_state['phenocluster__subset_values_options'] = list(pd.unique(st.session_state['phenocluster__clustering_adata'].obs[st.session_state["phenocluster__subset_col"]]))
                phenocluster__subset_vals = st.multiselect('Select a group for subsetting:', options = st.session_state['phenocluster__subset_values_options'], key='phenocluster__subset_vals_1')
                st.session_state["phenocluster__subset_vals"] = phenocluster__subset_vals 
                if st.button('Subset Data'):
                    phenocluster__subset_data(st.session_state['phenocluster__clustering_adata'],
                                            st.session_state["phenocluster__subset_col"],
                                            st.session_state["phenocluster__subset_vals"])
                
                            
            clusteringMethods = ['phenograph', 'scanpy', 'parc', 'utag']
            selected_clusteringMethod = st.selectbox('Select Clustering method:', clusteringMethods, 
                                                    key='clusteringMethods_dropdown') 

            # Update session state on every change
            st.session_state['phenocluster__cluster_method'] = selected_clusteringMethod

            # default widgets
            
            st.number_input(label = "Number of Principal Components", key='phenocluster__n_principal_components', step = 1, 
                            help='''Number of principal components to use for clustering.
                            If 0, Clustering will be performed on a numeric matrx (0 cannot be used for UTAG clustering)''')
            
            #st.session_state['phenocluster__n_neighbors_state']  = st.number_input(label = "K Nearest Neighbors", 
            #                        value=st.session_state['phenocluster__n_neighbors_state'])
            st.number_input(label = "K Nearest Neighbors", 
                                    key='phenocluster__n_neighbors_state', step = 1,
                                    help = '''The size of local neighborhood (in terms of number of neighboring data points) used for manifold approximation. 
                                    Larger values result in more global views of the manifold, while smaller values result in more local data being preserved. 
                                    In general values should be in the range 2 to 100''')
            
            st.number_input(label = "Clustering resolution", key='phenocluster__resolution', step = 0.1,format="%.1f",
                            help = '''A parameter value controlling the coarseness of the clustering. 
                            Higher values lead to more clusters''')
            
            st.number_input(label = "n_jobs", key='phenocluster__n_jobs', step=1,
                help = '''N threads to use''')
            
            st.number_input(label = "n_iterations", key='phenocluster__n_iterations', step=1,
                help = '''N iterations to use for leiden clustering''')
            
            if st.session_state['phenocluster__cluster_method'] == "phenograph":
                st.selectbox('Phenograph clustering algorithm:', ['louvain', 'leiden'], key='phenocluster__phenograph_clustering_algo')
                st.number_input(label = "Phenograph min cluster size", key='phenocluster__phenograph_min_cluster_size', step = 1,
                                help = '''
                                Cells that end up in a cluster smaller than min_cluster_size are considered
                                outliers and are assigned to -1 in the cluster labels
                                ''')
                st.selectbox('Distance metric:', ['euclidean', 'manhattan', 'correlation', 'cosine'], key='phenocluster__metric',
                            help='''Distance metric to define nearest neighbors.''')
                st.selectbox('Phenograph nn method:', ['kdtree', 'brute'], key='phenocluster__phenograph_nn_method',
                            help = '''Whether to use brute force or kdtree for nearest neighbor search.''')
                st.checkbox('Fast', key='phenocluster__fast', help = '''Use aproximate nearest neigbour search''')
            
            elif st.session_state['phenocluster__cluster_method'] == "scanpy":
                st.selectbox('Distance metric:', ['euclidean', 'manhattan', 'correlation', 'cosine'], key='phenocluster__metric',
                            help='''Distance metric to define nearest neighbors.''')
                st.checkbox('Fast', key='phenocluster__scanpy_fast', help = '''Use aproximate nearest neigbour search''')
                if st.session_state['phenocluster__scanpy_fast'] == True:
                    st.selectbox('Transformer:', ['Annoy', 'PNNDescent'], key='phenocluster__scanpy_transformer',
                            help = '''Transformer for the approximate nearest neigbours search''')
                else:
                    st.session_state["phenocluster__scanpy_transformer"] = None
            
            elif st.session_state['phenocluster__cluster_method'] == "parc":
                # make parc specific widgets
                st.number_input(label = "Parc dist std local", key='phenocluster__parc_dist_std_local', step = 1,
                                help = '''local pruning threshold: the number of standard deviations above the mean minkowski 
                                distance between neighbors of a given node. 
                                The higher the parameter, the more edges are retained.''')
                st.number_input(label = "Parc jac std global", key='phenocluster__parc_jac_std_global', step = 0.01,
                                help = '''Global level graph pruning. This threshold can also be set as the number of standard deviations below the network's 
                                mean-jaccard-weighted edges. 0.1-1 provide reasonable pruning. higher value means less pruning. 
                                e.g. a value of 0.15 means all edges that are above mean(edgeweight)-0.15*std(edge-weights) are retained.''')
                st.number_input(label = "Minimum cluster size to be considered a separate population",
                                key='phenocluster__parc_small_pop', step = 1,
                                help = '''Smallest cluster population to be considered a community.''')
                st.number_input(label = "Random seed", key='phenocluster__random_seed', step = 1,
                                help = '''enable reproducible Leiden clustering''')
                st.number_input(label = "HNSW exploration factor for construction", 
                                key='phenocluster__hnsw_param_ef_construction', step = 1,
                                help = '''Higher value increases accuracy of index construction. 
                                Even for several 100,000s of cells 150-200 is adequate''')
                st.checkbox('Fast', key='phenocluster__fast', help = '''Use aproximate nearest neigbour search''')
            elif st.session_state['phenocluster__cluster_method'] == "utag":
                # make utag specific widgets
                #st.selectbox('UTAG clustering method:', ['leiden', 'parc'], key='phenocluster__utag_clustering_method')
                st.number_input(label = "UTAG max dist", key='phenocluster__utag_max_dist', step = 1,
                                help = '''Threshold euclidean distance to determine whether a pair of cell is adjacent in graph structure. 
                                Recommended values are between 10 to 100 depending on magnification.''')
                st.checkbox('Fast', key='phenocluster__utag_fast', help = '''Use aproximate nearest neigbour search''')
                if st.session_state['phenocluster__utag_fast'] == True:
                    st.selectbox('Transformer:', ['Annoy', 'PNNDescent'], key='phenocluster__utag_transformer',
                            help = '''Transformer for the approximate nearest neigbours search''')
                else:
                    st.session_state["phenocluster__utag_transformer"] = None
            
            # add options if clustering has been run
            # add options if clustering has been run
            # if st.button('Run Clustering'):
            #     start_time = time.time()
            #     if st.session_state['phenocluster__cluster_method'] == "phenograph":
            #         with st.spinner('Wait for it...'):
            #             st.session_state['phenocluster__clustering_adata'] = RunPhenographClust(adata=st.session_state['phenocluster__clustering_adata'], 
            #                                                                                     n_neighbors=st.session_state['phenocluster__n_neighbors_state'],
            #                                                                                     clustering_algo=st.session_state['phenocluster__phenograph_clustering_algo'],
            #                                                                                     min_cluster_size=st.session_state['phenocluster__phenograph_min_cluster_size'],
            #                                                                                     primary_metric=st.session_state['phenocluster__metric'],
            #                                                                                     resolution_parameter=st.session_state['phenocluster__resolution'],
            #                                                                                     nn_method=st.session_state['phenocluster__phenograph_nn_method'],
            #                                                                                     random_seed=st.session_state['phenocluster__random_seed'],
            #                                                                                     n_principal_components=st.session_state['phenocluster__n_principal_components'],
            #                                                                                     n_jobs=st.session_state['phenocluster__n_jobs'],
            #                                                                                     n_iterations= st.session_state['phenocluster__n_iterations'],
            #                                                                                     fast=st.session_state["phenocluster__fast"]
            #                                                                                     )
            #     elif st.session_state['phenocluster__cluster_method'] == "scanpy":
            #         with st.spinner('Wait for it...'):
            #             st.session_state['phenocluster__clustering_adata'] = RunNeighbClust(adata=st.session_state['phenocluster__clustering_adata'], 
            #                                                                                 n_neighbors=st.session_state['phenocluster__n_neighbors_state'],
            #                                                                                 metric=st.session_state['phenocluster__metric'],
            #                                                                                 resolution=st.session_state['phenocluster__resolution'],
            #                                                                                 random_state=st.session_state['phenocluster__random_seed'],
            #                                                                                 n_principal_components=st.session_state['phenocluster__n_principal_components'],
            #                                                                                 n_jobs=st.session_state['phenocluster__n_jobs'],
            #                                                                                 n_iterations= st.session_state['phenocluster__n_iterations'],
            #                                                                                 fast=st.session_state["phenocluster__scanpy_fast"],
            #                                                                                 transformer = st.session_state["phenocluster__scanpy_transformer"]
            #                                                                                 )
            #     #st.session_state['phenocluster__clustering_adata'] = adata
            #     elif st.session_state['phenocluster__cluster_method'] == "parc":
            #         with st.spinner('Wait for it...'):                  
            #             st.session_state['phenocluster__clustering_adata'] = run_parc_clust(adata=st.session_state['phenocluster__clustering_adata'], 
            #                                                                                 n_neighbors=st.session_state['phenocluster__n_neighbors_state'],
            #                                                                                 dist_std_local=st.session_state['phenocluster__parc_dist_std_local'],
            #                                                                                 jac_std_global= st.session_state['phenocluster__parc_jac_std_global'],
            #                                                                                 small_pop=st.session_state['phenocluster__parc_small_pop'],
            #                                                                                 random_seed=st.session_state['phenocluster__random_seed'],
            #                                                                                 resolution_parameter=st.session_state['phenocluster__resolution'],
            #                                                                                 hnsw_param_ef_construction=st.session_state['phenocluster__hnsw_param_ef_construction'],
            #                                                                                 n_principal_components=st.session_state['phenocluster__n_principal_components'],
            #                                                                                 n_jobs=st.session_state['phenocluster__n_jobs'],
            #                                                                                 n_iterations= st.session_state['phenocluster__n_iterations'],
            #                                                                                 fast=st.session_state["phenocluster__fast"]
            #                                                                                 )
            #     elif st.session_state['phenocluster__cluster_method'] == "utag":
            #         #phenocluster__utag_resolutions = [st.session_state['phenocluster__resolution']]
            #         with st.spinner('Wait for it...'):
            #             st.session_state['phenocluster__clustering_adata'] = run_utag_clust(adata=st.session_state['phenocluster__clustering_adata'], 
            #                                                                                 n_neighbors=st.session_state['phenocluster__n_neighbors_state'], 
            #                                                                                 resolution=st.session_state['phenocluster__resolution'],
            #                                                                                 clustering_method=st.session_state['phenocluster__utag_clustering_method'],
            #                                                                                 max_dist=st.session_state['phenocluster__utag_max_dist'],
            #                                                                                 n_principal_components=st.session_state['phenocluster__n_principal_components'],
            #                                                                                 random_state=st.session_state['phenocluster__random_seed'],
            #                                                                                 n_jobs=st.session_state['phenocluster__n_jobs'],
            #                                                                                 n_iterations= st.session_state['phenocluster__n_iterations'],
            #                                                                                 fast=st.session_state["phenocluster__utag_fast"],
            #                                                                                 transformer = st.session_state["phenocluster__utag_transformer"]
            #                                                                                 )
                # save clustering result
                #st.session_state['phenocluster__clustering_adata'].write("input/clust_dat.h5ad")
                # end_time = time.time()
                # execution_time = end_time - start_time
                # rounded_time = round(execution_time, 2)
                # st.write('Execution time: ', rounded_time, 'seconds')       

            if st.session_state['phenocluster__cluster_method'] == "scanpy":
                analysis_framework.job_submission(
                    job_name="run_scanpy_clust",
                    inputs = dict(
                        adata=st.session_state['phenocluster__clustering_adata'], 
                        n_neighbors=st.session_state['phenocluster__n_neighbors_state'],
                        metric=st.session_state['phenocluster__metric'],
                        resolution=st.session_state['phenocluster__resolution'],
                        random_state=st.session_state['phenocluster__random_seed'],
                        n_principal_components=st.session_state['phenocluster__n_principal_components'],
                        n_jobs=st.session_state['phenocluster__n_jobs'],
                        n_iterations= st.session_state['phenocluster__n_iterations'],
                        fast=st.session_state["phenocluster__scanpy_fast"],
                        transformer = st.session_state["phenocluster__scanpy_transformer"]
                    ),
                    analysis_purpose = "run scanpy clustering",
                    st_key_prefix = "",
                )
                unsup_clust_key = "run_scanpy_clustering_results"
                if unsup_clust_key not in st.session_state:
                    st.warning("Unsupervised clustering results are not yet available.")
                    return
                
                st.session_state['phenocluster__clustering_adata'] = st.session_state[unsup_clust_key]['adata']

            elif st.session_state['phenocluster__cluster_method'] == "phenograph":
                analysis_framework.job_submission(
                    job_name="run_phenograph_clust",
                    inputs = dict(
                        adata=st.session_state['phenocluster__clustering_adata'], 
                        n_neighbors=st.session_state['phenocluster__n_neighbors_state'],
                        clustering_algo=st.session_state['phenocluster__phenograph_clustering_algo'],
                        min_cluster_size=st.session_state['phenocluster__phenograph_min_cluster_size'],
                        primary_metric=st.session_state['phenocluster__metric'],
                        resolution_parameter=st.session_state['phenocluster__resolution'],
                        nn_method=st.session_state['phenocluster__phenograph_nn_method'],
                        random_seed=st.session_state['phenocluster__random_seed'],
                        n_principal_components=st.session_state['phenocluster__n_principal_components'],
                        n_jobs=st.session_state['phenocluster__n_jobs'],
                        n_iterations= st.session_state['phenocluster__n_iterations'],
                        fast=st.session_state["phenocluster__fast"]
                    ),
                    analysis_purpose = "run phenograph clustering",
                    st_key_prefix = "",
                )
                unsup_clust_key = "run_phenograph_clustering_results"
                if unsup_clust_key not in st.session_state:
                    st.warning("Unsupervised clustering results are not yet available.")
                    return
                
                st.session_state['phenocluster__clustering_adata'] = st.session_state[unsup_clust_key]['adata']

            elif st.session_state['phenocluster__cluster_method'] == "parc":
                analysis_framework.job_submission(
                    job_name="run_parc_clust",
                    inputs = dict(
                        adata=st.session_state['phenocluster__clustering_adata'], 
                        n_neighbors=st.session_state['phenocluster__n_neighbors_state'],
                        dist_std_local=st.session_state['phenocluster__parc_dist_std_local'],
                        jac_std_global= st.session_state['phenocluster__parc_jac_std_global'],
                        small_pop=st.session_state['phenocluster__parc_small_pop'],
                        random_seed=st.session_state['phenocluster__random_seed'],
                        resolution_parameter=st.session_state['phenocluster__resolution'],
                        hnsw_param_ef_construction=st.session_state['phenocluster__hnsw_param_ef_construction'],
                        n_principal_components=st.session_state['phenocluster__n_principal_components'],
                        n_jobs=st.session_state['phenocluster__n_jobs'],
                        n_iterations= st.session_state['phenocluster__n_iterations'],
                        fast=st.session_state["phenocluster__fast"]
                    ),
                    analysis_purpose = "run parc clustering",
                    st_key_prefix = "",
                )
                unsup_clust_key = "run_parc_clustering_results"
                if unsup_clust_key not in st.session_state:
                    st.warning("Unsupervised clustering results are not yet available.")
                    return
                
                st.session_state['phenocluster__clustering_adata'] = st.session_state[unsup_clust_key]['adata']

            elif st.session_state['phenocluster__cluster_method'] == "utag":
                analysis_framework.job_submission(
                    job_name="run_utag_clust",
                    inputs = dict(
                        adata=st.session_state['phenocluster__clustering_adata'], 
                        n_neighbors=st.session_state['phenocluster__n_neighbors_state'], 
                        resolution=st.session_state['phenocluster__resolution'],
                        clustering_method=st.session_state['phenocluster__utag_clustering_method'],
                        max_dist=st.session_state['phenocluster__utag_max_dist'],
                        n_principal_components=st.session_state['phenocluster__n_principal_components'],
                        random_state=st.session_state['phenocluster__random_seed'],
                        n_jobs=st.session_state['phenocluster__n_jobs'],
                        n_iterations= st.session_state['phenocluster__n_iterations'],
                        fast=st.session_state["phenocluster__utag_fast"],
                        transformer = st.session_state["phenocluster__utag_transformer"]
                    ),
                    analysis_purpose = "run utag clustering",
                    st_key_prefix = "",
                )
                unsup_clust_key = "run_utag_clustering_results"
                if unsup_clust_key not in st.session_state:
                    st.warning("Unsupervised clustering results are not yet available.")
                    return
                
                st.session_state['phenocluster__clustering_adata'] = st.session_state[unsup_clust_key]['adata']
            
            # umap
            if 'Cluster' in st.session_state['phenocluster__clustering_adata'].obs.columns:
                
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
                st.session_state['phenocluster__umap_color_col'],
                phenocluster__col2
                ]
                        )
                
                st.button("Compute UMAP", on_click=phenocluster__scanpy_umap, args = [st.session_state['phenocluster__clustering_adata'],
                                                                                    st.session_state['phenocluster__n_neighbors_state'],
                                                                                    st.session_state['phenocluster__metric'],
                                                                                    st.session_state['phenocluster__n_principal_components']
                                                                                    ]
                        )
                if 'X_umap' in st.session_state['phenocluster__clustering_adata'].obsm.keys():
                    st.button('Plot UMAPs' , on_click=phenocluster__plotly_umaps, 
                            args = [st.session_state['phenocluster__clustering_adata'], 
                                    st.session_state['phenocluster__umap_cur_col'], 
                                    st.session_state['phenocluster__umap_cur_groups'],
                                    st.session_state['phenocluster__umap_color_col'],
                                    phenocluster__col2
                                    ]
                            )
                
                st.button('Add Clusters to Input Data' , on_click=phenocluster__add_clusters_to_input_df)
        
                

# Run the main function
if __name__ == '__main__':
    main()
    
# need to make differential expression on another page 
