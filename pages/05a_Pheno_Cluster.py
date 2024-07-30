# Import relevant libraries
import streamlit as st
import subprocess


try:
    import parc
except ImportError:
    subprocess.run("pip install parc", shell=True)
    try:
        import parc
    except ImportError:
       print("Failed to import parc.")
    
try:
    import annoy
except ImportError:
    subprocess.run("pip install annoy", shell=True)
    try:
        import annoy
    except ImportError:
       print("Failed to import annoy.")

try:
    import sklearn_ann
except ImportError:
    subprocess.run("pip install sklearn-ann", shell=True)
    try:
        import sklearn_ann
    except ImportError:
       print("Failed to import sklearn-ann.") 


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
import numpy as np
import scanpy.external as sce
import plotly.express as px
import time
from pynndescent import PyNNDescentTransformer
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

def z_score(x):
    """
    Scale (divide by standard deviation) and center (subtract mean) array-like objects.
    """
    return (x - x.min()) / (x.max() - x.min())

def sparse_matrix_dstack(
    matrices: tp.Sequence[scipy.sparse.csr_matrix],
) -> scipy.sparse.csr_matrix:
    """
    Diagonally stack sparse matrices.
    """
    import scipy
    from tqdm import tqdm

    n = sum([x.shape[0] for x in matrices])
    _res = list()
    i = 0
    for x in tqdm(matrices):
        v = scipy.sparse.csr_matrix((x.shape[0], n))
        v[:, i : i + x.shape[0]] = x
        _res.append(v)
        i += x.shape[0]
    return scipy.sparse.vstack(_res)

def utag(
    adata,
    channels_to_use = None,
    slide_key = "Slide",
    save_key: str = "UTAG Label",
    filter_by_variance: bool = False,
    max_dist: float = 20.0,
    normalization_mode: str = "l1_norm",
    keep_spatial_connectivity: bool = False,
    pca_kwargs: tp.Dict[str, tp.Any] = dict(n_comps=10),
    apply_umap: bool = False,
    umap_kwargs: tp.Dict[str, tp.Any] = dict(),
    apply_clustering: bool = True,
    clustering_method: tp.Sequence[str] = ["leiden", "parc", "kmeans"],
    resolutions: tp.Sequence[float] = [0.05, 0.1, 0.3, 1.0],
    leiden_kwargs: tp.Dict[str, tp.Any] = None,
    parc_kwargs: tp.Dict[str, tp.Any] = None,
    parallel: bool = True,
    processes: int = None,
):
    """
    Discover tissue architechture in single-cell imaging data
    by combining phenotypes and positional information of cells.

    Parameters
    ----------
    adata: AnnData
        AnnData object with spatial positioning of cells in obsm 'spatial' slot.
    channels_to_use: Optional[Sequence[str]]
        An optional sequence of strings used to subset variables to use.
        Default (None) is to use all variables.
    max_dist: float
        Maximum distance to cut edges within a graph.
        Should be adjusted depending on resolution of images.
        For imaging mass cytometry, where resolution is 1um, 20 often gives good results.
        Default is 20.
    slide_key: {str, None}
        Key of adata.obs containing information on the batch structure of the data.
        In general, for image data this will often be a variable indicating the image
        so image-specific effects are removed from data.
        Default is "Slide".
    save_key: str
        Key to be added to adata object holding the UTAG clusters.
        Depending on the values of `clustering_method` and `resolutions`,
        the final keys will be of the form: {save_key}_{method}_{resolution}".
        Default is "UTAG Label".
    filter_by_variance: bool
        Whether to filter vairiables by variance.
        Default is False, which keeps all variables.
    max_dist: float
        Recommended values are between 20 to 50 depending on magnification.
        Default is 20.
    normalization_mode: str
        Method to normalize adjacency matrix.
        Default is "l1_norm", any other value will not use normalization.
    keep_spatial_connectivity: bool
        Whether to keep sparse matrices of spatial connectivity and distance in the obsp attribute of the
        resulting anndata object. This could be useful in downstream applications.
        Default is not to (False).
    pca_kwargs: Dict[str, Any]
        Keyword arguments to be passed to scanpy.pp.pca for dimensionality reduction after message passing.
        Default is to pass n_comps=10, which uses 10 Principal Components.
    apply_umap: bool
        Whether to build a UMAP representation after message passing.
        Default is False.
    umap_kwargs: Dict[str, Any]
        Keyword arguments to be passed to scanpy.tl.umap for dimensionality reduction after message passing.
        Default is 10.0.
    apply_clustering: bool
        Whether to cluster the message passed matrix.
        Default is True.
    clustering_method: Sequence[str]
        Which clustering method(s) to use for clustering of the message passed matrix.
        Default is ["leiden", "parc"].
    resolutions: Sequence[float]
        What resolutions should the methods in `clustering_method` be run at.
        Default is [0.05, 0.1, 0.3, 1.0].
    leiden_kwargs: dict[str, Any]
        Keyword arguments to pass to scanpy.tl.leiden.
    parc_kwargs: dict[str, Any]
        Keyword arguments to pass to parc.PARC.
    parallel: bool
        Whether to run message passing part of algorithm in parallel.
        Will accelerate the process but consume more memory.
        Default is True.
    processes: int
        Number of processes to use in parallel.
        Default is to use all available (-1).

    Returns
    -------
    adata: AnnData
        AnnData object with UTAG domain predictions for each cell in adata.obs, column `save_key`.
    """
    ad = adata.copy()

    if channels_to_use:
        ad = ad[:, channels_to_use]

    if filter_by_variance:
        ad = low_variance_filter(ad)

    if isinstance(clustering_method, list):
        clustering_method = [m.upper() for m in clustering_method]
    elif isinstance(clustering_method, str):
        clustering_method = [clustering_method.upper()]
    else:
        print(
            "Invalid Clustering Method. Clustering Method Should Either be a string or a list"
        )
        return
    assert all(m in ["LEIDEN", "PARC", "KMEANS"] for m in clustering_method)

    if "PARC" in clustering_method:
        from parc import PARC  # early fail if not available
    if "KMEANS" in clustering_method:
        from sklearn.cluster import KMeans

    print("Applying UTAG Algorithm...")
    if slide_key:
        ads = [
            ad[ad.obs[slide_key] == slide].copy() for slide in ad.obs[slide_key].unique()
        ]
        ad_list = parmap.map(
            _parallel_message_pass,
            ads,
            radius=max_dist,
            coord_type="generic",
            set_diag=True,
            mode=normalization_mode,
            pm_pbar=True,
            pm_parallel=parallel,
            pm_processes=processes,
        )
        ad_result = anndata.concat(ad_list)
        if keep_spatial_connectivity:
            ad_result.obsp["spatial_connectivities"] = sparse_matrix_dstack(
                [x.obsp["spatial_connectivities"] for x in ad_list]
            )
            ad_result.obsp["spatial_distances"] = sparse_matrix_dstack(
                [x.obsp["spatial_distances"] for x in ad_list]
            )
    else:
        sq.gr.spatial_neighbors(ad, radius=max_dist, coord_type="generic", set_diag=True)
        ad_result = custom_message_passing(ad, mode=normalization_mode)

    if apply_clustering:
        if "n_comps" in pca_kwargs:
            if pca_kwargs["n_comps"] > ad_result.shape[1]:
                pca_kwargs["n_comps"] = ad_result.shape[1] - 1
                print(
                    f"Overwriding provided number of PCA dimensions to match number of features: {pca_kwargs['n_comps']}"
                )
        sc.tl.pca(ad_result, **pca_kwargs)
        sc.pp.neighbors(ad_result)

        if apply_umap:
            print("Running UMAP on Input Dataset...")
            sc.tl.umap(ad_result, **umap_kwargs)

        for resolution in tqdm(resolutions):

            res_key1 = save_key + "_leiden_" + str(resolution)
            res_key2 = save_key + "_parc_" + str(resolution)
            res_key3 = save_key + "_kmeans_" + str(resolution)
            if "LEIDEN" in clustering_method:
                print(f"Applying Leiden Clustering at Resolution: {resolution}...")
                kwargs = dict()
                kwargs.update(leiden_kwargs or {})
                sc.tl.leiden(
                    ad_result, resolution=resolution, key_added=res_key1, **kwargs
                )
                add_probabilities_to_centroid(ad_result, res_key1)

            if "PARC" in clustering_method:
                from parc import PARC

                print(f"Applying PARC Clustering at Resolution: {resolution}...")

                kwargs = dict(random_seed=1, small_pop=1000)
                kwargs.update(parc_kwargs or {})
                model = PARC(
                    ad_result.obsm["X_pca"],
                    neighbor_graph=ad_result.obsp["connectivities"],
                    resolution_parameter=resolution,
                    **kwargs,
                )
                model.run_PARC()
                ad_result.obs[res_key2] = pd.Categorical(model.labels)
                ad_result.obs[res_key2] = ad_result.obs[res_key2].astype("category")
                add_probabilities_to_centroid(ad_result, res_key2)

            if "KMEANS" in clustering_method:
                print(f"Applying K-means Clustering at Resolution: {resolution}...")
                k = int(np.ceil(resolution * 10))
                kmeans = KMeans(n_clusters=k, random_state=1).fit(ad_result.obsm["X_pca"])
                ad_result.obs[res_key3] = pd.Categorical(kmeans.labels_.astype(str))
                add_probabilities_to_centroid(ad_result, res_key3)

    return ad_result


def _parallel_message_pass(
    ad,
    radius: int,
    coord_type: str,
    set_diag: bool,
    mode: str,
):
    sq.gr.spatial_neighbors(ad, radius=radius, coord_type=coord_type, set_diag=set_diag)
    ad = custom_message_passing(ad, mode=mode)
    return ad


def custom_message_passing(adata, mode: str = "l1_norm"):
    # from scipy.linalg import sqrtm
    # import logging
    if mode == "l1_norm":
        A = adata.obsp["spatial_connectivities"]
        from sklearn.preprocessing import normalize
        affinity = normalize(A, axis=1, norm="l1")
    else:
        # Plain A_mod multiplication
        A = adata.obsp["spatial_connectivities"]
        affinity = A
    # logging.info(type(affinity))
    adata.X = affinity @ adata.X
    return adata


def low_variance_filter(adata):
    return adata[:, adata.var["std"] > adata.var["std"].median()]


def add_probabilities_to_centroid(
    adata, col: str, name_to_output: str = None
):
    from scipy.special import softmax

    if name_to_output is None:
        name_to_output = col + "_probabilities"

    mean = z_score(adata.to_df()).groupby(adata.obs[col]).mean()
    probs = softmax(adata.to_df() @ mean.T, axis=1)
    adata.obsm[name_to_output] = probs
    return adata


class AnnoyTransformer(TransformerChecksMixin, TransformerMixin, BaseEstimator):
    """Wrapper for using annoy.AnnoyIndex as sklearn's KNeighborsTransformer"""

    def __init__(self, n_neighbors=5, *, metric="euclidean", 
                 n_trees=10, search_k=-1, n_jobs=-1):
        self.n_neighbors = n_neighbors
        self.n_trees = n_trees
        self.search_k = search_k
        self.metric = metric
        self.n_jobs = n_jobs

    def fit(self, X, y=None):
        X = self._validate_data(X)
        self.n_samples_fit_ = X.shape[0]
        metric = self.metric if self.metric != "sqeuclidean" else "euclidean"
        self.annoy_ = annoy.AnnoyIndex(X.shape[1], metric=metric)
        for i, x in enumerate(X):
            self.annoy_.add_item(i, x.tolist())
        self.annoy_.build(self.n_trees, n_jobs = self.n_jobs)
        return self

    def transform(self, X):
        X = self._transform_checks(X, "annoy_")
        return self._transform(X)

    def fit_transform(self, X, y=None):
        return self.fit(X)._transform(X=None)

    def _transform(self, X):
        """As `transform`, but handles X is None for faster `fit_transform`."""

        n_samples_transform = self.n_samples_fit_ if X is None else X.shape[0]

        # For compatibility reasons, as each sample is considered as its own
        # neighbor, one extra neighbor will be computed.
        n_neighbors = self.n_neighbors + 1

        indices = np.empty((n_samples_transform, n_neighbors), dtype=int)
        distances = np.empty((n_samples_transform, n_neighbors))

        if X is None:
            for i in range(self.annoy_.get_n_items()):
                ind, dist = self.annoy_.get_nns_by_item(
                    i, n_neighbors, self.search_k, include_distances=True
                )

                indices[i], distances[i] = ind, dist
        else:
            for i, x in enumerate(X):
                indices[i], distances[i] = self.annoy_.get_nns_by_vector(
                    x.tolist(), n_neighbors, self.search_k, include_distances=True
                )

        if self.metric == "sqeuclidean":
            distances **= 2

        indptr = np.arange(0, n_samples_transform * n_neighbors + 1, n_neighbors)
        kneighbors_graph = csr_matrix(
            (distances.ravel(), indices.ravel(), indptr),
            shape=(n_samples_transform, self.n_samples_fit_),
        )

        return kneighbors_graph

    def _more_tags(self):
        return {
            "_xfail_checks": {"check_estimators_pickle": "Cannot pickle AnnoyIndex"},
            "requires_y": False,
        }


def phenocluster__make_adata(df, x_cols, meta_cols):
    mat = df[x_cols]
    meta = df[meta_cols]
    adata = ad.AnnData(mat)
    adata.obs = meta
    adata.layers["counts"] = adata.X.copy()
    #adata.write("input/clust_dat.h5ad")
    return adata

# scanpy clustering
def RunNeighbClust(adata, n_neighbors, metric, resolution, random_state, n_principal_components, 
                   n_jobs, n_iterations, fast, transformer):
    
    if fast == True:
        sc.pp.pca(adata, n_comps=n_principal_components)
        if transformer == "Annoy":
            sc.pp.neighbors(adata, transformer=AnnoyTransformer(n_neighbors=n_neighbors, metric=metric, n_jobs=n_jobs), 
                       n_pcs=n_principal_components, random_state=random_state)
        elif transformer == "PNNDescent":
            transformer = PyNNDescentTransformer(n_neighbors=n_neighbors, metric=metric, n_jobs=n_jobs, random_state=random_state)
            sc.pp.neighbors(adata, transformer=transformer, n_pcs=n_principal_components, random_state=random_state)
        else:
            sc.pp.neighbors(adata, n_neighbors=n_neighbors, metric=metric, n_pcs=n_principal_components, random_state=random_state)
        
    else:
        if n_principal_components > 0:
            sc.pp.pca(adata, n_comps=n_principal_components)
            sc.pp.neighbors(adata, n_neighbors=n_neighbors, metric=metric, n_pcs=n_principal_components, random_state=random_state)
            
    sc.tl.leiden(adata,resolution=resolution, random_state=random_state, n_iterations=n_iterations, flavor="igraph")
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
def run_utag_clust(adata, n_neighbors, resolution, clustering_method, max_dist, n_principal_components,
                   random_state, n_jobs, n_iterations, fast, transformer):
    #sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=0)
    #sc.tl.umap(adata)
    
    resolutions = [resolution]
    print(resolutions)
    adata.obsm['spatial'] = np.array(adata.obs[["Centroid X (µm)_(standardized)", "Centroid Y (µm)_(standardized)"]])
    
    if fast == True:
        utag_results = utag(adata,
        slide_key="Image ID_(standardized)",
        max_dist=max_dist,
        normalization_mode='l1_norm',
        apply_clustering=False,
        parallel = True,
        processes = n_jobs)
        
        sc.pp.pca(utag_results, n_comps=n_principal_components)
        print("start k graph")
        if transformer == "Annoy":
            sc.pp.neighbors(utag_results, transformer=AnnoyTransformer(n_neighbors=n_neighbors, n_jobs=n_jobs), 
                        n_pcs=n_principal_components, random_state=random_state)
        elif transformer == "PNNDescent":
            transformer = PyNNDescentTransformer(n_neighbors=n_neighbors, n_jobs=n_jobs, random_state=random_state)
            sc.pp.neighbors(utag_results, transformer=transformer, n_pcs=n_principal_components, random_state=random_state)
        else:
            sc.pp.neighbors(utag_results, n_neighbors=n_neighbors,n_pcs=n_principal_components, random_state=random_state)

        resolution_parameter = resolution
        sc.tl.leiden(utag_results,resolution=resolution_parameter, random_state=random_state,
                n_iterations=n_iterations, flavor="igraph")
        utag_results.obs['Cluster'] = utag_results.obs['leiden'].copy()
        adata.uns["leiden"] = utag_results.uns["leiden"].copy() 
        
    else:
        utag_results = utag(adata,
        slide_key="Image ID_(standardized)",
        max_dist=max_dist,
        normalization_mode='l1_norm',
        apply_clustering=True,
        clustering_method = "leiden", 
        resolutions = resolutions,
        leiden_kwargs={"n_iterations": n_iterations, "random_state": random_state},
        pca_kwargs = {"n_comps": n_principal_components},
        parallel = True,
        processes = n_jobs)
        
        curClusterCol = 'UTAG Label_leiden_'  + str(resolution)
        utag_results.obs['Cluster'] = utag_results.obs[curClusterCol].copy()

    cluster_list = list(utag_results.obs['Cluster'])
    print(pd.unique(cluster_list))    
    adata.obsp["distances"] = utag_results.obsp["distances"].copy()
    adata.obsp["connectivities"] = utag_results.obsp["connectivities"].copy()
    adata.obsm["X_pca"] = utag_results.obsm["X_pca"].copy()
    adata.uns["neighbors"] = utag_results.uns["neighbors"].copy()
    adata.varm["PCs"] = utag_results.varm["PCs"].copy()
    adata.obs["Cluster"] = cluster_list
    #utag_results.X = adata.X
    
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
        
    if 'phenocluster__n_jobs' not in st.session_state:
        st.session_state['phenocluster__n_jobs'] = 7
        
    if 'phenocluster__n_iterations' not in st.session_state:
        st.session_state['phenocluster__n_iterations'] = 5

    # phenograph options
    if 'phenocluster__n_neighbors_state' not in st.session_state:
        st.session_state['phenocluster__n_neighbors_state'] = 10
    
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
    
def phenocluster__add_clusters_to_input_df():
    if "phenocluster__phenotype_cluster_cols" in st.session_state:
        cur_df = st.session_state['input_dataset'].data
        cur_df = cur_df.drop(columns=st.session_state["phenocluster__phenotype_cluster_cols"])
        st.session_state['input_dataset'].data = cur_df
    print(pd.unique(st.session_state['phenocluster__clustering_adata'].obs['Cluster']))
    st.session_state['input_dataset'].data["Phenotype_Cluster"] = 'Phenotype ' + str(st.session_state['phenocluster__clustering_adata'].obs["Cluster"])
    print(st.session_state['input_dataset'].data["Phenotype_Cluster"])
    dummies = pd.get_dummies(st.session_state['phenocluster__clustering_adata'].obs["Cluster"], prefix='Phenotype Cluster').astype(int)
    dummies = dummies.replace({1: '+', 0: '-'})
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
    
    # set default values
    phenocluster__default_session_state()
    
    # make layout with columns    
    # options
    
    with phenocluster__col_0[0]:
        
        #select numeric 
        numeric_cols = st.multiselect('Select numeric columns for clustering:', options = st.session_state['input_dataset'].data.columns, 
                    key='phenocluster__X_cols')
        
        phenocluster__check_input_dat(input_dat=st.session_state['input_dataset'].data, numeric_cols=numeric_cols)
        
        meta_columns = st.multiselect('Select columns for metadata:', options = st.session_state['input_dataset'].data.columns, 
                key='phenocluster__meta_cols')
        
        items_to_add = ['Image ID_(standardized)','Centroid X (µm)_(standardized)', 'Centroid Y (µm)_(standardized)']
        
        #meta_columns = st.session_state['phenocluster__meta_cols']
        # Add the new items if they don't already exist in the list
        for item in items_to_add:
            if item not in meta_columns:
                meta_columns.append(item)
        
        if st.button('Submit columns', help = '''Confirm columns selections to create the AnnData object'''):
            st.session_state['phenocluster__clustering_adata'] = phenocluster__make_adata(st.session_state['input_dataset'].data, 
                                            numeric_cols,
                                            meta_columns)
    
    if 'phenocluster__clustering_adata' in st.session_state:
    
        with phenocluster__col1:
            
            # subset data
            st.checkbox('Subset Data:', key='phenocluster__subset_data', help = '''Subset data based on a variable''')
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
                # st.session_state['phenocluster__phenograph_k'] = st.number_input(label = "Phenograph k", 
                #                     value=st.session_state['phenocluster__phenograph_k'])
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
            
            elif st.session_state['phenocluster__cluster_method'] == "scanpy":
                st.selectbox('Distance metric:', ['euclidean', 'manhattan', 'correlation', 'cosine'], key='phenocluster__metric',
                             help='''Distance metric to define nearest neighbors.''')
                st.checkbox('Fast:', key='phenocluster__scanpy_fast', help = '''Use aproximate nearest neigbour search''')
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
            elif st.session_state['phenocluster__cluster_method'] == "utag":
                # make utag specific widgets
                #st.selectbox('UTAG clustering method:', ['leiden', 'parc'], key='phenocluster__utag_clustering_method')
                st.number_input(label = "UTAG max dist", key='phenocluster__utag_max_dist', step = 1,
                                help = '''Threshold euclidean distance to determine whether a pair of cell is adjacent in graph structure. 
                                Recommended values are between 10 to 100 depending on magnification.''')
                st.checkbox('Fast:', key='phenocluster__utag_fast', help = '''Use aproximate nearest neigbour search''')
                if st.session_state['phenocluster__utag_fast'] == True:
                    st.selectbox('Transformer:', ['Annoy', 'PNNDescent'], key='phenocluster__utag_transformer',
                             help = '''Transformer for the approximate nearest neigbours search''')
                else:
                    st.session_state["phenocluster__utag_transformer"] = None
            
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
                                                                                            n_principal_components=st.session_state['phenocluster__n_principal_components'],
                                                                                            n_jobs=st.session_state['phenocluster__n_jobs'],
                                                                                            n_iterations= st.session_state['phenocluster__n_iterations'],
                                                                                            fast=st.session_state["phenocluster__scanpy_fast"],
                                                                                            transformer = st.session_state["phenocluster__scanpy_transformer"]
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
                    #phenocluster__utag_resolutions = [st.session_state['phenocluster__resolution']]
                    with st.spinner('Wait for it...'):
                        st.session_state['phenocluster__clustering_adata'] = run_utag_clust(adata=st.session_state['phenocluster__clustering_adata'], 
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
                                                                                            )
                # save clustering result
                #st.session_state['phenocluster__clustering_adata'].write("input/clust_dat.h5ad")
                end_time = time.time()
                execution_time = end_time - start_time
                rounded_time = round(execution_time, 2)
                st.write('Execution time: ', rounded_time, 'seconds')       
                    
                
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
                st.session_state['phenocluster__umap_color_col']
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
                                      st.session_state['phenocluster__umap_color_col']]
                              )
                
                st.button('Add Clusters to Input Data' , on_click=phenocluster__add_clusters_to_input_df)
            
                
            
    

# Run the main function
if __name__ == '__main__':

    # Set page settings
    page_name = 'Unsupervised Phenotype Clustering'
    st.set_page_config(layout='wide', page_title=page_name)
    st.title(page_name)
    phenocluster__col_0 = st.columns(1)
    phenocluster__col1, phenocluster__col2 = st.columns([2, 6])
    
    # Run streamlit-dataframe-editor library initialization tasks at the top of the page
    st.session_state = sde.initialize_session_state(st.session_state)

    # Run Top of Page (TOP) functions
    st.session_state = top.top_of_page_reqs(st.session_state)

    # Call the main function
    main()

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)
    
# need to make differential expression on another page 
