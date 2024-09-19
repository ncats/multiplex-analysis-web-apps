# Import relevant libraries
import streamlit as st
import hnswlib
import parc
from parc import PARC
import annoy
import sklearn_ann
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
import leidenalg
import igraph as ig
from scipy import stats
from igraph.community import _community_leiden
community_leiden = _community_leiden

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

class PARC_2:
    def __init__(self, data, true_label=None, dist_std_local=3, jac_std_global='median', keep_all_local_dist='auto',
                 too_big_factor=0.4, small_pop=10, jac_weighted_edges=True, knn=30, n_iter_leiden=5, random_seed=42,
                 num_threads=-1, distance='l2', time_smallpop=15, partition_type="ModularityVP",
                 resolution_parameter=1.0,
                 knn_struct=None, neighbor_graph=None, hnsw_param_ef_construction=150):
        # higher dist_std_local means more edges are kept
        # highter jac_std_global means more edges are kept
        if keep_all_local_dist == 'auto':
            if data.shape[0] > 300000:
                keep_all_local_dist = True  # skips local pruning to increase speed
            else:
                keep_all_local_dist = False
        if resolution_parameter != 1:
            partition_type = "RBVP"  # Reichardt and Bornholdt’s Potts model. Note that this is the same as ModularityVertexPartition when setting 𝛾 = 1 and normalising by 2m
        self.data = data
        self.true_label = true_label
        self.dist_std_local = dist_std_local  # similar to the jac_std_global parameter. avoid setting local and global pruning to both be below 0.5 as this is very aggresive pruning.
        self.jac_std_global = jac_std_global  # 0.15 is also a recommended value performing empirically similar to 'median'. Generally values between 0-1.5 are reasonable.
        self.keep_all_local_dist = keep_all_local_dist  # decides whether or not to do local pruning. default is 'auto' which omits LOCAL pruning for samples >300,000 cells.
        self.too_big_factor = too_big_factor  # if a cluster exceeds this share of the entire cell population, then the PARC will be run on the large cluster. at 0.4 it does not come into play
        self.small_pop = small_pop  # smallest cluster population to be considered a community
        self.jac_weighted_edges = jac_weighted_edges  # boolean. whether to partition using weighted graph
        self.knn = knn
        self.n_iter_leiden = n_iter_leiden  # the default is 5 in PARC
        self.random_seed = random_seed  # enable reproducible Leiden clustering
        self.num_threads = num_threads  # number of threads used in KNN search/construction
        self.distance = distance  # Euclidean distance 'l2' by default; other options 'ip' and 'cosine'
        self.time_smallpop = time_smallpop  # number of seconds trying to check an outlier
        self.partition_type = partition_type  # default is the simple ModularityVertexPartition where resolution_parameter =1. In order to change resolution_parameter, we switch to RBConfigurationVP
        self.resolution_parameter = resolution_parameter  # defaults to 1. expose this parameter in leidenalg
        self.knn_struct = knn_struct  # the hnsw index of the KNN graph on which we perform queries
        self.neighbor_graph = neighbor_graph  # CSR affinity matrix for pre-computed nearest neighbors
        self.hnsw_param_ef_construction = hnsw_param_ef_construction  # set at 150. higher value increases accuracy of index construction. Even for several 100,000s of cells 150-200 is adequate

    def make_knn_struct(self, too_big=False, big_cluster=None):
        if self.knn > 190: print('consider using a lower K_in for KNN graph construction')
        ef_query = max(100, self.knn + 1)  # ef always should be >K. higher ef, more accurate query
        if too_big == False:
            num_dims = self.data.shape[1]
            n_elements = self.data.shape[0]
            p = hnswlib.Index(space=self.distance, dim=num_dims)  # default to Euclidean distance
            p.set_num_threads(self.num_threads)  # allow user to set threads used in KNN construction
            if n_elements < 10000:
                ef_query = min(n_elements - 10, 500)
                ef_construction = ef_query
            else:
                ef_construction = self.hnsw_param_ef_construction
            if (num_dims > 30) & (n_elements <= 50000):
                p.init_index(max_elements=n_elements, ef_construction=ef_construction,
                             M=48)  ## good for scRNA seq where dimensionality is high
            else:
                p.init_index(max_elements=n_elements, ef_construction=ef_construction, M=24)  # 30
            p.add_items(self.data)
        if too_big == True:
            num_dims = big_cluster.shape[1]
            n_elements = big_cluster.shape[0]
            p = hnswlib.Index(space='l2', dim=num_dims)
            p.init_index(max_elements=n_elements, ef_construction=200, M=30)
            p.add_items(big_cluster)
        p.set_ef(ef_query)  # ef should always be > k

        return p

    def knngraph_full(self):  # , neighbor_array, distance_array):
        k_umap = 15
        t0 = time.time()
        # neighbors in array are not listed in in any order of proximity
        self.knn_struct.set_ef(k_umap + 1)
        neighbor_array, distance_array = self.knn_struct.knn_query(self.data, k=k_umap)

        row_list = []
        n_neighbors = neighbor_array.shape[1]
        n_cells = neighbor_array.shape[0]

        row_list.extend(list(np.transpose(np.ones((n_neighbors, n_cells)) * range(0, n_cells)).flatten()))

        row_min = np.min(distance_array, axis=1)
        row_sigma = np.std(distance_array, axis=1)

        distance_array = (distance_array - row_min[:, np.newaxis]) / row_sigma[:, np.newaxis]

        col_list = neighbor_array.flatten().tolist()
        distance_array = distance_array.flatten()
        distance_array = np.sqrt(distance_array)
        distance_array = distance_array * -1

        weight_list = np.exp(distance_array)

        threshold = np.mean(weight_list) + 2 * np.std(weight_list)

        weight_list[weight_list >= threshold] = threshold

        weight_list = weight_list.tolist()

        graph = csr_matrix((np.array(weight_list), (np.array(row_list), np.array(col_list))),
                           shape=(n_cells, n_cells))

        graph_transpose = graph.T
        prod_matrix = graph.multiply(graph_transpose)

        graph = graph_transpose + graph - prod_matrix
        return graph

    def make_csrmatrix_noselfloop(self, neighbor_array, distance_array):
        # neighbor array not listed in in any order of proximity
        row_list = []
        col_list = []
        weight_list = []

        n_neighbors = neighbor_array.shape[1]
        n_cells = neighbor_array.shape[0]
        rowi = 0
        discard_count = 0
        if self.keep_all_local_dist == False:  # locally prune based on (squared) l2 distance

            print('commencing local pruning based on Euclidean distance metric at',
                  self.dist_std_local, 's.dev above mean')
            distance_array = distance_array + 0.1
            for row in neighbor_array:
                distlist = distance_array[rowi, :]
                to_keep = np.where(distlist < np.mean(distlist) + self.dist_std_local * np.std(distlist))[0]  # 0*std
                updated_nn_ind = row[np.ix_(to_keep)]
                updated_nn_weights = distlist[np.ix_(to_keep)]
                discard_count = discard_count + (n_neighbors - len(to_keep))

                for ik in range(len(updated_nn_ind)):
                    if rowi != row[ik]:  # remove self-loops
                        row_list.append(rowi)
                        col_list.append(updated_nn_ind[ik])
                        dist = np.sqrt(updated_nn_weights[ik])
                        weight_list.append(1 / (dist + 0.1))

                rowi = rowi + 1

        if self.keep_all_local_dist == True:  # dont prune based on distance
            row_list.extend(list(np.transpose(np.ones((n_neighbors, n_cells)) * range(0, n_cells)).flatten()))
            col_list = neighbor_array.flatten().tolist()
            weight_list = (1. / (distance_array.flatten() + 0.1)).tolist()

        csr_graph = csr_matrix((np.array(weight_list), (np.array(row_list), np.array(col_list))),
                               shape=(n_cells, n_cells))
        return csr_graph

    def func_mode(self, ll):  # return MODE of list
        # If multiple items are maximal, the function returns the first one encountered.
        return max(set(ll), key=ll.count)

    def run_toobig_subPARC(self, X_data, jac_std_toobig=0.3,
                           jac_weighted_edges=True):
        n_elements = X_data.shape[0]
        hnsw = self.make_knn_struct(too_big=True, big_cluster=X_data)
        if n_elements <= 10: print('consider increasing the too_big_factor')
        if n_elements > self.knn:
            knnbig = self.knn
        else:
            knnbig = int(max(5, 0.2 * n_elements))

        neighbor_array, distance_array = hnsw.knn_query(X_data, k=knnbig)
        # print('shapes of neigh and dist array', neighbor_array.shape, distance_array.shape)
        csr_array = self.make_csrmatrix_noselfloop(neighbor_array, distance_array)
        sources, targets = csr_array.nonzero()

        #mask = np.zeros(len(sources), dtype=bool)

        #mask |= (csr_array.data > (        np.mean(csr_array.data) + np.std(csr_array.data) * 5))  # weights are set as 1/dist. so smaller distance means stronger edge

        #csr_array.data[mask] = 0
        #csr_array.eliminate_zeros()
        #sources, targets = csr_array.nonzero()
        edgelist = list(zip(sources.tolist(), targets.tolist()))
        edgelist_copy = edgelist.copy()
        G = ig.Graph(edgelist, edge_attrs={'weight': csr_array.data.tolist()})
        sim_list = G.similarity_jaccard(pairs=edgelist_copy)  # list of jaccard weights
        new_edgelist = []
        sim_list_array = np.asarray(sim_list)
        if jac_std_toobig == 'median':
            threshold = np.median(sim_list)
        else:
            threshold = np.mean(sim_list) - jac_std_toobig * np.std(sim_list)
        print('jac threshold %.3f' % threshold)
        print('jac std %.3f' % np.std(sim_list))
        print('jac mean %.3f' % np.mean(sim_list))
        strong_locs = np.where(sim_list_array > threshold)[0]
        for ii in strong_locs: new_edgelist.append(edgelist_copy[ii])
        sim_list_new = list(sim_list_array[strong_locs])

        if jac_weighted_edges == True:
            G_sim = ig.Graph(n=n_elements, edges=list(new_edgelist), edge_attrs={'weight': sim_list_new})
        else:
            G_sim = ig.Graph(n=n_elements, edges=list(new_edgelist))
        G_sim.simplify(combine_edges='sum')
        import random
        random.seed(self.random_seed)
        if jac_weighted_edges == True:
            if self.partition_type == 'ModularityVP':
                partition = leidenalg.find_partition(G_sim, leidenalg.ModularityVertexPartition, weights='weight',
                                                     n_iterations=self.n_iter_leiden, seed=self.random_seed)
                print('partition type MVP')
            else:
                # partition = leidenalg.find_partition(G_sim, leidenalg.RBConfigurationVertexPartition, weights='weight',
                #                                      n_iterations=self.n_iter_leiden, seed=self.random_seed,
                #                                      resolution_parameter=self.resolution_parameter)
                print("custom partition")
                print(self.resolution_parameter)
                print(self.n_iter_leiden)
                partition = G_sim.community_leiden(objective_function='modularity', weights='weight', n_iterations=self.n_iter_leiden, 
                                                   resolution=self.resolution_parameter)
                print('partition type RBC')
        else:
            if self.partition_type == 'ModularityVP':
                print('partition type MVP')
                partition = leidenalg.find_partition(G_sim, leidenalg.ModularityVertexPartition,
                                                     n_iterations=self.n_iter_leiden, seed=self.random_seed)
            else:
                # print('partition type RBC')
                # partition = leidenalg.find_partition(G_sim, leidenalg.RBConfigurationVertexPartition,
                #                                      n_iterations=self.n_iter_leiden, seed=self.random_seed,
                #                                      resolution_parameter=self.resolution_parameter)
                print("custom partition")
                print(self.resolution_parameter)
                print(self.n_iter_leiden)
                partition = G_sim.community_leiden(objective_function='modularity', n_iterations=self.n_iter_leiden, 
                                                   resolution=self.resolution_parameter)
        # print('Q= %.2f' % partition.quality())
        PARC_labels_leiden = np.asarray(partition.membership)
        PARC_labels_leiden = np.reshape(PARC_labels_leiden, (n_elements, 1))
        small_pop_list = []
        small_cluster_list = []
        small_pop_exist = False
        dummy, PARC_labels_leiden = np.unique(list(PARC_labels_leiden.flatten()), return_inverse=True)
        for cluster in set(PARC_labels_leiden):
            population = len(np.where(PARC_labels_leiden == cluster)[0])
            if population < 10:
                small_pop_exist = True
                small_pop_list.append(list(np.where(PARC_labels_leiden == cluster)[0]))
                small_cluster_list.append(cluster)

        for small_cluster in small_pop_list:
            for single_cell in small_cluster:
                old_neighbors = neighbor_array[single_cell]
                group_of_old_neighbors = PARC_labels_leiden[old_neighbors]
                group_of_old_neighbors = list(group_of_old_neighbors.flatten())
                available_neighbours = set(group_of_old_neighbors) - set(small_cluster_list)
                if len(available_neighbours) > 0:
                    available_neighbours_list = [value for value in group_of_old_neighbors if
                                                 value in list(available_neighbours)]
                    best_group = max(available_neighbours_list, key=available_neighbours_list.count)
                    PARC_labels_leiden[single_cell] = best_group

        time_smallpop_start = time.time()
        print('handling fragments')
        while (small_pop_exist) == True & (time.time() - time_smallpop_start < self.time_smallpop):
            small_pop_list = []
            small_pop_exist = False
            for cluster in set(list(PARC_labels_leiden.flatten())):
                population = len(np.where(PARC_labels_leiden == cluster)[0])
                if population < 10:
                    small_pop_exist = True

                    small_pop_list.append(np.where(PARC_labels_leiden == cluster)[0])
            for small_cluster in small_pop_list:
                for single_cell in small_cluster:
                    old_neighbors = neighbor_array[single_cell, :]
                    group_of_old_neighbors = PARC_labels_leiden[old_neighbors]
                    group_of_old_neighbors = list(group_of_old_neighbors.flatten())
                    best_group = max(set(group_of_old_neighbors), key=group_of_old_neighbors.count)
                    PARC_labels_leiden[single_cell] = best_group

        dummy, PARC_labels_leiden = np.unique(list(PARC_labels_leiden.flatten()), return_inverse=True)

        return PARC_labels_leiden

    def run_subPARC(self):

        X_data = self.data
        too_big_factor = self.too_big_factor
        small_pop = self.small_pop
        jac_std_global = self.jac_std_global
        jac_weighted_edges = self.jac_weighted_edges
        knn = self.knn
        n_elements = X_data.shape[0]

        if self.neighbor_graph is not None:
            csr_array = self.neighbor_graph
            neighbor_array = np.split(csr_array.indices, csr_array.indptr)[1:-1]
        else:
            if self.knn_struct is None:
                print('knn struct was not available, so making one')
                self.knn_struct = self.make_knn_struct()
            else:
                print('knn struct already exists')
            neighbor_array, distance_array = self.knn_struct.knn_query(X_data, k=knn)
            csr_array = self.make_csrmatrix_noselfloop(neighbor_array, distance_array)

        sources, targets = csr_array.nonzero()

        edgelist = list(zip(sources, targets))

        edgelist_copy = edgelist.copy()

        G = ig.Graph(edgelist, edge_attrs={'weight': csr_array.data.tolist()})
        # print('average degree of prejacard graph is %.1f'% (np.mean(G.degree())))
        # print('computing Jaccard metric')
        sim_list = G.similarity_jaccard(pairs=edgelist_copy)

        print('commencing global pruning')

        sim_list_array = np.asarray(sim_list)
        edge_list_copy_array = np.asarray(edgelist_copy)

        if jac_std_global == 'median':
            threshold = np.median(sim_list)
        else:
            threshold = np.mean(sim_list) - jac_std_global * np.std(sim_list)
        strong_locs = np.where(sim_list_array > threshold)[0]
        # print('Share of edges kept after Global Pruning %.2f' % (len(strong_locs) / len(sim_list)), '%')
        new_edgelist = list(edge_list_copy_array[strong_locs])
        sim_list_new = list(sim_list_array[strong_locs])

        G_sim = ig.Graph(n=n_elements, edges=list(new_edgelist), edge_attrs={'weight': sim_list_new})
        # print('average degree of graph is %.1f' % (np.mean(G_sim.degree())))
        G_sim.simplify(combine_edges='sum')  # "first"
        # print('average degree of SIMPLE graph is %.1f' % (np.mean(G_sim.degree())))
        print('commencing community detection')
        import random
        random.seed(self.random_seed)
        if jac_weighted_edges == True:
            start_leiden = time.time()
            if self.partition_type == 'ModularityVP':
                print('partition type MVP')
                partition = leidenalg.find_partition(G_sim, leidenalg.ModularityVertexPartition, weights='weight',
                                                     n_iterations=self.n_iter_leiden, seed=self.random_seed)
            else:
                print('partition type RBC')
                # partition = leidenalg.find_partition(G_sim, leidenalg.RBConfigurationVertexPartition, weights='weight',
                #                                      n_iterations=self.n_iter_leiden, seed=self.random_seed,
                #                                      resolution_parameter=self.resolution_parameter)
                print("custom partition")
                print(self.resolution_parameter)
                print(self.n_iter_leiden)
                partition = G_sim.community_leiden(objective_function='modularity',weights='weight', n_iterations=self.n_iter_leiden, 
                                                   resolution=self.resolution_parameter)
            # print(time.time() - start_leiden)
        else:
            start_leiden = time.time()
            if self.partition_type == 'ModularityVP':
                partition = leidenalg.find_partition(G_sim, leidenalg.ModularityVertexPartition,
                                                     n_iterations=self.n_iter_leiden, seed=self.random_seed)
                print('partition type MVP')
            else:
                # partition = leidenalg.find_partition(G_sim, leidenalg.RBConfigurationVertexPartition,
                #                                      n_iterations=self.n_iter_leiden, seed=self.random_seed,
                #                                      resolution_parameter=self.resolution_parameter)
                print("custom partition")
                print(self.resolution_parameter)
                print(self.n_iter_leiden)
                partition = G_sim.community_leiden(objective_function='modularity', n_iterations=self.n_iter_leiden, 
                                                   resolution=self.resolution_parameter)
                
                #print('partition type RBC')
            # print(time.time() - start_leiden)
        time_end_PARC = time.time()
        # print('Q= %.1f' % (partition.quality()))
        PARC_labels_leiden = np.asarray(partition.membership)
        PARC_labels_leiden = np.reshape(PARC_labels_leiden, (n_elements, 1))

        too_big = False

        # print('labels found after Leiden', set(list(PARC_labels_leiden.T)[0])) will have some outlier clusters that need to be added to a cluster if a cluster has members that are KNN

        cluster_i_loc = np.where(PARC_labels_leiden == 0)[
            0]  # the 0th cluster is the largest one. so if cluster 0 is not too big, then the others wont be too big either
        pop_i = len(cluster_i_loc)
        if pop_i > too_big_factor * n_elements:  # 0.4
            too_big = True
            cluster_big_loc = cluster_i_loc
            list_pop_too_bigs = [pop_i]
            cluster_too_big = 0

        while too_big == True:

            X_data_big = X_data[cluster_big_loc, :]
            PARC_labels_leiden_big = self.run_toobig_subPARC(X_data_big)
            # print('set of new big labels ', set(PARC_labels_leiden_big.flatten()))
            PARC_labels_leiden_big = PARC_labels_leiden_big + 100000
            # print('set of new big labels +1000 ', set(list(PARC_labels_leiden_big.flatten())))
            pop_list = []

            for item in set(list(PARC_labels_leiden_big.flatten())):
                pop_list.append([item, list(PARC_labels_leiden_big.flatten()).count(item)])
            print('pop of big clusters', pop_list)
            jj = 0
            print('shape PARC_labels_leiden', PARC_labels_leiden.shape)
            for j in cluster_big_loc:
                PARC_labels_leiden[j] = PARC_labels_leiden_big[jj]
                jj = jj + 1
            dummy, PARC_labels_leiden = np.unique(list(PARC_labels_leiden.flatten()), return_inverse=True)
            print('new set of labels ', set(PARC_labels_leiden))
            too_big = False
            set_PARC_labels_leiden = set(PARC_labels_leiden)

            PARC_labels_leiden = np.asarray(PARC_labels_leiden)
            for cluster_ii in set_PARC_labels_leiden:
                cluster_ii_loc = np.where(PARC_labels_leiden == cluster_ii)[0]
                pop_ii = len(cluster_ii_loc)
                not_yet_expanded = pop_ii not in list_pop_too_bigs
                if pop_ii > too_big_factor * n_elements and not_yet_expanded == True:
                    too_big = True
                    print('cluster', cluster_ii, 'is too big and has population', pop_ii)
                    cluster_big_loc = cluster_ii_loc
                    cluster_big = cluster_ii
                    big_pop = pop_ii
            if too_big == True:
                list_pop_too_bigs.append(big_pop)
                print('cluster', cluster_big, 'is too big with population', big_pop, '. It will be expanded')
        dummy, PARC_labels_leiden = np.unique(list(PARC_labels_leiden.flatten()), return_inverse=True)
        small_pop_list = []
        small_cluster_list = []
        small_pop_exist = False

        for cluster in set(PARC_labels_leiden):
            population = len(np.where(PARC_labels_leiden == cluster)[0])

            if population < small_pop:  # 10
                small_pop_exist = True

                small_pop_list.append(list(np.where(PARC_labels_leiden == cluster)[0]))
                small_cluster_list.append(cluster)

        for small_cluster in small_pop_list:

            for single_cell in small_cluster:
                old_neighbors = neighbor_array[single_cell]
                group_of_old_neighbors = PARC_labels_leiden[old_neighbors]
                group_of_old_neighbors = list(group_of_old_neighbors.flatten())
                available_neighbours = set(group_of_old_neighbors) - set(small_cluster_list)
                if len(available_neighbours) > 0:
                    available_neighbours_list = [value for value in group_of_old_neighbors if
                                                 value in list(available_neighbours)]
                    best_group = max(available_neighbours_list, key=available_neighbours_list.count)
                    PARC_labels_leiden[single_cell] = best_group
        time_smallpop_start = time.time()
        while (small_pop_exist == True) & ((time.time() - time_smallpop_start) < self.time_smallpop):
            small_pop_list = []
            small_pop_exist = False
            for cluster in set(list(PARC_labels_leiden.flatten())):
                population = len(np.where(PARC_labels_leiden == cluster)[0])
                if population < small_pop:
                    small_pop_exist = True
                    print(cluster, ' has small population of', population, )
                    small_pop_list.append(np.where(PARC_labels_leiden == cluster)[0])
            for small_cluster in small_pop_list:
                for single_cell in small_cluster:
                    old_neighbors = neighbor_array[single_cell, :]
                    group_of_old_neighbors = PARC_labels_leiden[old_neighbors]
                    group_of_old_neighbors = list(group_of_old_neighbors.flatten())
                    best_group = max(set(group_of_old_neighbors), key=group_of_old_neighbors.count)
                    PARC_labels_leiden[single_cell] = best_group

        dummy, PARC_labels_leiden = np.unique(list(PARC_labels_leiden.flatten()), return_inverse=True)
        PARC_labels_leiden = list(PARC_labels_leiden.flatten())
        # print('final labels allocation', set(PARC_labels_leiden))
        pop_list = []
        for item in set(PARC_labels_leiden):
            pop_list.append((item, PARC_labels_leiden.count(item)))
        print('list of cluster labels and populations', len(pop_list), pop_list)

        self.labels = PARC_labels_leiden  # list
        return

    def accuracy(self, onevsall=1):

        true_labels = self.true_label
        Index_dict = {}
        PARC_labels = self.labels
        N = len(PARC_labels)
        n_cancer = list(true_labels).count(onevsall)
        n_pbmc = N - n_cancer

        for k in range(N):
            Index_dict.setdefault(PARC_labels[k], []).append(true_labels[k])
        num_groups = len(Index_dict)
        sorted_keys = list(sorted(Index_dict.keys()))
        error_count = []
        pbmc_labels = []
        thp1_labels = []
        fp, fn, tp, tn, precision, recall, f1_score = 0, 0, 0, 0, 0, 0, 0

        for kk in sorted_keys:
            vals = [t for t in Index_dict[kk]]
            majority_val = self.func_mode(vals)
            if majority_val == onevsall: print('cluster', kk, ' has majority', onevsall, 'with population', len(vals))
            if kk == -1:
                len_unknown = len(vals)
                print('len unknown', len_unknown)
            if (majority_val == onevsall) and (kk != -1):
                thp1_labels.append(kk)
                fp = fp + len([e for e in vals if e != onevsall])
                tp = tp + len([e for e in vals if e == onevsall])
                list_error = [e for e in vals if e != majority_val]
                e_count = len(list_error)
                error_count.append(e_count)
            elif (majority_val != onevsall) and (kk != -1):
                pbmc_labels.append(kk)
                tn = tn + len([e for e in vals if e != onevsall])
                fn = fn + len([e for e in vals if e == onevsall])
                error_count.append(len([e for e in vals if e != majority_val]))

        predict_class_array = np.array(PARC_labels)
        PARC_labels_array = np.array(PARC_labels)
        number_clusters_for_target = len(thp1_labels)
        for cancer_class in thp1_labels:
            predict_class_array[PARC_labels_array == cancer_class] = 1
        for benign_class in pbmc_labels:
            predict_class_array[PARC_labels_array == benign_class] = 0
        predict_class_array.reshape((predict_class_array.shape[0], -1))
        error_rate = sum(error_count) / N
        n_target = tp + fn
        tnr = tn / n_pbmc
        fnr = fn / n_cancer
        tpr = tp / n_cancer
        fpr = fp / n_pbmc

        if tp != 0 or fn != 0: recall = tp / (tp + fn)  # ability to find all positives
        if tp != 0 or fp != 0: precision = tp / (tp + fp)  # ability to not misclassify negatives as positives
        if precision != 0 or recall != 0:
            f1_score = precision * recall * 2 / (precision + recall)
        majority_truth_labels = np.empty((len(true_labels), 1), dtype=object)

        for cluster_i in set(PARC_labels):
            cluster_i_loc = np.where(np.asarray(PARC_labels) == cluster_i)[0]
            true_labels = np.asarray(true_labels)
            majority_truth = self.func_mode(list(true_labels[cluster_i_loc]))
            majority_truth_labels[cluster_i_loc] = majority_truth

        majority_truth_labels = list(majority_truth_labels.flatten())
        accuracy_val = [error_rate, f1_score, tnr, fnr, tpr, fpr, precision,
                        recall, num_groups, n_target]

        return accuracy_val, predict_class_array, majority_truth_labels, number_clusters_for_target

    def run_PARC(self):
        print('input data has shape', self.data.shape[0], '(samples) x', self.data.shape[1], '(features)')
        if self.true_label is None:
            self.true_label = [1] * self.data.shape[0]
        list_roc = []

        time_start_total = time.time()

        time_start_knn = time.time()

        time_end_knn_struct = time.time() - time_start_knn
        # Query dataset, k - number of closest elements (returns 2 numpy arrays)
        self.run_subPARC()
        run_time = time.time() - time_start_total
        print('time elapsed {:.1f} seconds'.format(run_time))

        targets = list(set(self.true_label))
        N = len(list(self.true_label))
        self.f1_accumulated = 0
        self.f1_mean = 0
        self.stats_df = pd.DataFrame({'jac_std_global': [self.jac_std_global], 'dist_std_local': [self.dist_std_local],
                                      'runtime(s)': [run_time]})
        self.majority_truth_labels = []
        if len(targets) > 1:
            f1_accumulated = 0
            f1_acc_noweighting = 0
            for onevsall_val in targets:
                print('target is', onevsall_val)
                vals_roc, predict_class_array, majority_truth_labels, numclusters_targetval = self.accuracy(
                    onevsall=onevsall_val)
                f1_current = vals_roc[1]
                print('target', onevsall_val, 'has f1-score of %.2f' % (f1_current * 100))
                f1_accumulated = f1_accumulated + f1_current * (list(self.true_label).count(onevsall_val)) / N
                f1_acc_noweighting = f1_acc_noweighting + f1_current

                list_roc.append(
                    [self.jac_std_global, self.dist_std_local, onevsall_val] + vals_roc + [numclusters_targetval] + [
                        run_time])

            f1_mean = f1_acc_noweighting / len(targets)
            print("f1-score (unweighted) mean %.2f" % (f1_mean * 100), '%')
            print('f1-score weighted (by population) %.2f' % (f1_accumulated * 100), '%')

            df_accuracy = pd.DataFrame(list_roc,
                                       columns=['jac_std_global', 'dist_std_local', 'onevsall-target', 'error rate',
                                                'f1-score', 'tnr', 'fnr',
                                                'tpr', 'fpr', 'precision', 'recall', 'num_groups',
                                                'population of target', 'num clusters', 'clustering runtime'])

            self.f1_accumulated = f1_accumulated
            self.f1_mean = f1_mean
            self.stats_df = df_accuracy
            self.majority_truth_labels = majority_truth_labels
        return

    def run_umap_hnsw(self, X_input, graph, n_components=2, alpha: float = 1.0, negative_sample_rate: int = 5,
                      gamma: float = 1.0, spread=1.0, min_dist=0.1, init_pos='spectral', random_state=1, ):

        from umap.umap_ import find_ab_params, simplicial_set_embedding
        import matplotlib.pyplot as plt

        a, b = find_ab_params(spread, min_dist)
        print('a,b, spread, dist', a, b, spread, min_dist)
        t0 = time.time()
        X_umap = simplicial_set_embedding(data=X_input, graph=graph, n_components=n_components, initial_alpha=alpha,
                                          a=a, b=b, n_epochs=0, metric_kwds={}, gamma=gamma,
                                          negative_sample_rate=negative_sample_rate, init=init_pos,
                                          random_state=np.random.RandomState(random_state), metric='euclidean',
                                          verbose=1)
        return X_umap

def knngraph_full_2(self):  # , neighbor_array, distance_array):
        k_umap = 15
        t0 = time.time()
        # neighbors in array are not listed in in any order of proximity
        self.knn_struct.set_ef(k_umap + 1)
        neighbor_array, distance_array = self.knn_struct.knn_query(self.data, k=k_umap)

        row_list = []
        n_neighbors = neighbor_array.shape[1]
        n_cells = neighbor_array.shape[0]

        row_list.extend(list(np.transpose(np.ones((n_neighbors, n_cells)) * range(0, n_cells)).flatten()))

        row_min = np.min(distance_array, axis=1)
        row_sigma = np.std(distance_array, axis=1)

        distance_array = (distance_array - row_min[:, np.newaxis]) / row_sigma[:, np.newaxis]

        col_list = neighbor_array.flatten().tolist()
        distance_array = distance_array.flatten()
        distance_array = np.sqrt(distance_array)
        distance_array = distance_array * -1

        weight_list = np.exp(distance_array)

        threshold = np.mean(weight_list) + 2 * np.std(weight_list)

        weight_list[weight_list >= threshold] = threshold

        weight_list = weight_list.tolist()

        graph = csr_matrix((np.array(weight_list), (np.array(row_list), np.array(col_list))),
                           shape=(n_cells, n_cells))

        graph_transpose = graph.T
        prod_matrix = graph.multiply(graph_transpose)

        #graph = graph_transpose + graph - prod_matrix
        return graph , prod_matrix

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
        else:
            sc.pp.neighbors(adata, n_neighbors=n_neighbors, metric=metric, n_pcs=0, random_state=random_state)
            
    sc.tl.leiden(adata,resolution=resolution, random_state=random_state, n_iterations=n_iterations, flavor="igraph")
    adata.obs['Cluster'] = adata.obs['leiden']
    #sc.tl.umap(adata)
    
    adata.obsm['spatial'] = np.array(adata.obs[["Centroid X (µm)_(standardized)", "Centroid Y (µm)_(standardized)"]])
    return adata

# phenograph clustering
def RunPhenographClust(adata, n_neighbors, clustering_algo, min_cluster_size, 
                       primary_metric, resolution_parameter, nn_method, random_seed, n_principal_components,
                       n_jobs, n_iterations, fast):
    print(adata.shape)
    if fast == True:
        print("fast phenograph selected")
        sc.pp.pca(adata, n_comps=n_principal_components)
        print("PCA done")
        p1 = PARC(adata.obsm["X_pca"],  keep_all_local_dist=True, num_threads=n_jobs)  # without labels
        print("Parc object created")
        p1.knn_struct = p1.make_knn_struct()
        print("Parc knn struct created")
        graph_parc_2, mat_parc_2 = knngraph_full_2(p1)
        print("Parc graph prepared")
        communities, graph_phen, Q = phenograph.cluster(graph_parc_2, clustering_algo=None, k=n_neighbors, 
                                                min_cluster_size=min_cluster_size, primary_metric=primary_metric, 
                                                resolution_parameter=resolution_parameter, nn_method=nn_method,
                                                seed=random_seed, n_iterations=n_iterations, n_jobs=n_jobs)
        adata.obsp["pheno_jaccard_ig"] = graph_phen.tocsr()
        print("start_leiden")
        sc.tl.leiden(adata,resolution=resolution_parameter, random_state=random_seed,
                    n_iterations=n_iterations, flavor="igraph", obsp="pheno_jaccard_ig")
        adata.obs['Cluster'] = adata.obs['leiden'].astype(str)
        print("Fast phenograph clustering done")
    else:
        if n_principal_components == 0:
            communities, graph, Q = phenograph.cluster(adata.X, clustering_algo=clustering_algo, k=n_neighbors, 
                                                    min_cluster_size=min_cluster_size, primary_metric=primary_metric, 
                                                    resolution_parameter=resolution_parameter, nn_method=nn_method,
                                                    seed=random_seed, n_iterations=n_iterations, n_jobs=n_jobs)
        else:
            sc.pp.pca(adata, n_comps=n_principal_components)
            communities, graph, Q = phenograph.cluster(adata.obsm['X_pca'], clustering_algo=clustering_algo, k=n_neighbors, 
                                                    min_cluster_size=min_cluster_size, primary_metric=primary_metric, 
                                                    resolution_parameter=resolution_parameter, nn_method=nn_method,
                                                    seed=random_seed, n_iterations=n_iterations, n_jobs=n_jobs)
        adata.obs['Cluster'] = communities
        adata.obs['Cluster'] = adata.obs['Cluster'].astype(str)
        print("Regular phenograph clustering done")
    #sc.tl.umap(adata)
    adata.obsm['spatial'] = np.array(adata.obs[["Centroid X (µm)_(standardized)", "Centroid Y (µm)_(standardized)"]])
    return adata

# parc clustering
def run_parc_clust(adata, n_neighbors, dist_std_local, jac_std_global, small_pop, 
                   random_seed, resolution_parameter, hnsw_param_ef_construction, 
                   n_principal_components, n_iterations, n_jobs, fast):
    if fast == True:
        sc.pp.pca(adata, n_comps=n_principal_components)
        parc_results = PARC_2(adata.obsm["X_pca"], dist_std_local=dist_std_local, jac_std_global=jac_std_global, 
                            small_pop=small_pop, random_seed=random_seed, knn=n_neighbors,
                            resolution_parameter=resolution_parameter, 
                            hnsw_param_ef_construction=hnsw_param_ef_construction,
                            partition_type="RBConfigurationVP",
                            n_iter_leiden=n_iterations, num_threads=n_jobs)  # without labels
        parc_results.run_PARC()
        adata.obs['Cluster'] = parc_results.labels
        adata.obs['Cluster'] = adata.obs['Cluster'].astype(str)
    else: 
        if n_principal_components == 0:
            parc_results = parc.PARC(adata.X, dist_std_local=dist_std_local, jac_std_global=jac_std_global, 
                                    small_pop=small_pop, random_seed=random_seed, knn=n_neighbors,
                                    resolution_parameter=resolution_parameter, 
                                    hnsw_param_ef_construction=hnsw_param_ef_construction,
                                    partition_type="RBConfigurationVP",
                                    n_iter_leiden=n_iterations, num_threads=n_jobs)
        else:
            sc.pp.pca(adata, n_comps=n_principal_components)
            parc_results = parc.PARC(adata.obsm['X_pca'], dist_std_local=dist_std_local, jac_std_global=jac_std_global, 
                                    small_pop=small_pop, random_seed=random_seed, knn=n_neighbors,
                                    resolution_parameter=resolution_parameter, 
                                    hnsw_param_ef_construction=hnsw_param_ef_construction,
                                    partition_type="RBConfigurationVP",
                                    n_iter_leiden=n_iterations, num_threads=n_jobs)
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
                )
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
    
    try:
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
    except:
        st.warning("container issue")
        
    try:
            
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
                if st.button('Run Clustering'):
                    start_time = time.time()
                    if st.session_state['phenocluster__cluster_method'] == "phenograph":
                        with st.spinner('Wait for it...'):
                            st.session_state['phenocluster__clustering_adata'] = RunPhenographClust(adata=st.session_state['phenocluster__clustering_adata'], 
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
                                                                                                n_principal_components=st.session_state['phenocluster__n_principal_components'],
                                                                                                n_jobs=st.session_state['phenocluster__n_jobs'],
                                                                                                n_iterations= st.session_state['phenocluster__n_iterations'],
                                                                                                fast=st.session_state["phenocluster__fast"]
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
            
                
    except:
        st.warning("container issue")
                
                
            
    

# Run the main function
if __name__ == '__main__':

    # Set page settings
    page_name = 'Unsupervised Phenotype Clustering'
    #st.set_page_config(layout='wide', page_title=page_name)
    st.title(page_name)
    
    # Run streamlit-dataframe-editor library initialization tasks at the top of the page
    st.session_state = sde.initialize_session_state(st.session_state)

    # Run Top of Page (TOP) functions
    st.session_state = top.top_of_page_reqs(st.session_state)

    # Call the main function
    main()

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)
    
# need to make differential expression on another page 
