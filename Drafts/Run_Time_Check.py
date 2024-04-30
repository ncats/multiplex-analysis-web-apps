import pandas as pd 
import scanpy as sc
import phenograph
import parc
from utag import utag
import os 
import numpy as np
import time
import anndata as ad

times_list = []

def timer(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        execution_time = end_time - start_time
        rounded_time = round(execution_time, 2)
        print(f'Execution time: {rounded_time} seconds')
        times_list.append(rounded_time)
        return result
    return wrapper

os.chdir("/Users/bombina2/github/multiplex-analysis-web-apps/")

# phenograph clustering
@timer
def RunPhenographClust(adata, n_neighbors, clustering_algo, min_cluster_size, primary_metric, resolution_parameter, nn_method):
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=0)
    communities, graph, Q = phenograph.cluster(adata.X, clustering_algo=clustering_algo, k=n_neighbors, 
                                               min_cluster_size=min_cluster_size, primary_metric=primary_metric, 
                                               resolution_parameter=resolution_parameter, nn_method=nn_method,
                                               seed=42, n_iterations=-1)
    adata.obs['Cluster'] = communities
    adata.obs['Cluster'] = adata.obs['Cluster'].astype(str)
    sc.tl.umap(adata)
    adata.obsm['spatial'] = np.array(adata.obs[["Centroid X (µm)_(standardized)", "Centroid Y (µm)_(standardized)"]])
    return adata

# scanpy clustering
@timer
def RunNeighbClust(adata, n_neighbors, metric, resolution, random_state):
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, metric=metric, n_pcs=0)
    sc.tl.leiden(adata,resolution=resolution, random_state=random_state, n_iterations=-1, flavor="igraph")
    adata.obs['Cluster'] = adata.obs['leiden']
    sc.tl.umap(adata)
    adata.obsm['spatial'] = np.array(adata.obs[["Centroid X (µm)_(standardized)", "Centroid Y (µm)_(standardized)"]])
    return adata

# parc clustering
@timer
def run_parc_clust(adata, n_neighbors, dist_std_local, jac_std_global, small_pop, random_seed, resolution_parameter, hnsw_param_ef_construction):
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=0)
    parc_results = parc.PARC(adata.X, dist_std_local=dist_std_local, jac_std_global=jac_std_global, 
                             small_pop=small_pop, random_seed=random_seed, knn=n_neighbors,
                             resolution_parameter=resolution_parameter, 
                             hnsw_param_ef_construction=hnsw_param_ef_construction,
                             partition_type="RBConfigurationVP",
                             n_iter_leiden=-1)
    parc_results.run_PARC()
    adata.obs['Cluster'] = parc_results.labels
    adata.obs['Cluster'] = adata.obs['Cluster'].astype(str)
    sc.tl.umap(adata)
    adata.obsm['spatial'] = np.array(adata.obs[["Centroid X (µm)_(standardized)", "Centroid Y (µm)_(standardized)"]])
    return adata

# utag clustering
# need to make image selection based on the variable
@timer
def run_utag_clust(adata, n_neighbors, resolutions, clustering_method, max_dist):
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=0)
    sc.tl.umap(adata)
    adata.obsm['spatial'] = np.array(adata.obs[["Centroid X (µm)_(standardized)", "Centroid Y (µm)_(standardized)"]])
    utag_results = utag(adata,
        slide_key="Image ID_(standardized)",
        max_dist=max_dist,
        normalization_mode='l1_norm',
        apply_clustering=True,
        clustering_method = clustering_method, 
        resolutions = resolutions,
        leiden_kwargs={"n_iterations": -1, "random_state": 42}
    )
     
    curClusterCol = 'UTAG Label_leiden_'  + str(resolutions[0])
    utag_results.obs['Cluster'] = utag_results.obs[curClusterCol]
    adata.obs['Cluster'] = utag_results.obs[curClusterCol]
        
    curClusterCol = 'UTAG Label_leiden_'  + str(resolutions[0])
    utag_results.obs['Cluster'] = utag_results.obs[curClusterCol]
    adata.obs['Cluster'] = utag_results.obs[curClusterCol]
    return utag_results

## test with more than a million cells 
# df = pd.read_csv("input/measurementsthymus-exported.csv")
# list(df.columns)
# meta = df.iloc[:, :4]
# list(meta.columns)
# mat = df.iloc[:, 4:]
# list(mat.columns)
# meta = meta.rename(columns={"Image": "Image ID_(standardized)", 
#                             "Centroid X µm": "Centroid X (µm)_(standardized)", 
#                             "Centroid Y µm": "Centroid Y (µm)_(standardized)"})

# adata = ad.AnnData(mat)
# adata.obs = meta
# adata.layers["counts"] = adata.X.copy()

# test phenograph
# adata = RunPhenographClust(adata, n_neighbors=10, clustering_algo="leiden", 
#                            min_cluster_size=10, primary_metric="euclidean", 
#                            resolution_parameter=1, nn_method="kdtree")

# # test scanpy
# adata = RunNeighbClust(adata, n_neighbors=10, metric="euclidean", 
#                        resolution=1, random_state=42)

# # test parc
# adata = run_parc_clust(adata, n_neighbors=10, 
#                        dist_std_local=3, jac_std_global=0.15, 
#                        small_pop=50, random_seed=42, 
#                        resolution_parameter=1, hnsw_param_ef_construction=150)

# # test utag
# adata = run_utag_clust(adata, n_neighbors=10, 
#                        resolutions=[1], clustering_method="leiden",
#                        max_dist=20)

# with open('../times.list', 'w') as f:
#     # Write each item on a new line
#     f.writelines(f'{item}\n' for item in times_list)

# test second dataset
#del adata

df = pd.read_csv("input/measurementsEpCAMLy51MHCII-exported.csv")
list(df.columns)
meta = df.iloc[:, :4]
list(meta.columns)
mat = df.iloc[:, 4:]
list(mat.columns)
meta = meta.rename(columns={"Image": "Image ID_(standardized)", 
                             "Centroid X µm": "Centroid X (µm)_(standardized)", 
                             "Centroid Y µm": "Centroid Y (µm)_(standardized)"})

adata = ad.AnnData(mat)
adata.obs = meta
adata.layers["counts"] = adata.X.copy()

# run tests
times_list = []
#test phenograph
adata = RunPhenographClust(adata, n_neighbors=10, clustering_algo="leiden", 
                           min_cluster_size=10, primary_metric="euclidean", 
                           resolution_parameter=1, nn_method="kdtree")

# test scanpy
adata = RunNeighbClust(adata, n_neighbors=10, metric="euclidean", 
                       resolution=1, random_state=42)

# test parc
adata = run_parc_clust(adata, n_neighbors=10, 
                       dist_std_local=3, jac_std_global=0.15, 
                       small_pop=50, random_seed=42, 
                       resolution_parameter=1, hnsw_param_ef_construction=150)

# test utag
adata = run_utag_clust(adata, n_neighbors=10, 
                       resolutions=[1], clustering_method="leiden",
                       max_dist=20)

# save times
with open('../times.list', 'w') as f:
    # Write each item on a new line
    f.writelines(f'{item}\n' for item in times_list)
    
#