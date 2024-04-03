import pandas as pd 
import scanpy as sc


adata = sc.read("input/clust_dat.h5ad")

def find_markers(adata, umap_cur_col):
    sc.tl.rank_genes_groups(adata, umap_cur_col, method="wilcoxon", layer="counts")
    

find_markers(adata, "Cluster")

adata.uns["rank_genes_groups"].keys()