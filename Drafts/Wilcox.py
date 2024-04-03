import pandas as pd 
import scanpy as sc


adata = sc.read("input/clust_dat.h5ad")

def find_markers(adata, umap_cur_col):
    sc.tl.rank_genes_groups(adata, umap_cur_col, method="wilcoxon", layer="counts")
    

find_markers(adata, "Cluster")


result = adata.uns["rank_genes_groups"]
groups = result["names"].dtype.names
difExprDf = pd.DataFrame(
    {
        group + "_" + key[:1]: result[key][group]
        for group in groups
        for key in ["names", "pvals"]
    }
)

geneNames = pd.DataFrame(adata.uns["rank_genes_groups"]["names"])
genePvals = pd.DataFrame(adata.uns["rank_genes_groups"]["pvals"])

len(pd.unique(adata.obs["Cluster"]))
len(geneNames.columns)