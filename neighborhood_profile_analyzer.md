---
nav_order: 8
---

# Neighborhood profiles analyzer

## More to come!

### UMAP Differences Analyzer
After completing the UMAP decomposition and clustering analysis, the user may now take a look at the down-stream figures generated as a result of these analyses. While there are not many levers and knobs to change the data implicitly here, the user can generate different figures.

1. Before starting to view these Clustering Differences, you must complete at least the UMAP processing seen on the previous page. To experience the full offering of the Clustering Differences page, you must also complete the Clustering step on the previous page. There are warnings on the page to help you remember what needs to be completed in order to see each figure.

2. The Figures that are available for viewing:  
`Full 2D UMAP`  
`2D UMAP filtered by lineage and features`  
`Difference UMAP scaled by features`  

#### Full 2D UMAP
This is purely the outcome of UMAP decomposition colored by density of cells. This is meant to be a template to compare all additional UMAPs, both within this study design as well as for future reproduceability. 

#### 2D UMAP filtered by lineage and features
This is a copy of the UMAP on the left, but with additional filtering options.

#### Differnce UMAP scaled by features
These are again further copies of the 2D UMAP. This give you options to check differences between feature conditions. At this time, this works best with Boolean data (True/False), but also works with range data (x>0, x<0). 

### Cluster Analyzer
After completing the UMAP decomposition and clustering analysis, the user may now take a look at the down-stream figures generated as a result of these analyses. The Cluster Analyzer page contains two figures generated from the upstream data analysis:

`Cluster/Phenotype Heatmap`  
`Incidence Lineplot`  

These figures have been created to investigate the composition of phenotypes of cells in assigned clusters, and the feature expression of cells in assigned clusters. Once a given figure is generated, you can change the name of the output in the text below each and add it as a figure to be exported in the Data Input and Output Page. The following sections are some general information about each figure:

#### Cluster/Phenotype Heatmap
The heatmap offers a view of the number of each phenotyped cell located within each cluster. It offers three nromalization options for viewing the heatmap:

1. `No Norm`: No normalization is applied to the heatmap. The relative colors for each cell is scaled for all cells in all phenotypes in all clusters. If you were to sum the numbers shown in the grid, they would sum to the total number of cells fit to the spatial-umap model.
2. `Norm within Clusters`: The grid values are decimal values of the number of cells within a cluster assigned to a given phenotype. In this schema, the relative color of the grid is based on the within-cluster distribution.
3. `Norm within Phenotypes`: The grid values are decimal values of the number of cells within a phenotype assigned to a given cluster. In this schema, the relative color of the grid is based on the within-phenotype distribution. 

#### Incidence Lineplot
The incidence lineplot details how the cells within each cluster differ in their expression of the data features recorded alongside the cell positions and marker values. These features range from boolean values (True/False), continuous values (-1, 0, 1), and string values('time0'). There are two selection boxes to augment the indicence line plot, and a radio button to select the type of comparison to perform. They are the following:

`Feature Selectbox`: Features that can be considered for the Incidence lineplot.

* Cell Counts: The number of cells assigned to a given cluster
* HYPOXIC, NORMOXIC, NucArea, RelOrientation, etc: Columns from your dataset by which you want to directly compare TWO conditions. At this time, this works best with Boolean data (True/False), but also works with range data (x>0, x<0). Once a feature is selected, the incidence plot no longer shows a pure count, but instead a comparison of the two conditions within the feature.

`Phenotype Selectionbox`: The phenotype the cells being plotted. The options shown are:

* All Phenotypes: Shows all cells irrespective of phenotype
* VIM+, ECAD+, VIM+ECAD+, Other, etc...: Shows only the cells that express for the specifically chosen phenotype (created during the Phenotyping stage of the workflow).

`DisplayAs Radio Button`: How the values of the Feature selectbox should be displayed. This radio button is disabled for the Cell Counts condition, but is enabled for any other Feature selection. The options to be displayed are:

* Count Differences: The value shown on the y-axis is the difference between the number of cells in a cluster in the Y>0 condition subtracted from the number of cells in that cluster in the Y<0 condition.
* Percentages: The value shown on the y-axis is the percentage of cells that match a feature condition in that given cluster. The Sum of the values across a given cluster would be to 100%.
* Ratios: The value shown on the y-axis is the ratio of r1/r0 where r1 is the precentage of cells that match the feature of condition shown on y>0 in that cluster, and r0 is the percentage of cells that match the feature of the condition show on y<0 in that cluster.