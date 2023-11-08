# Neighborhood Profiles Web App 
## NCATS-NCI-DMAP

## Welcome
Welcome to the Neighborhood Profiles presented by NCATS-NCI-DMAP. This is your one stop resource for data exploration, algorithm tuning, and figure generation. The below is a guide to using this app with suggested workflow, step-by-step instructions, and FAQs

## Neighborhoood Profiles workflow
Based on the conversations with our working groups, the following is the intended workflow for using this application.  
1. Import/Export Data
2. Phenotype Assignment
3. UMAP and Clustering Analysis
4. Review and Explore Differences in UMAP patterns
5. Review and Explore Differences in Cluster compositions.

Towards these end, the app has been organized into a series of pages that perform each of these actions. The instructions below will be outlined for each of these given pages.

### Phenotyping Instructions
The first page to start your analysis pipeline is the phenotyping page. This is where you will load your data, view your different feature conditions, and select your phenotyping. 

### 2. Neighborhood Profiles (UMAP) Instructions
Once data is loaded and phenotyped appropriately, further analyses can commence.   
1. Perform Cell Counts/Areas: Clicking on this button starts the process of calculating the density of the cells in a given neighborhood profile.
2. Perform UMAP: Clicking on this button starts the process of performing a 2-dimensional spatial UMAP on the dataset.
3. Perform Clustering: Clicking this button 

### 3. UMAP Differences ANalyzer
After completing the UMAP decomposition and clustering analysis, the user may now take a look at the down-stream figures generated as a result of these analyses. While there are not many levers and knobs to change the data implicitly here, the user can generate different figures.
1. Before starting to view these Clustering Differences, you must complete at least the UMAP processing seen on the previous page. To experience the full offering of the Clustering Differences page, you must also complete the Clustering step on the previous page. There are warnings on the page to help you remember what needs to be completed in order to see each figure.
2. The Figures that are available for viewing:  
    1. Full 2D UMAP
    2. 2D UMAP filtered by lineage and features 
    3. Different UMAP scaled by features 

### 4. Clusters Analyzer
After completing the UMAP decomposition and clustering analysis, the user may now take a look at the down-stream figures generated as a result of these analyses. The Cluster Analyzer page contains two figures generated from the upstream data analysis:
1. `Cluster/Phenotype Heatmap`  
2. `Incidence Lineplot`  

These figures have been created to investigate the composition of phenotypes of cells in assigned clustera, and the feature expression of cells in assigned s. Once a given figure is generated, you can change the name of the output in the text below each and add it as a figure to be exported in the `Data Input and Output` Page. The following sections are some general information about each figure:  

#### Cluster/Phenotype Heatmap
The heatmap offers a view of the number of each phenotyped cell located within each cluster. It offers three nromalization options for viewing the heatmap:

1. `No Norm`: No normalization is applied to the heatmap. The relative colors for each cell is scaled for all cells in all phenotypes in all clusters. If you were to sum the numbers shown in the grid, they would sum to the total number of cells fit to the spatial-umap model. 
2. `Norm within Clusters`: The grid values are decimal values of the number of cells within a cluster assigned to a given phenotype. In this schema, the relative color of the grid is based on the within- 
3. `Norm within Phenotypes`: The grid values are decimal values

#### Incidence Lineplot
The incidence lineplot details how the cells within each cluster differ in their expression of the data features recorded alongside the cell positions and marker values. These features range from boolean values (True/False), continuous values (-1, 0, 1), and string values('time0'). There are two selection boxes to augment the indicence line plot, and a radio button to select the type of comparison to perform. They are the following:

`Feature Selectbox`: Features that can be considered for the Incidence lineplot.
- Cell Counts: The number of cells assigned to a given cluster
- HYPOXIC, NORMOXIC, NucArea, RelOrientation, etc: Any other feature that specified to Dante/Andrew as one that is worth showing. YOU MUST tell us which ones you want and we will set it up for you. 

`Phenotype Selectionbox`: The phenotype the cells being plotted. The options shown are:
- All Phenotypes: Shows all cells irrespective of phenotype
- VIM+, ECAD+, VIM+ECAD+, Other, etc...: The other phenotypes that have been selected in the Phenotyping stage of the workflow.

`DisplayAs Radio Button`: How the values of the Feature selectbox should be displayed. This radio button is disabled for the Cell Counts condition, but is enabled for any other Feature selection. The options to be displayed are:
- Count Differences: The value shown on the y-axis is the difference between the number of cells in a cluster in the Y>0 condition subtracted from the number of cells in that cluster in the Y<0 condition.
- Percentages: The value shown on the y-axis is the percentage of cells that match a feature condition in that given cluster. If you were to sum all the values across the clusters, they would sum to 100%.  
- Ratios: The value shown on the y-axis is the ratio of r1/r0 where r1 is the precentage of cells that match the feature of condition shown on y>0 in that cluster, and r0 is the percentage of cells that match the feature of the condition show on y<0 in that cluster.

## FAQs
Q: How do I add more features to view in the filtering step?  
A: Ask a member of NCATS-NCI-DMAP(Dante or Andrew) to add that column name to the dropdown menus  

Q: What do if I cant find the data I want in the drop-down menus?  
A: The easiest way to add new data for access is to add it to the following NIDAP dataset:
   
Q: Where is the data/figure I exported to NIDAP?  
Q: How do I know if my data is in the right format for use in this app?  
Q: Can I load a previously generated phenotyping summary file?  