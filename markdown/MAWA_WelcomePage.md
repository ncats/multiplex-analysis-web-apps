## Multiplex Analysis Web Apps (MAWA)

### NCI CBIIT

## Welcome

Welcome to the Multiple Analysis Web Apps (MAWA) presented by NCATS-NCI-DMAP. This is your one stop resource for data exploration, algorithm tuning, and figure generation. The below is a guide to using this app with suggested workflow, step-by-step instructions, and FAQs

## Available Apps

1. File Handling
1. Multiaxial Gating
1. Phenotyping
1. Spatial Interaction Tool
1. Neighborhood Profiles

### 1. File Handling

As with most projects, the first step in starting your analysis workflow is with importing the data you intend to use. This is argueabely the most important step and one where the most issues may arise. If you any questions at all about importing and exporting data to not hesitate in seeking help from a DMAP team member (Dante and Andrew).

#### NIDAP Infrastructure

NIDAP manages files using a product called Compass. Compass is akin to other file management systems like Windows File Explorer and Apple Finder. Data (.csv, .txt, .png, etc files) are stored in a DATASET, which can be thought of as a file folder. It is the goal of this app that as data is made available to be processed, it is stored in a NIDAP DATASET. Similarly, as results, figures, and updated datatables are generated, these data objects will be placed in a NIDAP DATASET for later sharing, downloading, and storing.

#### 1a. Data Import and Export

1. Select any and all files that you want to import into MAWA from the the left-hand side of the screen in the section titled: **Available input data on NIDAP**.
2. Once those files are selected, click on the button in the middle of the screen that reads: **Load selected (at left) Input Data**.
3. Once the files have finished loading, make sure all selected files are visible on the right-hand side of the screen in the section titled: **Input data available to the tool.**

#### 1B. Datafile Unification

1. On Step 1 of this page, select any and all files that you wish to combine in the table on the left titled: **Select Datafiles**
2. Click the button that says **Combine selected files into single dataframe**
3. Once the files are done combining, review the table at the very bottom of the screen to identify that your data was loaded correctly.
4. Skip Step 2 in this section
5. On Step 3 of Datafile Unification, select the column labeled as `ShortName` for your Image identification. Click the **Assign Images** button to continue
6. On Step 4 of DU, simply click the button that says Assign ROIs and move on to the next step
7. On Step 5 of DU, keep the toggle which reads: **Select number of columns that specify one coordinate axis** to One Column. For the x-coordinate drop down select 'CentroidX' and for the y-coordinate select 'CentroidY'. Keep the micron conversion factor = 1um/pixel. Finalize this step by clicking the button titled: **Assign Coordinates**
8. For the first part of step 6 of the DU, expand the collapsed container and in the first drop down box select the option titled 'One Column for all Phenotypes'. After the page updates, in the second drop down box, select the column titled: `Pheno` or `Pheno_InEx`. Finish this section by selecting the button titled **Select Phenotypes**
9. For the second part of the Step 6 of the DU, once the phenotypes have been extracted, you should see a green box that says, `# phenotypes extracted`. In a table below that, you will see the names of the  phenotypes. You have the option now to rename them to something else, should you choose. I would recommend removing the '+' sign from the phenotype labels in the right-hand column of the renamed phenotypes (There are 5 phenotypes that require this). Once this is done, click on the button titled **Assign Phenotypes**
10. For Step 7 of the DU, save your version of the unified dataset with a memorable name. Keep in mind that it will be prefaced with 'mawa-unified_dataset_'. Click the button titled: 'Save dataframe to csv'.
11. Return to the Data Import and Export Tab and look for your titled dataframe in the dropdown select box in the center of the screen titled: 'Save MAWA-unified datafile to NIDAP'. Once your dataframe is selected, click the button titled: ‘Save Selected (above) MAWA unified datafile to NIDAP’. This may take some time (~1min) to complete. Once done, you will have a permeant version of the unified datafile to use in the future. Load it into the MAWA input space anytime you want to just like you would the base files.

#### 1C. Open File

1. If you have just completed the datafile unification process you will see the toggle at the top of the screen titled **Load Dataset from Datafile Unifier** is ON. Feel free to keep it ON.  
2. Go ahead and click the button below that reads: **Load the selected input dataset**. This may take a moment to complete (~1min). When it has you will see a sample of the dataset. Feel free to review it.  
3. If in the future, you have loaded a previously created mawa-unified dataset, once it is loaded into MAWA memory, you can move directly to this screen and load it using the drop down select box.  

### 2. Multiaxial Gating

### 3. Phenotyping

The second page to start your analysis pipeline is the phenotyping page. This is where you will load your data, view your different feature conditions, and select your phenotyping. There are two primary steps in performing phenotyping.

**Step 0**: Click the button at the top of the screen titled: Load Data. It may take a few minutes to complete (<3 min, working to improve this).
**Step 1**: Once the data is loaded, select a phenotyping method in the top right. The different phenotyping methods are described below. Once a methods has been selected, click the button titled **Apply Phenotyping Methond**. The app currently offers three different phenotyping methods. They are:
    `Species`: The phenotype is set as a direct copy of the species name (species name short to be more precise). The species name is the combination of markers that the cell type is positive for.
    `Markers`: The phenotype is set as one of the given Marker that the cell is positive for. If the cell is positive for more than one Marker, the cell entry in the dataset is duplicated and represented by each positive marker value being studied.
    `Custom`: The phenotype is set as a value of your own choise. Once the Custom phenotyping is selected, the Phenotype Summary Table on the right side of the page becomes editable.

### 4. Spatial Interaction Tool

### 5. Neighborhood Profiles (UMAP) Workflow

#### 5A. Neighborhood Profiles

1. Expand the collapsed container labeled: *Neighborhood Profiles Settings*. In this container make sure that `Number of CPUs` is set to 7 and `Calculate unique areas` is set to OFF. For the middle number boxes, set Percentage of cells to Subset for Fitting Step equal to 50 and Percentage of cells to Subset for Transforming Step equal to 50. Make sure the toggle titled Subset data transformed by UMAP is ON, and the toggle titled Load pre-generated UMAP is set to OFF.

2. Click the button titled **Perform Cell Density Analysis**. This should complete in under 5 min
3. Click the button titled **Perform UMAP Analysis**. This should complete in 10-15 min
4. Once the UMAP is complete, you will see options for performing the clustering. There are two ways to perform clustering.
    1. Perform Clustering on UMAP Density Difference = OFF: This will not perform any difference metrics, but instead will perform clustering on the whole UMAP distribution. No distinction is made based on a feature of the data. Try selecting a random cluster number between 1-10. Then click the button titled Perform Clustering Analysis. An elbow plot will appear to allow you to adjust the number of clusters to a number of your choosing. Each time you adjust the number of clusters, you will need to resubmit the Perform Clustering Analysis button. This should be completed in roughly 3 min (we are working to improve the timing). 
    2. Perform Clustering on UMAP Density Difference = ON:  This will allow you to perform individual clustering steps on regions of the UMAP which include cells of a specific feature condition. For example, how do differences between large nuclei and small nuclei cell contribute to the distribution of the UMAP? Select a column from the dropdown select box that is a numeric value (like area) or has exactly two unique values (For example: TRUE/FALSE). If you attempt to choose a categorical or string feature that has only 1 unique value or more than 2 unique values, you cannot perform the difference UMAP clustering. Choose any number of clusters for the FALSE and TRUE clustering to start off. For now, ignore the box titled: `Cutoff Percentage`. Rerun the clustering by hitting the box titled Perform Clustering Analysis. Once it has completed, elbow plots will appear under the cluster values for your investigation. You will also see many figures appear as well. Anytime you want to adjust the column being observed, or the number of clusters to use, you will need to resubmit the Perform Clustering Analysis button. This should be completed in roughly 3 min (we are working to improve the timing).
5. Once clustering is complete peruse the figures, as well as moving on to sections of Neighborhood Profiles//UMAP Differences and Neighborhood Profiles//Clusters Analyzer.

#### 5b. UMAP Differences Analyzer

After completing the UMAP decomposition and clustering analysis, the user may now take a look at the down-stream figures generated as a result of these analyses. While there are not many levers and knobs to change the data implicitly here, the user can generate different figures.

1. Before starting to view these Clustering Differences, you must complete at least the UMAP processing seen on the previous page. To experience the full offering of the Clustering Differences page, you must also complete the Clustering step on the previous page. There are warnings on the page to help you remember what needs to be completed in order to see each figure.
2. The Figures that are available for viewing:  
    1. Full 2D UMAP
    2. 2D UMAP filtered by lineage and features
    3. Different UMAP scaled by features

#### 5c. Clusters Analyzer

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