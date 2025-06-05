---
title: Neighborhood Profiles Workflow
nav_order: 8
layout: default
---

# Neighborhood Profiles Workflow

!['Neighborhood Profiles'](/assets/images/NeiPro_Splash.png)

## Introduction and Theory

Neighborhood Profiles identifies the types of cells that often cluster with one another, and find patterns of these types of clusterings occurring more often in certain tissue types or health-conditions. This goal of these analyses is to help identify the specific neighborhoods present in a given dataset, and characterize their makeup based on the phenotypes present.

## Workflow

The Neighborhood Profiles Workflow can be used once your data has been properly loaded and phenotyped. Of particular importance is that the your data contains X/Y coordinates for each cell, and each cell is categorized as a phenotype in the phenotyping page. The sections below will show examples of settings and visualizations that can be created in the Neighborhood Profiles Workflow. For these examples, we used a sample dataset included with MAWA called **Multiplex_Image_Sample.csv** and selected three Markers (ECAD, HistoneH3, and CD45RO) for our phenotyping step. When combined, these three markers create eight unique phenotypes with the *Species phenotyping* method in Thresholded intensities.

### Neighborhood Profiles

This page is the starting point for running the Neighborhood Profiles. On the top of this page you will see options for running each step of the analysis. At the bottom of the page are placeholder spaces for two figures.

In the top panel, there are the following buttons

* Perform Cell Counts/Areas Analysis
* Perform UMAP Analysis
* Perform Clustering Analysis

There is also a collapsable container for the Neighborhood Profiles Settings. This container allows you to set the parameters for the Neighborhood Profiles analysis. These settings include:

* Number of CPUs
* Calculate Unique Areas
* Area Filter Ratio
* Percentage of cells to subset for UMAP fitting step
* Subset data transformed for UMAP
* Percentage of cells to subset for UMAP transforming step
* Load pre-generated UMAP model

The bottom figure panels are labeled

* Clusters Plot
* Neighborhood Profiles

#### Instructions for Use

1. Start by checking the message text that can be seen in the middle of the screen. If you have not completed the previous phenotyping step, this message will appear as the following.

!['Neighborhood Profiles Step 0'](./assets/images/NeiPro_Step0.png)

2. If you have finished the phenotying step, the middle message text should look like the following. This means the app is ready to be used for processing Neighborhood Profiles.

!['Neighborhood Profiles Step 1'](./assets/images/NeiPro_Step1.png)

3. Begin the Cell Counts/Area Analysis by clicking the button at the top of the page. This process will take varying amounts on time depending on the size of your dataset, the size of your images, and the number of phenotypes you have selected in the phenotyping step. For a dataset of 48k cells, and 8 phenotypes, this process takes approximately 5min. When this step has completed correctly, you will see a message in the middle of the screen that reads as the following:

!['Neighborhood Profiles Step 2'](./assets/images/NeiPro_Step2.png)

4. Next, begin the UMAP Analysis by clicking the button at the top of the page. Running the UMAP decomposition will take varying amounts on time depending on the size of your dataset, the size of your images, and the number of phenotypes you have selected in the phenotyping step. For a dataset of 48k cells, and 8 phenotypes, this process takes approximately 1 min. When this step has completed correctly, you will see a message in the middle of the screen that reads as the following:

!['Neighborhood Profiles Step 3'](./assets/images/NeiPro_Step3.png)

5. Finally, we can begin the process of clustering the results of the UMAP. At the lower part of the analysis section, select a number of clusters to use for the k-means clustering algorithm. Once a number is selected, click on the button to *Perform Clustering Analysis*. This step is the fastest, and depending on the number of cells in your datasets, this should take under 1 min to run. When its complete, the user should see be able to see the scatter below populated and colored by cluster, and the Neighborhood Profiles line plots populated and drawn for the phenotypes included in each cluster.

6. In the clusters figure, observe that any time you can swap back and forth between plot the colors of the scatter plot by the Cluster label or by the Phenotype label. You can also progress through all the images that are included in your dataset and see how the clusters have partitioned individual tissues samples

!['Neighborhood Profiles Clusters Plot'](./assets/images/NeiPro_ClustersPlot.png)

7. In the Neighborhood Profiles Line plot, the user can observe the make up (profiles) of each of the clusters (neighborhoods) created by the K-means clustering algorithm. These line plots show the density measurement of the number of cells for a given phenotype (counts /mm2) within different annuli surrounding the cells of that given cluster. The user can swap between different clusters by selecting them from the drop down menu. All figures are scaled to fit on the same axes, and as such any one phenotype from any one cluster might contribute more to the overall scale of these line plots. As is often case, there are many cells (often the majority) that are assigned to the 'Other' phenotype. When this other category dwarfs the other phenotypes, it might be helpful to hide the 'Other' phenotype from the figure. This can be done in the Options menu seen in the second image below. This action can be similarly done for the 'No Cluster' cluster if that has been created in your workflow. As you test differences in your UMAP results as a factor of the k-means cluster size, your Neighborhood Profiles line plots will change, giving the user an opportunity to tune and test their assumptions of the data.

!['Neighborhood Profiles Line Plot'](./assets/images/NeiPro_NeighPlot.png)

!['Neighborhood Profiles Line Plot Options'](./assets/images/NeiPro_NeighPlot_options.png)

### UMAP Differences

UMAP Differences if the next step in your Neighborhood Profiles analysis pipeline. After the UMAP decomposition is completed, and the clustering has completed, this page will offer further down-stream analsyses to consider. While there are not many levers and knobs to change the data implicitly here, the user can generate different figures based other outcome variables from the input dataset, and the phenotypes that were defined earlier in the MAWA workflow.

!['UMAP Differences Page'](./assets/images/UMDiff_Sample.png)

Before starting to view these Clustering Differences, you must complete at least the UMAP processing seen on the previous page. To experience the full offering on this page, the user will also want to complete the Clustering step on the Neighborhood Profiles page. There are warnings on the page to help you remember what needs to be completed in order to see each figure.

The Figures that are available for viewing are:

1. `Full 2D density UMAP`  
1. `2D UMAP density filtered by Features and Phenotypes`  
1. `Difference UMAP filtered by Features and Phenotypes (Both density and clustered)`  

#### Full 2D UMAP

This is purely the outcome of UMAP decomposition colored by density of cells. This is meant to be a template to compare to, after the other UMAP permutations are generated, both within this study design as well as for future reproducibility.

#### 2D UMAP filtered by lineage and features

This is a copy of the UMAP on the left, but with additional filtering options. Specifically, this allows you adjust which cells from your dataset contribute to which parts of the UMAP. Filtering by a specific phenotype will show the parts of the UMAP that include that phenotype. If instead the user filters by a Feature of the dataset(AlivePatient or Nucleus size for examples), then again the UMAP will filter by the data that for that feature. At any time, the user can swap between viewing the density UMAP and the clustered UMAP to evaluate different information. Additionally if at any time you want to look at Markers vs Phenotypes, that option is also available.

#### Differnce UMAP scaled by features

These are again further copies of the 2D UMAP. These are displayed by default both as denisty and clustered, side by side. These difference UMAPs aim to show you This give you options to check differences between feature conditions. At this time, this works best with Boolean data (True/False), but also works with range data (x>0, x<0).

### Clusters Analyzer

The final step in the Neighborhood Profiles workflow is the clustering analysis. Again, this page will be most useful after all three analysis steps from the [Neighborhood Profiles](#neighborhood-profiles) page are completed. The Cluster Analyzer page contains two figures generated from the Neighborhood Profiles analysis:

1. [Phenotype/Cluster Heatmap](#phenotypecluster-heatmap)
1. [Incidence Figure](#incidence-figure)

!['Clusters Analyzer Page'](./assets/images/clust_analyzer_main.png)

These figures have been created to investigate the composition of phenotypes of cells in assigned clusters, and the feature expression of cells in assigned clusters. Each figure can be customized further using the options available in the interface. Each figure can be exported for use in other applications by right-clicking on the image and clicking 'save as'.

IMPORTANT: All of the values displayed in the figures on this page are measured from the cells used in the UMAP processing step (Step 2 in the [Neighborhood Profiles](#neighborhood-profiles) page). It is likely the case that the cells used in the UMAP step are a subset of the full dataset, and as such, the values shown in these figures may seem smaller than what the full dataset describes. If you want the UMAP model to be applied to all cells in your dataset, you can change the Neighborhood Profiles Settings to transform the full dataset. Using the full dataset will take longer to run, and will not necessarily improve the results of the clustering analysis.

The sections below will show examples of settings and visualizations that can be created in the Clusters Analyzer page. For these examples, we used a sample dataset included with MAWA called **Multiplex_Image_Sample.csv** and selected three Markers (ECAD, HistoneH3, and CD45RO) for our phenotyping step. When combined, these three markers create eight unique phenotypes with the *Species phenotyping* method in Thresholded intensities.

#### Phenotype/Cluster Heatmap

The heatmap offers a view of the number of each phenotyped cell located within each cluster. The heatmap can be modified using a toggle switch for heatmap normalization. These widgets have the following properties.

`Normalization Toggle`: This toggle switch allows you to change the normalization method applied to the heatmap. The normalization options are as follows:

1. No Norm (default): No normalization is applied to the heatmap. The relative colors for each cell is scaled for all cells in all phenotypes and in all clusters. The sum of each number shown in the grid corresponds to the total number of cells transformed by the UMAP model.

![Heatmap No Normalization](./assets/images/clust_analyzer_heatmap1.png)

2. Norm within Clusters: The grid values are decimal values of the number of cells within a cluster assigned to a given phenotype. In this schema, the relative color of the grid is based on the within-cluster distribution. The sum of the numbers in each row sum to 1.

![Heatmap Norm within Clusters](./assets/images/clust_analyzer_heatmap2.png)

3. Norm within Phenotypes: The grid values are decimal values of the number of cells within a phenotype assigned to a given cluster. In this schema, the relative color of the grid is based on the within-phenotype distribution. The sum of the numbers in each column sum to 1.

![Heatmap Norm within Phenotypes](./assets/images/clust_analyzer_heatmap3.png)

#### Incidence Figure

The incidence figure is one way to represent the counts of the cells present in each cluster. In this example I choose the standard form of clustering analysis with 5 clusters. When it is first loaded, it looks like the following (Figure 1):

The incidence figure details how the cells within each cluster differ in their expression of the data features recorded alongside the cell positions and marker values. These features range from boolean values (True/False), continuous values (-1, 0, 1), and string values('time0'). There are two selection boxes to augment the indicence figure, and a radio button to select the type of comparison to perform. They are the following:

`Feature Select Box`: Features that can be considered for the Incidence figure.

* Cell Counts: The number of cells assigned to a given cluster (Default)
* Other features from your datasets: Columns from your dataset by which you want to directly compare TWO conditions. At this time, this works best with Boolean data (True/False), but also works with range data (x>0, x<0). Once a feature is selected, the incidence plot no longer shows a pure count, but instead a comparison of the two conditions within the feature.

`Phenotype Select Box`: The phenotype the cells being plotted. The options shown are:

* All Phenotypes: Shows all cells irrespective of phenotype (Default)
* VIM+, ECAD+, VIM+ECAD+, Other, etc...: Shows only the cells that express for the specifically chosen phenotype (created during the Phenotyping stage of the workflow).

`Display-as Radio Button`: How the values of the Feature select box should be displayed. This radio button is disabled for the Cell Counts condition, but is enabled for any other Feature selection. For each of the options shown below, there are equations detailing how the values for each condition are calculated. For equations 1-4, *d* represents the full datasets, and the *Condition* is the Feature and Value combinaton being condisdered. Each equation is also considered at each cluster in the dataset. Therefore when a feature is selected, the dataset will be split as follows

$$
\begin{aligned}
d_{\text{cond0}} = d \subset Condition0 \quad\quad\quad \text{Equation 1a}\\
d_{\text{cond1}} = d \subset Condition1 \quad\quad\quad \text{Equation 1b}\\
\end{aligned}
$$

The options to be displayed are:

* Count Differences: The value shown on the y-axis is the difference between the quantity of cells in a cluster in the Y>0 condition subtracted from the quantity of cells in that cluster in the Y<0 condition.

$$
\begin{aligned}
d\_count_{\text{clust, cond0}} &= |d_{\text{cond0}}| \hspace{4.5cm} \text{Equation 2a} \\
d\_count_{\text{clust, cond1}} &= |d_{\text{cond1}}| \hspace{4.5cm} \text{Equation 2b} \\
d\_diff_{\text{clust}} &= d\_count_{\text{clust, cond1}} - d\_count_{\text{clust, cond0}} \hspace{1.35cm} \text{Equation 2c}
\end{aligned}
$$

* Percentages: The value shown on the y-axis is the percentage of cells that match a feature condition in that given cluster. The sum of the values across a given cluster would be to 100%.

$$
\begin{aligned}

\end{aligned}
$$

* Ratios: The value shown on the y-axis is the ratio of r1/r0 where r1 is the percentage of cells that match the feature of condition shown on y>0 in that cluster, and r0 is the percentage of cells that match the feature of the condition show on y<0 in that cluster.

`Show Raw Counts Toggle`: This toggle switch is enabled for any Feature selection other than Cell Counts. When a feature is selected from the `Feature Select Box`, one of the three Display-As analyses are performed to show a line plot of the Incidence for that feature across clusters. This line plot is a very simply comparison of two numbers, and it can hide the magnitude of the raw cell counts. By turning this toggle to true, a bar plot is overlayed 

#### Example workflow for Incidence Figure
There are two drop-down menus for Feature and Phenotype, the defaults of which are Cell Counts and All Phenotypes respectively. Below the drop-down tables is a set of radio buttons and a toggle switch, the default states of which are disabled. Following these widgets is the Incidence Figure. On the figure, the specific clusters are measured on the horizontal axis. The vertical axis is the incidence being described. In its default state (Cell Counts), the bar plot shows the number of cells present in each cluster, irrespective of phenotype. To be specific, these are the counts of cells that were clustered. The clustered cells are not always your full dataset; it depends on how you perform the UMAP step before the clustering. See the Neighborhood Profiles Settings section for how to sample your data for UMAP. In this example, we sampled about ~7300 total cells. This bar plot, shows the distribution of where those 7300 cells fall in each cluster.

!['Cell Counts by Cluster'](./assets/images/clust_analyzer_inci_fig1.png)

It looks like Cluster 3 had least number of cells, and Cluster 5 had the most. At its core, this is the starting point for using this figure. Here we can begin to understand the presence of cells within each cluster. We can next start to drill down into the data a bit further and look at how cells of a specific phenotype are represented in each cluster. For example, to investigate how CD45RO+ HistoneH3+ is represented in each cluster, you can select that phenotype from the right menu and the figure will adjust to look something like Figure 2. Here the scale is much smaller, and that the bar plot distribution changes. This is now a subset of the data seen in the previous figure, only CD45RO+ HistoneH3+ cells are shown instead of all cells. It also appears that the great majority of CD45RO+ HistoneH3+ cells are found in Cluster 3 and not many at all in the other clusters.

!['CD45RO+ HistoneH3+ Cell Counts by Cluster'](./assets/images/clust_analyzer_inci_fig2.png)

Now this first example focuses on pure Cell Counts and does not consider specific feature differences. This figure was built to also show some differences for specific Features or measurements of the data other than phenotype. In this sample dataset, there is a feature called DNA1. The values of DNA1 are all numerical and were a continuous measurement. The left side drop-down menu contains a list of all available features in the dataset as well as the default entry of *Cell Counts*.  When a Feature other than *Cell Counts* is selected, the dataset will split into two halves based on the values in the dataset. Selecting a Feature will also enable the radio buttons and the toggle switch.

If the selected Feature contains exactly two unique values, then the dataset will be evenly split between the values, and the figure will look similar to Figure 3. If the feature column contains has more than 2 unique values and the values are numerical, then the median value of the range will be found, and the data will be split evenly around the median. If the feature column only has 1 or fewer unique values, or if the data has more than 2 unique values and is a string value, MAWA will tell the user that the data cannot be easily split and comparison on this Feature is innapropriate.

Once the feature is split into two parts, the figure will display differences between the parts. To that end, the first dataset of the feature will be drawn to the top of the figure (Upper, above the horizontal axis), and the second half will be drawn to the bottom of the figure (Lower, below the horizontal axis). 

When the first radio button is selected, the figure draws the visualization as *Counts Differences*. This line graph is the difference of the quantiy of cells that match the Upper condition subtracted from the quantity of cells that match the Lower condition (Eq 1). For a given cluster, if the drawn line falls above the horizontal axis, then there are more values in the Upper condition than there are in the Lower condition, and vice versa. As you might expect, you can drill down further into the data by selecting different Phenotype subsets in combination with the Feature selection.

When the second radio button is selected, the figure draws the visualization as a *Percentage*. Presently, this *Percentage* focuses only on the Upper condition of the Feature. This calculates the percentage of Upper value cells in each cluster that across all Upper value cells. This means that the percentages seen in this line plot (in our example, five percentages) should sum to 100% (Figure 4). Now there are some cases, where for one reason or another, there are 0 cells that match a condition which can make this sort of calculation troublesome. For this reason, the actual calculation that is performed adds a 1 to the numerator and denominator in equations 2a and 2b.

Finally, when the third radio button is selected, the figure draws the visualization as a *Ratio*. This *Ratio* returns to a comparison between the Upper and Lower conditions of the selected Feature (Figure 5). Specifically, each value on the visualization is the ratio of the percentages generated from equations 2a and 2b (equation 3). The aim of this *Ratio* is to illustrate the magnitude of how many more (or fewer) cells appear for a condition in the Upper condition compared to the Lower condition.
