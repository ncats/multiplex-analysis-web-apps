---
nav_order: 7
---

# Spatial Interaction Tool (SIT)

## Introduction and Theory


## Workflow
The Spatial Interaction Tool (SIT) can be used once your dataset has been loaded, and phenotypes have been assigned to each cell in your dataset. The first place to visit for your SIT workflow is the Tool Parameter Selection page. Once settings are set in the Tool Parameter Selection 

### Tool Parameter Selection

This is the first page in the SIT workflow and the most important when it comes to identifying 

!['SIT Tool Parameter Selection'](.\assets\images\SIT_ToolParameterSelection.png)

### Run SIT Workflow

This page is the second page to visit when setting up SIT. Once the settings are set on the previous page, this workflow page will allow the user to identify which parts of the analysis to run.

!['SIT Run Workflow Page](.\assets\images\SIT_RunWorkflow.png)

From within this page you can see the following settings:

Tool Components for running Analysis

* Instantiate TIME Class
* Plot ROIs
* Calculate P values
* Check Metrics, impose plotting settings and convert to numpy format
* Plot desnity heatmaps per ROI
* Plot ROI outlines individually on the whole slides
* Average density P values over ROIS for each slide
* Average Density P values over TOIs for each annotation Region Type
* Plot density P values for each ROI over slide spatial plot

Job Execution Parameters:

* Should we use multiple logical CPUS to speed up calculations?
* Select number of threads for calculations:

### Display Individual ROI Heatmaps

### Display Average Heatmaps

### Display Average Heatmaps per Annotations

### Display ROI P Values Overlaid on Slides
