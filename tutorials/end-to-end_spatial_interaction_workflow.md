---
layout: default
parent: Tutorials
nav_order: 1
---

# End-to-End Spatial Interaction Workflow
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

This two-part tutorial outlines some recommended steps for using MAWA to (1) perform manual phenotyping on a CSV file containing multiplex intensities for one or more images and (2) obtain measures of the degrees to which cells of the defined phenotypes are interacting. This could be used for correlating with response variables to detect biomarkers. For example, one could say that in image A, for which a patient survived a cancer, their M2 macrophages interacted heavily with CD8 T cells, whereas in image B, corresponding to a patient who died, their M2 macrophages generally repelled away from CD8 T cells.

## Part I: Manual Phenotyping Using Multiaxial Gating

### Step 1: Load a sample multiplex datafile into MAWA's Code Workspace

When you start MAWA (by clicking on the link or file sent to you by the admins) you should see this:

<img src="end-to-end_spatial_interaction_workflow/welcome_page.png" alt="Welcome page" width="1024"/>

Click on the "Data Import and Export" page,

<img src="end-to-end_spatial_interaction_workflow/welcome_page-highlighted_data_import_and_export.png" alt="Highlighted Data Import and Export page" width="400"/>

loading a page that looks like this:

<img src="end-to-end_spatial_interaction_workflow/data_import_and_export_page.png" alt="Data Import and Export page" width="1024"/>

Load the sample dataset by clicking the checkbox next to "Multiplexed_Image_Sample.csv" in the "Available input data on NIDAP" section and clicking the button "Load selected (at left) input data",

<img src="end-to-end_spatial_interaction_workflow/data_import_and_export_page-highlighted_data_loading_options.png" alt="Highlighted data loading options" width="1024"/>

resulting in the file being loaded from NIDAP's main filesystem to the Code Workspace:

<img src="end-to-end_spatial_interaction_workflow/data_import_and_export_page-loaded_datafile.png" alt="Loaded datafile" width="1024"/>

This datafile is now ready to be processed by the MAWA app suite in the Code Workspace.

### Step 2: Load the sample multiplex datafile into the multiaxial gater app

Click on the "Multiaxial Gating" page,

<img src="end-to-end_spatial_interaction_workflow/data_import_and_export_page-highlighted_multiaxial_gating.png" alt="Highlighted Multiaxial Gating page" width="400"/>

loading a page that looks like this:

<img src="end-to-end_spatial_interaction_workflow/multiaxial_gating_page.png" alt="Multiaxial Gating page" width="1024"/>

Note that the file you just imported into the app, "Multiplexed_Image_Sample.csv", is availble in the "Filename" field. Since the coordinates in this file are already in microns, change the "x-y coordinate units (microns)" value to 1 by selecting what's in the box and typing "1", and if you wish, turn on the switch labeled "Perform batch normalization", which transforms the values in each numeric column in each image to a Z score. Then click the "Load data" button,

<img src="end-to-end_spatial_interaction_workflow/multiaxial_gating_page-dataset_settings.png" alt="Dataset settings" width="400"/>

loading a page that looks like this:

<img src="end-to-end_spatial_interaction_workflow/multiaxial_gating_page-loaded_dataset.png" alt="Multiaxial Gating page with loaded dataset" width="1024"/>

---

### Step 3: Build a phenotype by gating on multiple intensity channels

Select a column corresponding to an intensity channel in the "Column for filtering" dropdown, say, CD16:

<img src="end-to-end_spatial_interaction_workflow/multiaxial_gating_page-cd16_selected.png" alt="CD16 selected" width="400"/>

☝️ ***Tip:** As for all dropdown menus, you can easily search for a column name by clicking on the dropdown and immediately typing in the name (case-insensitive).*

A histogram (smoothed by a kernel density estimate) of the distribution of Z scores over the CD16 column for all images (since "All images" is selected in the "Image selection" field) is displayed in red. To select only a subset of intensities, limit the selected range by adjusting either end of the double-ended slider in the "Selected value range" field. For example, to limit the selection to cells with CD16 intensity values that are at least 1.5 standard deviations above the mean, move the left end of the slider to about 1.5, and click on the "Add column filter to current phenotype" button,

<img src="end-to-end_spatial_interaction_workflow/multiaxial_gating_page-cd16_range_modified.png" alt="CD16 range modified" width="400"/>

leading to this:

<img src="end-to-end_spatial_interaction_workflow/multiaxial_gating_page-cd16_filtered.png" alt="CD16 filtered" width="800"/>

You'll see the "Column filter" you just created gets added to the "Current phenotype" that you are in the process of defining. The CD16 column in the "Column filter" section disappears as an option because it has already been added to the current phenotype, resulting in the "Column for filtering" selection to reset to the first numeric column in the dataset.

Note that you can refine ranges at any point by double clicking in the "Minimum value" or "Maximum value" columns in the "Current phenotype" table and editing the bound using the keyboard. For example, here we can refine "1.51" to "1.5", hitting Enter to finalize the edit:

<img src="end-to-end_spatial_interaction_workflow/multiaxial_gating_page-cd16_refined_minimum_bound.png" alt="CD16 refined minimum bound" width="400"/>

As a cell phenotype can be defined by thresholding multiple fields, we can filter on additional columns in the dataset to define a phenotype. (E.g., one may want to set a lower bound on the DAPI channel in order to be confident that the object in the dataset has a sufficient nucleus to be a real cell.) Let's do this by adding a filter for low PDL1 expression by selecting PDL1 in the "Column for filtering" dropdown and setting the maximum bound to 0:

<img src="end-to-end_spatial_interaction_workflow/multiaxial_gating_page-pdl1_with_maximum_bound.png" alt="PDL1 with maximum bound" width="400"/>

Hit the "Add column filter to current phenotype" button and update the maximum bound value if you didn't get it quite right with the slider:

<img src="end-to-end_spatial_interaction_workflow/multiaxial_gating_page-first_phenotype_defined.png" alt="First phenotype defined" width="600"/>

Now enter a name for the phenotype in the "Phenotype name" box and click the "Add phenotype to assignments table" button,

<img src="end-to-end_spatial_interaction_workflow/multiaxial_gating_page-adding_first_phenotype_to_gating_table.png" alt="Adding first phenotype to gating table" width="600"/>

leading to this:

<img src="end-to-end_spatial_interaction_workflow/multiaxial_gating_page-first_phenotype_in_gating_table.png" alt="First phenotype in gating table" width="600"/>

You have now defined the "CD16 high PDL1 low" phenotype as can be seen in the "Phenotype assignments" section. Note you can further modify your phenotype criteria by editing the values in the table, as done previously.

### Step 4: Add additional phenotypes

Now add two more phenotypes in the same manner, adding one or more column filters to the table in the "Current phenotype" section, and when you're done defining each phenotype, assigning a name to the phenotype and pressing the "Add phenotype to assignments table" button to transfer it to the full gating table in the "Phenotype assignments" section. Define the following phenotypes:

<img src="end-to-end_spatial_interaction_workflow/multiaxial_gating_page-full_phenotype_gating_table.png" alt="Full phenotype gating table" width="800"/>

☝️ ***Tip:** You can view a large image of the phenotype assignments table (gating table) by clicking on the "maximize" icon that appears at the top right of the table when hoving your mouse over it:*

<img src="end-to-end_spatial_interaction_workflow/multiaxial_gating_page-full_phenotype_gating_table-small.png" alt="Full phenotype gating table - small" width="600"/>

### Step 5: Generate the defined phenotypes, adding them to the dataset

Click the "Append phenotype assignments to the dataset" button to generate the phenotypes specified in the "Phenotype assignments" table:

<img src="end-to-end_spatial_interaction_workflow/multiaxial_gating_page-add_phenotype_assignments_to_dataset.png" alt="Add phenotypes to the dataset" width="600"/>

Note that the sample of the dataset that's always displayed in the "New dataset" section now has three new phenotype columns appended to the end:

<img src="end-to-end_spatial_interaction_workflow/multiaxial_gating_page-augmented_dataset_sample.png" alt="Augmented dataset sample" width="600"/>

### Step 6: Visualize the new phenotypes

Select an image and phenotype to plot in the "Image to plot" and "Phenotype to plot" fields, respectively, and click on the "Plot the selected phenotype in the selected image" button to display a corresponding scatterplot:

<img src="end-to-end_spatial_interaction_workflow/multiaxial_gating_page-visualize_new_phenotypes.png" alt="Plot of the selected phenotype in the selected image" width="600"/>

☝️ ***Tip:** The displayed scatterplot is fully interactive (zoom, pan, etc.). Double click on it to reset to the original view.*

### Step 7: Resolve any phenotype conflicts

Generally, the phenotypes you define will be mutually exlusive because it is useful for a cell to be described by one and only one phenotype. If a single cell is described by more than one phenotype, you should refine your phenotype definitions until the cells in your dataset are uniquely defined by them. However, there will often be cases in which this cannot be done precisely.

For example, say that you defined the "X low" phenotype as having the Z score range [-10, 0.1] in the "X" column and the "X high" phenotype as having the Z score range [0.0, 15]. You would then expect some cells to be both "X low" positive and "X high" positive. This is generally a bit nonsensical and refers to the small population of cells having "X" column Z scores in beween 0.0 and 0.1.

Instead of refining your phenotype definition, you may wish to simply say "Define all 'X low' positive and 'X high' positive cells as 'X high' only (i.e., not 'X low' as well)." This can be addressed in this step of the workflow.

In the three phenotype definitions specified above, we do not have any cells that fall into this category. But, we still did not ensure that all our phenotype definitions were mutually exclusive. While the first two phenotypes "CD16 high PDL1 low" and "CD16 high PDL1 high" *are* mutually exclusive, we made no effort to ensure that the third phenotype "Ecad mid" is different from the first two. Indeed, 7.5% of the cells in the dataset are both "Ecad mid" and either of the first two phenotypes. This can also be addressed in this step of the workflow, as will now be demonstrated.

☝️ ***Tip:** There may be situations in which you want cells to be labeled with multiple phenotypes, and this is a highly supported feature: Simply perform subsequent phenotyping in the "Phenotyping" app using the "species" phenotyping method. Instead, below you will use the "Custom" phenotyping method to re-assign such compound phenotypes.*

Click on the "Phenotyping" page,

<img src="end-to-end_spatial_interaction_workflow/multiaxial_gating_page-highlighted_phenotyping.png" alt="Highlighted Phenotyping page" width="400"/>

and on the page that loads, click on the "Load Multiaxial Gating Data" button:

<img src="end-to-end_spatial_interaction_workflow/phenotyping_page-load_multiaxial_gating_button.png" alt="Load Multiaxial Gata Data button" width="400"/>

In the "Choose a Phenotyping Method" box, select "Custom" and click the "Apply Phenotyping Method" button:

<img src="end-to-end_spatial_interaction_workflow/phenotyping_page-set_custom_phenotyping_method.png" alt="Select Custom phenotyping method" width="1024"/>

The table in the "Phenotype Assignments" section now looks like:

<img src="end-to-end_spatial_interaction_workflow/phenotyping_page-phenotype_assignments_table_unassigned.png" alt="Phenotype assignments table unassigned" width="800"/>

Double click in each row in the "phenotype" column to potentially re-assign compound phenotypes with your keyboard. For example:

<img src="end-to-end_spatial_interaction_workflow/phenotyping_page-phenotype_assignments_table_assigned.png" alt="Phenotype assignments table unassigned" width="800"/>

☝️ ***Tip:** There is no need to assign a "phenotype" to every row in the "Phenotype Assignments" table. The downstream spatial analysis apps will (in different ways) exclude phenotypes labeled "unassigned", as the Spatial Interaction Tool does when it loads the settings from the Phenotyping app and "Custom" phenotyping is selected.*

Note that the "Phenotype Summary" table below and the "Phenotype Plot" to the left get dynamically updated as phenotype assignments are made.
