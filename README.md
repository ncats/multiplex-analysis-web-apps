# Multiplex Analysis Web Apps

Multiplex Analysis Web Apps (MAWA) is a collection of tools developed by data scientists at NIH/NCI/CBIIT to support researchers utilizing multiplex imaging techniques. This suite of apps employs a Streamlit interface and for NIH users can be accessed using NIDAP free of charge.

## Local installation and execution (tested on Ubuntu)

Remember, at NIH, no installation is required; NIDAP should be used instead. NIH users should contact the authors below for more information.

```bash
conda update -n base -c conda-forge conda
conda env create -f environment-2024-02-08.yml
conda activate mawa-2024-02-08
streamlit run Multiplex_Analysis_Web_Apps.py
```

## Using MAWA

Each page or app within MAWA is located in the `pages` folder in this repo. The pages are as follows:

1. 01_data_import_and_export.py
1. multiaxial_gating.py
1. 02_phenotyping.py
1. 03a_Tool_parameter_selection.py
1. 03b_Run_workflow.py
1. 03c_Display_individual_ROI_heatmaps.py
1. 03d_Display_average_heatmaps.py
1. 03e_Display_average_heatmaps_per_annotation.py
1. 03f_Display_ROI_P_values_overlaid_on_slides.py
1. 04a_Neighborhood Profiles.py
1. 04b_UMAP Analyzer.py
1. 04c_Clusters Analyzer.py

## Authors

[Andrew Weisman](mailto:andrew.weisman@nih.gov)  
[Dante Smith](mailto:dante.smith@nih.gov)
