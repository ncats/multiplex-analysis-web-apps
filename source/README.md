# Multiplex Analysis Web Apps

Multiplex Analysis Web Apps (MAWA) is a collection of tools developed by data scientists at NIH/NCI/CBIIT to support researchers utilizing multiplex imaging techniques. This suite of apps employs a Streamlit interface and for NIH users can be accessed using NIDAP free of charge.

## Local installation and execution (working on both Ubuntu and Windows as of 2024-02-08)

Remember, at NIH, no installation is required; NIDAP should be used instead. NIH users should contact the authors below for more information.

```bash
conda update -n base -c conda-forge conda
conda env create -f environment-2024-02-08.yml
conda activate mawa-2024-02-08
streamlit run Multiplex_Analysis_Web_Apps.py
```

## Using MAWA

Each page or app within MAWA is located in the `pages` folder in this repo. The pages are as follows:

1. data_import_and_export.py
1. datafile_format_unifier.py
1. open_file.py
1. robust_scatter_plotter.py
1. thresholded_phenotyping.py
1. multiaxial_gating.py
1. adaptive_phenotyping.py
1. Pheno_Cluster_a.py
1. Pheno_Cluster_b.py
1. Tool_parameter_selection.py
1. Run_workflow.py
1. Display_individual_ROI_heatmaps.py
1. Display_average_heatmaps.py
1. Display_average_heatmaps_per_annotation.py
1. Display_ROI_P_values_overlaid_on_slides.py
1. Neighborhood_Profiles.py
1. UMAP_Analyzer.py
1. Clusters_Analyzer.py
1. radial_bins_plots.py
1. radial_profiles_analysis.py
1. preprocessing.py
1. memory_analyzer.py
1. results_transfer.py

## Authors

[Andrew Weisman](mailto:andrew.weisman@nih.gov)  
[Dante Smith](mailto:dante.smith@nih.gov)  
[Andrei Bombin](mailto:andrei.bombin@nih.gov)
