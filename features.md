---
nav_order: 2
title: Feature Summary
---

# Feature Summary

General:

* Intuitive, point-and-click graphical user interface
* Employment of multiple CPUs
* Aimed for users of all computational skill levels

Phenotyping:

* Multiaxial gating of numeric or categorical data, with optional batch normalization preprocessing
* Multiple methods for phenotyping on pre-thresholded markers: "species", "marker", or "custom"

Spatial analysis - Spatial Interaction Tool:

* "Patching" of slides into regions of interest
* Multiple methods for calculating the degree of interaction between cell phenotypes: Poisson (radius), permutation (radius), or permutation (k-nearest neighbors)
* Averaging of the interactions over whole slides or individual tissue types (tumor, stroma, necrosis, etc.)

Spatial analysis - Neighborhood Profiler:

* UMAP decomposition of counts cells of different phenotypes at multiple distances from each cell
* Clustering of the UMAP components to determine a discrete number of "neighborhood profiles" in the dataset
* Exploratory visualizations to inspect and compare neighborhood makeups and to investigate potential correlations with other features or response variables in the dataset
