## Status

We are currently in the process of porting MAWA to scale efficiently to large datasets (~100M input rows).

We have completed the hard parts (high-performance library implementation and corresponding infrastructure) and are porting prior functionalities as priority dictates.

Everything that has been ported over so far is located in the "High-performance workflow" section on the left sidebar (topmost section). Note this includes:

* Improved app session management (saving/loading app sessions).
* Loading your dataset from a single unified input file.
* "Marker" phenotyping.
* Complete neighborhood profiles / spatial UMAP functionality. We have replaced some standard features (e.g., clustering on the UMAPs) with more effective, fully interactive analysis.

Next up on our list to port to high performance:

* "Species" phenotyping.
* Datafile unification.
* Custom color selection throughout the high-performance workflow.
* Robust scatter plots.

Please let us know if anything urgent is not on this list; we will prioritize accordingly!

## Suggested workflow

1. Use the standard datafile unifier (Datafile Unification page in the File Handling section) to create a unified datafile.
1. Save it back to the server as usual. Ensure it appears in the new "Load unified input file" page in the "High-performance workflow" section at left.
1. If you see it there, we want to stay in the high-performance workflow and start fresh. So, we'd suggest: (1) save your current session using the "Manage sessions" page and (2) press the "🧹 Reset app" button near the bottom of the left sidebar. This will clear the memory of the app session, and from here on out, as long as you stay in the "High-performance workflow" section, the app will only save to memory what is absolutely necessary.
1. Load in your unified datafile using the new "Load unified input file" page, and then step through the pages in the "High-performance workflow" section one-by-one. This will get you a complete neighborhood profiles / spatial UMAP workflow using "marker" phenotyping on the input datafile.

## Other tasks on our plate

1. Video demo of full high-performance workflow.
1. Enable loading of Ana's old archive.
1. Fix coloring on scatterplot for Lisa.
1. Remove lines from scatter plotter per Leandro's 11/25/25 email.
1. Anything important that we missed?

Please reach out to [Andrew Weisman](mailto:andrew.weisman@nih.gov) or [Andrei Bombin](mailto:andrei.bombin@nih.gov) with any questions, suggestions, or comments!
