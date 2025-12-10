# Import relevant libraries.
import streamlit as st
from fast_neighborhood_profiles import main as fnp_main
import framework.analysis_framework as analysis_framework
import framework.utils as framework_utils
import streamlit_dataframe_editor as sde
import pandas as pd

# Define constants.
ST_KEY_PREFIX = "run_spatial_umap.py__"
ST_KEY_PREFIX_PHENOTYPE = "phenotype.py__"


# Get a color map for True/False values.
@st.cache_data()
def get_true_false_color_map():
    return fnp_main.get_true_false_color_map()


# Define the main function.
def main():

    # Ensure the phenotyped lazyframe is ready for usage.
    if not ("LAZYFRAMES" in st.session_state and "phenotyped" in st.session_state["LAZYFRAMES"]):
        st.warning("Please perform phenotyping (at left).")
        return

    # In the first of two main columns...
    main_columns = st.columns(2)
    with main_columns[0]:
        st.header("Analysis parameters")

        # Allow the user to reset the algorithm parameters to defaults.
        if st.button("Reset defaults"):
            widget_keys = ["de_min_coords", "dist_bin_um_list", "custom_areas", "area_downsample", "area_threshold", "keep_images_with_too_little_data", "n", "train_sample_frac", "test_sample_frac", "set_seed_for_train_test_split", "set_cpu_pool_size", "cpu_pool_size"]
            for key in widget_keys:
                del st.session_state[ST_KEY_PREFIX + key]

        # Set whether to de-min the coordinates.
        key = ST_KEY_PREFIX + "de_min_coords"
        st.session_state.setdefault(key, True)
        de_min_coords = st.checkbox("Shift each image's coordinates to origin (highly recommended!)", key=key, help="Giraldo et. al. did not do this, but we highly recommend this for stability of the algorithm and minimal dropped cells/images.")

        # Set distance bins.
        key = ST_KEY_PREFIX + "de_dist_bin_um"
        if key not in st.session_state:
            st.session_state[key] = sde.DataframeEditor(df_name=ST_KEY_PREFIX + "df_dist_bin_um", default_df_contents=pd.DataFrame({"dist_bin_um": [25, 50, 100, 150, 200]}))
        st.write("Distance bins (µm):")
        st.session_state[ST_KEY_PREFIX + "de_dist_bin_um"].dataframe_editor(reset_data_editor_button_text='Reset bins to defaults')
        key = ST_KEY_PREFIX + "dist_bin_um_list"
        st.session_state[key] = st.session_state[ST_KEY_PREFIX + "de_dist_bin_um"].reconstruct_edited_dataframe()["dist_bin_um"].to_list()
        dist_bin_um_list = st.session_state[key]

        # Set whether to use custom areas.
        key = ST_KEY_PREFIX + "custom_areas"
        st.session_state.setdefault(key, False)
        custom_areas = st.checkbox("Calculate a custom area for each cell", key=key, help="Giraldo et. al. calculated custom cell areas, but this can discard significant tissue when not analyzing TMA cores.")

        # Set area downsample.
        key = ST_KEY_PREFIX + "area_downsample"
        st.session_state.setdefault(key, 0.2)
        area_downsample = st.number_input("Area downsample:", min_value=0.0001, max_value=1.0, step=0.05, key=key, disabled=not st.session_state[ST_KEY_PREFIX + "custom_areas"])

        # Set area threshold.
        key = ST_KEY_PREFIX + "area_threshold"
        st.session_state.setdefault(key, 0.8)
        area_threshold = st.number_input("Area threshold:", min_value=0.0001, max_value=1.0, step=0.05, key=key, disabled=not st.session_state[ST_KEY_PREFIX + "custom_areas"])

        # Set whether to keep images with too little data.
        key = ST_KEY_PREFIX + "keep_images_with_too_little_data"
        st.session_state.setdefault(key, True)
        keep_images_with_too_little_data = st.checkbox("Keep images with too little data", key=key, help="Giraldo et. al. did not do this, i.e., they dropped entire images with too little non-filtered-out data.")

        # Set param related to the minimum number of cells per image before the image is discarded for UMAP (min=int((train_sample_frac+test_sample_frac)*n)).
        key = ST_KEY_PREFIX + "n"
        st.session_state.setdefault(key, 2500)
        n = st.number_input("Param related to min. # of cells per image before the image is discarded for UMAP (`min = int((train_sample_frac + test_sample_frac) * n)`):", min_value=1, step=1, key=key)

        # Set UMAP train sample fraction.
        key = ST_KEY_PREFIX + "train_sample_frac"
        st.session_state.setdefault(key, 1.0)
        train_sample_frac = st.number_input("UMAP training sample fraction:", min_value=0.0001, max_value=1.0, step=0.01, key=key)

        # Set UMAP test sample fraction.
        key = ST_KEY_PREFIX + "test_sample_frac"
        st.session_state.setdefault(key, 1.0)
        test_sample_frac = st.number_input("UMAP inference sample fraction:", min_value=0.0001, max_value=1.0, step=0.01, key=key)

        # Set seed for test-train split.
        key = ST_KEY_PREFIX + "set_seed_for_train_test_split"
        st.session_state.setdefault(key, False)
        set_seed_for_train_test_split = st.checkbox("Set seed for UMAP data splitting", key=key, help="Giraldo et. al. set a seed for reproducibility, but conclusions should be independent of sampling, so we do not recommend this.")
        if set_seed_for_train_test_split:
            seed_for_train_test_split = 54321
        else:
            seed_for_train_test_split = None

        # Whether to specify the number of threads to use.
        key = ST_KEY_PREFIX + "set_cpu_pool_size"
        st.session_state.setdefault(key, False)
        set_cpu_pool_size = st.checkbox("Specify number of CPU threads to use", key=key, help="If unchecked, algorithm will use all available CPUs.")

        # Set number of threads to use.
        key = ST_KEY_PREFIX + "cpu_pool_size"
        st.session_state.setdefault(key, 4)
        cpu_pool_size = st.number_input("Number of CPU threads to use:", min_value=1, step=1, key=key, disabled=not st.session_state[ST_KEY_PREFIX + "set_cpu_pool_size"])
        if not set_cpu_pool_size:
            cpu_pool_size = None

    # Get shortcuts to some variables.
    sumap_cell_file_format = "parquet"
    key = ST_KEY_PREFIX + "spatial_UMAP_results"
    unique_labels = st.session_state[ST_KEY_PREFIX_PHENOTYPE + "unique_labels"]

    # Assemble the inputs (less the polars dataframe, to be part of the job preprocessing) to the spatial UMAP analysis.
    # We need all lazyframes so we can rebuild what we need inside the job worker, and we will parse down from that, once we've regenerated what's needed, inside the spatial UMAP wrapper. Then we'll input only what we need (the single input lazyframe, plus parameters) into the core spatial UMAP function.
    lazyframes = {k: st.session_state["LAZYFRAMES"][k] for k in ["unified_input_file", "phenotyped"] if k in st.session_state["LAZYFRAMES"]}
    inputs = dict(LAZYFRAMES=lazyframes, unique_labels=unique_labels, dist_bin_um_list=dist_bin_um_list, area_downsample=area_downsample, um_per_px=1, cpu_pool_size=cpu_pool_size, subdir="spatial_umap", counts_method="andrew", area_threshold=area_threshold, custom_areas=custom_areas, seed_for_train_test_split=seed_for_train_test_split, n=n, keep_images_with_too_little_data=keep_images_with_too_little_data, train_sample_frac=train_sample_frac, test_sample_frac=test_sample_frac, de_min_coords=de_min_coords, mp_start_method='forkserver')

    # Allow the user to run the spatial UMAP analysis asynchronously.
    analysis_framework.job_submission(
        job_name="spatial_umap",
        inputs=inputs,
        analysis_purpose="spatial UMAP",
        st_key_prefix=ST_KEY_PREFIX,
        # preprocess={
        #     "function": fnp_main.format_lazyframe,
        #     "args": dict(lf=st.session_state["LAZYFRAMES"]["phenotyped"]["lf"], sample_size=None, sample_seed=42),
        # }
    )

    # Ensure the job results are available in the session state.
    if key not in st.session_state:
        st.info("Please press the button above to generate spatial UMAP results.")
        return
    
    # Get a shortcut to the spatial UMAP results.
    spatial_umap = st.session_state[key]["spatial_umap"]

    # If the spatial UMAP job just completed, save the results to a lazyframe and store it in the session state.
    if "JOB_JUST_COMPLETED" in st.session_state and st.session_state["JOB_JUST_COMPLETED"] == "spatial_umap":
        params = dict(handle="sumap_cells", file_format=sumap_cell_file_format, index_column_name="sumap_cell_index")
        lf = fnp_main.save_and_load_pandas_df_to_lf(spatial_umap.cells, **params, topdir=framework_utils.session_dir())
        st.session_state["LAZYFRAMES"]["sumap_cells"] = {
            "lf": lf,
            "function_metadata": {"module_name": "fast_neighborhood_profiles.main", "qualpath": "save_and_load_pandas_df_to_lf"},
            "input_dataset": {"type": "pandas_df", "keys": (key, "spatial_umap", "cells")},
            "params": params,
            }
        del st.session_state["JOB_JUST_COMPLETED"]

    # Get a shortcut to the cells lazyframe.
    lf = st.session_state["LAZYFRAMES"]["sumap_cells"]["lf"]

    # In the second of two main columns...
    with main_columns[1]:
        st.header("Check images for which (if any) cells were dropped")

        # Write a note.
        st.write("`area_filter==True` means the cell was not filtered out by any custom areas calculation (if custom areas were used).")

        # Allow for image selection and plot it colored by area filter.
        image_colname = "TMA_core_id"
        unique_image_ids = st.session_state[ST_KEY_PREFIX_PHENOTYPE + "unique_image_ids"]
        with st.container(horizontal=True, vertical_alignment="bottom"):
            if ST_KEY_PREFIX + "selected_image_to_plot" in st.session_state and st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"] not in unique_image_ids:
                del st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"]
            st.session_state.setdefault(ST_KEY_PREFIX + "selected_image_to_plot", unique_image_ids[0])
            selected_image_to_plot = st.selectbox("Select image to plot:", options=unique_image_ids, key=ST_KEY_PREFIX + "selected_image_to_plot")
            st.button("Previous", on_click=lambda: st.session_state.update({ST_KEY_PREFIX + "selected_image_to_plot": unique_image_ids[max(0, unique_image_ids.index(st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"]) - 1)]}), disabled=(st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"] == unique_image_ids[0]))
            st.button("Next", on_click=lambda: st.session_state.update({ST_KEY_PREFIX + "selected_image_to_plot": unique_image_ids[min(len(unique_image_ids) - 1, unique_image_ids.index(st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"]) + 1)]}), disabled=(st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"] == unique_image_ids[-1]))
        marker_size = st.slider("Marker size:", min_value=2, max_value=10, value=5)
        st.plotly_chart(fnp_main.plot_image_from_frame(lf, image_colname=image_colname, selected_images=[selected_image_to_plot], marker_size=marker_size, xcol="Xcor", ycol="Ycor", color_col="area_filter", color_map=get_true_false_color_map(), custom_columns=["input_index", "sumap_cell_index"]))


# Run the main function if this script is executed.
if __name__ == "__main__":
    main()
