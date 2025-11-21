# Import relevant libraries.
import streamlit as st
from fast_neighborhood_profiles import main as fnp_main
import framework.utils as framework_utils
import os
import polars as pl

ST_KEY_PREFIX = "run_spatial_umap.py__"
ST_KEY_PREFIX_PHENOTYPE = "phenotype.py__"


# Format the lazyframe, optionally sample it, and convert to a polars dataframe.
def format_lazyframe(lf, sample_size=None, sample_seed=42):

    # Load in cells and patient data.
    lf = (
        lf
        .rename({"Image ID_(standardized)": "TMA_core_id", "Centroid X (µm)_(standardized)": "Xcor", "Centroid Y (µm)_(standardized)": "Ycor", "label": "Lineage"})
        .select(pl.col(["TMA_core_id", "Xcor", "Ycor", "Lineage"]))
        )

    # Load in cells and patient data. Sampling will aid in faster testing and development. The sorting after the sampling is crucial to ensure consistent ordering.
    if sample_size is None:
        pldf = (
            lf
            .sort(by="TMA_core_id")
            .collect()
            )
    else:
        pldf = (
            lf
            .collect()
            .sample(n=sample_size, seed=sample_seed)
            .sort(by="TMA_core_id")
            )
        
    return pldf


# Define the main function.
def main():

    # Ensure the phenotyped lazyframe is ready for usage.
    if not (
        ("LAZYFRAMES" in st.session_state)
        and ("marker_phenotyping" in st.session_state["LAZYFRAMES"])
        and (os.path.exists(os.path.join(framework_utils.session_dir(), st.session_state["LAZYFRAMES"]["unified_input_file"]["input_params"]["local_filepath"])))
        ):
        st.warning("Please perform phenotyping first (at left).")
        return

    # Get the main lazyframe from session state.
    lf_phenotyped = st.session_state["LAZYFRAMES"]["marker_phenotyping"]["lf"]

    # Allow the user to reset the algorithm parameters to defaults.
    if st.button("Reset defaults"):
        widget_keys = ["de_min_coords", "dist_bin_um_list", "custom_areas", "area_downsample", "area_threshold", "keep_images_with_too_little_data", "n", "train_sample_frac", "test_sample_frac", "set_seed_for_train_test_split", "set_cpu_pool_size", "cpu_pool_size"]
        for key in widget_keys:
            del st.session_state[ST_KEY_PREFIX + key]

    settings_columns = st.columns(2)
    with settings_columns[0]:

        # Set whether to de-min the coordinates.
        key = ST_KEY_PREFIX + "de_min_coords"
        st.session_state.setdefault(key, True)
        de_min_coords = st.checkbox("Shift each image's coordinates to origin (highly recommended!)", key=key, help="Giraldo et. al. did not do this, but we highly recommend this for stability of the algorithm and minimal dropped cells/images.")

        # Set (in the future) distance bins.
        key = ST_KEY_PREFIX + "dist_bin_um_list"
        st.session_state.setdefault(key, [25, 50, 100, 150, 200])
        dist_bin_um_list = st.session_state[key]
        st.write(f"Distance bins (µm) (not editable *yet*!): `{dist_bin_um_list}`")

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

    with settings_columns[1]:

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

    sumap_cell_file_format = "parquet"
    if st.button("Run spatial UMAP"):
        with st.spinner("Running spatial UMAP..."):

            pldf_phenotyped = format_lazyframe(lf_phenotyped, sample_size=None, sample_seed=42)  # Not making these two parameters editable as haven't used for a while.
            unique_labels = st.session_state[ST_KEY_PREFIX_PHENOTYPE + "unique_labels"]
            topdir = framework_utils.session_dir()
            subdir = os.path.join("output", "spatial_umap")

            spatial_umap, _ = fnp_main.generate_umap(pldf_phenotyped, unique_labels, dist_bin_um_list=dist_bin_um_list, area_downsample=area_downsample, um_per_px=1, cpu_pool_size=cpu_pool_size, topdir=topdir, subdir=subdir, counts_method="andrew", area_threshold=area_threshold, custom_areas=custom_areas, seed_for_train_test_split=seed_for_train_test_split, n=n, keep_images_with_too_little_data=keep_images_with_too_little_data, train_sample_frac=train_sample_frac, test_sample_frac=test_sample_frac, de_min_coords=de_min_coords)

            st.session_state[ST_KEY_PREFIX + "spatial_umap"] = spatial_umap

            sumap_cells_filepath = os.path.join(framework_utils.session_dir(), "input", f"sumap_cells.{sumap_cell_file_format}")
            if os.path.exists(sumap_cells_filepath):
                os.remove(sumap_cells_filepath)

    if ST_KEY_PREFIX + "spatial_umap" not in st.session_state:
        st.info("Please run spatial UMAP first.")
        return

    st.success("Spatial UMAP results are ready.")

    spatial_umap = st.session_state[ST_KEY_PREFIX + "spatial_umap"]

    if not os.path.exists(os.path.join(framework_utils.session_dir(), "input", f"sumap_cells.{sumap_cell_file_format}")):
        with st.spinner("Saving spatial UMAP results to file..."):
            fnp_main.save_pandas_df_to_file(spatial_umap.cells, handle="sumap_cells", file_format=sumap_cell_file_format, topdir=framework_utils.session_dir(), subdir="input")


if __name__ == "__main__":
    main()
