# Import relevant libraries.
import streamlit as st
from fast_neighborhood_profiles import main as fnp_main
import polars as pl
import numpy as np
from functools import partial

# Define session state key prefixes.
ST_KEY_PREFIX = "study_results.py__"
ST_KEY_PREFIX_PHENOTYPE = "phenotype.py__"
ST_KEY_PREFIX_SUMAP = "run_spatial_umap.py__"


# Obtain the main indices from a selection on one of the scatter plots.
def get_selected_indices(selected_handle):
    both_handles = {"umap", "real_space"}
    other_handle = (both_handles - {selected_handle}).pop()
    selection = st.session_state[ST_KEY_PREFIX + f"{selected_handle}_plot__do_not_persist"]
    if "selection" in selection and "points" in selection["selection"] and selection["selection"]["points"]:
        points_list = selection["selection"]["points"]
        indices = [point["customdata"][4] for point in points_list]  # Note this means that if the "index" column is added to the plot data when calling main.plot_image_from_frame(), it must be the very first custom_column, i.e., at position 4 (0-based indexing).
        st.session_state[ST_KEY_PREFIX + "selected_indices_for_" + other_handle] = indices
        st.session_state[ST_KEY_PREFIX + "selected_indices_for_neighborhood_profile"] = indices
    else:
        st.session_state[ST_KEY_PREFIX + "selected_indices_for_" + other_handle] = []
        st.session_state[ST_KEY_PREFIX + "selected_indices_for_neighborhood_profile"] = []


# Main function.
def main():

    # Ensure the lazyframe is ready for usage.
    if not ("LAZYFRAMES" in st.session_state and "sumap_cells" in st.session_state["LAZYFRAMES"]):
        st.warning("Please run the spatial UMAP analysis (at left).")
        return

    # Get the main lazyframe from session state.
    lf = st.session_state["LAZYFRAMES"]["sumap_cells"]["lf"]

    # Add a row index to the lazyframe. These are indices *after* potentially dropping entire images. They are consistent with spatial_umap.cells and spatial_umap.density.
    lf_indexed = lf.with_row_index(name="index")  

    # Grab values we'll need downstream.
    image_colname = "TMA_core_id"
    unique_image_ids = st.session_state[ST_KEY_PREFIX_PHENOTYPE + "unique_image_ids"]
    phenotype_color_map = st.session_state[ST_KEY_PREFIX_PHENOTYPE + "phenotype_color_map"]
    unique_labels = st.session_state[ST_KEY_PREFIX_PHENOTYPE + "unique_labels"]  # should be correct, i.e. spatial_umap.species, i.e. lf_phenotyped.select(pl.col("label").unique().sort()).collect().to_series().to_list()
    dist_bin_um_list = st.session_state[ST_KEY_PREFIX_SUMAP + "dist_bin_um_list"]  # should be correct
    spatial_umap = st.session_state[ST_KEY_PREFIX_SUMAP + "spatial_umap"]
    selected_indices_for_umap = []
    if ST_KEY_PREFIX + "selected_indices_for_umap" in st.session_state and st.session_state[ST_KEY_PREFIX + "selected_indices_for_umap"]:
        selected_indices_for_umap = st.session_state[ST_KEY_PREFIX + "selected_indices_for_umap"]
    selected_indices_for_real_space = []
    if ST_KEY_PREFIX + "selected_indices_for_real_space" in st.session_state and st.session_state[ST_KEY_PREFIX + "selected_indices_for_real_space"]:
        selected_indices_for_real_space = st.session_state[ST_KEY_PREFIX + "selected_indices_for_real_space"]
    selected_indices_for_neighborhood_profile = []
    if ST_KEY_PREFIX + "selected_indices_for_neighborhood_profile" in st.session_state and st.session_state[ST_KEY_PREFIX + "selected_indices_for_neighborhood_profile"]:
        selected_indices_for_neighborhood_profile = st.session_state[ST_KEY_PREFIX + "selected_indices_for_neighborhood_profile"]

    # In the first of two columns...
    main_columns = st.columns(2)
    with main_columns[0]:
    
        # Allow the user to select which images to plot.
        st.session_state.setdefault(ST_KEY_PREFIX + "selected_images_to_plot", unique_image_ids)
        selected_images_to_plot = st.multiselect("Select images whose UMAP to plot:", options=unique_image_ids, key=ST_KEY_PREFIX + "selected_images_to_plot")

        # Allow the user to select marker size.
        st.session_state.setdefault(ST_KEY_PREFIX + "marker_size_umap", 3)
        marker_size_umap = st.slider("Marker size:", min_value=2, max_value=10, key=ST_KEY_PREFIX + "marker_size_umap")

        # Write the number of selected points in the UMAP. Remember it says _for_real_space even though the selection is done on the UMAP because it's the selection of points on the UMAP that will be highlighted *for* the real space plot.
        with st.container(horizontal=True):
            st.write(f"Number of selected points in UMAP: {len(selected_indices_for_real_space):_}")
            st.button("Clear selection", on_click=lambda: st.session_state.update({ST_KEY_PREFIX + "selected_indices_for_real_space": []}), key=ST_KEY_PREFIX + "clear_umap_selection_button__do_not_persist")

        # Plot the UMAP with selectable points.
        fig = fnp_main.plot_image_from_frame(lf_indexed, image_colname=image_colname, selected_images=selected_images_to_plot, marker_size=marker_size_umap, xcol="umap_1", ycol="umap_2", color_col="Lineage", custom_columns=["index"], color_map=phenotype_color_map, highlight_indices=selected_indices_for_umap)
        fig.update_layout(uirevision="static")  # this doesn't seem to be honored; investigate in the future
        st.plotly_chart(fig, on_select=partial(get_selected_indices, selected_handle="umap"), selection_mode=("points", "box", "lasso"), key=ST_KEY_PREFIX + "umap_plot__do_not_persist")

    # In the second of two columns...
    with main_columns[1]:
        if selected_images_to_plot:
            with st.container(horizontal=True, vertical_alignment="bottom"):

                # Image selection drop-down.
                if ST_KEY_PREFIX + "selected_image_to_plot" in st.session_state and st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"] not in selected_images_to_plot:
                    del st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"]
                st.session_state.setdefault(ST_KEY_PREFIX + "selected_image_to_plot", selected_images_to_plot[0])
                selected_image_to_plot = st.selectbox("Select image to plot:", options=selected_images_to_plot, key=ST_KEY_PREFIX + "selected_image_to_plot")

                # Previous button.
                st.button("Previous", on_click=lambda: st.session_state.update({ST_KEY_PREFIX + "selected_image_to_plot": selected_images_to_plot[max(0, selected_images_to_plot.index(st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"]) - 1)]}), disabled=(st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"] == selected_images_to_plot[0]))

                # Next button.
                st.button("Next", on_click=lambda: st.session_state.update({ST_KEY_PREFIX + "selected_image_to_plot": selected_images_to_plot[min(len(selected_images_to_plot) - 1, selected_images_to_plot.index(st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"]) + 1)]}), disabled=(st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"] == selected_images_to_plot[-1]))

            # Allow the user to select marker size.
            st.session_state.setdefault(ST_KEY_PREFIX + "marker_size_real_space", 3)
            marker_size_real_space = st.slider("Marker size:", min_value=2, max_value=10, key=ST_KEY_PREFIX + "marker_size_real_space")

            # Write the number of selected points in real space. Remember it says _for_umap even though the selection is done on real space because it's the selection of points in real space that will be highlighted *for* the UMAP plot.
            with st.container(horizontal=True):
                st.write(f"Number of selected points in real space: {len(selected_indices_for_umap):_}")
                st.button("Clear selection", on_click=lambda: st.session_state.update({ST_KEY_PREFIX + "selected_indices_for_umap": []}), key=ST_KEY_PREFIX + "clear_real_space_selection_button__do_not_persist")

            # Give the user the option to only plot real space points that were used for UMAP inference.
            st.session_state.setdefault(ST_KEY_PREFIX + "display_only_real_space_coords_with_umap_coords", False)
            display_only_real_space_coords_with_umap_coords = st.checkbox("Display only real space coords with UMAP coords", key=ST_KEY_PREFIX + "display_only_real_space_coords_with_umap_coords")
            if display_only_real_space_coords_with_umap_coords:
                lf_indexed = lf_indexed.filter(pl.col("umap_test"))

            # Plot the real space with selectable points.
            fig = fnp_main.plot_image_from_frame(lf_indexed, image_colname="TMA_core_id", xcol="Xcor", ycol="Ycor", color_col="Lineage", selected_images=[selected_image_to_plot], marker_size=marker_size_real_space, highlight_indices=selected_indices_for_real_space, custom_columns=["index"], color_map=phenotype_color_map)
            fig.update_layout(uirevision="static")  # this doesn't seem to be honored; investigate in the future
            st.plotly_chart(fig, on_select=partial(get_selected_indices, selected_handle="real_space"), selection_mode=("points", "box", "lasso"), key=ST_KEY_PREFIX + "real_space_plot__do_not_persist")

            # Display a note about selecting points.
            if not display_only_real_space_coords_with_umap_coords:
                st.write("Keep in mind that not every point in real space was used for UMAP inference. So while selecting points in UMAP space will render the same number of selections in real space (over all the images), selecting points in real space will often render fewer selections in UMAP space. However, selecting points in real space still allows you to faithfully see their neighborhood profiles below.")

    # If there are selected points...
    if selected_indices_for_neighborhood_profile:

        # Write the number of selected points for the neighborhood profile plot.
        with st.container(horizontal=True):
            st.write(f"Last number of selected points for neighborhood profile: {len(selected_indices_for_neighborhood_profile):_}")
            st.button("Clear selection", on_click=lambda: st.session_state.update({ST_KEY_PREFIX + "selected_indices_for_neighborhood_profile": []}), key=ST_KEY_PREFIX + "clear_neighborhood_profile_selection_button__do_not_persist")

        # Obtain from it the mean density for the selected points.
        density = spatial_umap.density[selected_indices_for_neighborhood_profile, :, :]

        # Plot the neighborhood profiles.
        # fig = fnp_main.line_plot_with_series(density.mean(axis=0, dtype=np.float32), dist_bin_um_list, unique_labels, axis_0_name="Distance bin (µm)", axis_1_name="Phenotype", value_name="Mean density", color_map=phenotype_color_map)
        fig = fnp_main.violin_plot_with_series(density, dist_bin_um_list, unique_labels, axis_1_name="Distance bin (µm)", axis_2_name="Phenotype", value_name="Density", color_map=phenotype_color_map)
        st.plotly_chart(fig)


# Run the main function if this script is executed.
if __name__ == "__main__":
    main()
