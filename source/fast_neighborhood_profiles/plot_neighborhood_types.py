# Import relevant libraries.
import streamlit as st
from fast_neighborhood_profiles import main as fnp_main
import polars as pl
from functools import partial

# Define session state key prefixes.
ST_KEY_PREFIX = "plot_neighborhood_types.py__"
ST_KEY_PREFIX_PHENOTYPE = "phenotype.py__"
ST_KEY_PREFIX_SUMAP = "run_spatial_umap.py__"
ST_KEY_PREFIX_ASSIGN = "assign_neighborhood_types.py__"


# Obtain the main indices from a selection on one of the scatter plots. Activate just the "other" plot and the neighborhood profile plot.
def get_neighborhood_type(selected_handle):
    lf = st.session_state["LAZYFRAMES"]["neighborhood_types"]["lf"]
    selection = st.session_state[ST_KEY_PREFIX + f"{selected_handle}_plot__do_not_persist"]
    st.session_state[ST_KEY_PREFIX + "selected_indices_for_neighborhood_profile"] = []
    if "selection" in selection and "points" in selection["selection"] and selection["selection"]["points"]:
        points_list = selection["selection"]["points"]
        indices = [point["customdata"][4] for point in points_list]  # Note this means that if the "index" column is added to the plot data when calling main.plot_image_from_frame(), it must be the very first custom_column, i.e., at position 4 (0-based indexing).
        selected_neighborhood_types = lf.filter(pl.col("sumap_cell_index").is_in(indices)).select(pl.col("neighborhood_type").unique()).collect().to_series().to_list()
        if len(selected_neighborhood_types) == 1:
            df = st.session_state[ST_KEY_PREFIX_ASSIGN + "df_reconstructed_selections"]
            st.session_state[ST_KEY_PREFIX + "selected_indices_for_neighborhood_profile"] = df.loc[df["label"] == selected_neighborhood_types[0], "sumap_cell_indices"].values[0]


# Main function.
def main():

    # Ensure the lazyframe is ready for usage.
    if not ("LAZYFRAMES" in st.session_state and "neighborhood_types" in st.session_state["LAZYFRAMES"]):
        st.warning("Please assign the neighborhood types (at left).")
        return

    # Get the main lazyframe from session state.
    lf = st.session_state["LAZYFRAMES"]["neighborhood_types"]["lf"]

    # Grab values we'll need downstream.
    image_colname = "TMA_core_id"
    unique_image_ids = st.session_state[ST_KEY_PREFIX_PHENOTYPE + "unique_image_ids"]
    phenotype_color_map = st.session_state[ST_KEY_PREFIX_PHENOTYPE + "phenotype_color_map"]
    neighborhood_type_color_map = st.session_state[ST_KEY_PREFIX_ASSIGN + "neighborhood_type_color_map"]
    unique_phenotypes = st.session_state[ST_KEY_PREFIX_PHENOTYPE + "unique_labels"]  # should be correct, i.e. spatial_umap.species, i.e. lf_phenotyped.select(pl.col("label").unique().sort()).collect().to_series().to_list()
    dist_bin_um_list = st.session_state[ST_KEY_PREFIX_SUMAP + "dist_bin_um_list"]  # should be correct
    spatial_umap = st.session_state[ST_KEY_PREFIX_SUMAP + "spatial_UMAP_results"]["spatial_umap"]
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

        # Allow the user to select whether to color by phenotype or neighborhood types.
        st.session_state.setdefault(ST_KEY_PREFIX + "umap_color_by", "Neighborhood type")
        umap_color_by = st.radio("Color UMAP by:", options=["Phenotype", "Neighborhood type"], key=ST_KEY_PREFIX + "umap_color_by", horizontal=True)

        # Plot the UMAP with selectable points.
        color_col_mapping = {"Phenotype": "Lineage", "Neighborhood type": "neighborhood_type"}
        color_map_mapping = {"Phenotype": phenotype_color_map, "Neighborhood type": neighborhood_type_color_map}
        fig = fnp_main.plot_image_from_frame(lf, image_colname=image_colname, selected_images=selected_images_to_plot, marker_size=marker_size_umap, xcol="umap_1", ycol="umap_2", color_col=color_col_mapping[umap_color_by], custom_columns=["sumap_cell_index", "input_index"], color_map=color_map_mapping[umap_color_by])
        fig.update_layout(uirevision="static")  # this doesn't seem to be honored; investigate in the future
        st.plotly_chart(fig, on_select=partial(get_neighborhood_type, selected_handle="umap"), selection_mode=("points"), key=ST_KEY_PREFIX + "umap_plot__do_not_persist")

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

            # Give the user the option to only plot real space points that were used for UMAP inference.
            st.session_state.setdefault(ST_KEY_PREFIX + "display_only_real_space_coords_with_umap_coords", False)
            display_only_real_space_coords_with_umap_coords = st.checkbox("Display only real space coords with UMAP coords", key=ST_KEY_PREFIX + "display_only_real_space_coords_with_umap_coords")
            if display_only_real_space_coords_with_umap_coords:
                lf = lf.filter(pl.col("umap_test"))

            # Allow the user to select whether to color by phenotype or neighborhood types.
            st.session_state.setdefault(ST_KEY_PREFIX + "real_space_color_by", "Neighborhood type")
            real_space_color_by = st.radio("Color real space by:", options=["Phenotype", "Neighborhood type"], key=ST_KEY_PREFIX + "real_space_color_by", horizontal=True)

            # Allow user to plot rectangles faithful to the object sizes, if possible.
            st.session_state.setdefault(ST_KEY_PREFIX + "plot_faithful_object_sizes", False)
            plot_faithful_object_sizes = st.checkbox("Plot faithful object sizes (if available)", key=ST_KEY_PREFIX + "plot_faithful_object_sizes")

            # Allow the user to select marker size.
            st.session_state.setdefault(ST_KEY_PREFIX + "marker_size_real_space", 3)
            marker_size_real_space = st.slider("Marker size:", min_value=2, max_value=10, key=ST_KEY_PREFIX + "marker_size_real_space", disabled=plot_faithful_object_sizes)

            # Plot the real space with selectable points.
            color_col_mapping = {"Phenotype": "Lineage", "Neighborhood type": "neighborhood_type"}
            color_map_mapping = {"Phenotype": phenotype_color_map, "Neighborhood type": neighborhood_type_color_map}
            frame_with_faithful_columns = st.session_state["LAZYFRAMES"]["unified_input_file"]["lf"] if plot_faithful_object_sizes else None
            fig = fnp_main.plot_image_from_frame(lf, image_colname="TMA_core_id", xcol="Xcor", ycol="Ycor", color_col=color_col_mapping[real_space_color_by], selected_images=[selected_image_to_plot], marker_size=marker_size_real_space, custom_columns=["sumap_cell_index", "input_index"], color_map=color_map_mapping[real_space_color_by], plot_faithful_object_sizes=plot_faithful_object_sizes, frame_with_faithful_columns=frame_with_faithful_columns)
            fig.update_layout(uirevision="static")  # this doesn't seem to be honored; investigate in the future
            st.plotly_chart(fig, on_select=partial(get_neighborhood_type, selected_handle="real_space"), selection_mode="points", key=ST_KEY_PREFIX + "real_space_plot__do_not_persist")

    # If there are selected points...
    if selected_indices_for_neighborhood_profile:

        # Write the number of selected points for the neighborhood profile plot.
        with st.container(horizontal=True):
            st.write(f"Last number of selected points for neighborhood profile: {len(selected_indices_for_neighborhood_profile):_}")
            st.button("Clear selection", on_click=lambda: st.session_state.update({ST_KEY_PREFIX + "selected_indices_for_neighborhood_profile": []}), key=ST_KEY_PREFIX + "clear_neighborhood_profile_selection_button__do_not_persist")

        # Allow the user to select the neighborhood profile plot type.
        st.session_state.setdefault(ST_KEY_PREFIX + "neighborhood_profile_plot_type", "line")
        neighborhood_profile_plot_type = st.radio("Select plot type:", options=["line", "box", "violin"], key=ST_KEY_PREFIX + "neighborhood_profile_plot_type")

        # Grab the density for the selected indices for all distance bins and all phenotypes.
        density_counts_per_sq_mm = spatial_umap.density[selected_indices_for_neighborhood_profile, :, :] * 1e6

        # Plot the neighborhood profiles.
        fig, extra_return_info = fnp_main.plot_neighborhood_profile(density_counts_per_sq_mm, neighborhood_profile_plot_type, dist_bin_um_list, unique_phenotypes, axis_1_name="Distance bin (µm)", axis_2_name="Phenotype", value_name="Density (count/mm²)", color_map=phenotype_color_map)
        st.plotly_chart(fig)
        if extra_return_info:
            st.write(extra_return_info)


# Run the main function if this script is executed.
if __name__ == "__main__":
    main()
