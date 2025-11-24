# Import relevant libraries.
import streamlit as st
from fast_neighborhood_profiles import main as fnp_main
import polars as pl
import plotly.express as px
import numpy as np
from functools import partial

# Define session state key prefixes.
ST_KEY_PREFIX = "study_results.py__"
ST_KEY_PREFIX_PHENOTYPE = "phenotype.py__"
ST_KEY_PREFIX_SUMAP = "run_spatial_umap.py__"


def get_selected_indices(selected_handle):
    both_handles = {"umap", "real_space"}
    other_handle = (both_handles - {selected_handle}).pop()
    selection = st.session_state[ST_KEY_PREFIX + f"{selected_handle}_plot__do_not_persist"]
    if "selection" in selection and "points" in selection["selection"] and selection["selection"]["points"]:
        points_list = selection["selection"]["points"]
        indices = [point["customdata"][4] for point in points_list]  # Note this means that if the "index" column is added to the plot data when calling main.plot_image_from_frame(), it must be the very first custom_column, i.e., at position 4 (0-based indexing).
        st.session_state[ST_KEY_PREFIX + "selected_indices_for_" + other_handle] = indices
    else:
        st.session_state[ST_KEY_PREFIX + "selected_indices_for_" + other_handle] = []


# Function to create a line plot with multiple series.
def line_plot_with_series(data, labels_axis_0, labels_axis_1, axis_0_name="Distance bin (µm)", axis_1_name="Phenotype", value_name="Mean density", color_map=None):
    # Note the axes here refer to the axes of data, not the plot axes.

    # Create a DataFrame from the data.
    data = pl.DataFrame(data, schema=labels_axis_1).with_columns(
        pl.Series(axis_0_name, labels_axis_0)
    )

    # Unpivot the DataFrame for plotting.
    data_unpivoted = data.unpivot(
        index=axis_0_name,
        on=labels_axis_1,
        variable_name=axis_1_name,
        value_name=value_name,
    )

    # Define color map if not provided.
    if not color_map:
        colors = px.colors.qualitative.Plotly
        unique_labels = sorted(labels_axis_1)
        color_map = {label: colors[i % len(colors)] for i, label in enumerate(unique_labels)}

    # Create the line plot.
    fig = px.line(
        data_unpivoted,
        x=axis_0_name,
        y=value_name,
        color=axis_1_name,
        color_discrete_map=color_map,
        markers=True,
    )

    # Update layout for clarity.
    fig.update_layout(
        xaxis_title=axis_0_name,
        yaxis_title=value_name,
        legend_title=axis_1_name,
    )

    # Return the figure.
    return fig


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

    # Grab values we'll need soon.
    image_colname = "TMA_core_id"
    unique_image_ids = st.session_state[ST_KEY_PREFIX_PHENOTYPE + "unique_image_ids"]
    phenotype_color_map = st.session_state[ST_KEY_PREFIX_PHENOTYPE + "phenotype_color_map"]
    unique_labels = st.session_state[ST_KEY_PREFIX_PHENOTYPE + "unique_labels"]  # should be correct, i.e. spatial_umap.species, i.e. lf_phenotyped.select(pl.col("label").unique().sort()).collect().to_series().to_list()
    dist_bin_um_list = st.session_state[ST_KEY_PREFIX_SUMAP + "dist_bin_um_list"]  # should be correct
    spatial_umap = st.session_state[ST_KEY_PREFIX_SUMAP + "spatial_umap"]

    # In the first of two columns...
    main_columns = st.columns(2)
    with main_columns[0]:
    
        # Allow the user to select which images to plot.
        st.session_state.setdefault(ST_KEY_PREFIX + "selected_images_to_plot", unique_image_ids)
        selected_images_to_plot = st.multiselect("Select images whose UMAP to plot:", options=unique_image_ids, key=ST_KEY_PREFIX + "selected_images_to_plot")

        # Allow the user to select marker size.
        st.session_state.setdefault(ST_KEY_PREFIX + "marker_size_umap", 3)
        marker_size_umap = st.slider("Marker size:", min_value=2, max_value=10, key=ST_KEY_PREFIX + "marker_size_umap")

        if ST_KEY_PREFIX + "selected_indices_for_umap" in st.session_state and st.session_state[ST_KEY_PREFIX + "selected_indices_for_umap"]:
            selected_indices = st.session_state[ST_KEY_PREFIX + "selected_indices_for_umap"]
        else:
            selected_indices = []
        st.write(f"Number of selected points in real space: {len(selected_indices):_}")

        # Plot the UMAP with selectable points.
        fig = fnp_main.plot_image_from_frame(lf_indexed, image_colname=image_colname, selected_images=selected_images_to_plot, marker_size=marker_size_umap, xcol="umap_1", ycol="umap_2", color_col="Lineage", custom_columns=["index"], color_map=phenotype_color_map, highlight_indices=selected_indices)
        fig.update_layout(uirevision="static")  # this doesn't seem to be honored; investigate in the future
        st.plotly_chart(fig, on_select=partial(get_selected_indices, selected_handle="umap"), selection_mode=("points", "box", "lasso"), key=ST_KEY_PREFIX + "umap_plot__do_not_persist")

    # In the second of two columns...
    with main_columns[1]:
        if selected_images_to_plot:
            with st.container(horizontal=True, vertical_alignment="bottom"):

                if ST_KEY_PREFIX + "selected_image_to_plot" in st.session_state and st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"] not in selected_images_to_plot:
                    del st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"]
                st.session_state.setdefault(ST_KEY_PREFIX + "selected_image_to_plot", selected_images_to_plot[0])
                selected_image_to_plot = st.selectbox("Select image to plot:", options=selected_images_to_plot, key=ST_KEY_PREFIX + "selected_image_to_plot")

                st.button("Previous", on_click=lambda: st.session_state.update({ST_KEY_PREFIX + "selected_image_to_plot": selected_images_to_plot[max(0, selected_images_to_plot.index(st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"]) - 1)]}), disabled=(st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"] == selected_images_to_plot[0]))

                st.button("Next", on_click=lambda: st.session_state.update({ST_KEY_PREFIX + "selected_image_to_plot": selected_images_to_plot[min(len(selected_images_to_plot) - 1, selected_images_to_plot.index(st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"]) + 1)]}), disabled=(st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"] == selected_images_to_plot[-1]))

            st.session_state.setdefault(ST_KEY_PREFIX + "marker_size_real_space", 3)
            marker_size_real_space = st.slider("Marker size:", min_value=2, max_value=10, key=ST_KEY_PREFIX + "marker_size_real_space")

            if ST_KEY_PREFIX + "selected_indices_for_real_space" in st.session_state and st.session_state[ST_KEY_PREFIX + "selected_indices_for_real_space"]:
                selected_indices = st.session_state[ST_KEY_PREFIX + "selected_indices_for_real_space"]
            else:
                selected_indices = []
            st.write(f"Number of selected points in UMAP: {len(selected_indices):_}")

            st.session_state.setdefault(ST_KEY_PREFIX + "display_only_real_space_coords_with_umap_coords", False)
            display_only_real_space_coords_with_umap_coords = st.checkbox("Display only real space coords with UMAP coords", key=ST_KEY_PREFIX + "display_only_real_space_coords_with_umap_coords")
            if display_only_real_space_coords_with_umap_coords:
                lf_indexed = lf_indexed.filter(pl.col("umap_test"))

            fig = fnp_main.plot_image_from_frame(lf_indexed, image_colname="TMA_core_id", xcol="Xcor", ycol="Ycor", color_col="Lineage", selected_images=[selected_image_to_plot], marker_size=marker_size_real_space, highlight_indices=selected_indices, custom_columns=["index"], color_map=phenotype_color_map)
            fig.update_layout(uirevision="static")  # this doesn't seem to be honored; investigate in the future
            st.plotly_chart(fig, on_select=partial(get_selected_indices, selected_handle="real_space"), selection_mode=("points", "box", "lasso"), key=ST_KEY_PREFIX + "real_space_plot__do_not_persist")
            if not display_only_real_space_coords_with_umap_coords:
                st.write("Keep in mind that not every point in real space was used for UMAP inference. So while selecting points in UMAP space will render the same number of selections in real space (over all the images), selecting points in real space will often render fewer selections in UMAP space. However, selecting points in real space still allows you to faithfully see their neighborhood profiles below.")

    # # If there are selected points...
    # if selected_indices:

    #     # Obtain from it the mean density for the selected points.
    #     density = spatial_umap.density[selected_indices, :, :]
    #     density_mean = density.mean(axis=0, dtype=np.float32)

    #     # Plot the neighborhood profiles.
    #     fig = line_plot_with_series(density_mean, dist_bin_um_list, unique_labels, axis_0_name="Distance bin (µm)", axis_1_name="Phenotype", value_name="Mean density", color_map=phenotype_color_map)
    #     st.plotly_chart(fig)


# Run the main function if this script is executed.
if __name__ == "__main__":
    main()
