# Import relevant libraries.
import streamlit as st
from fast_neighborhood_profiles import main as fnp_main
import polars as pl
import plotly.express as px
import numpy as np

# Define session state key prefixes.
ST_KEY_PREFIX = "study_results.py__"
ST_KEY_PREFIX_PHENOTYPE = "phenotype.py__"
ST_KEY_PREFIX_SUMAP = "run_spatial_umap.py__"


# Function to create a line plot with multiple series.
def line_plot_with_series(data, labels_axis_0, labels_axis_1, axis_0_name="Distance bin (µm)", axis_1_name="Phenotype", value_name="Mean density"):
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

    # Create the line plot.
    fig = px.line(
        data_unpivoted,
        x=axis_0_name,
        y=value_name,
        color=axis_1_name,
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

    # In the first of two columns...
    main_columns = st.columns(2)
    with main_columns[0]:
    
        # Store some values we'll need soon.
        image_colname = "TMA_core_id"
        unique_image_ids = st.session_state[ST_KEY_PREFIX_PHENOTYPE + "unique_image_ids"]

        # Allow the user to select which images to plot.
        st.session_state.setdefault(ST_KEY_PREFIX + "selected_images_to_plot", unique_image_ids)
        selected_images_to_plot = st.multiselect("Select images whose UMAP to plot:", options=unique_image_ids, key=ST_KEY_PREFIX + "selected_images_to_plot")

        # Allow the user to select marker size.
        marker_size = st.slider("Marker size:", min_value=1, max_value=10, value=1)

        # Plot the UMAP with selectable points.
        selection = st.plotly_chart(fnp_main.plot_image_from_frame(lf_indexed, image_colname=image_colname, selected_images=selected_images_to_plot, marker_size=marker_size, xcol="umap_1", ycol="umap_2", color_col="Lineage", custom_columns=["index"]), on_select="rerun", selection_mode=("points", "box", "lasso"))

    # In the second of two columns...
    with main_columns[1]:

        # If there are selected points...
        if "selection" in selection and "points" in selection["selection"] and selection["selection"]["points"]:

            # Obtain the indices that align with the "index" column.
            points_list = selection["selection"]["points"]
            indices = [point["customdata"][2] for point in points_list]
        
            # Grab the spatial UMAP object from session state.
            spatial_umap = st.session_state[ST_KEY_PREFIX_SUMAP + "spatial_umap"]

            # Obtain from it the mean density for the selected points.
            density = spatial_umap.density[indices, :, :]
            density_mean = density.mean(axis=0, dtype=np.float32)

            # Get the labels for the columns in the mean density array.
            dist_bin_um_list = st.session_state[ST_KEY_PREFIX_SUMAP + "dist_bin_um_list"]  # should be correct
            unique_labels = st.session_state[ST_KEY_PREFIX_PHENOTYPE + "unique_labels"]  # should be correct, i.e. spatial_umap.species, i.e. lf_phenotyped.select(pl.col("label").unique().sort()).collect().to_series().to_list()

            # Plot the neighborhood profiles.
            fig = line_plot_with_series(density_mean, dist_bin_um_list, unique_labels, axis_0_name="Distance bin (µm)", axis_1_name="Phenotype", value_name="Mean density")
            st.plotly_chart(fig)


# Run the main function if this script is executed.
if __name__ == "__main__":
    main()
