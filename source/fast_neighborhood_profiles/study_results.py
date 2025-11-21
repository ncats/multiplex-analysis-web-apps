import streamlit as st
from fast_neighborhood_profiles import main as fnp_main
import polars as pl
import plotly.express as px

ST_KEY_PREFIX = "study_results.py__"
ST_KEY_PREFIX_PHENOTYPE = "phenotype.py__"
ST_KEY_PREFIX_SUMAP = "run_spatial_umap.py__"


def main():

    # Ensure the phenotyped lazyframe is ready for usage.
    if not ("LAZYFRAMES" in st.session_state and "sumap_cells" in st.session_state["LAZYFRAMES"]):
        st.warning("Please run the spatial UMAP analysis (at left).")
        return

    # Get the main lazyframe from session state.
    lf = st.session_state["LAZYFRAMES"]["sumap_cells"]["lf"]

    # sumap_cells_df = 

    lf_indexed = lf.with_row_index(name="index1")  # These are indices *after* potentially dropping entire images. They are consistent with spatial_umap.cells and spatial_umap.density.

    main_columns = st.columns(2)

    with main_columns[0]:
    
        # key = ST_KEY_PREFIX + "umap_plot"
        image_colname = "TMA_core_id"
        unique_image_ids = st.session_state[ST_KEY_PREFIX_PHENOTYPE + "unique_image_ids"]
        with st.container(horizontal=True, vertical_alignment="bottom"):
            st.session_state.setdefault(ST_KEY_PREFIX + "selected_images_to_plot", unique_image_ids)
            selected_images_to_plot = st.multiselect("Select images whose UMAP to plot:", options=unique_image_ids, key=ST_KEY_PREFIX + "selected_images_to_plot")
        marker_size = st.slider("Marker size:", min_value=1, max_value=10, value=1)
        selection = st.plotly_chart(fnp_main.plot_image_from_frame(lf_indexed, image_colname=image_colname, selected_images=selected_images_to_plot, marker_size=marker_size, xcol="umap_1", ycol="umap_2", color_col="Lineage", existing_index_columns=["index1"]), on_select="rerun", selection_mode=("points", "box", "lasso"))

        # st.write(lf_indexed.filter(pl.col(image_colname).is_in(selected_images_to_plot)).describe())

    with main_columns[1]:

        # st.write(selection)

        if "selection" in selection and "points" in selection["selection"] and selection["selection"]["points"]:
            points_list = selection["selection"]["points"]
            indices = [point["customdata"][2] for point in points_list]
            # st.write(indices)
            spatial_umap = st.session_state[ST_KEY_PREFIX_SUMAP + "spatial_umap"]
            # st.write(spatial_umap.cells.iloc[indices])
            # st.write(spatial_umap.density.shape)
            # st.write(type(indices))
            # import numpy as np
            # density = spatial_umap.density[indices, :, :]
            # st.write(density.reshape((len(indices), -1)))
            # st.write(density.mean(axis=0))

            density = spatial_umap.density[indices, :, :]
            density_mean = density.mean(axis=0)
            distance_bins = st.session_state[ST_KEY_PREFIX_SUMAP + "dist_bin_um_list"]
            phenotype_labels = st.session_state[ST_KEY_PREFIX_PHENOTYPE + "unique_labels"]

            density_df = pl.DataFrame(density_mean, schema=phenotype_labels).with_columns(
                pl.Series("distance_um", distance_bins)
            )
            density_long = density_df.unpivot(
                index="distance_um",
                on=phenotype_labels,
                variable_name="phenotype",
                value_name="density",
            )

            fig = px.line(
                density_long.to_pandas(),
                x="distance_um",
                y="density",
                color="phenotype",
                markers=True,
                template="plotly_white",
            )
            fig.update_layout(
                xaxis_title="Distance bin (µm)",
                yaxis_title="Mean density",
                legend_title="Phenotype",
            )

            st.plotly_chart(fig, use_container_width=True)

            # selected_point_indices = selection["selection"]["point_indices"]
            # st.write(f"Number of selected points: {len(selected_point_indices)}")
            # st.write("First 10 selected points:")
            # st.write(lf.filter(pl.arange(0, pl.count()).is_in(selected_point_indices)).collect())
            # st.write(lf[selected_point_indices].head(10).collect())

            # st.write(lf_indexed.with_row_index().filter(pl.col("index").is_in(selected_point_indices)).collect())

    # st.write(lf_indexed.head().collect())


    # lf = pl.scan_csv("data.csv").with_row_count("row_nr")
    # indices = [0, 5, 10]
    # subset = lf.filter(pl.col("row_nr").is_in(indices))
    # result = subset.collect()


if __name__ == "__main__":
    main()
