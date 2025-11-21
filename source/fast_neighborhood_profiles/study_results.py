import streamlit as st
from fast_neighborhood_profiles import main as fnp_main
import polars as pl

ST_KEY_PREFIX = "study_results.py__"
ST_KEY_PREFIX_PHENOTYPE = "phenotype.py__"


def main():

    # Ensure the phenotyped lazyframe is ready for usage.
    if not ("LAZYFRAMES" in st.session_state and "sumap_cells" in st.session_state["LAZYFRAMES"]):
        st.warning("Please run the spatial UMAP analysis (at left).")
        return

    # Get the main lazyframe from session state.
    lf = st.session_state["LAZYFRAMES"]["sumap_cells"]["lf"]

    main_columns = st.columns(2)

    with main_columns[0]:
    
        # key = ST_KEY_PREFIX + "umap_plot"
        image_colname = "TMA_core_id"
        unique_image_ids = st.session_state[ST_KEY_PREFIX_PHENOTYPE + "unique_image_ids"]
        with st.container(horizontal=True, vertical_alignment="bottom"):
            st.session_state.setdefault(ST_KEY_PREFIX + "selected_images_to_plot", unique_image_ids)
            selected_images_to_plot = st.multiselect("Select images whose UMAP to plot:", options=unique_image_ids, key=ST_KEY_PREFIX + "selected_images_to_plot")
        marker_size = st.slider("Marker size:", min_value=1, max_value=10, value=1)
        selection = st.plotly_chart(fnp_main.plot_image_from_frame(lf, image_colname=image_colname, selected_images=selected_images_to_plot, marker_size=marker_size, xcol="umap_1", ycol="umap_2", color_col="Lineage"), on_select="rerun", selection_mode=("points", "box", "lasso"))

    with main_columns[1]:

        st.write(selection)

        if "selection" in selection and "point_indices" in selection["selection"] and selection["selection"]["point_indices"]:
            selected_point_indices = selection["selection"]["point_indices"]
            st.write(f"Number of selected points: {len(selected_point_indices)}")
            st.write("First 10 selected points:")
            # st.write(lf.filter(pl.arange(0, pl.count()).is_in(selected_point_indices)).collect())
            # st.write(lf[selected_point_indices].head(10).collect())

            st.write(lf.with_row_index().filter(pl.col("index").is_in(selected_point_indices)).collect())


    # lf = pl.scan_csv("data.csv").with_row_count("row_nr")
    # indices = [0, 5, 10]
    # subset = lf.filter(pl.col("row_nr").is_in(indices))
    # result = subset.collect()


if __name__ == "__main__":
    main()
