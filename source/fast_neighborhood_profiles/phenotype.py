# Import relevant libraries.
import streamlit as st
from fast_neighborhood_profiles import main as fnp_main
import polars as pl
import plotly.express as px

# Define constants.
ST_KEY_PREFIX = "phenotype.py__"
ST_KEY_PREFIX_LOAD = "load_unified_input_file.py__"


# Define the main function.
def main():

    # Ensure the main lazyframe is ready for usage.
    if not ("LAZYFRAMES" in st.session_state and "unified_input_file" in st.session_state["LAZYFRAMES"]):
        st.warning("Please load a unified input file (at left).")
        return

    # Get the main lazyframe from session state.
    lf = st.session_state["LAZYFRAMES"]["unified_input_file"]["lf"]

    # Create two main columns on the page.
    main_columns = st.columns([1/3, 2/3], border=1)
    with main_columns[0]:

        # Optionally add a suffix to exclude when detecting marker columns.
        key = ST_KEY_PREFIX + "exclusion_suffix"
        st.session_state.setdefault(key, "")
        exclusion_suffix = st.text_input("Exclusion suffix for marker columns:", key=key)

        # Button to get the marker columns.
        key = ST_KEY_PREFIX + "marker_columns"
        if st.button("Get marker columns"):
            st.session_state[key] = fnp_main.get_marker_columns(lf, exclusion_suffix=exclusion_suffix)

        # Ensure the marker columns are in session state.
        if key not in st.session_state:
            st.info("Please press the button above to obtain the marker columns.")
            return
        
        # Display the marker columns.
        marker_columns = st.session_state[key]
        st.write(f"Found {len(marker_columns)} marker columns:")
        st.write(marker_columns)

        # Allow the user to perform phenotyping.
        if st.button("Perform marker phenotyping"):
            params = {"marker_columns": marker_columns}
            lf_phenotyped = fnp_main.perform_marker_phenotyping_on_lazyframe(lf, **params)
            st.session_state["LAZYFRAMES"]["marker_phenotyping"] = {
                "lf": lf_phenotyped,
                "function": fnp_main.perform_marker_phenotyping_on_lazyframe,
                "input_dataset": {"type": "lf", "keys": ("unified_input_file",)},
                "params": params,
                "extras": None,
            }
            st.session_state[ST_KEY_PREFIX + "num_phenotyped_rows"] = lf_phenotyped.select(pl.len()).collect().item()
            st.session_state[ST_KEY_PREFIX + "unique_labels"] = lf_phenotyped.select(pl.col("label").unique().sort()).collect().to_series().to_list()
            st.session_state[ST_KEY_PREFIX + "unique_image_ids"] = lf_phenotyped.select(pl.col("Image ID_(standardized)").unique().sort()).collect().to_series().to_list()
            colors = px.colors.qualitative.Plotly
            st.session_state[ST_KEY_PREFIX + "phenotype_color_map"] = {label: colors[i % len(colors)] for i, label in enumerate(st.session_state[ST_KEY_PREFIX + "unique_labels"])}


        # Ensure the phenotyped lazyframe is in session state.
        if "marker_phenotyping" not in st.session_state["LAZYFRAMES"]:
            st.info("Please press the button above to perform marker phenotyping.")
            return
        
        # Display the number of rows in the phenotyped lazyframe.
        lf_phenotyped = st.session_state["LAZYFRAMES"]["marker_phenotyping"]["lf"]
        num_phenotyped_rows = st.session_state[ST_KEY_PREFIX + "num_phenotyped_rows"]
        unique_labels = st.session_state[ST_KEY_PREFIX + "unique_labels"]
        unique_image_ids = st.session_state[ST_KEY_PREFIX + "unique_image_ids"]
        information = f'''
        :small_orange_diamond: # of phenotyped rows: `{num_phenotyped_rows:_}`  
        :small_orange_diamond: Unique labels: `{unique_labels}`  
        :small_orange_diamond: # of unique images: `{len(unique_image_ids)}`  
        '''
        st.markdown(information)

    # Plot the phenotyped data.
    with main_columns[1]:
        image_colname = "Image ID_(standardized)"
        with st.container(horizontal=True, vertical_alignment="bottom"):
            if ST_KEY_PREFIX + "selected_image_to_plot" in st.session_state and st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"] not in unique_image_ids:
                del st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"]
            st.session_state.setdefault(ST_KEY_PREFIX + "selected_image_to_plot", unique_image_ids[0])
            selected_image_to_plot = st.selectbox("Select image to plot:", options=unique_image_ids, key=ST_KEY_PREFIX + "selected_image_to_plot")
            st.button("Previous", on_click=lambda: st.session_state.update({ST_KEY_PREFIX + "selected_image_to_plot": unique_image_ids[max(0, unique_image_ids.index(st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"]) - 1)]}), disabled=(st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"] == unique_image_ids[0]))
            st.button("Next", on_click=lambda: st.session_state.update({ST_KEY_PREFIX + "selected_image_to_plot": unique_image_ids[min(len(unique_image_ids) - 1, unique_image_ids.index(st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"]) + 1)]}), disabled=(st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"] == unique_image_ids[-1]))
        marker_size = st.slider("Marker size:", min_value=2, max_value=10, value=3)
        st.plotly_chart(fnp_main.plot_image_from_frame(lf_phenotyped, image_colname=image_colname, selected_images=[selected_image_to_plot], marker_size=marker_size, color_map=st.session_state[ST_KEY_PREFIX + "phenotype_color_map"], custom_columns=["input_index"]))


# Run the main function if this script is executed.
if __name__ == "__main__":
    main()
