# Import relevant libraries.
import streamlit as st
from fast_neighborhood_profiles import main as fnp_main
import os
import framework.utils as framework_utils
import polars as pl

# Define constant.
ST_KEY_PREFIX = "phenotype.py__"
ST_KEY_PREFIX_LOAD = "load_unified_input_file.py__"


# Get the marker column names from the lazyframe.
def get_marker_columns(lf, exclusion_suffix=""):
    if exclusion_suffix:
        marker_columns = [column for column in lf.collect_schema().names() if column.startswith("Phenotype_(standardized) ") and not column.endswith(exclusion_suffix)]
    else:
        marker_columns = [column for column in lf.collect_schema().names() if column.startswith("Phenotype_(standardized) ")]
    return sorted(marker_columns)


# Define the main function.
def main():

    # Ensure the main lazyframe is ready for usage.
    if not (
        ("LAZYFRAMES" in st.session_state)
        and ("unified_input_file" in st.session_state["LAZYFRAMES"])
        and (os.path.exists(os.path.join(framework_utils.session_dir(), st.session_state["LAZYFRAMES"]["unified_input_file"]["input_params"]["local_filepath"])))
        ):
        st.warning("Please load a unified input file first (at left).")
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
            st.session_state[key] = get_marker_columns(lf, exclusion_suffix=exclusion_suffix)

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
            lf_phenotyped = fnp_main.perform_marker_phenotyping_on_lazyframe(lf, marker_columns)
            st.session_state["LAZYFRAMES"]["marker_phenotyping"] = {
                "lf": lf_phenotyped,
                "input_params": {"input_key": "unified_input_file", "function": fnp_main.perform_marker_phenotyping_on_lazyframe, "inputs": {"marker_columns": marker_columns}},
            }
            st.session_state[ST_KEY_PREFIX + "num_phenotyped_rows"] = lf_phenotyped.select(pl.len()).collect().item()
            st.session_state[ST_KEY_PREFIX + "unique_labels"] = lf_phenotyped.select(pl.col("label").unique().sort()).collect().to_series().to_list()
            st.session_state[ST_KEY_PREFIX + "unique_image_ids"] = lf_phenotyped.select(pl.col("Image ID_(standardized)").unique().sort()).collect().to_series().to_list()

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

    # Temporarily write something that access the lazyframe so we can test the framework.
    with main_columns[1]:
        image_colname = "Image ID_(standardized)"
        with st.container(horizontal=True, vertical_alignment="bottom"):
            st.session_state.setdefault(ST_KEY_PREFIX + "selected_image_to_plot", unique_image_ids[0])
            selected_image_to_plot = st.selectbox("Select image to plot:", options=unique_image_ids, key=ST_KEY_PREFIX + "selected_image_to_plot")
            st.button("Previous", on_click=lambda: st.session_state.update({ST_KEY_PREFIX + "selected_image_to_plot": unique_image_ids[max(0, unique_image_ids.index(st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"]) - 1)]}), disabled=(st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"] == unique_image_ids[0]))
            st.button("Next", on_click=lambda: st.session_state.update({ST_KEY_PREFIX + "selected_image_to_plot": unique_image_ids[min(len(unique_image_ids) - 1, unique_image_ids.index(st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"]) + 1)]}), disabled=(st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"] == unique_image_ids[-1]))
        marker_size = st.slider("Marker size:", min_value=1, max_value=10, value=1)
        st.plotly_chart(fnp_main.plot_image_from_frame(lf_phenotyped, image_colname=image_colname, selected_images=[selected_image_to_plot], marker_size=marker_size))


# Run the main function if this script is executed.
if __name__ == "__main__":
    main()
