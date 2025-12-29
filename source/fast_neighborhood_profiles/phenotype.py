# Import relevant libraries.
import streamlit as st
from fast_neighborhood_profiles import main as fnp_main
import streamlit_dataframe_editor as sde
import polars as pl
import framework.utils as framework_utils

# Define constants.
ST_KEY_PREFIX = "phenotype.py__"
ST_KEY_PREFIX_LOAD = "load_unified_input_file.py__"


# GUI interface for marker phenotyping.
def marker_phenotyping(lf, marker_columns_with_prefix):

    # Allow the user to perform phenotyping.
    if st.button("Perform marker phenotyping"):

        # Generate and save generation metadata for the phenotyped lazyframe.
        params = {"marker_columns_with_prefix": marker_columns_with_prefix}
        lf_phenotyped = fnp_main.perform_marker_phenotyping_on_lazyframe(lf, **params)
        st.session_state["LAZYFRAMES"]["phenotyped"] = {
            "lf": lf_phenotyped,
            "function_metadata": {"module_name": "fast_neighborhood_profiles.main", "qualpath": "perform_marker_phenotyping_on_lazyframe"},
            "input_dataset": {"type": "lf", "keys": ("unified_input_file",)},
            "params": params,
        }

        # Obtain information related to the phenotyped data.
        metadata = fnp_main.get_phenotyped_metadata(lf_phenotyped)

        # Save some metadata about the phenotypes.
        st.session_state[ST_KEY_PREFIX + "num_phenotyped_rows"] = metadata["num_phenotyped_rows"]
        st.session_state[ST_KEY_PREFIX + "unique_labels"] = metadata["unique_labels"]
        st.session_state[ST_KEY_PREFIX + "unique_image_ids"] = metadata["unique_image_ids"]
        st.session_state[ST_KEY_PREFIX + "phenotype_color_map"] = metadata["phenotype_color_map"]

        # Store a random string in the session state to indicate that phenotyping has been performed.
        st.session_state[ST_KEY_PREFIX + "phenotyping_random_str"] = framework_utils.get_unique_id()


# GUI interface for species phenotyping.
def species_phenotyping(lf, marker_columns_with_prefix):

    # Allow the user to detect the species in the dataset.
    if st.button("Detect species in dataset"):

        # Define marker columns without prefix.
        marker_columns = [col.removeprefix("Phenotype_(standardized) ") for col in marker_columns_with_prefix]

        # Add the species column to the main lazyframe.
        lf = fnp_main.obtain_species_column_from_markers_columns(lf, marker_columns_with_prefix, marker_columns)

        # Create a table that is the base for the assignments table.
        lf_species_counts = fnp_main.create_species_assignments_table(lf)

        # Create a new dataframe editor for the species assignments.
        if ST_KEY_PREFIX + "de_species_assignments" in st.session_state:
            del st.session_state[ST_KEY_PREFIX + "de_species_assignments"]
        st.session_state[ST_KEY_PREFIX + "de_species_assignments"] = sde.DataframeEditor(df_name=ST_KEY_PREFIX + "df_species_assignments", default_df_contents=lf_species_counts.collect(engine="streaming").to_pandas())

    # Ensure the species assignments data editor is in session state before we can display it.
    if ST_KEY_PREFIX + "de_species_assignments" not in st.session_state:
        st.info("Please press the button above to detect species in the dataset.")
        return

    # Allow the user to modify species assignments.
    st.write("Edit species names (and/or deleted unwanted species):")
    st.session_state[ST_KEY_PREFIX + "de_species_assignments"].dataframe_editor(reset_data_editor_button_text='Reset to default names', disabled=["Species", "Count in dataset"])

    # Allow the user to perform species phenotyping.
    if st.button("Perform species phenotyping"):

        # Generate and save generation metadata for the phenotyped lazyframe.
        df_species_assignments = st.session_state[ST_KEY_PREFIX + "de_species_assignments"].reconstruct_edited_dataframe()
        params = dict(marker_columns_with_prefix=marker_columns_with_prefix, df_species_assignments=df_species_assignments)
        lf_phenotyped = fnp_main.perform_species_phenotyping_on_lazyframe(lf, **params)
        st.session_state["LAZYFRAMES"]["phenotyped"] = {
            "lf": lf_phenotyped,
            "function_metadata": {"module_name": "fast_neighborhood_profiles.main", "qualpath": "perform_species_phenotyping_on_lazyframe"},
            "input_dataset": {"type": "lf", "keys": ("unified_input_file",)},
            "params": params,
        }

        # Obtain information related to the phenotyped data.
        metadata = fnp_main.get_phenotyped_metadata(lf_phenotyped)

        # Save some metadata about the phenotypes.
        st.session_state[ST_KEY_PREFIX + "num_phenotyped_rows"] = metadata["num_phenotyped_rows"]
        st.session_state[ST_KEY_PREFIX + "unique_labels"] = metadata["unique_labels"]
        st.session_state[ST_KEY_PREFIX + "unique_image_ids"] = metadata["unique_image_ids"]
        st.session_state[ST_KEY_PREFIX + "phenotype_color_map"] = metadata["phenotype_color_map"]

        # Store a random string in the session state to indicate that phenotyping has been performed.
        st.session_state[ST_KEY_PREFIX + "phenotyping_random_str"] = framework_utils.get_unique_id()


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
        st.header("Options")

        # Button to extract the marker columns.
        if st.button("Get marker column options"):
            marker_column_options, marker_column_options_with_prefix = fnp_main.get_marker_columns(lf, prefix="Phenotype_(standardized) ")
            st.session_state[ST_KEY_PREFIX + "marker_column_options"] = marker_column_options
            st.session_state[ST_KEY_PREFIX + "marker_column_options_with_prefix"] = marker_column_options_with_prefix
            if ST_KEY_PREFIX + "marker_columns" in st.session_state:
                del st.session_state[ST_KEY_PREFIX + "marker_columns"]

        # Ensure the marker columns are in session state.
        if ST_KEY_PREFIX + "marker_column_options" not in st.session_state:
            st.info("Please press the button above to obtain the marker column options.")
            return
        
        # Get a shortcut to the marker column options.
        marker_column_options = st.session_state[ST_KEY_PREFIX + "marker_column_options"]
        marker_column_options_with_prefix = st.session_state[ST_KEY_PREFIX + "marker_column_options_with_prefix"]

        # Allow the user to select the marker columns they want to use.
        st.session_state.setdefault(ST_KEY_PREFIX + "marker_columns", marker_column_options)
        marker_columns = st.multiselect("Select marker columns to use for phenotyping:", options=marker_column_options, key=ST_KEY_PREFIX + "marker_columns")
        marker_columns = [x for x in marker_column_options if x in marker_columns]  # Maintain decreasing frequency order.
        marker_columns_with_prefix = [f"Phenotype_(standardized) {x}" for x in marker_columns]

        # Create tabs for the two different phenotyping methods.
        marker_tab, species_tab = st.tabs(["Marker phenotyping", "Species phenotyping"])

        # For marker phenotyping...
        with marker_tab:
            marker_phenotyping(lf, marker_columns_with_prefix)

        # For species phenotyping...
        with species_tab:
            species_phenotyping(lf, marker_columns_with_prefix)

        # Ensure the phenotyped lazyframe is in session state.
        if "phenotyped" not in st.session_state["LAZYFRAMES"]:
            st.info("Please perform phenotyping above.")
            return
        
        # Display the number of rows in the phenotyped lazyframe.
        lf_phenotyped = st.session_state["LAZYFRAMES"]["phenotyped"]["lf"]
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
        st.header("Result")
        image_colname = "Image ID_(standardized)"

        more_options_columns = st.columns(2)
        with more_options_columns[0]:

            # Make grid lines optional.
            st.session_state.setdefault(ST_KEY_PREFIX + 'show_grid_lines', True)
            show_grid_lines = st.checkbox("Show grid lines", key=ST_KEY_PREFIX + 'show_grid_lines')

            st.session_state.setdefault(ST_KEY_PREFIX + 'use_coordinate_mins_and_maxs', False)
            use_coordinate_mins_and_maxs = st.checkbox("Use coordinate mins and maxs for faithful plotting of object sizes, if possible", key=ST_KEY_PREFIX + 'use_coordinate_mins_and_maxs')

            # Allow the user to set the marker size if not using coordinate mins and maxs, which don't use markers.
            # Note if we plot faithful object sizes, the reported units on the plot are likely wrong (they should be pixels probably).
            if not use_coordinate_mins_and_maxs:
                st.session_state.setdefault(ST_KEY_PREFIX + 'marker_size', 5)
                marker_size = st.slider('Marker size:', min_value=1, max_value=50, step=1, key=ST_KEY_PREFIX + 'marker_size')
            else:
                marker_size = None

            # Allow the user to set export options.
            st.session_state.setdefault(ST_KEY_PREFIX + 'print_width_in', 7.0)
            st.session_state.setdefault(ST_KEY_PREFIX + 'print_height_in', 4.5)
            st.session_state.setdefault(ST_KEY_PREFIX + 'target_dpi', 600)
            st.session_state.setdefault(ST_KEY_PREFIX + 'figure_name', "figure")
            figure_name = st.text_input('Figure name:', key=ST_KEY_PREFIX + 'figure_name')

        with more_options_columns[1]:
            print_width_in = st.number_input('Print width (inches):', min_value=1.0, max_value=20.0, step=0.1, key=ST_KEY_PREFIX + 'print_width_in')
            print_height_in = st.number_input('Print height (inches):', min_value=1.0, max_value=20.0, step=0.1, key=ST_KEY_PREFIX + 'print_height_in')
            target_dpi = st.number_input('Target DPI:', min_value=72, max_value=1200, step=1, key=ST_KEY_PREFIX + 'target_dpi')
            final_figure_name = f"{figure_name}_{print_width_in}in_x_{print_height_in}in_{target_dpi}dpi"

        with st.container(horizontal=True, vertical_alignment="bottom"):
            if ST_KEY_PREFIX + "selected_image_to_plot" in st.session_state and st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"] not in unique_image_ids:
                del st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"]
            st.session_state.setdefault(ST_KEY_PREFIX + "selected_image_to_plot", unique_image_ids[0])
            selected_image_to_plot = st.selectbox("Select image to plot:", options=unique_image_ids, key=ST_KEY_PREFIX + "selected_image_to_plot")
            st.button("Previous", on_click=lambda: st.session_state.update({ST_KEY_PREFIX + "selected_image_to_plot": unique_image_ids[max(0, unique_image_ids.index(st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"]) - 1)]}), disabled=(st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"] == unique_image_ids[0]))
            st.button("Next", on_click=lambda: st.session_state.update({ST_KEY_PREFIX + "selected_image_to_plot": unique_image_ids[min(len(unique_image_ids) - 1, unique_image_ids.index(st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"]) + 1)]}), disabled=(st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"] == unique_image_ids[-1]))

        fig = fnp_main.plot_image_from_frame(lf_phenotyped, image_colname=image_colname, selected_images=[selected_image_to_plot], marker_size=marker_size, color_map=st.session_state[ST_KEY_PREFIX + "phenotype_color_map"], custom_columns=["input_index"], plot_faithful_object_sizes=use_coordinate_mins_and_maxs, frame_with_faithful_columns=lf_phenotyped)

        # Optionally show the grid lines.
        fig.update_xaxes(showgrid=show_grid_lines)
        fig.update_yaxes(showgrid=show_grid_lines)

        # Derived pixel size for export
        export_width_px = int(print_width_in * target_dpi)   # e.g., 7.0 * 600 = 4200 px
        export_height_px = int(print_height_in * target_dpi) # e.g., 4.5 * 600 = 2700 px

        # ---- Configure the built-in "Download as PNG" button ----
        config = {
            "toImageButtonOptions": {
                "format": "png",
                "filename": final_figure_name,
                "height": export_height_px,
                "width": export_width_px,
                "scale": 1,  # keep 1 since width/height already encode 600 DPI
            },
        }

        # Plot the plotly chart in Streamlit
        st.plotly_chart(fig, config=config)

        value_counts_columns = st.columns(2)
        with value_counts_columns[0]:
            st.subheader("Full dataset counts")
            @st.cache_data(show_spinner="Computing full dataset counts...", show_time=True)
            def full_dataset_counts(phenotyping_random_str):
                framework_utils.multiprint(f"Computing full dataset counts for phenotyping random str: {phenotyping_random_str}", (print,))
                return lf_phenotyped.group_by("label").agg(pl.count().alias("Count in dataset")).sort("Count in dataset", descending=True).collect(engine="streaming")
            st.write(full_dataset_counts(st.session_state[ST_KEY_PREFIX + "phenotyping_random_str"]))
        with value_counts_columns[1]:
            st.subheader("Selected image counts")
            @st.cache_data(show_spinner="Computing image counts...", show_time=True)
            def image_counts(selected_image_to_plot, phenotyping_random_str):
                framework_utils.multiprint(f"Computing image counts for image {selected_image_to_plot} and phenotyping random str: {phenotyping_random_str}", (print,))
                return lf_phenotyped.filter(pl.col(image_colname) == selected_image_to_plot).group_by("label").agg(pl.count().alias(f"Count in {selected_image_to_plot}")).sort(f"Count in {selected_image_to_plot}", descending=True).collect(engine="streaming")
            st.write(image_counts(selected_image_to_plot, st.session_state[ST_KEY_PREFIX + "phenotyping_random_str"]))


# Run the main function if this script is executed.
if __name__ == "__main__":
    main()
