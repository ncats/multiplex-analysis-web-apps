# Import relevant libraries.
import streamlit as st
from fast_neighborhood_profiles import main as fnp_main
import polars as pl
from functools import partial
import streamlit_dataframe_editor as sde
import pandas as pd
import framework.utils as framework_utils

# Define session state key prefixes.
ST_KEY_PREFIX = "delete_cells.py__"
ST_KEY_PREFIX_PHENOTYPE = "phenotype.py__"


# Activate highlights.
def activate_selection_group():
    selections_table = st.session_state[ST_KEY_PREFIX + "selections_table__do_not_persist"]
    rows = selections_table["selection"]["rows"]
    if rows:
        if len(rows) > 1:
            framework_utils.multiprint("Somehow multiple rows are selected, which is unexpected.", (print, st.warning))
            return
        df = st.session_state[ST_KEY_PREFIX + "de_selections"].reconstruct_edited_dataframe()
        input_indices = df.iloc[rows[0]]["input_indices"]  # input_indices for the selected selection group.
        st.session_state[ST_KEY_PREFIX + "selected_indices"] = input_indices


# Obtain the selected indices from the scatter plot.
def get_selected_indices():
    selection = st.session_state[ST_KEY_PREFIX + f"scatter_plot__do_not_persist"]
    if "selection" in selection and "points" in selection["selection"] and selection["selection"]["points"]:
        points_list = selection["selection"]["points"]
        indices = [point["customdata"][4] for point in points_list]  # Note this means that if the "input_index" column is added to the plot data when calling main.plot_image_from_frame(), it must be the very first custom_column, i.e., at position 4 (0-based indexing) since there are four required columns in front of it.
        st.session_state[ST_KEY_PREFIX + "selected_indices"] = indices
    else:
        st.session_state[ST_KEY_PREFIX + "selected_indices"] = []


# Main function.
def main():

    # Ensure the phenotyped lazyframe is ready for usage.
    if not ("LAZYFRAMES" in st.session_state and "phenotyped" in st.session_state["LAZYFRAMES"]):
        st.warning("Please perform phenotyping (at left).")
        return

    # Get the main lazyframe from session state.
    lf = st.session_state["LAZYFRAMES"]["phenotyped"]["lf"]

    unique_image_ids = st.session_state[ST_KEY_PREFIX_PHENOTYPE + "unique_image_ids"]
    selected_indices = []
    if ST_KEY_PREFIX + "selected_indices" in st.session_state and st.session_state[ST_KEY_PREFIX + "selected_indices"]:
        selected_indices = st.session_state[ST_KEY_PREFIX + "selected_indices"]
    image_colname = "Image ID_(standardized)"
    xcol = "Centroid X (µm)_(standardized)"
    ycol = "Centroid Y (µm)_(standardized)"
    color_col = "label"
    key = ST_KEY_PREFIX + "de_selections"
    if key not in st.session_state:
        st.session_state[key] = sde.DataframeEditor(df_name=ST_KEY_PREFIX + "df_selections", default_df_contents=pd.DataFrame(columns=["label", "number_of_cells", "input_indices", "color"]))

    st.write(lf.head().collect(engine="streaming"))

    with st.container(horizontal=True, vertical_alignment="bottom"):
        selected_image_to_plot = st.selectbox("Select image to plot:", options=unique_image_ids, key=ST_KEY_PREFIX + "selected_image_to_plot")
        st.button("Previous", on_click=lambda: st.session_state.update({ST_KEY_PREFIX + "selected_image_to_plot": unique_image_ids[max(0, unique_image_ids.index(st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"]) - 1)]}), disabled=(st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"] == unique_image_ids[0]))
        st.button("Next", on_click=lambda: st.session_state.update({ST_KEY_PREFIX + "selected_image_to_plot": unique_image_ids[min(len(unique_image_ids) - 1, unique_image_ids.index(st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"]) + 1)]}), disabled=(st.session_state[ST_KEY_PREFIX + "selected_image_to_plot"] == unique_image_ids[-1]))

    # Allow the user to select marker size.
    st.session_state.setdefault(ST_KEY_PREFIX + "marker_size", 5)
    marker_size = st.slider("Marker size:", min_value=2, max_value=10, key=ST_KEY_PREFIX + "marker_size")

    # Write the number of selected points.
    with st.container(horizontal=True):
        st.write(f"Number of selected points: {len(selected_indices):_}")
        st.button("Clear selection", on_click=lambda: st.session_state.update({ST_KEY_PREFIX + "selected_indices": []}), key=ST_KEY_PREFIX + "clear_selection_button__do_not_persist")

    # Plot the real space with selectable points.
    fig = fnp_main.plot_image_from_frame(lf, image_colname=image_colname, xcol=xcol, ycol=ycol, color_col=color_col, selected_images=[selected_image_to_plot], marker_size=marker_size, custom_columns=["input_index"], color_map=st.session_state[ST_KEY_PREFIX_PHENOTYPE + "phenotype_color_map"], highlight_index_col="input_index", highlight_indices=selected_indices, sort_index_col="input_index")  # this highlight_indices=selected_indices might be a bit weird
    fig.update_layout(uirevision="static")  # this doesn't seem to be honored; investigate in the future... actually, maybe it is?
    st.plotly_chart(fig, on_select=get_selected_indices, selection_mode=("points", "box", "lasso"), key=ST_KEY_PREFIX + "scatter_plot__do_not_persist")

    # If there are selected points...
    if selected_indices:

        # Allow user to choose a label for the selected cells.
        key = ST_KEY_PREFIX + "selected_label"
        st.session_state.setdefault(key, "")
        selected_label = st.text_input("Enter label for selected cells (can edit later)", key=key)

        # Allow user to pick color of selected cells.
        key = ST_KEY_PREFIX + "selected_color"
        st.session_state.setdefault(key, "#FF0000")  # FF0000 is red
        selected_color = st.color_picker("Select color for selected cells (can edit later)", key=key)

        # Allow user to add the selected cells to a selections dataframe.
        if st.button("Add selected cells to selections table"):
            df = st.session_state[ST_KEY_PREFIX + "de_selections"].reconstruct_edited_dataframe()
            new_row = {
                "label": selected_label if selected_label else f"Selection {len(df) + 1}",
                "number_of_cells": len(selected_indices),
                "input_indices": selected_indices,
                "color": selected_color,
            }
            df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)
            st.session_state[ST_KEY_PREFIX + "de_selections"].update_editor_contents(new_df_contents=df)

    # Plot the editable and selectable tables side-by-side.
    st.write("Select a row in the neighborhood types table below to highlight above.")
    selections_table_columns = st.columns(2)
    with selections_table_columns[0]:
        st.session_state[ST_KEY_PREFIX + "de_selections"].dataframe_editor(reset_data_editor_button_text='Reset selections', disabled=["number_of_cells", "input_indices"])
    with selections_table_columns[1]:
        st.dataframe(st.session_state[ST_KEY_PREFIX + "de_selections"].reconstruct_edited_dataframe(), on_select=activate_selection_group, key=ST_KEY_PREFIX + "selections_table__do_not_persist", selection_mode="single-row")



    # # Allow user to register the selected neighborhood types.
    # if st.button("Register selected neighborhood types"):
    #     df = st.session_state[ST_KEY_PREFIX + "de_selections"].reconstruct_edited_dataframe()
    #     params = dict(updates_pd=df, keep=keep_strategy, missing_label_value=missing_label_value)
    #     lf_neighborhoods = fnp_main.add_new_label_column(lf=lf, **params)
    #     st.session_state["LAZYFRAMES"]["neighborhood_types"] = {
    #         "lf": lf_neighborhoods,
    #         "function_metadata": {"module_name": "fast_neighborhood_profiles.main", "qualpath": "add_new_label_column"},
    #         "input_dataset": {"type": "lf", "keys": ("sumap_cells",)},
    #         "params": params,
    #     }
    #     color_map = dict(zip(df["label"], df["color"]))
    #     color_map[missing_label_value] = "#808080"
    #     st.session_state[ST_KEY_PREFIX + "neighborhood_type_color_map"] = color_map
    #     st.session_state[ST_KEY_PREFIX + "unique_neighborhood_types"] = list(set(df["label"].to_list() + [missing_label_value]))
    #     st.session_state[ST_KEY_PREFIX + "df_reconstructed_selections"] = df


# Run the main function if this script is executed.
if __name__ == "__main__":
    main()
