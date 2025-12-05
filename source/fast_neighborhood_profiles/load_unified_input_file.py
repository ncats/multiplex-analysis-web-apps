# Import relevant libraries.
import streamlit as st
import polars as pl
import framework.utils as framework_utils
from fast_neighborhood_profiles import main as fnp_main

# Define constant.
ST_KEY_PREFIX = "load_unified_input_file.py__"


# Get the list of objects in the bucket.
@st.cache_data()
def get_objects_list(upload_location):
    return fnp_main.get_objects_list(upload_location)


# Store information about possible upload locations.
@st.cache_data()
def get_location_settings():
    return fnp_main.get_location_settings()


# Define the main function.
def main():

    # Show the current contents of the selected upload location using a selectable dataframe.
    with st.columns(2)[0]:
        upload_location = "Available input files"
        objects_list = get_objects_list(upload_location)
        unified_datafile_mapping = {fullname.removeprefix("mawa-unified_datafile-").removesuffix(".csv.zip").removesuffix(".csv.gz"): fullname for fullname in objects_list if fullname.startswith("mawa-unified_datafile-") and fullname.endswith((".csv.zip", ".csv.gz"))}
        objects_list = unified_datafile_mapping.keys()
        column_heading = "Unified input file"
        key = "current_contents_table__do_not_persist"
        if objects_list:
            df = pl.DataFrame({column_heading: objects_list})
            st.dataframe(df, on_select="rerun", key=key, selection_mode="single-row")
            st.write(f"{len(objects_list)} unified input file(s) found.")
        else:
            st.write("No unified input files found.")
        st.button("Refresh file list", on_click=get_objects_list.clear)

        # If some files are selected...
        if key in st.session_state:
            rows = st.session_state[key]["selection"]["rows"]
            if rows:
                
                # Get a list of the selected shortnames (short versions of the filenames).
                selected_filenames = df[rows][column_heading].to_list()

                # Allow the user to select the intermediate file format.
                available_file_formats = ["parquet (recommended)", "arrow", "csv"]
                intermediate_file_format = st.selectbox("Select intermediate file format:", options=available_file_formats, index=available_file_formats.index("parquet (recommended)"))

                # Load the lazyframe from the selected row.
                if st.button(f"Load unified input file"):
                    object_filename = unified_datafile_mapping[selected_filenames[0]]
                    file_format = "parquet" if intermediate_file_format == "parquet (recommended)" else intermediate_file_format
                    db_schema = get_location_settings()[upload_location]["db_schema"]
                    bucket_name = get_location_settings()[upload_location]["bucket_name"]
                    params = dict(file_format=file_format, db_schema=db_schema, bucket_name=bucket_name, object_filename=object_filename)
                    with st.spinner("Loading file..."):
                        lf = fnp_main.load_unified_input_file_data(**params, topdir=framework_utils.session_dir())
                    st.session_state["LAZYFRAMES"] = {}  # Clear existing lazyframes.
                    st.session_state["LAZYFRAMES"]["unified_input_file"] = {
                        "lf": lf,
                        "function_metadata": {"module_name": "fast_neighborhood_profiles.main", "qualpath": "load_unified_input_file_data"},
                        "input_dataset": None,
                        "params": params,
                    }

    # If there's lazyframe information in the session state...
    if not ("LAZYFRAMES" in st.session_state and "unified_input_file" in st.session_state["LAZYFRAMES"]):
        st.info("Please load a unified input file (whether above or from a session archive) to see its details here.")
        return

    # Get information about the lazyframe from the metadata in the session state.
    file_format = st.session_state["LAZYFRAMES"]["unified_input_file"]["params"]["file_format"]
    db_schema = st.session_state["LAZYFRAMES"]["unified_input_file"]["params"]["db_schema"]
    bucket_name = st.session_state["LAZYFRAMES"]["unified_input_file"]["params"]["bucket_name"]
    object_filename = st.session_state["LAZYFRAMES"]["unified_input_file"]["params"]["object_filename"]

    # Get the lazyframe from the session state now.
    lf = st.session_state["LAZYFRAMES"]["unified_input_file"]["lf"]

    # At this point the lazyframe must be working, so display information about it.
    information = f'''
    Properties:

    :small_orange_diamond: File format: `{file_format}`  
    :small_orange_diamond: database.schema: `{db_schema}`  
    :small_orange_diamond: Bucket name: `{bucket_name}`  
    :small_orange_diamond: Object filename: `{object_filename}`  
    :small_orange_diamond: Number of rows: `{lf.select(pl.len()).collect().item():_}`  
    :small_orange_diamond: Number of columns: `{len(lf.collect_schema())}`  
    '''
    st.markdown(information)

    # Show a sample of 100 rows from the lazyframe.
    st.write(lf.collect().sample(100).sort(pl.col("Image ID_(standardized)")))
    st.button("Resample dataset")


# Run the main function if this script is executed.
if __name__ == "__main__":
    main()
