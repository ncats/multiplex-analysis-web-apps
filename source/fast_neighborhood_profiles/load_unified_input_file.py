# Import relevant libraries.
import streamlit as st
import polars as pl
import framework.platform_abstraction as pa
import framework.utils as framework_utils
import os
import zipfile
from fast_neighborhood_profiles import main as fnp_main

# Define constant.
ST_KEY_PREFIX = "load_unified_input_file.py__"


# Load the lazyframe from the specified file using an intermediate file.
def load_unified_input_file_data(file_format, db_schema, bucket_name, object_filename):
    with st.spinner("Loading file..."):

        # Shortcuts, the first for generalizability.
        full_filenames = [object_filename]
        input_dir = os.path.join(framework_utils.session_dir(), "input")

        # Download the file from the server.
        pa.download_objects_parallel(
            object_names=full_filenames,
            db_schema=db_schema,
            bucket_name=bucket_name,
            dest_dir=input_dir,
        )

        # Unzip the downloaded file.
        unzipped_paths = []
        for full_filename in full_filenames:
            with zipfile.ZipFile(os.path.join(input_dir, full_filename), 'r') as zip_ref:
                zip_ref.extractall(input_dir)
            base_name = full_filename.removesuffix(".zip")
            os.remove(os.path.join(input_dir, full_filename))
            unzipped_paths.append(os.path.join(input_dir, base_name))

        # Generate an intermediate file from which to load the lazyframe.
        for unzipped_path in unzipped_paths:
            filename = os.path.basename(unzipped_path)
            filepath = fnp_main.subset_csv_to_file(csv_filename=filename, handle="unified_input_file", topdir=framework_utils.session_dir(), subdir="input", file_format=file_format)
            os.remove(unzipped_path)

        # Load the lazyframe.
        lf = fnp_main.get_lf("unified_input_file", topdir=framework_utils.session_dir(), subdir="input", file_format=file_format)

        # Save the lazyframe and its metadata to the session state.
        return {"lf": lf, "input_params": {"file_format": file_format, "db_schema": db_schema, "bucket_name": bucket_name, "object_filename": full_filenames[0], "local_filepath": filepath.removeprefix(framework_utils.session_dir() + os.sep)}}


# Get the list of objects in the bucket.
@st.cache_data()
def get_objects_list(upload_location):
    objects_list = pa.list_objects_in_bucket(
        db_schema=get_location_settings()[upload_location]["db_schema"],
        bucket_name=get_location_settings()[upload_location]["bucket_name"],
    )
    return objects_list


# Store information about possible upload locations.
@st.cache_data()
def get_location_settings():
    return {
        "Available input files": {
            "bucket_name": pa.DATA_OBJECTS_BUCKET_NAME,
            "db_schema": f"{pa.get_user_group(pa.get_current_username())}_group_db.curated_schema",
        },
    }


# Define the main function.
def main():

    # Show the current contents of the selected upload location using a selectable dataframe.
    with st.columns(2)[0]:
        upload_location = "Available input files"
        objects_list = get_objects_list(upload_location)
        unified_datafile_mapping = {fullname.removeprefix("mawa-unified_datafile-").removesuffix(".csv.zip"): fullname for fullname in objects_list if fullname.startswith("mawa-unified_datafile-") and fullname.endswith(".csv.zip")}
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
                    st.session_state["LAZYFRAMES"] = {}  # Clear existing lazyframes.
                    st.session_state["LAZYFRAMES"]["unified_input_file"] = load_unified_input_file_data(file_format, db_schema, bucket_name, object_filename)

    # If there's lazyframe information in the session state...
    if "LAZYFRAMES" in st.session_state and "unified_input_file" in st.session_state["LAZYFRAMES"]:

        # Get information about the lazyframe from the metadata in the session state.
        file_format = st.session_state["LAZYFRAMES"]["unified_input_file"]["input_params"]["file_format"]
        db_schema = st.session_state["LAZYFRAMES"]["unified_input_file"]["input_params"]["db_schema"]
        bucket_name = st.session_state["LAZYFRAMES"]["unified_input_file"]["input_params"]["bucket_name"]
        object_filename = st.session_state["LAZYFRAMES"]["unified_input_file"]["input_params"]["object_filename"]
        filepath = os.path.join(framework_utils.session_dir(), st.session_state["LAZYFRAMES"]["unified_input_file"]["input_params"]["local_filepath"])

        # If the intermediate file doesn't actually exist, we know we have to load the intermediate file and define a lazyframe to point to it.
        if not os.path.exists(filepath):
            # Load the unified input file into a lazyframe.
            st.session_state["LAZYFRAMES"]["unified_input_file"] = load_unified_input_file_data(file_format=file_format, db_schema=db_schema, bucket_name=bucket_name, object_filename=object_filename)

            # Load any other lazyframes.
            other_lazyframes = [lf_name for lf_name in st.session_state["LAZYFRAMES"] if lf_name != "unified_input_file"]
            for lf_name in other_lazyframes:
                input_params = st.session_state["LAZYFRAMES"][lf_name]["input_params"]
                st.session_state["LAZYFRAMES"][lf_name]["lf"] = input_params["function"](st.session_state["LAZYFRAMES"][input_params["input_key"]]["lf"], **input_params["inputs"])

        # Get the lazyframe from the session state now that it's certainly up-to-date using either loading method (choosing a row or reading metadata from the session state).
        lf = st.session_state["LAZYFRAMES"]["unified_input_file"]["lf"]

        # At this point the lazyframe must be working, so display information about it.
        information = f'''
        Properties:

        :small_orange_diamond: File format: `{file_format}`  
        :small_orange_diamond: database.schema: `{db_schema}`  
        :small_orange_diamond: Bucket name: `{bucket_name}`  
        :small_orange_diamond: Object filename: `{object_filename}`  
        :small_orange_diamond: Filepath: `{filepath}`  
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
