# Import relevant libraries.
import streamlit as st
import polars as pl
import framework.platform_abstraction as pa
import framework.utils as framework_utils
import os
import zipfile
from fast_neighborhood_profiles import main as fnp_main

ST_KEY_PREFIX = "load_unified_input_file.py__"


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

    # Show the current contents of the selected upload location.
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
            selected_filenames = df[rows][column_heading].to_list()

            available_file_formats = ["parquet (recommended)", "arrow", "csv"]
            intermediate_file_format = st.selectbox("Select intermediate file format:", options=available_file_formats, index=available_file_formats.index("parquet (recommended)"))

            # Allow the user to download them from the server to the machine where the app is running.
            if st.button(f"Load unified input file"):
                with st.spinner("Loading file..."):
                    full_filenames = [unified_datafile_mapping[short_filename] for short_filename in selected_filenames]
                    input_dir = os.path.join(framework_utils.session_dir(), "input")

                    pa.download_objects_parallel(
                        object_names=full_filenames,
                        db_schema=get_location_settings()[upload_location]["db_schema"],
                        bucket_name=get_location_settings()[upload_location]["bucket_name"],
                        dest_dir=input_dir,
                    )

                    unzipped_paths = []
                    for full_filename in full_filenames:
                        with zipfile.ZipFile(os.path.join(input_dir, full_filename), 'r') as zip_ref:
                            zip_ref.extractall(input_dir)
                        base_name = full_filename.removesuffix(".zip")
                        os.remove(os.path.join(input_dir, full_filename))
                        unzipped_paths.append(os.path.join(input_dir, base_name))

                    file_format = "parquet" if intermediate_file_format == "parquet (recommended)" else intermediate_file_format
                    for unzipped_path in unzipped_paths:
                        filename = os.path.basename(unzipped_path)
                        filepath = fnp_main.subset_csv_to_file(csv_filename=filename, handle="unified_input_file", topdir=framework_utils.session_dir(), subdir="input", file_format=file_format)
                        os.remove(unzipped_path)

                    lf = fnp_main.get_lf("unified_input_file", topdir=framework_utils.session_dir(), subdir="input", file_format=file_format)
                    st.session_state[ST_KEY_PREFIX + "unified_input_file"] = {"lf": lf, "file_format": file_format, "db_schema": get_location_settings()["Available input files"]["db_schema"], "bucket_name": get_location_settings()["Available input files"]["bucket_name"], "object_filename": full_filenames[0], "local_filepath": filepath}


# Run the main function if this script is executed.
if __name__ == "__main__":
    main()
