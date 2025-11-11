# Import relevant libraries.
import streamlit as st
import polars as pl
import framework.platform_abstraction as pa
import tempfile
import os
import io
import zipfile


# Store information about possible upload locations.
@st.cache_data()
def get_location_settings():
    return {
        "Relatively tidy data (e.g., MAWA input files)": {
            "bucket_name": pa.DATA_OBJECTS_BUCKET_NAME,
            "db_schema": f"{pa.get_user_group(pa.get_current_username())}_group_db.curated_schema",
        },
        "NIDAP-formatted MAWA archives": {
            "bucket_name": pa.OLD_ARCHIVES_BUCKET_NAME,
            "db_schema": f"{pa.get_user_group(pa.get_current_username())}_group_db.mawa_schema",
        },
    }


# Define the main function.
def main():

    # Write some information.
    st.write(f"Running as user: **{pa.get_current_username()}** in group: **{pa.get_user_group(pa.get_current_username())}**")

    # Have the user select to where they want to upload data.
    upload_location_options = ["Relatively tidy data (e.g., MAWA input files)", "NIDAP-formatted MAWA archives"]
    upload_location = st.selectbox("Select upload location:", options=upload_location_options)

    # Show the current contents of the selected upload location.
    st.write(f"Current contents of this upload location:")
    objects_list = pa.list_objects_in_bucket(
        db_schema=get_location_settings()[upload_location]["db_schema"],
        bucket_name=get_location_settings()[upload_location]["bucket_name"],
    )
    key = "current_contents_table__do_not_persist"
    if objects_list:
        df = pl.DataFrame({"Filename": objects_list})
        st.dataframe(df, on_select="rerun", key=key)
        st.write(f"{len(objects_list)} file(s) found in this location.")
    else:
        st.write("No files found in this location.")

    # If some files are selected...
    if key in st.session_state:
        rows = st.session_state[key]["selection"]["rows"]
        if rows:
            selected_filenames = df[rows]["Filename"].to_list()

            # Allow the user to download them from the server to the machine where the app is running.
            gunzip_if_gz = st.checkbox("Decompress .gz files upon download", value=True)
            if st.button(f"Download {len(selected_filenames)} file(s) to app"):
                with st.spinner("Downloading files..."):
                    with tempfile.TemporaryDirectory() as tempdir:
                        pa.download_objects_parallel(
                            object_names=selected_filenames,
                            db_schema=get_location_settings()[upload_location]["db_schema"],
                            bucket_name=get_location_settings()[upload_location]["bucket_name"],
                            dest_dir=tempdir,
                            gunzip_if_gz=gunzip_if_gz,
                        )

                        # Normalize whatever filenames actually exist after potential gunzip.
                        downloaded_paths = []
                        for original in selected_filenames:
                            if gunzip_if_gz and original.lower().endswith(".gz"):
                                base_name = os.path.splitext(original)[0]
                            else:
                                base_name = original
                            downloaded_paths.append(os.path.join(tempdir, base_name))

                        # Create a ZIP file in memory from the downloaded files.
                        zip_buf = io.BytesIO()
                        with zipfile.ZipFile(zip_buf, "w", compression=zipfile.ZIP_DEFLATED) as z:
                            for p in downloaded_paths:
                                z.write(p, arcname=os.path.basename(p))
                        zip_buf.seek(0)

                        # Store in session state for download button outside spinner.
                        st.session_state["zip_buffer"] = zip_buf
                        st.session_state["num_files_in_buffer"] = len(selected_filenames)
                        st.session_state["zip_filename"] = f"{get_location_settings()[upload_location]['bucket_name']}_files.zip"

    # If some data have been downloaded, offer download button.
    if "zip_buffer" in st.session_state:
        st.download_button(
            label=f"Download {st.session_state['num_files_in_buffer']} file(s) to computer",
            data=st.session_state["zip_buffer"],
            file_name=st.session_state["zip_filename"],
            mime="application/zip",
        )


# Run the main function if this script is executed.
if __name__ == "__main__":
    main()
