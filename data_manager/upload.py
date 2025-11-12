# Import relevant libraries.
import streamlit as st
import polars as pl
import framework.platform_abstraction as pa
import zipfile
import tempfile
import os


def count_files_in_zip(uploaded_zip_file, exclude_dirs=True):
    """
    Return the number of members in the uploaded ZIP.
    exclude_dirs=True counts only regular file entries (skips directory markers).
    """
    # Ensure pointer at start (Streamlit UploadedFile is a BytesIO-like object)
    uploaded_zip_file.seek(0)
    with zipfile.ZipFile(uploaded_zip_file) as zf:
        if exclude_dirs:
            # Python 3.11+: use zi.is_dir(); for older versions fallback to name check.
            return sum(
                1 for zi in zf.infolist()
                if not (getattr(zi, "is_dir", lambda: zi.filename.endswith("/"))())
            )
        else:
            return len(zf.infolist())


# Define a function to upload files.
def upload_files_zip(zip_file, upload_location, compress_files_upon_upload, overwrite_existing_files):

    # Rewind because count_files_in_zip() advanced the pointer.
    zip_file.seek(0)
    with st.spinner("Unzipping and uploading files..."):
        with tempfile.TemporaryDirectory() as tmpdirname:
            with zipfile.ZipFile(zip_file) as zf:
                for member in zf.infolist():
                    # Mitigate zip-slip
                    target_path = os.path.realpath(os.path.join(tmpdirname, member.filename))
                    if not target_path.startswith(os.path.realpath(tmpdirname) + os.sep):
                        st.error(f"Unsafe path in zip: {member.filename}")
                        return
                zf.extractall(path=tmpdirname)

            uploaded_files = []
            for root, _, files in os.walk(tmpdirname):
                for file in files:
                    uploaded_files.append(os.path.join(root, file))

            results = pa.upload_objects_parallel(
                file_paths=uploaded_files,
                db_schema=get_location_settings()[upload_location]["db_schema"],
                bucket_name=get_location_settings()[upload_location]["bucket_name"],
                gzip_if_possible=compress_files_upon_upload,
                overwrite=overwrite_existing_files,
            )
    if results:
        st.success(f"Successfully uploaded {len(uploaded_files)} file(s).")
    else:
        st.error("An error occurred while uploading files.")


# Define a function to upload files.
def upload_files_multiple(uploaded_files, upload_location, compress_files_upon_upload, overwrite_existing_files):
    with st.spinner(f"Uploading files..."):
        results = pa.upload_objects_parallel(
            file_paths=uploaded_files,
            db_schema=get_location_settings()[upload_location]["db_schema"],
            bucket_name=get_location_settings()[upload_location]["bucket_name"],
            gzip_if_possible=compress_files_upon_upload,
            overwrite=overwrite_existing_files,
            )
    if results:
        st.success(f"Successfully uploaded {len(uploaded_files)} file(s).")
    else:
        st.error("An error occurred while uploading files.")


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

    st.write("**⚠️ Note**: If you are having trouble uploading files, please try uploading less data or uploading a zip of the data. If problems persist (likely due to large file sizes), just contact us and we will upload the data for you! (We also have the option to increase the compute power of the Data Manger node.)")

    # Have the user select to where they want to upload data.
    upload_location_options = ["Relatively tidy data (e.g., MAWA input files)", "NIDAP-formatted MAWA archives"]
    upload_location = st.selectbox("Select upload location:", options=upload_location_options)

    # Show the current contents of the selected upload location.
    st.write(f"Current contents of this upload location:")
    objects_list = pa.list_objects_in_bucket(
        db_schema=get_location_settings()[upload_location]["db_schema"],
        bucket_name=get_location_settings()[upload_location]["bucket_name"],
    )
    if objects_list:
        df = pl.DataFrame({"Filename": objects_list})
        st.dataframe(df)
        st.write(f"{len(objects_list)} file(s) found in this location.")
    else:
        st.write("No files found in this location.")

    st.write("**Note:** Please don't upload more than ~6GB *total* of files at a time. And please be judicious about space.")

    upload_columns = st.columns(2)

    with upload_columns[0]:
        st.header("Multi-file upload")

        # Create a file uploader.
        if "uploader_key_multiple" not in st.session_state:
            st.session_state["uploader_key_multiple"] = 0
        uploaded_files = st.file_uploader("Upload multiple files at once:", accept_multiple_files=True, key=f"uploader_multiple-{st.session_state['uploader_key_multiple']}__do_not_persist")

        # Optionally push the uploaded files to the server.
        if uploaded_files:
            st.write("**Files won't actually be saved until you push them to the server.**")
            key = "multiple_compress_checkbox"
            if key not in st.session_state:
                st.session_state[key] = True
            compress_files_upon_upload = st.checkbox("Compress files upon upload (recommended)", key=key)
            key = "multiple_overwrite_checkbox"
            if key not in st.session_state:
                st.session_state[key] = False
            overwrite_existing_files = st.checkbox("Overwrite existing files with the same name", key=key)
            st.button(f"Upload {len(uploaded_files)} file(s) to server", on_click=upload_files_multiple, kwargs={"uploaded_files": uploaded_files, "upload_location": upload_location, "compress_files_upon_upload": compress_files_upon_upload, "overwrite_existing_files": overwrite_existing_files})

    with upload_columns[1]:
        st.header("Single zip upload")

        # Create a zip file uploader.
        if "uploader_key_zip" not in st.session_state:
            st.session_state["uploader_key_zip"] = 0
        zip_file = st.file_uploader("Upload a single zip file containing multiple files:", type=["zip"], accept_multiple_files=False, key=f"uploader_zip-{st.session_state['uploader_key_zip']}__do_not_persist")

        # Optionally push the uploaded zip file to the server.
        if zip_file:
            st.write("**Files won't actually be saved until you push them to the server.**")
            key = "zip_compress_checkbox"
            if key not in st.session_state:
                st.session_state[key] = True
            compress_files_upon_upload = st.checkbox("Compress files upon upload (recommended)", key=key)
            key = "zip_overwrite_checkbox"
            if key not in st.session_state:
                st.session_state[key] = False
            overwrite_existing_files = st.checkbox("Overwrite existing files with the same name", key=key)
            st.button(f"Upload {count_files_in_zip(zip_file)} file(s) to server", on_click=upload_files_zip, kwargs={"zip_file": zip_file, "upload_location": upload_location, "compress_files_upon_upload": compress_files_upon_upload, "overwrite_existing_files": overwrite_existing_files})


# Run the main function if this script is executed.
if __name__ == "__main__":
    main()
