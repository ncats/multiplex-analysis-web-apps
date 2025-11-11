# Import relevant libraries.
import streamlit as st
import polars as pl
import framework.platform_abstraction as pa


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
    if objects_list:
        df = pl.DataFrame({"Filename": objects_list})
        st.dataframe(df)
        st.write(f"{len(objects_list)} file(s) found in this location.")
    else:
        st.write("No files found in this location.")

    # Create a file uploader.
    st.write("**Note:** You can upload multiple files at once. Please don't upload more than ~6GB *total* of files at a time. And please be judicious about space.")
    uploaded_files = st.file_uploader("Upload files to container:", accept_multiple_files=True)

    # Optionally push the uploaded files to the server.
    if uploaded_files:
        st.write("**Files won't actually be saved until you push them to the server.**")
        compress_files_upon_upload = st.checkbox("Compress files upon upload", value=True)
        overwrite_existing_files = st.checkbox("Overwrite existing files with the same name", value=False)
        if st.button(f"Upload {len(uploaded_files)} file(s) to server"):
            with st.spinner(f"Uploading files..."):
                results = pa.upload_objects_parallel(
                    file_paths=uploaded_files,
                    db_schema=get_location_settings()[upload_location]["db_schema"],
                    bucket_name=get_location_settings()[upload_location]["bucket_name"],
                    gzip_if_possible=compress_files_upon_upload,
                    overwrite=overwrite_existing_files,
                    )
                if results:
                    st.rerun()


# Run the main function if this script is executed.
if __name__ == "__main__":
    main()
