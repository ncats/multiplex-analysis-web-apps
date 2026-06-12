# Import relevant libraries.
import streamlit as st
import polars as pl
import framework.platform_abstraction as pa


# Define a function to delete files.
def delete_files(upload_location, selected_filenames):
    with st.spinner("Deleting files..."):
        success = pa.delete_objects(
            bucket_name=get_location_settings()[upload_location]["bucket_name"],
            object_names=selected_filenames,
            db_schema=get_location_settings()[upload_location]["db_schema"],
        )
    if success:
        st.success(f"Successfully deleted {len(selected_filenames)} file(s).")
    else:
        st.error("An error occurred while deleting files.")


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

    # Write a warning that deletions are permanent.
    st.warning("⚠️ **Warning:** Deletions are permanent and cannot be undone. Please proceed with caution and ensure you have backups of any important data before deleting files.")
    
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
    if objects_list and (key in st.session_state):
        rows = st.session_state[key]["selection"]["rows"]
        if rows:
            selected_filenames = df[rows]["Filename"].to_list()

            # Show a button to delete the selected files.
            st.button(f"⚠️ Delete {len(selected_filenames)} file(s) from server", on_click=delete_files, kwargs={"upload_location": upload_location, "selected_filenames": selected_filenames}, help="This is permanent; please ensure you have backups.", type="primary")


# Run the main function if this script is executed.
if __name__ == "__main__":
    main()
