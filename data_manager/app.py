# Import relevant libraries.
import streamlit as st
import framework.startup as startup
import framework.platform_abstraction as pa
import upload
import download
import delete
import home
import os

ST_KEY_PREFIX = "app.py__"


# Bump the key so Streamlit creates a brand-new widget with no value.
def clear_data_from_memory():
    if "uploader_key_multiple" in st.session_state:
        st.session_state["uploader_key_multiple"] += 1
    if "uploader_key_zip" in st.session_state:
        st.session_state["uploader_key_zip"] += 1
    if "zip_buffer" in st.session_state:
        del st.session_state["zip_buffer"]
        del st.session_state["num_files_in_buffer"]
        del st.session_state["zip_filename"]


@st.cache_data()
def get_app_title():
    return os.getenv("APP_TITLE")


# Define the main function.
def main():

    # Run one-time initialization.
    key = ST_KEY_PREFIX + "app_initialized"
    if key not in st.session_state:
        startup.initialize()
        st.session_state[key] = True

    # Define the pages for the navigation bar.
    pg = st.navigation(
        {
            'Menu':
                [
                    st.Page(home.main, title="🏠 Home", url_path='home'),
                    st.Page(upload.main, title="☁️ Upload", url_path='upload'),
                    st.Page(download.main, title="📥 Download", url_path='download'),
                    st.Page(delete.main, title="🔥 Delete", url_path='delete'),
                ],
        }
    )

    # For widget persistence between pages, we need always copy the session state to itself.
    for key in st.session_state:
        if not key.endswith('__do_not_persist'):  # Could add things like: "(not key.endswith('_button'))".
            st.session_state[key] = st.session_state[key]

    # This is needed for the st.dataframe_editor() class (https://github.com/andrew-weisman/streamlit-dataframe-editor) but is also useful for seeing where we are and where we've been.
    st.session_state['current_page_name'] = pg.url_path if pg.url_path != '' else 'Home'
    if 'previous_page_name' not in st.session_state:
        st.session_state['previous_page_name'] = st.session_state['current_page_name']

    # Sidebar organization
    with st.sidebar:

        # Allow user to shut down entire app cleanly.
        with st.container(horizontal=True):
            st.button("🔄 Refresh page", help="If you want to refresh the page, press this button, *not* your browser's refresh button.")
            st.button("🧹 Clear data from memory", on_click=clear_data_from_memory)
            if st.button("🛑 Shut down app"):
                pa.shut_down_app()


    # Display the title of the page.
    st.title(get_app_title() + " - " + pg.title)

    # Display the page.
    pg.run()

    # Update the previous page location.
    st.session_state['previous_page_name'] = st.session_state['current_page_name']


# Run the main function.
if __name__ == "__main__":
    main()
