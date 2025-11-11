# Import relevant libraries.
import streamlit as st
import framework.startup as startup
import framework.platform_abstraction as pa
import data_uploader
import data_downloader

ST_KEY_PREFIX = "app.py__"


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
                    st.Page(data_uploader.main, title="☁️ Data Uploader", url_path='data_uploader'),
                    st.Page(data_downloader.main, title="📥 Data Downloader", url_path='data_downloader'),
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
            if st.button("🛑 Shut down app"):
                pa.shut_down_app()

    # Display the title of the page.
    st.title(pg.title)

    # Display the page.
    pg.run()

    # Update the previous page location.
    st.session_state['previous_page_name'] = st.session_state['current_page_name']


# Run the main function.
if __name__ == "__main__":
    main()
