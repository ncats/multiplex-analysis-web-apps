# Import relevant libraries
import streamlit as st
import app_top_of_page as top
import streamlit_dataframe_editor as sde


def main():
    """
    Main function for the page.
    """

    st.write('Insert your code here.')


# Run the main function
if __name__ == '__main__':

    # Set page settings
    page_name = 'Your Page Name Here'
    st.set_page_config(layout='wide', page_title=page_name)
    st.title(page_name)

    # Run streamlit-dataframe-editor library initialization tasks at the top of the page
    st.session_state = sde.initialize_session_state(st.session_state)

    # Run Top of Page (TOP) functions
    st.session_state = top.top_of_page_reqs(st.session_state)

    # Call the main function
    main()

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)
