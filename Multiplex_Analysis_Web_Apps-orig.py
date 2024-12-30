'''
Top level Streamlit Application for MAWA
'''

# Import relevant libraries
import os
import streamlit as st

# Import relevant libraries
import nidap_dashboard_lib as ndl   # Useful functions for dashboards connected to NIDAP
import app_top_of_page as top
import streamlit_dataframe_editor as sde

def main():
    '''
    Define the single, main function
    '''

    # Set a wide layout
    st.set_page_config(layout="wide")

    input_path = './input'
    if not os.path.exists(input_path):
        os.makedirs(input_path)

    output_path = './output'
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Run streamlit-dataframe-editor library initialization tasks at the top of the page
    st.session_state = sde.initialize_session_state(st.session_state)

    # Run Top of Page (TOP) functions
    st.session_state = top.top_of_page_reqs(st.session_state)

    # Markdown text
    intro_markdown = ndl.read_markdown_file('markdown/MAWA_WelcomePage.md')
    st.markdown(intro_markdown, unsafe_allow_html=True)

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)

# Call the main function
if __name__ == '__main__':
    main()
