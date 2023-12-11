# Import relevant libraries
import streamlit as st
import platform_io
from streamlit_javascript import st_javascript
from st_pages import show_pages_from_config, add_indentation
from streamlit_extras.app_logo import add_logo
import streamlit_dataframe_editor as sde

# Import relevant libraries
import nidap_dashboard_lib as ndl   # Useful functions for dashboards connected to NIDAP

# Function to determine whether we're on NIDAP
def platform_is_nidap():
    import numpy as np
    import subprocess
    return np.any(['nidap.nih.gov' in x for x in subprocess.run('conda config --show channels', shell=True, capture_output=True).stdout.decode().split('\n')[1:-1]])

# Define the single, main function
def main():

    # Set a wide layout
    st.set_page_config(layout="wide")

    # Run streamlit-dataframe-editor library initialization tasks at the top of the page
    st.session_state = sde.initialize_session_state(st.session_state)

    # Apply pages order and indentation
    add_indentation()
    show_pages_from_config()

    # Get link for app to be opened in new tab and add link to docs
    url = st_javascript("await fetch('').then(r => window.parent.location.href)")  # this single line is causing a repeatable "KeyError: 'platform'" that seems to not have any ill effects. This is fixed by putting this line after the "if 'platform' not in st.session_state:" block and "platform = st.session_state['platform']" line
    # Sidebar organization
    with st.sidebar:
        st.write('**:book: [Documentation](https://ncats.github.io/multiplex-analysis-web-apps)**')
        st.write('*[Open app in new tab]({}) (may only work in Chrome or Edge, not Firefox)*'.format(url))

    # Add logo to page
    add_logo('app_images/mawa_logo-width315.png', height=150)

    # Initialize the platform object
    if 'platform' not in st.session_state:
        st.session_state['platform'] = platform_io.Platform(platform=('nidap' if platform_is_nidap() else 'local'))
    
    intro_markdown = ndl.read_markdown_file('markdown/MAWA_WelcomePage.md')
    st.markdown(intro_markdown, unsafe_allow_html=True)

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)

# Call the main function
if __name__ == '__main__':
    main()
