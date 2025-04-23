'''
Set of methods that organize the functions and processes
that must be run first at the top of every Streamlit page
'''
import subprocess
import numpy as np
import platform_io
import streamlit as st
# from st_pages import show_pages_from_config, add_indentation
from streamlit_extras.app_logo import add_logo
import streamlit_session_state_management
import streamlit_utils

# Import relevant libraries
import nidap_dashboard_lib as ndl   # Useful functions for dashboards connected to NIDAP

def top_of_page_reqs(session_state):
    '''
    Top of the page requirements. These are the commands
    that are run at the top of the page that help maintain 
    '''

    # Add Logo
    add_logo('app_images/mawa_logo-width315.png', height=250)

    # Determine whether this is the first time the app has been run
    if 'app_has_been_run_at_least_once' not in session_state:
        session_state['app_has_been_run_at_least_once'] = True
        first_app_run = True
    else:
        first_app_run = False

    # Apply pages order and indentation
    if first_app_run:
        show_pages_from_config()  # this is slow so only do it once
    add_indentation()

    # Run session state management in the sidebar
    streamlit_session_state_management.execute(first_app_run)

    # Initalize session_state values for streamlit processing
    if 'init' not in session_state:
        session_state = ndl.init_session_state(session_state)

    # Sidebar organization
    with st.sidebar:
        st.write('**:book: [Documentation](https://ncats.github.io/multiplex-analysis-web-apps/)**')
        with st.expander('Advanced:'):
            benchmark_button = True
            if benchmark_button:
                st.button('Record Benchmarking', on_click = session_state.bc.save_run_to_csv)
            if st.button('Calculate memory used by Python session'):
                streamlit_utils.write_python_session_memory_usage()

    # Check the platform
    session_state = check_for_platform(session_state)

    tooltip_style = """
        <style>
        div[data-baseweb="tooltip"] {
        width: 250px;
        }
        </style>
    """
    st.markdown(tooltip_style,unsafe_allow_html=True)

    return session_state

def platform_is_nidap():
    '''
    Check if the Streamlit application is operating on NIDAP
    '''
    return np.any(['foundry-artifacts' in x for x in subprocess.run('conda config --show channels', shell=True, capture_output=True).stdout.decode().split('\n')[1:-1]])

def check_for_platform(session_state):
    '''
    Set the platform parameters based on the platform the Streamlit app is running on
    '''
    # Initialize the platform object
    if 'platform' not in session_state:
        session_state['platform'] = platform_io.Platform(platform=('nidap' if platform_is_nidap() else 'local'))
    return session_state
