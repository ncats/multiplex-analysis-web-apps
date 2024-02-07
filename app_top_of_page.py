'''
Set of methods that organize the functions and processes
that must be run first at the top of every Streamlit page
'''
import subprocess
import numpy as np
import platform_io
import streamlit as st
from st_pages import show_pages_from_config, add_indentation
from streamlit_extras.app_logo import add_logo

def top_of_page_reqs(session_state):
    '''
    Top of the page requirements. These are the commands
    that are run at the top of the page that help maintain 
    '''

    # # Remove key values from session_state that should not persist
    # for key, val in session_state.items():
    #     if (not key.endswith('__do_not_persist')) and (not key.startswith('FormSubmitter:')):
    #         session_state[key] = val

    # Add Logo
    add_logo('app_images/mawa_logo-width315.png', height=150)

    # Apply pages order and indentation
    add_indentation()
    show_pages_from_config()

    # Sidebar organization
    with st.sidebar:
        st.write('**:book: [Documentation](https://ncats.github.io/multiplex-analysis-web-apps/)**')

    # Check the platform
    session_state = check_for_platform(session_state)

    return session_state

def platform_is_nidap():
    '''
    Check if the Streamlit application is operating on NIDAP
    '''
    return np.any(['nidap.nih.gov' in x for x in subprocess.run('conda config --show channels', shell=True, capture_output=True).stdout.decode().split('\n')[1:-1]])

def check_for_platform(session_state):
    '''
    Set the platform parameters based on the platform the Streamlit app is running on
    '''
    # Initialize the platform object
    if 'platform' not in session_state:
        session_state['platform'] = platform_io.Platform(platform=('nidap' if platform_is_nidap() else 'local'))
    return session_state
