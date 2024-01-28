import numpy as np
import platform_io
import subprocess
import streamlit as st
from st_pages import show_pages_from_config, add_indentation
from streamlit_extras.app_logo import add_logo

def top_of_page_reqs(session_state, newtab_flag = False):
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
        st.write('**:book: [Documentation](https://nidap.nih.gov/workspace/module/view/latest/ri.workshop.main.module.6b54d691-4d02-4f77-8acd-734f42e3eee0)**')
        if newtab_flag:
            # Get link for app to be opened in new tab and add link to docs
            # This line is causing a repeatable "KeyError: 'platform'" that seems 
            # not to have any ill effects. This is fixed by putting this line after the 
            # "if 'platform' not in st.session_state:" block and "platform = #st.session_state['platform']" line
            # url = st_javascript("await fetch('').then(r => window.parent.location.href)")  
            # st.write('*[Open app in new tab]({}) (may only work in Chrome or Edge, not Firefox)*'.format(url))
            pass


    # Check the platform
    session_state = check_for_platform(session_state)

    return session_state

def platform_is_nidap():
    return np.any(['nidap.nih.gov' in x for x in subprocess.run('conda config --show channels', shell=True, capture_output=True).stdout.decode().split('\n')[1:-1]])

def check_for_platform(session_state):
    # Initialize the platform object
    if 'platform' not in session_state:
        session_state['platform'] = platform_io.Platform(platform=('nidap' if platform_is_nidap() else 'local'))
    return session_state