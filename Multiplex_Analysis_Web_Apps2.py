import streamlit as st
import os
from streamlit_extras.app_logo import add_logo
import streamlit_session_state_management
import nidap_dashboard_lib as ndl   # Useful functions for dashboards connected to NIDAP
import streamlit_utils
import numpy as np
import subprocess
import platform_io


def welcome_page():
    # Markdown text
    intro_markdown = ndl.read_markdown_file('markdown/MAWA_WelcomePage.md')
    st.markdown(intro_markdown, unsafe_allow_html=True)


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


def main():

    # Set a wide layout
    st.set_page_config(layout="wide")

    # Use the new st.naviation()/st.Page() API to create a multi-page app
    pg = st.navigation({
        'Home':
            [st.Page(welcome_page, title="Home", url_path='home')],
        # 'first section':
        #     [st.Page(one.main, title="Home", url_path='home'),
        #      st.Page(two.main, title="Second page", url_path='two')],
        # 'second section':
        #     [st.Page(three.main, title="Third page", url_path='three')]
        })

    # Ensure the input/output directories exist
    input_path = './input'
    if not os.path.exists(input_path):
        os.makedirs(input_path)
    output_path = './output'
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # For widget persistence, we need always copy the session state to itself, being careful with widgets that cannot be persisted, like st.data_editor() (where we use the "__do_not_persist" suffix to avoid persisting it)
    for key in st.session_state.keys():
        if (not key.endswith('__do_not_persist')) and (not key.startswith('FormSubmitter:')):
            st.session_state[key] = st.session_state[key]

    # This is needed for the st.dataframe_editor() class (https://github.com/andrew-weisman/streamlit-dataframe-editor) but is also useful for seeing where we are and where we've been
    st.session_state['current_page_name'] = pg.url_path if pg.url_path != '' else 'Home'
    if 'previous_page_name' not in st.session_state:
        st.session_state['previous_page_name'] = st.session_state['current_page_name']

    # Add logo to sidebar
    add_logo('app_images/mawa_logo-width315.png', height=250)

    # Determine whether this is the first time the app has been run
    if 'app_has_been_run_at_least_once' not in st.session_state:
        st.session_state['app_has_been_run_at_least_once'] = True
        first_app_run = True
    else:
        first_app_run = False

    # Run session state management in the sidebar
    streamlit_session_state_management.execute(first_app_run)

    # Initalize session_state values for streamlit processing
    if 'init' not in st.session_state:
        st.session_state = ndl.init_session_state(st.session_state)

    # Sidebar organization
    with st.sidebar:
        st.write('**:book: [Documentation](https://ncats.github.io/multiplex-analysis-web-apps/)**')
        with st.expander('Advanced:'):
            benchmark_button = True
            if benchmark_button:
                st.button('Record Benchmarking', on_click = st.session_state.bc.save_run_to_csv)
            if st.button('Calculate memory used by Python session'):
                streamlit_utils.write_python_session_memory_usage()

    # Check the platform
    st.session_state = check_for_platform(st.session_state)

    # Format tooltips
    tooltip_style = """
        <style>
        div[data-baseweb="tooltip"] {
        width: 250px;
        }
        </style>
    """
    st.markdown(tooltip_style,unsafe_allow_html=True)

    # On every page, display its title
    st.title(pg.title)

    # Render the select page
    pg.run()

    # Update the previous page location
    st.session_state['previous_page_name'] = st.session_state['current_page_name']


# Needed for rendering pages which use multiprocessing (https://docs.python.org/3/library/multiprocessing.html#the-spawn-and-forkserver-start-methods)
if __name__ == '__main__':
    main()
