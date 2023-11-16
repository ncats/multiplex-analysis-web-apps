# Import relevant libraries
import streamlit as st
import platform_io
from streamlit_javascript import st_javascript
from st_pages import show_pages_from_config, add_indentation
from streamlit_extras.app_logo import add_logo
import streamlit_utils

# Function to determine whether we're on NIDAP
def platform_is_nidap():
    import numpy as np
    import subprocess
    return np.any(['nidap.nih.gov' in x for x in subprocess.run('conda config --show channels', shell=True, capture_output=True).stdout.decode().split('\n')[1:-1]])

# Define the single, main function
def main():

    # Set a wide layout
    st.set_page_config(layout="wide")

    # Restore previous session state values, including from other pages; see https://discuss.streamlit.io/t/simultaneous-multipage-widget-state-persistence-data-editors-with-identical-contents-and-multiprocessing-capability/52554 for more information
    for key, val in st.session_state.items():
        if (not key.endswith('__do_not_persist')) and (not key.startswith('FormSubmitter:')):
            st.session_state[key] = val

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

    # Display page heading
    st.title('Data import/export')

    # Initialize the platform object
    if 'platform' not in st.session_state:
        st.session_state['platform'] = platform_io.Platform(platform=('nidap' if platform_is_nidap() else 'local'))

    # Store a copy (not a link) of the platform object for clarity below
    platform = st.session_state['platform']

    # Define widget default values
    streamlit_utils.assign_default_values_in_session_state('basename_suffix_for_new_results_archive', 'sample_project_name_panel_12')
    streamlit_utils.assign_default_values_in_session_state('basename_suffix_for_new_local_archive_dir', 'sample_project_name_panel_12')

    # One of two main rows, one for Input and one for Results (i.e., output)
    st.header('Input')

    # Split into three non-evenly spaced columns
    cols = st.columns(spec=[0.4, 0.2, 0.4])

    # In the first column...
    with cols[0]:
        platform.display_available_inputs_df()
        platform.add_refresh_available_inputs_button()

    # In the second column...
    with cols[1]:
        platform.load_selected_inputs()

    # In the third column...
    with cols[2]:
        platform.display_local_inputs_df()
        platform.add_delete_local_inputs_button()

    # Second of two main rows
    st.header('Results')

    # Split into three non-evenly spaced columns
    cols = st.columns(spec=[0.4, 0.2, 0.4])

    # In the first column...
    with cols[0]:
        platform.display_archives_df()
        platform.add_delete_archives_button()
        platform.add_refresh_archives_button()

    # In the second column...
    with cols[1]:
        platform.load_selected_archive()
        st.divider()
        platform.save_results_to_archive()

    # In the third column...
    with cols[2]:
        platform.display_local_results_df()
        platform.add_delete_local_results_button()

    # Use first third of page for additional functionality
    st.divider()
    cols = st.columns(3)
    with cols[0]:
        with st.expander('Additional functionality:'):
            platform.write_settings_to_local_results()
            platform.write_environment_to_local_results()
            platform.add_create_empty_archive_button()
            platform.add_delete_empty_archives_button()

    # Since the platform object above was a local copy of that in Streamlit, and it may have been modified, save it back to Streamlit
    st.session_state['platform'] = platform

# Call the main function
if __name__ == '__main__':
    main()
