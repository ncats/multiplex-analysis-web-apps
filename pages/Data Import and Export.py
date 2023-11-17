'''
This is the python script which produces the PHENOTYPING PAGE
'''
import streamlit as st
from st_pages import show_pages_from_config, add_indentation
from streamlit_extras.app_logo import add_logo

# Import relevant libraries
import streamlit_utils

def main():
    '''
    Main function for running the page
    '''

    # Use the whole page width
    st.set_page_config(page_title="Data Import and Export",
                       layout="wide")

    # Remove key values from session_state that should not persist
    for key, val in st.session_state.items():
        if (not key.endswith('__do_not_persist')) and (not key.startswith('FormSubmitter:')):
            st.session_state[key] = val

    # Apply pages order and indentation
    add_indentation()
    show_pages_from_config()

    # Sidebar organization
    with st.sidebar:
        st.write('**:book: [Documentation](https://ncats.github.io/multiplex-analysis-web-apps)**')

    # Add logo to page
    add_logo('app_images/mawa_logo-width315.png', height=150)

    st.title('Data Import and Export')

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

if __name__ == '__main__':
    main()
