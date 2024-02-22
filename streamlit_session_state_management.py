# Import relevant libraries
import streamlit as st
import pickle
import os
from datetime import datetime
import os
import app_top_of_page as top

def save_session_state(saved_streamlit_session_states_dir, saved_streamlit_session_state_prefix, saved_streamlit_session_state_key):
    """
    Save the session state to a pickle file.

    The session state is saved to a pickle file in the `saved_streamlit_session_states_dir`
    directory. The session state is saved as a dictionary and each key-value pair is saved
    individually to the pickle file. The pickle filename is generated with the save date
    and time.

    Args:
        saved_streamlit_session_states_dir (str): The directory in which to save the session state pickle files
        saved_streamlit_session_state_prefix (str): The prefix to add to the pickle filename
        saved_streamlit_session_state_key (str): The key to exclude from the session state when saving

    Returns:
        None
    """
    
    # Create the output directory for saving session state if it doesn't exist
    os.makedirs(saved_streamlit_session_states_dir, exist_ok=True)

    # Generate the pickle filename with the save date and time
    now = datetime.now()
    filename = os.path.join(saved_streamlit_session_states_dir, saved_streamlit_session_state_prefix + now.strftime('date%Y-%m-%d_time%H-%M-%S') + '.pkl')

    # Save each key-value pair in the session_state to the pickle file
    session_dict = {key: value for key, value in st.session_state.items() if key != saved_streamlit_session_state_key}
    with open(filename, 'wb') as f:
        pickle.dump(session_dict, f)
    
    # Output a success message
    st.success('State saved to ' + filename)

def load_session_state(saved_streamlit_session_states_dir, saved_streamlit_session_state_prefix, saved_streamlit_session_state_key, selected_session=None):
    """
    Load the session state from a pickle file.

    The session state is loaded from a pickle file in the `saved_streamlit_session_states_dir`
    directory. The session state is loaded as a dictionary and each key-value pair is loaded
    individually into the session state.

    Args:
        saved_streamlit_session_states_dir (str): The directory from which to load the session state pickle files
        saved_streamlit_session_state_prefix (str): The prefix to add to the pickle filename
        saved_streamlit_session_state_key (str): The key to exclude from the session state when loading
        selected_session (str, optional): The selected session to load. If not provided, the selectbox key will be used.
    
    Returns:
        None
    """

    # Get the selected session basename to load
    if selected_session is None:
        selected_session = st.session_state[saved_streamlit_session_state_key]

    # If no session file was explicitly input and if no session files exist, do nothing; otherwise, load one of these selected sessions (if not a manually input one [nominally the most recent], then the one selected in the session state file selection dropdown)
    if selected_session is not None:

        # Generate the filename based on the selected session
        filename = os.path.join(saved_streamlit_session_states_dir, saved_streamlit_session_state_prefix + selected_session + '.pkl')

        # Load the state (as a dictionary) from the pickle file
        with open(filename, 'rb') as f:
            session_dict = pickle.load(f)

        # Delete every key in the current session state except for the selected session
        for key in st.session_state.keys():
            if key != saved_streamlit_session_state_key:
                del st.session_state[key]

        # Load each key-value pair individually into session_state
        for key, value in session_dict.items():
            st.session_state[key] = value

        # Output a success message
        st.success('State loaded from ' + filename)

    # If no session state files exist, output a warning
    else:
        st.warning('No session state files exist so none were loaded')

def reset_session_state(saved_streamlit_session_state_key):
    """
    Reset the session state.

    This function deletes every key in the session state except for the selected session.

    Args:
        saved_streamlit_session_state_key (str): The key to exclude from the session state when resetting

    Returns:
        None
    """

    # Delete every key in the session state except for the selected session
    for key in st.session_state.keys():
        if key != saved_streamlit_session_state_key:
            del st.session_state[key]

    # Output a success message
    st.success('Session state reset')

def app_session_management(saved_streamlit_session_states_dir, saved_streamlit_session_state_prefix, saved_streamlit_session_state_key):
    """
    Create the sidebar for app session management.

    This function creates the sidebar for app session management, including the 'Save',
    'Load', and 'Reset' buttons, as well as the selectbox to select the session to load.

    Args:
        saved_streamlit_session_states_dir (str): The directory in which to save the session state pickle files
        saved_streamlit_session_state_prefix (str): The prefix to add to the pickle filename
        saved_streamlit_session_state_key (str): The key to exclude from the session state when saving, loading, and resetting

    Returns:
        str: The name of the most recent session state file
    """

    # Create the output directory for saving session state if it doesn't exist
    os.makedirs(saved_streamlit_session_states_dir, exist_ok=True)

    # Delete all symbolic links in the saved_streamlit_session_states_dir directory
    for f in os.listdir(saved_streamlit_session_states_dir):
        if os.path.islink(os.path.join(saved_streamlit_session_states_dir, f)):
            os.unlink(os.path.join(saved_streamlit_session_states_dir, f))
    
    # Check if the right type of pickle file exists in the "output" directory, and if so, create a symbolic link to it from the saved_streamlit_session_states_dir directory
    session_state_files_in_output_dir = []
    if os.path.exists('output'):
        for f in os.listdir('output'):
            if f.endswith('.pkl') and f.startswith(saved_streamlit_session_state_prefix):
                symlink_path = os.path.join(saved_streamlit_session_states_dir, f)
                os.symlink(os.path.join('..', 'output', f), symlink_path)
                session_state_files_in_output_dir.append(f)

    # Get the list of pickle files in the saved session state directory (unsorted)
    files = [f for f in os.listdir(saved_streamlit_session_states_dir) if (f.endswith('.pkl') and f.startswith(saved_streamlit_session_state_prefix))]

    # Name the relevant section in the sidebar
    st.sidebar.subheader('App Session Management')

    # Create the sidebar with the 'Save', 'Load', and 'Reset' buttons
    col1, col2, col3 = st.sidebar.columns(3)
    col1.button('Save', on_click=save_session_state, use_container_width=True, kwargs={'saved_streamlit_session_states_dir': saved_streamlit_session_states_dir, 'saved_streamlit_session_state_prefix': saved_streamlit_session_state_prefix, 'saved_streamlit_session_state_key': saved_streamlit_session_state_key})
    col2.button('Load', on_click=load_session_state, use_container_width=True, disabled=not files, kwargs={'saved_streamlit_session_states_dir': saved_streamlit_session_states_dir, 'saved_streamlit_session_state_prefix': saved_streamlit_session_state_prefix, 'saved_streamlit_session_state_key': saved_streamlit_session_state_key})
    col3.button('Reset', on_click=reset_session_state, use_container_width=True, kwargs={'saved_streamlit_session_state_key': saved_streamlit_session_state_key})

    # Create the dropdown to select the checkpoint to load
    session_basenames_in_reverse_order = sorted([f.removeprefix(saved_streamlit_session_state_prefix).removesuffix('.pkl') for f in files], reverse=True)
    st.sidebar.selectbox('Session to load:', session_basenames_in_reverse_order, key=saved_streamlit_session_state_key)

    # Write to screen the session state files in the output directory
    if session_state_files_in_output_dir:
        st.sidebar.write('Session state files in the "output" directory: `{}`'.format([f.removeprefix(saved_streamlit_session_state_prefix).removesuffix('.pkl') for f in session_state_files_in_output_dir]))

    # Return the most recent session state file, or None
    return session_basenames_in_reverse_order[0] if session_basenames_in_reverse_order else None

def execute(first_app_run):
    """
    Execute the session state management functions.

    This function executes the session state management functions, including the creation of the sidebar for app session
    management and the initialization of the session state with an existing session
    if the app is first run or if Streamlit has been restarted.

    Args:
        first_app_run (bool): Whether the app has been run before

    Returns:
        None
    """

    # Set some main variables
    saved_streamlit_session_states_dir = 'saved_streamlit_session_states'
    saved_streamlit_session_state_prefix = 'streamlit_session_state-'
    saved_streamlit_session_state_key = 'session_selection'

    # Create the sidebar content for app session management
    most_recent_session_management_file = app_session_management(saved_streamlit_session_states_dir, saved_streamlit_session_state_prefix, saved_streamlit_session_state_key)

    # Initialize the session state with an existing session if the app is first run or if Streamlit has been restarted
    if first_app_run:
        load_session_state(saved_streamlit_session_states_dir, saved_streamlit_session_state_prefix, saved_streamlit_session_state_key, selected_session=most_recent_session_management_file)

def main():
    """
    Main function for the Streamlit app.

    The following is a sample of how to use this module.

    Args:
        None

    Returns:
        None
    """

    # Run Top of Page (TOP) functions
    st.session_state = top.top_of_page_reqs(st.session_state)

    # Initialize the selected_option if it doesn't exist
    if 'selected_option' not in st.session_state:
        st.session_state['selected_option'] = 'Option 1'

    # Sample widget that's part of the session state
    st.selectbox('Select an option', ['Option 1', 'Option 2', 'Option 3'], key='selected_option')

# Run the main function
if __name__ == '__main__':
    main()
