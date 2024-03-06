# Import relevant libraries
import streamlit as st
import pickle
import os
from datetime import datetime
import os
import app_top_of_page as top
import streamlit_dataframe_editor as sde

def save_session_state(saved_streamlit_session_states_dir, saved_streamlit_session_state_prefix='streamlit_session_state-', saved_streamlit_session_state_key='session_selection'):
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

    # Create a dictionary of most items in the session state
    session_dict = {}
    keys_to_exclude = []
    for key, value in st.session_state.items():
        if (not key.endswith('__do_not_persist')) and (not key.startswith('FormSubmitter:')) and (key != saved_streamlit_session_state_key):
            if isinstance(value, sde.DataframeEditor):
                print(f'Saving components for dataframe editor {key}')
                dataframe_editor_components = {
                    'df_name': value.df_name,
                    'default_df_contents': value.default_df_contents,
                    'edited_dataframe': value.reconstruct_edited_dataframe()
                }
                dataframe_editor_components_name = 'dataframe_editor_components__' + key
                session_dict[dataframe_editor_components_name] = dataframe_editor_components
                keys_to_exclude.append(value.df_name)
                keys_to_exclude.append(value.df_name + '_changes_dict')
                keys_to_exclude.append(value.df_name + '_key')
            else:
                print(f'Saving {key} of type {type(value)}')
                session_dict[key] = value

    # For each key to exclude, delete it from session_dict, i.e., the thing being saved to disk. This should allow a new key to be assigned to the st.data_editor() object within any dataframe editor objects that are later initialized in the load_session_state() function. Otherwise, the key will always be the same so upon loading a session state, the st.data_editor() object will not be redrawn. The solution to that is to force the st.data_editor() object to be redrawn by forcing its key to change.
    for key in keys_to_exclude:
        if key in session_dict:
            print(f'Not actually saving {key} of type {type(session_dict[key])}')
            del session_dict[key]

    # Save the dictionary to the pickle file. Note this no longer randomly crashes with "PicklingError: Can't pickle <class 'streamlit_dataframe_editor.DataframeEditor'>: it's not the same object as streamlit_dataframe_editor.DataframeEditor" because we're no longer saving the DataframeEditor object itself, but rather the initialization data and the current contents. Note that using dill did also solve the problem, which if we were to use dill, we could try saving the entire session at once (instead of individual objects) and also thereby include difficult items such as st.form objects. NOTE: If we start getting the error again, try either using dill or probably better yet, excluding other custom object types from being saved in the first place, e.g., class 'platform_io.Platform'. Such exclusion would be done in keys_to_exclude as above.
    with open(filename, 'wb') as f:
        pickle.dump(session_dict, f)
    
    # Output a success message
    st.success('State saved to ' + filename)

def load_session_state(saved_streamlit_session_states_dir, saved_streamlit_session_state_prefix='streamlit_session_state-', saved_streamlit_session_state_key='session_selection', selected_session=None):
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

        # Delete every key in the current session state except for the selected session
        for key in st.session_state.keys():
            if key != saved_streamlit_session_state_key:
                del st.session_state[key]

        # Load the state (as a dictionary) from the pickle file
        with open(filename, 'rb') as f:
            session_dict = pickle.load(f)

        # Load each key-value pair individually into session_state
        for key, value in session_dict.items():
            if key.startswith('dataframe_editor_components__'):
                print(f'Initializing dataframe editor {key.removeprefix("dataframe_editor_components__")}')
                dataframe_editor_key = key.removeprefix('dataframe_editor_components__')
                dataframe_editor_components = value
                # curr_dataframe_editor = sde.DataframeEditor(df_name=dataframe_editor_components['df_name'], default_df_contents=dataframe_editor_components['default_df_contents'])
                # curr_dataframe_editor.update_editor_contents(new_df_contents=dataframe_editor_components['edited_dataframe'], reset_key=True)
                # st.session_state[dataframe_editor_key] = curr_dataframe_editor
                st.session_state[dataframe_editor_key] = sde.DataframeEditor(df_name=dataframe_editor_components['df_name'], default_df_contents=dataframe_editor_components['default_df_contents'])
                st.session_state[dataframe_editor_key].update_editor_contents(new_df_contents=dataframe_editor_components['edited_dataframe'], reset_key=True)  # no real point of setting reset_key=True here since that's the default, but keeping it as a reminder that it's something we can modify
            else:
                print(f'Loading {key} of type {type(value)}')
                st.session_state[key] = value

        # Output a success message
        st.success('State loaded from ' + filename)

    # If no session state files exist, output a warning
    else:
        st.warning('No session state files exist so none were loaded')

def reset_session_state(saved_streamlit_session_state_key='session_selection', success_message=True):
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

    # Flag that the reset button has been hit
    st.session_state['reset_button_hit'] = True

    # Output a success message
    if success_message:
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
        load_most_recent_session = False
        if 'reset_button_hit' in st.session_state:
            if not st.session_state['reset_button_hit']:
                load_most_recent_session = True
            else:
                st.session_state['reset_button_hit'] = False
        else:
            load_most_recent_session = True
        if load_most_recent_session:
            load_session_state(saved_streamlit_session_states_dir, saved_streamlit_session_state_prefix, saved_streamlit_session_state_key, selected_session=most_recent_session_management_file)

def main():
    """
    Main function for the Streamlit app.

    The following is a sample of how to use this module.
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
