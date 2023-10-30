import streamlit as st
from streamlit_extras.app_logo import add_logo

# Define a function to easily assign default values to items in the session state
def assign_default_values_in_session_state(key, value, args=None):
    if key not in st.session_state:
        if callable(value):
            if args is None:
                st.session_state[key] = value()
            else:
                st.session_state[key] = value(args)
        else:
            st.session_state[key] = value

# Function to go to a previous option using a button
def previous_option(options, selected_option='selected_image_id', button_name='image_button'):
    curr_index = options.index(st.session_state[selected_option])
    if curr_index != 0:
        curr_index = curr_index - 1
    st.session_state[selected_option] = options[curr_index]
    update_button_activities(options, selected_option=selected_option, button_name=button_name)

# Function to go to a next option using a button
def next_option(options, selected_option='selected_image_id', button_name='image_button'):
    curr_index = options.index(st.session_state[selected_option])
    if curr_index != (len(options) - 1):
        curr_index = curr_index + 1
    st.session_state[selected_option] = options[curr_index]
    update_button_activities(options, selected_option=selected_option, button_name=button_name)

# Function to determine whether previous or next buttons should be disabled
def update_button_activities(options, selected_option='selected_image_id', button_name='image_button'):
    curr_index = options.index(st.session_state[selected_option])
    if curr_index == 0:
        st.session_state['disable_previous_{}'.format(button_name)] = True
    else:
        st.session_state['disable_previous_{}'.format(button_name)] = False
    if curr_index == (len(options) - 1):
        st.session_state['disable_next_{}'.format(button_name)] = True
    else:
        st.session_state['disable_next_{}'.format(button_name)] = False

def propagate_settings_to_session_state(settings):
    for key0 in settings:
        for key1 in settings[key0]:
            current_setting = settings[key0][key1]
            combined_key = 'settings__' + key0 + '__' + key1
            st.session_state[combined_key] = current_setting

def get_current_settings():

    # Initialize the dictionary holding the actual settings dictionary to be used in the workflow
    settings = dict()
    settings['input_datafile'], settings['phenotyping'], settings['analysis'], settings['annotation'], settings['plotting'] = dict(), dict(), dict(), dict(), dict()

    # Assign current values of the settings
    for key in st.session_state:
        if key.startswith('settings__'):
            key0 = key.split('__')[1]
            key1 = key.split('__')[2]
            settings[key0][key1] = st.session_state[key]

    # Return the nested dictionary of settings
    return settings

def streamlit_write_function(message):
    import streamlit as st
    # st.write(message)
    st.error(message, icon="🚨")

def update_session_state_keys(transformation_dict):
    for key in transformation_dict:
        st.session_state[key] = transformation_dict[key]
