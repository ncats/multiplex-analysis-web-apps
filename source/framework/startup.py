# These should all be one-time startup operations.

import streamlit as st
import yaml
import framework.utils as utils
import framework.platform_abstraction as pa
import nidap_dashboard_lib as ndl   # Useful functions for dashboards connected to NIDAP

ST_KEY_PREFIX = "startup.py__"
SETTINGS_FILENAME = "settings.yaml"


def initialize():

    # Store the app settings in the session state if not already present.
    # I should probably get rid of this eventually in favor of environment variables.
    key = ST_KEY_PREFIX + "app_settings"
    with open(SETTINGS_FILENAME, 'r') as f:
        st.session_state[key] = yaml.safe_load(f)

    # Set page configuration.
    st.set_page_config(
        page_title=st.session_state[key]['general']['app_title'],
        layout='wide'
        )

    # Set up database, object storage, and session directory.
    pa.set_up_database()
    pa.set_up_object_storage()

    # Generate a unique session ID.
    app_session_id = utils.get_unique_id()
    st.session_state[ST_KEY_PREFIX + "app_session_id"] = app_session_id

    # Create the session directory.
    utils.session_dir()  # This creates the session directory if it doesn't already exist.

    # Get the current username.
    current_username = pa.get_current_username()

    # Create app session entry.
    pa.log_app_session((app_session_id, current_username, pa.get_user_group(current_username), pa.get_frontend_image_id()))

    # Dante's session state initialization.
    ndl.init_session_state(st.session_state)
