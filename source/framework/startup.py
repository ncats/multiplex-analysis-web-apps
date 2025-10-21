# These should all be one-time startup operations.
# Note that when the sidebar "Reset app" button gets pressed, everything in this script should be considered for resetting, as in manage_sessions.reset_session_state().
# We should also consider whether things get appropriately reset if e.g. the user refreshes the page.
# Also when the user loads a session state, as in manage_sessions.load_session_state().

import streamlit as st
import framework.utils as utils
import framework.platform_abstraction as pa
import os
import platform_io

ST_KEY_PREFIX = "startup.py__"
APP_TITLE = os.getenv("APP_TITLE")


def initialize():

    # Set page configuration.
    st.set_page_config(
        page_title=APP_TITLE,
        layout='wide'
        )

    # Set up database, object storage, and session directory.
    pa.set_up_database()
    pa.set_up_object_storage()

    # Generate a unique session ID.
    app_session_id = utils.get_unique_id()
    st.session_state[ST_KEY_PREFIX + "app_session_id"] = app_session_id

    # Create the session directory directly (avoid calling utils.session_dir() here to prevent circular dependency).
    session_dir = f"/tmp/{utils._app_title_simple()}/app_session_data/{app_session_id}"
    os.makedirs(session_dir, exist_ok=True)

    # Get the current username.
    current_username = pa.get_current_username()

    # Create app session entry.
    pa.log_app_session((app_session_id, current_username, pa.get_user_group(current_username), pa.get_frontend_image_id()))

    # Ensure the input/output directories exist
    os.makedirs(os.path.join(session_dir, "input"), exist_ok=True)
    os.makedirs(os.path.join(session_dir, "output"), exist_ok=True)

    # Initialize the platform object.
    st.session_state['platform'] = platform_io.Platform(platform=os.getenv("APP_PLATFORM"))
