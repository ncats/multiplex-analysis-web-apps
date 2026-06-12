# These should all be one-time startup operations.
# We should also consider whether things get appropriately reset if e.g. the user refreshes the page.

import streamlit as st
import framework.platform_abstraction as pa
import os

APP_TITLE = os.getenv("APP_TITLE")


def initialize():

    # Set page configuration.
    st.set_page_config(
        page_title=APP_TITLE,
        )

    # Set up object storage and database.
    pa.set_up_minio()
    pa.set_up_postgresql()
