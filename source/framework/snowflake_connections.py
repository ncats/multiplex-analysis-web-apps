import os
# from snowflake.core import Root  # Snowflake SDK: pip install snowflake
from snowflake.snowpark import Session  # Snowpark: pip install snowflake-snowpark-python
# import snowflake.connector  # Python connector: pip install snowflake-connector-python
import streamlit as st
import atexit
import tomllib
import pathlib


def _load_local_settings(settings_if_on_local: dict, file_path: pathlib.Path | None = None) -> dict:
    """
    Load a Snowflake connection profile from ~/.snowflake/connections.toml.

    :param file_path: Optional explicit path to the TOML file.
    :return: Dict containing keys like 'account', 'user', etc.
    :raises FileNotFoundError, KeyError, ValueError
    """
    # Default to ~/.snowflake/connections.toml
    file_path = file_path or pathlib.Path("~/.snowflake/connections.toml").expanduser()

    if not file_path.is_file():
        raise FileNotFoundError(f"Connections file not found: {file_path}")

    with file_path.open("rb") as f:
        data = tomllib.load(f)

    local_profile = settings_if_on_local["local_profile"]

    if local_profile not in data:
        raise KeyError(f"Profile [{local_profile}] not found in {file_path}.")

    settings = data[local_profile]

    # Minimal validation (customize as needed)
    required = ["account", "user", "password"]
    missing = [k for k in required if k not in settings or not str(settings[k]).strip()]
    if missing:
        raise ValueError(f"Missing required keys in [{local_profile}]: {', '.join(missing)}")
    
    del settings_if_on_local["local_profile"]
    settings.update(settings_if_on_local)  # Merge any additional settings

    return settings


def _read_spcs_token() -> str:
    with open("/snowflake/session/token", "r", encoding="utf-8") as f:
        return f.read()


def _snowflake_conn_params(settings_if_on_local: dict = {}) -> dict:
    # Minimal auto-detect: if the SPCS token file exists, use OAuth token mode.
    if os.path.exists("/snowflake/session/token"):
        return {
            "account": os.getenv("SNOWFLAKE_ACCOUNT"),
            "host": os.getenv("SNOWFLAKE_HOST"),            # required with the token
            "authenticator": "oauth",
            "token": _read_spcs_token(),  # newer Python connector versions also support token_file_path
            "warehouse": os.getenv("SNOWFLAKE_WAREHOUSE"),  # you set this; SPCS doesn't inject it
            "database": os.getenv("SNOWFLAKE_DATABASE"),
            "schema": os.getenv("SNOWFLAKE_SCHEMA"),
        }
    else:
        return _load_local_settings(settings_if_on_local)  # you need account, user, and password at minimum


@st.cache_resource()
def _create_snowpark_session(settings_if_on_local: dict = {}) -> Session:
    session = Session.builder.configs(_snowflake_conn_params(settings_if_on_local)).create()
    atexit.register(lambda: session.close())
    return session


def get_snowpark_session(settings_if_on_local: dict = {}) -> Session:
    session = _create_snowpark_session(settings_if_on_local)
    try:
        session.sql("select 1").collect()
    except Exception as e:
        print(f"Snowpark session interrupted. Creating a new one now. Error: {e}.")
        _create_snowpark_session.clear()
        session = _create_snowpark_session(settings_if_on_local)
    return session


# @st.cache_resource()
# def _create_connector_connection() -> snowflake.connector.SnowflakeConnection:
#     conn = snowflake.connector.connect(**_snowflake_conn_params())
#     atexit.register(lambda: conn.close())
#     return conn


# def get_connector_connection() -> snowflake.connector.SnowflakeConnection:
#     conn = _create_connector_connection()
#     try:
#         with conn.cursor() as cur:
#             cur.execute("select 1").fetchall()
#     except Exception as e:
#         print(f"Connector connection interrupted. Creating a new one now. Error: {e}.")
#         _create_connector_connection.clear()
#         conn = _create_connector_connection()
#     return conn


# def get_root(conn=None) -> Root:
#     if conn is None:
#         conn = get_connector_connection()
#     return Root(conn)
