import os
# from snowflake.core import Root  # Snowflake SDK: pip install snowflake
from snowflake.snowpark import Session  # Snowpark: pip install snowflake-snowpark-python
# import snowflake.connector  # Python connector: pip install snowflake-connector-python
import streamlit as st
import atexit


def _read_spcs_token() -> str:
    with open("/snowflake/session/token", "r", encoding="utf-8") as f:
        return f.read()


def _snowflake_conn_params() -> dict:
    return {
        "account": os.getenv("SNOWFLAKE_ACCOUNT"),
        "host": os.getenv("SNOWFLAKE_HOST"),            # required with the token
        "authenticator": "oauth",
        "token": _read_spcs_token(),  # newer Python connector versions also support token_file_path
        "warehouse": os.getenv("SNOWFLAKE_WAREHOUSE"),  # you set this; SPCS doesn't inject it
        "database": os.getenv("SNOWFLAKE_DATABASE"),
        "schema": os.getenv("SNOWFLAKE_SCHEMA"),
    }


@st.cache_resource()
def _create_snowpark_session() -> Session:
    session = Session.builder.configs(_snowflake_conn_params()).create()
    atexit.register(lambda: session.close())
    return session


def get_snowpark_session() -> Session:
    session = _create_snowpark_session()
    try:
        session.sql("select 1").collect()
    except Exception as e:
        print(f"Snowpark session interrupted. Creating a new one now. Error: {e}.")
        _create_snowpark_session.clear()
        session = _create_snowpark_session()
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
