import streamlit as st
import os
import io
import psycopg2.pool  # Probably change to version 3 (psycopg) soon!
import minio
import getpass
import polars as pl
import requests
import atexit
import framework.utils as framework_utils
import framework.analysis_framework as analysis_framework
import framework.snowflake_connections as snowflake_connections
import framework.snowflake_orchestrator as snowflake_orchestrator
from concurrent.futures import ThreadPoolExecutor, as_completed
import time
import hashlib
import mimetypes


ST_KEY_PREFIX_STARTUP = "startup.py__"


#### 1. DATABASE FUNCTIONALITY ####################################################################


DB_URL_GROUP = os.getenv('DB_URL_GROUP')
DB_URL_COMMON = os.getenv('DB_URL_COMMON')
APP_NAME = os.getenv('APP_NAME')


@st.cache_resource()
def get_connection_pool(db_url):
    if framework_utils.platform() == "local":
        # Note we could use atexit to gracefully close the db connection pool. Note that nothing is needed for minio as shutdown is already clean.
        try:
            pool = psycopg2.pool.ThreadedConnectionPool(
                minconn=1,
                maxconn=10,
                dsn=db_url
            )
            atexit.register(lambda: pool.closeall())
            return pool
        except Exception as e:
            st.error(f"Failed to create database pool: {e}")
            return None
    elif framework_utils.platform() == "snowflake":
        pass  # Like everything left as "pass", this is not needed on Snowflake.


# Note we could set this and the following up using @contextmanager as in 8/28/25 chat with GH Copilot, but keeping out for now for simplicity.
def get_database_connection(db_url):
    if framework_utils.platform() == "local":
        try:
            db_pool = get_connection_pool(db_url)
            if db_pool:
                return db_pool.getconn()
        except Exception as e:
            st.error(f"Failed to get database connection: {e}")
            return None
    elif framework_utils.platform() == "snowflake":
        pass


def return_database_connection(conn, db_url):
    if framework_utils.platform() == "local":
        try:
            db_pool = get_connection_pool(db_url)
            if db_pool and conn:
                db_pool.putconn(conn)
            return True
        except Exception as e:
            st.error(f"Failed to return database connection: {e}")
            return False
    elif framework_utils.platform() == "snowflake":
        pass

# Helper to avoid repeating rollback logic in postgresql.
def _rollback_and_return(conn, db_url):
    if conn:
        try:
            conn.rollback()
        except Exception:
            pass
        return_database_connection(conn, db_url)


# Create four tables for the app. Note it should be largely consistent with what's in 01_set_up_non_user_objects.sql for now, and later on we should probably have this function, if even still necessary, just run setup.sql so we don't have to maintain this logic in two places.
@st.cache_data()
def set_up_postgresql():
    if framework_utils.platform() == "local":
        try:
            conn_common = get_database_connection(DB_URL_COMMON)
            with conn_common.cursor() as cur:
                # Create schema if it doesn't exist
                cur.execute(f"""
                    CREATE SCHEMA IF NOT EXISTS admin_schema
                """)
                
                # Check if user_groups_table table exists and is empty
                cur.execute(f"""
                    SELECT EXISTS (
                        SELECT FROM information_schema.tables
                        WHERE table_schema = 'admin_schema'
                        AND table_name = 'user_groups_table'
                    )
                """)
                table_exists = cur.fetchone()[0]
                # Create user_groups_table table.
                cur.execute(f"""
                    CREATE TABLE IF NOT EXISTS admin_schema.user_groups_table (
                        id SERIAL PRIMARY KEY,
                        username VARCHAR(255) UNIQUE NOT NULL,
                        user_group VARCHAR(255),
                        user_added_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                        who_added VARCHAR(255),
                        user_email VARCHAR(255)
                    )
                """)
                # Only insert initial users if table was just created (didn't exist before)
                if not table_exists:
                    initial_users = [
                        ('andrew', 'dmap', 'andrew', 'andrew@example.com'),
                        ('jessica', 'ABC Lab', 'andrew', 'jessica@example.com'),
                        ('Tessa', 'dmap', 'andrew', 'tessa@example.com')
                    ]
                    cur.executemany(f"""
                        INSERT INTO admin_schema.user_groups_table (username, user_group, who_added, user_email)
                        VALUES (%s, %s, %s, %s)
                    """, initial_users)
            conn_common.commit()
            return_database_connection(conn_common, DB_URL_COMMON)
                
            conn_group = get_database_connection(DB_URL_GROUP)
            with conn_group.cursor() as cur:
                # Create schema if it doesn't exist
                cur.execute(f"""
                    CREATE SCHEMA IF NOT EXISTS {APP_NAME}_schema
                """)
                
                # Create app_sessions_table table.
                cur.execute(f"""
                    CREATE TABLE IF NOT EXISTS {APP_NAME}_schema.app_sessions_table (
                        id SERIAL PRIMARY KEY,
                        app_session_id VARCHAR(255) UNIQUE NOT NULL,
                        username VARCHAR(255),
                        user_group VARCHAR(255),
                        startup_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                        explicit_shutdown_time TIMESTAMP,
                        container_image_id VARCHAR(255)
                    )
                """)

                # Create archives_table table.
                cur.execute(f"""
                    CREATE TABLE IF NOT EXISTS {APP_NAME}_schema.archives_table (
                        id SERIAL PRIMARY KEY,
                        creator VARCHAR(255),
                        user_group VARCHAR(255),
                        archive_description TEXT,
                        current_git_commit VARCHAR(255),
                        container_image_id VARCHAR(255),
                        archive_id VARCHAR(255) UNIQUE NOT NULL,
                        app_session_id VARCHAR(255),
                        creation_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                    )
                """)

                # Create jobs_table table.
                cur.execute(f"""
                    CREATE TABLE IF NOT EXISTS {APP_NAME}_schema.jobs_table (
                        id SERIAL PRIMARY KEY,
                        job_id VARCHAR(255) UNIQUE NOT NULL,
                        job_name VARCHAR(255),
                        job_status VARCHAR(255),
                        submitter VARCHAR(255),
                        submitter_group VARCHAR(255),
                        app_session_id VARCHAR(255),
                        worker_image_id VARCHAR(255),
                        submission_time TIMESTAMP,
                        start_time TIMESTAMP,
                        completion_time TIMESTAMP,
                        failure_time TIMESTAMP
                    )
                """)
            conn_group.commit()
            return_database_connection(conn_group, DB_URL_GROUP)
            return True
        except Exception as e:
            st.error(f"Failed to set up databases: {e}")
            _rollback_and_return(conn_common, DB_URL_COMMON)
            _rollback_and_return(conn_group, DB_URL_GROUP)
            return False
    elif framework_utils.platform() == "snowflake":
        pass


def write_archive_database_data(row_tuple):
    if framework_utils.platform() == "local":
        try:
            conn = get_database_connection(DB_URL_GROUP)
            with conn.cursor() as cur:
                cur.execute(f"""
                    INSERT INTO {APP_NAME}_schema.archives_table (creator, user_group, archive_description, current_git_commit, container_image_id, archive_id, app_session_id)
                    VALUES (%s, %s, %s, %s, %s, %s, %s)
                """, row_tuple)
            conn.commit()
            return_database_connection(conn, DB_URL_GROUP)
            return True
        except Exception as e:
            st.error(f"Failed to write archive database data: {e}")
            _rollback_and_return(conn, DB_URL_GROUP)
            return False
    elif framework_utils.platform() == "snowflake":
        try:
            session = snowflake_connections.get_snowpark_session()
            session.sql(f"""
                INSERT INTO {get_user_group(get_current_username())}_group_db.{APP_NAME}_schema.archives_table (creator, user_group, archive_description, current_git_commit, container_image_id, archive_id, app_session_id)
                VALUES (?, ?, ?, ?, ?, ?, ?)
            """, row_tuple).collect()
            return True
        except Exception as e:
            st.error(f"Failed to write archive database data: {e}")
            return False


def log_app_session(row_tuple):
    if framework_utils.platform() == "local":
        try:
            conn = get_database_connection(DB_URL_GROUP)
            with conn.cursor() as cur:
                cur.execute(f"""
                    INSERT INTO {APP_NAME}_schema.app_sessions_table (app_session_id, username, user_group, container_image_id)
                    VALUES (%s, %s, %s, %s)
                """, row_tuple)
            conn.commit()
            return_database_connection(conn, DB_URL_GROUP)
            return True
        except Exception as e:
            st.error(f"Failed to log app session: {e}")
            _rollback_and_return(conn, DB_URL_GROUP)
            return False
    elif framework_utils.platform() == "snowflake":
        try:
            session = snowflake_connections.get_snowpark_session()
            session.sql(f"""
                INSERT INTO {get_user_group(get_current_username())}_group_db.{APP_NAME}_schema.app_sessions_table (app_session_id, username, user_group, container_image_id)
                VALUES (?, ?, ?, ?)
            """, row_tuple).collect()
            return True
        except Exception as e:
            st.error(f"Failed to log app session: {e}")
            return False


def set_app_session_shutdown_time(app_session_id):
    if framework_utils.platform() == "local":
        try:
            conn = get_database_connection(DB_URL_GROUP)
            with conn.cursor() as cur:
                cur.execute(f"""
                    UPDATE {APP_NAME}_schema.app_sessions_table
                    SET explicit_shutdown_time = %s
                    WHERE app_session_id = %s
                """, (framework_utils.get_timestamp(), app_session_id))
            conn.commit()
            return_database_connection(conn, DB_URL_GROUP)
            return True
        except Exception as e:
            st.error(f"Failed to set app session shutdown time: {e}")
            _rollback_and_return(conn, DB_URL_GROUP)
            return False
    elif framework_utils.platform() == "snowflake":
        try:
            session = snowflake_connections.get_snowpark_session()
            session.sql(f"""
                UPDATE {get_user_group(get_current_username())}_group_db.{APP_NAME}_schema.app_sessions_table
                SET explicit_shutdown_time = CURRENT_TIMESTAMP()
                WHERE app_session_id = ?
            """, (app_session_id,)).collect()
            return True
        except Exception as e:
            st.error(f"Failed to set app session shutdown time: {e}")
            return False


@st.cache_data()
def get_user_group(username):
    if framework_utils.platform() == "local":
        try:
            conn = get_database_connection(DB_URL_COMMON)
            with conn.cursor() as cur:
                cur.execute(f"""
                    SELECT user_group
                    FROM admin_schema.user_groups_table
                    WHERE username = %s
                """, (username,))
                user_group = cur.fetchone()
            return_database_connection(conn, DB_URL_COMMON)
            return user_group[0] if user_group else None
        except Exception as e:
            st.error(f"Failed to retrieve user group: {e}")
            if conn:
                return_database_connection(conn, DB_URL_COMMON)
            return None
    elif framework_utils.platform() == "snowflake":
        try:
            session = snowflake_connections.get_snowpark_session()
            result = session.sql(f"""
                SELECT user_group
                FROM common_db.admin_schema.user_groups_table
                WHERE username = ?
            """, (username,)).collect()
            return result[0]["USER_GROUP"] if result else None
        except Exception as e:
            st.error(f"Failed to retrieve user group: {e}")
            return None


@st.cache_data()
def get_user_groups_table_data():
    if framework_utils.platform() == "local":
        try:
            conn = get_database_connection(DB_URL_COMMON)
            with conn.cursor() as cur:
                cur.execute(f"""
                    SELECT username, user_group, user_added_time, who_added, user_email
                    FROM admin_schema.user_groups_table
                """)
                rows = cur.fetchall()
            return_database_connection(conn, DB_URL_COMMON)
            df = pl.DataFrame(rows, schema=["username", "user_group", "user_added_time", "who_added", "user_email"], strict=False, orient="row")
            return df
        except Exception as e:
            st.error(f"Failed to retrieve user groups table data: {e}")
            if conn:
                return_database_connection(conn, DB_URL_COMMON)
            return None
    elif framework_utils.platform() == "snowflake":
        try:
            session = snowflake_connections.get_snowpark_session()
            rows = session.sql(f"""
                SELECT username, user_group, user_added_time, who_added, user_email
                FROM common_db.admin_schema.user_groups_table
            """).collect()
            df = pl.DataFrame(rows, schema=["username", "user_group", "user_added_time", "who_added", "user_email"], strict=False, orient="row")
            return df
        except Exception as e:
            st.error(f"Failed to retrieve user groups table data: {e}")
            return None


@st.cache_data()
def get_app_sessions_table_data():
    if framework_utils.platform() == "local":
        try:
            conn = get_database_connection(DB_URL_GROUP)
            with conn.cursor() as cur:
                cur.execute(f"""
                    SELECT app_session_id, username, user_group, startup_time, explicit_shutdown_time, container_image_id
                    FROM {APP_NAME}_schema.app_sessions_table
                    ORDER BY startup_time DESC
                """)
                rows = cur.fetchall()
            return_database_connection(conn, DB_URL_GROUP)
            df = pl.DataFrame(rows, schema=["app_session_id", "username", "user_group", "startup_time", "explicit_shutdown_time", "container_image_id"], strict=False, orient="row")
            return df
        except Exception as e:
            st.error(f"Failed to retrieve app sessions table data: {e}")
            if conn:
                return_database_connection(conn, DB_URL_GROUP)
            return None
    elif framework_utils.platform() == "snowflake":
        try:
            session = snowflake_connections.get_snowpark_session()
            rows = session.sql(f"""
                SELECT app_session_id, username, user_group, startup_time, explicit_shutdown_time, container_image_id
                FROM {get_user_group(get_current_username())}_group_db.{APP_NAME}_schema.app_sessions_table
                ORDER BY startup_time DESC
            """).collect()
            df = pl.DataFrame(rows, schema=["app_session_id", "username", "user_group", "startup_time", "explicit_shutdown_time", "container_image_id"], strict=False, orient="row")
            return df
        except Exception as e:
            st.error(f"Failed to retrieve app sessions table data: {e}")
            return None


@st.cache_data()
def get_archives_table_data():
    if framework_utils.platform() == "local":
        try:
            conn = get_database_connection(DB_URL_GROUP)
            with conn.cursor() as cur:
                cur.execute(f"""
                    SELECT creator, user_group, archive_description, current_git_commit, container_image_id, archive_id, app_session_id, creation_time
                    FROM {APP_NAME}_schema.archives_table
                    ORDER BY creation_time DESC
                """)
                rows = cur.fetchall()
            return_database_connection(conn, DB_URL_GROUP)
            df = pl.DataFrame(rows, schema=["creator", "user_group", "archive_description", "current_git_commit", "container_image_id", "archive_id", "app_session_id", "creation_time"], strict=False, orient="row")
            return df
        except Exception as e:
            st.error(f"Failed to retrieve archives_table table data: {e}")
            if conn:
                return_database_connection(conn, DB_URL_GROUP)
            return None
    elif framework_utils.platform() == "snowflake":
        try:
            session = snowflake_connections.get_snowpark_session()
            rows = session.sql(f"""
                SELECT creator, user_group, archive_description, current_git_commit, container_image_id, archive_id, app_session_id, creation_time
                FROM {get_user_group(get_current_username())}_group_db.{APP_NAME}_schema.archives_table
                ORDER BY creation_time DESC
            """).collect()
            df = pl.DataFrame(rows, schema=["creator", "user_group", "archive_description", "current_git_commit", "container_image_id", "archive_id", "app_session_id", "creation_time"], strict=False, orient="row")
            return df
        except Exception as e:
            st.error(f"Failed to retrieve archives_table table data: {e}")
            return None


@st.cache_data()
def get_jobs_table_data():
    if framework_utils.platform() == "local":
        try:
            conn = get_database_connection(DB_URL_GROUP)
            with conn.cursor() as cur:
                columns = "job_id, job_name, job_status, submitter, submitter_group, app_session_id, worker_image_id, submission_time, start_time, completion_time, failure_time"
                cur.execute(f"""
                    SELECT {columns}
                    FROM {APP_NAME}_schema.jobs_table
                    ORDER BY submission_time DESC NULLS LAST
                """)
                rows = cur.fetchall()
            return_database_connection(conn, DB_URL_GROUP)
            df = pl.DataFrame(rows, schema=["job_id", "job_name", "job_status", "submitter", "submitter_group", "app_session_id", "worker_image_id", "submission_time", "start_time", "completion_time", "failure_time"], strict=False, orient="row")
            return df
        except Exception as e:
            st.error(f"Failed to retrieve jobs_table table data: {e}")
            if conn:
                return_database_connection(conn, DB_URL_GROUP)
            return None
    elif framework_utils.platform() == "snowflake":
        try:
            session = snowflake_connections.get_snowpark_session()
            columns = "job_id, job_name, job_status, submitter, submitter_group, app_session_id, worker_image_id, submission_time, start_time, completion_time, failure_time"
            rows = session.sql(f"""
                SELECT {columns}
                FROM {get_user_group(get_current_username())}_group_db.{APP_NAME}_schema.jobs_table
                ORDER BY submission_time DESC NULLS LAST
            """).collect()
            df = pl.DataFrame(rows, schema=["job_id", "job_name", "job_status", "submitter", "submitter_group", "app_session_id", "worker_image_id", "submission_time", "start_time", "completion_time", "failure_time"], strict=False, orient="row")
            return df
        except Exception as e:
            st.error(f"Failed to retrieve jobs_table table data: {e}")
            return None


@st.cache_data()
def get_available_archives():
    if framework_utils.platform() == "local":
        try:
            conn = get_database_connection(DB_URL_GROUP)
            with conn.cursor() as cur:
                cur.execute(f"""
                    SELECT creator, creation_time, archive_description, archive_id, app_session_id
                    FROM {APP_NAME}_schema.archives_table
                    ORDER BY creation_time DESC
                """)
                rows = cur.fetchall()
            return_database_connection(conn, DB_URL_GROUP)
            df = pl.DataFrame(rows, schema=["Creator", "Creation time", "Archive description", "Archive ID", "App session ID"], strict=False, orient="row")
            return df
        except Exception as e:
            st.error(f"Failed to retrieve available archives_table: {e}")
            if conn:
                return_database_connection(conn, DB_URL_GROUP)
            return None
    elif framework_utils.platform() == "snowflake":
        try:
            session = snowflake_connections.get_snowpark_session()
            rows = session.sql(f"""
                SELECT creator, creation_time, archive_description, archive_id, app_session_id
                FROM {get_user_group(get_current_username())}_group_db.{APP_NAME}_schema.archives_table
                ORDER BY creation_time DESC
            """).collect()
            df = pl.DataFrame(rows, schema=["Creator", "Creation time", "Archive description", "Archive ID", "App session ID"], strict=False, orient="row")
            return df
        except Exception as e:
            st.error(f"Failed to retrieve available archives_table: {e}")
            return None


def log_job(row_tuple):
    if framework_utils.platform() == "local":
        try:
            conn = get_database_connection(DB_URL_GROUP)
            with conn.cursor() as cur:
                cur.execute(f"""
                    INSERT INTO {APP_NAME}_schema.jobs_table (job_id, job_name, submitter, submitter_group, app_session_id)
                    VALUES (%s, %s, %s, %s, %s)
                """, row_tuple)
            conn.commit()
            return_database_connection(conn, DB_URL_GROUP)
            return True
        except Exception as e:
            st.error(f"Failed to log job: {e}. It's possible that job with ID {row_tuple[0]} already exists (unique constraint violated), which would indicate a job ID generation bug.")
            _rollback_and_return(conn, DB_URL_GROUP)
            return False
    elif framework_utils.platform() == "snowflake":
        try:
            session = snowflake_connections.get_snowpark_session()
            session.sql(f"""
                INSERT INTO {get_user_group(get_current_username())}_group_db.{APP_NAME}_schema.jobs_table (job_id, job_name, submitter, submitter_group, app_session_id)
                VALUES (?, ?, ?, ?, ?)
            """, row_tuple).collect()
            return True
        except Exception as e:
            st.error(f"Failed to log job: {e}. It's possible that job with ID {row_tuple[0]} already exists (unique constraint violated), which would indicate a job ID generation bug.")
            return False


def update_job_status(job_id, new_status, time_column):
    if framework_utils.platform() == "local":
        try:
            conn = get_database_connection(DB_URL_GROUP)
            with conn.cursor() as cur:
                cur.execute(f"""
                    UPDATE {APP_NAME}_schema.jobs_table
                    SET job_status = %s, {time_column} = %s
                    WHERE job_id = %s
                """, (new_status, framework_utils.get_timestamp(), job_id))
            conn.commit()
            return_database_connection(conn, DB_URL_GROUP)
            return True
        except Exception as e:
            st.error(f"Failed to update status of job {job_id} to {new_status} and update {time_column}: {e}")
            _rollback_and_return(conn, DB_URL_GROUP)
            return False
    elif framework_utils.platform() == "snowflake":
        try:
            session = snowflake_connections.get_snowpark_session()
            session.sql(f"""
                UPDATE {get_user_group(get_current_username())}_group_db.{APP_NAME}_schema.jobs_table
                SET job_status = ?, {time_column} = CURRENT_TIMESTAMP()
                WHERE job_id = ?
            """, (new_status, job_id)).collect()
            return True
        except Exception as e:
            st.error(f"Failed to update status of job {job_id} to {new_status} and update {time_column}: {e}")
            return False


def get_job_status(job_id):
    if framework_utils.platform() == "local":
        try:
            conn = get_database_connection(DB_URL_GROUP)
            with conn.cursor() as cur:
                cur.execute(f"""
                    SELECT job_status
                    FROM {APP_NAME}_schema.jobs_table
                    WHERE job_id = %s
                """, (job_id,))
                job_status = cur.fetchone()
            return_database_connection(conn, DB_URL_GROUP)
            return job_status[0] if job_status else None
        except Exception as e:
            st.error(f"Failed to retrieve status of job {job_id}: {e}")
            if conn:
                return_database_connection(conn, DB_URL_GROUP)
            return None
    elif framework_utils.platform() == "snowflake":
        try:
            session = snowflake_connections.get_snowpark_session()
            result = session.sql(f"""
                SELECT job_status
                FROM {get_user_group(get_current_username())}_group_db.{APP_NAME}_schema.jobs_table
                WHERE job_id = ?
            """, (job_id,)).collect()
            return result[0]["JOB_STATUS"] if result else None
        except Exception as e:
            st.error(f"Failed to retrieve status of job {job_id}: {e}")
            return None


@st.cache_data()
def get_job_function_name(job_id):
    if framework_utils.platform() == "local":
        try:
            conn = get_database_connection(DB_URL_GROUP)
            with conn.cursor() as cur:
                cur.execute(f"""
                    SELECT job_name
                    FROM {APP_NAME}_schema.jobs_table
                    WHERE job_id = %s
                """, (job_id,))
                job_name = cur.fetchone()
            return_database_connection(conn, DB_URL_GROUP)
            return job_name[0] if job_name else None
        except Exception as e:
            st.error(f"Failed to retrieve function name of job {job_id}: {e}")
            if conn:
                return_database_connection(conn, DB_URL_GROUP)
            return None
    elif framework_utils.platform() == "snowflake":
        try:
            session = snowflake_connections.get_snowpark_session()
            result = session.sql(f"""
                SELECT job_name
                FROM {get_user_group(get_current_username())}_group_db.{APP_NAME}_schema.jobs_table
                WHERE job_id = ?
            """, (job_id,)).collect()
            return result[0]["JOB_NAME"] if result else None
        except Exception as e:
            st.error(f"Failed to retrieve function name of job {job_id}: {e}")
            return None


def set_worker_image_id(job_id, worker_image_id):
    if framework_utils.platform() == "local":
        try:
            conn = get_database_connection(DB_URL_GROUP)
            with conn.cursor() as cur:
                cur.execute(f"""
                    UPDATE {APP_NAME}_schema.jobs_table
                    SET worker_image_id = %s
                    WHERE job_id = %s
                """, (worker_image_id, job_id))
            conn.commit()
            return_database_connection(conn, DB_URL_GROUP)
            return True
        except Exception as e:
            st.error(f"Failed to set worker image ID for job {job_id} to {worker_image_id}: {e}")
            _rollback_and_return(conn, DB_URL_GROUP)
            return False
    elif framework_utils.platform() == "snowflake":
        try:
            session = snowflake_connections.get_snowpark_session()
            session.sql(f"""
                UPDATE {get_user_group(get_current_username())}_group_db.{APP_NAME}_schema.jobs_table
                SET worker_image_id = ?
                WHERE job_id = ?
            """, (worker_image_id, job_id)).collect()
            return True
        except Exception as e:
            st.error(f"Failed to set worker image ID for job {job_id} to {worker_image_id}: {e}")
            return False


def record_explicit_shutdown_time(app_session_id):
    if framework_utils.platform() == "local":
        try:
            conn = get_database_connection(DB_URL_GROUP)
            with conn.cursor() as cur:
                cur.execute(f"""
                    UPDATE {APP_NAME}_schema.app_sessions_table
                    SET explicit_shutdown_time = %s
                    WHERE app_session_id = %s
                """, (framework_utils.get_timestamp(), app_session_id))
            conn.commit()
            return_database_connection(conn, DB_URL_GROUP)
            return True
        except Exception as e:
            st.error(f"Failed to record explicit shutdown time for app session {app_session_id}: {e}")
            _rollback_and_return(conn, DB_URL_GROUP)
            return False
    elif framework_utils.platform() == "snowflake":
        try:
            session = snowflake_connections.get_snowpark_session()
            session.sql(f"""
                UPDATE {get_user_group(get_current_username())}_group_db.{APP_NAME}_schema.app_sessions_table
                SET explicit_shutdown_time = CURRENT_TIMESTAMP()
                WHERE app_session_id = ?
            """, (app_session_id,)).collect()
            return True
        except Exception as e:
            st.error(f"Failed to record explicit shutdown time for app session {app_session_id}: {e}")
            return False


#### 2. OBJECT STORAGE FUNCTIONALITY ##############################################################


MINIO_ENDPOINT = os.getenv('MINIO_ENDPOINT')
MINIO_ACCESS_KEY = os.getenv('MINIO_ACCESS_KEY')
MINIO_SECRET_KEY = os.getenv('MINIO_SECRET_KEY')
ARCHIVES_BUCKET_NAME = os.getenv('ARCHIVES_BUCKET_NAME')
OLD_ARCHIVES_BUCKET_NAME = os.getenv('OLD_ARCHIVES_BUCKET_NAME')
JOB_INPUTS_BUCKET_NAME = os.getenv('JOB_INPUTS_BUCKET_NAME')
JOB_OUTPUTS_BUCKET_NAME = os.getenv('JOB_OUTPUTS_BUCKET_NAME')
DATA_OBJECTS_BUCKET_NAME = os.getenv('DATA_OBJECTS_BUCKET_NAME')


@st.cache_resource()
def get_object_storage_client():
    if framework_utils.platform() == "local":
        try:
            return minio.Minio(
                MINIO_ENDPOINT,
                access_key=MINIO_ACCESS_KEY,
                secret_key=MINIO_SECRET_KEY,
                secure=False  # Set to True for HTTPS
            )
        except Exception as e:
            st.error(f"Failed to create MinIO client: {e}")
            return None
    elif framework_utils.platform() == "snowflake":
        pass


@st.cache_data()
def set_up_minio():
    if framework_utils.platform() == "local":
        try:
            client = get_object_storage_client()
            if not client.bucket_exists(ARCHIVES_BUCKET_NAME):
                client.make_bucket(ARCHIVES_BUCKET_NAME)
            if not client.bucket_exists(JOB_INPUTS_BUCKET_NAME):
                client.make_bucket(JOB_INPUTS_BUCKET_NAME)
            if not client.bucket_exists(JOB_OUTPUTS_BUCKET_NAME):
                client.make_bucket(JOB_OUTPUTS_BUCKET_NAME)
            if not client.bucket_exists(DATA_OBJECTS_BUCKET_NAME):
                client.make_bucket(DATA_OBJECTS_BUCKET_NAME)
            if not client.bucket_exists(OLD_ARCHIVES_BUCKET_NAME):
                client.make_bucket(OLD_ARCHIVES_BUCKET_NAME)
            return True
        except Exception as e:
            st.error(f"Failed to set up object storage: {e}")
            return False
    elif framework_utils.platform() == "snowflake":
        pass


def upload_zip_object_data(bucket_name, zip_name, zip_buffer, db_schema: str = None):
    if framework_utils.platform() == "local":
        try:
            client = get_object_storage_client()

            # Get the size of the zip buffer
            zip_buffer.seek(0, io.SEEK_END)
            zip_size = zip_buffer.tell()
            zip_buffer.seek(0)

            client.put_object(
                bucket_name=bucket_name,
                object_name=f"{zip_name}.zip",
                data=zip_buffer,
                length=zip_size,
                content_type='application/zip'
            )
            return True
        except Exception as e:
            st.error(f"Failed to write {zip_name}.zip to bucket {bucket_name}: {e}")
            return False
    elif framework_utils.platform() == "snowflake":
        try:
            session = snowflake_connections.get_snowpark_session()
            zip_buffer.seek(0)
            if db_schema is None:
                db_schema = f"{get_user_group(get_current_username())}_group_db.{APP_NAME}_schema"
            results = session.file.put_stream(
                input_stream=zip_buffer,
                stage_location=f"@{db_schema}.{bucket_name}_stage/{zip_name}.zip",
                auto_compress=False,
            )
            return results
        except Exception as e:
            st.error(f"Failed to write {zip_name}.zip to bucket {bucket_name}: {e}")
            return None


# This could potentially be a lot of data, so we don't want to cache it using st.cache_data().
def download_zip_object_data(bucket_name, zip_name, db_schema: str = None):
    if framework_utils.platform() == "local":
        response = None
        try:
            client = get_object_storage_client()
            object_name = f"{zip_name}.zip"
            response = client.get_object(bucket_name, object_name)
            # NOTE: For very large objects consider streaming in chunks instead of reading all at once.
            data = response.read()
            response.close()
            return io.BytesIO(data)
        except Exception as e:
            st.error(f"Failed to download object data: {e}")
            if response:
                response.close()
            return None
    elif framework_utils.platform() == "snowflake":
        try:
            session = snowflake_connections.get_snowpark_session()
            if db_schema is None:
                db_schema = f"{get_user_group(get_current_username())}_group_db.{APP_NAME}_schema"
            stage_path = f"@{db_schema}.{bucket_name}_stage/{zip_name}.zip"
            bytes_io = session.file.get_stream(stage_location=stage_path)
            bytes_io.seek(0)
            return bytes_io
        except Exception as e:
            st.error(f"Failed to download {zip_name}.zip from bucket {bucket_name}: {e}")
            return None
        

def list_objects_in_bucket(bucket_name: str, db_schema: str = None):
    if framework_utils.platform() == "local":
        try:
            client = get_object_storage_client()
            objects = client.list_objects(bucket_name)
            object_list = [obj.object_name for obj in objects]
            return object_list
        except Exception as e:
            st.error(f"Failed to list objects in {bucket_name} bucket: {e}")
            return None
    elif framework_utils.platform() == "snowflake":
        try:
            session = snowflake_connections.get_snowpark_session()
            if db_schema is None:
                user_group = get_user_group(get_current_username())
                db_schema = f"{user_group}_group_db.curated_schema"
            stage_location = f"@{db_schema}.{bucket_name}_stage"
            files = session.file.list(stage_location=stage_location)
            object_list = [file['name'] for file in files]
            return object_list
        except Exception as e:
            st.error(f"Failed to list objects in {bucket_name} stage in database.schema {db_schema}: {e}")
            return None


def download_objects_parallel(
    bucket_name: str,
    object_names: list[str],
    dest_dir: str,
    db_schema: str = None,
    max_workers: int = 8,
    chunk_size: int = 1024 * 1024,
    retries: int = 3,
    backoff_base: float = 0.3,
    verify_etag: bool = False
):
    """
    Download arbitrary objects from MinIO in parallel.
    object_names: list of object keys exactly as stored (may include extensions or paths).
    dest_dir: local directory root to place downloaded files (object key subpaths preserved).
    Returns dict {object_name: {'status': 'ok', 'path': local_path} or {'status': 'error', 'error': Exception}}.
    """
    if framework_utils.platform() == "local":
        try:
            client = get_object_storage_client()

            os.makedirs(dest_dir, exist_ok=True)

            def download_one(obj_name: str):
                attempt = 0
                while attempt < retries:
                    response = None
                    try:
                        response = client.get_object(bucket_name, obj_name)
                        local_path = os.path.join(dest_dir, obj_name)
                        parent = os.path.dirname(local_path)
                        if parent:
                            os.makedirs(parent, exist_ok=True)
                        with open(local_path, "wb") as f:
                            while True:
                                data = response.read(chunk_size)
                                if not data:
                                    break
                                f.write(data)
                        response.close()

                        if verify_etag:
                            # NOTE: For multipart uploads the ETag is not a simple MD5; skip strict verify in that case.
                            stat = client.stat_object(bucket_name, obj_name)
                            etag = getattr(stat, "etag", None)
                            if etag and "-" not in etag:  # Heuristic: multipart ETags contain '-'
                                h = hashlib.md5()
                                with open(local_path, "rb") as f:
                                    for block in iter(lambda: f.read(chunk_size), b""):
                                        h.update(block)
                                if h.hexdigest() != etag.lower():
                                    raise ValueError(f"ETag mismatch for {obj_name}: expected {etag}, got {h.hexdigest()}")
                        return obj_name, {'status': 'ok', 'path': local_path}
                    except Exception as e:
                        if response:
                            try:
                                response.close()
                            except Exception:
                                pass
                        attempt += 1
                        if attempt < retries:
                            time.sleep(backoff_base * (2 ** (attempt - 1)))
                        else:
                            return obj_name, {'status': 'error', 'error': e}

            results: dict[str, dict] = {}
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                future_map = {executor.submit(download_one, name): name for name in object_names}
                for fut in as_completed(future_map):
                    name = future_map[fut]
                    try:
                        obj_name, res = fut.result()
                        results[obj_name] = res
                    except Exception as e:
                        results[name] = {'status': 'error', 'error': e}

            failures = [k for k, v in results.items() if v['status'] == 'error']
            if failures:
                st.warning(f"{len(failures)} downloads failed.")
            else:
                st.success(f"Downloaded {len(results)} objects.")
            return results
        except Exception as e:
            st.error(f"Failed to download objects from MinIO bucket: {e}")
            return None
    elif framework_utils.platform() == "snowflake":
        try:
            session = snowflake_connections.get_snowpark_session()
            os.makedirs(dest_dir, exist_ok=True)

            if db_schema is None:
                user_group = get_user_group(get_current_username())
                db_schema = f"{user_group}_group_db.curated_schema"
            stage_prefix = f"@{db_schema}.{bucket_name}_stage"

            def download_one(obj_name: str):
                attempt = 0
                while attempt < retries:
                    stream = None
                    try:
                        stage_location = f"{stage_prefix}/{obj_name}"
                        stream = session.file.get_stream(stage_location=stage_location)
                        local_path = os.path.join(dest_dir, obj_name)
                        parent = os.path.dirname(local_path)
                        if parent:
                            os.makedirs(parent, exist_ok=True)
                        with open(local_path, "wb") as f:
                            while True:
                                chunk = stream.read(chunk_size)
                                if not chunk:
                                    break
                                f.write(chunk)
                        stream.close()

                        if verify_etag:
                            # Use MD5 from stage listing when available.
                            dir_part = os.path.dirname(obj_name)
                            list_location = stage_prefix + ("/" + dir_part if dir_part else "")
                            files = session.file.list(stage_location=list_location)
                            meta = next((m for m in files if m.get("name") == os.path.basename(obj_name)), None)
                            remote_md5 = meta.get("md5") if meta else None
                            if remote_md5:
                                h = hashlib.md5()
                                with open(local_path, "rb") as f:
                                    for block in iter(lambda: f.read(chunk_size), b""):
                                        h.update(block)
                                if h.hexdigest().lower() != remote_md5.lower():
                                    raise ValueError(f"MD5 mismatch for {obj_name}: expected {remote_md5}, got {h.hexdigest()}")
                        return obj_name, {'status': 'ok', 'path': local_path}
                    except Exception as e:
                        attempt += 1
                        if stream:
                            try:
                                stream.close()
                            except Exception:
                                pass
                        if attempt < retries:
                            time.sleep(backoff_base * (2 ** (attempt - 1)))
                        else:
                            return obj_name, {'status': 'error', 'error': e}

            results: dict[str, dict] = {}
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                future_map = {executor.submit(download_one, name): name for name in object_names}
                for fut in as_completed(future_map):
                    name = future_map[fut]
                    try:
                        obj_name, res = fut.result()
                        results[obj_name] = res
                    except Exception as e:
                        results[name] = {'status': 'error', 'error': e}

            failures = [k for k, v in results.items() if v['status'] == 'error']
            if failures:
                st.warning(f"{len(failures)} downloads failed.")
            else:
                st.success(f"Downloaded {len(results)} objects.")
            return results
        except Exception as e:
            st.error(f"Failed to download objects from Snowflake stage: {e}")
            return None


def upload_objects_parallel(
    bucket_name: str,
    file_paths: list[str],
    base_dir: str | None = None,
    db_schema: str = None,
    max_workers: int = 8,
    retries: int = 3,
    backoff_base: float = 0.3,
    guess_content_type: bool = True,
    verify_etag: bool = False
):
    """
    Upload arbitrary local files to object storage in parallel.

    file_paths: list of absolute or relative local file paths.
    base_dir: if provided, object names are path-relative to this directory; else filename only.
    Returns dict {object_name: {'status': 'ok', 'size': bytes} or {'status': 'error', 'error': Exception}}.
    """

    def compute_object_name(path: str):
        if base_dir:
            b = os.path.abspath(base_dir)
            p = os.path.abspath(path)
            rel = os.path.relpath(p, b)
            if rel.startswith(".."):
                raise ValueError(f"{path} outside base_dir {base_dir}")
            return rel.replace("\\", "/")
        return os.path.basename(path)

    if framework_utils.platform() == "local":
        try:
            client = get_object_storage_client()

            def upload_one(local_path: str):
                object_name = compute_object_name(local_path)
                attempt = 0
                while attempt < retries:
                    f = None
                    try:
                        size = os.path.getsize(local_path)
                        f = open(local_path, "rb")
                        content_type = None
                        if guess_content_type:
                            content_type = mimetypes.guess_type(local_path)[0] or "application/octet-stream"
                        client.put_object(
                            bucket_name=bucket_name,
                            object_name=object_name,
                            data=f,
                            length=size,
                            content_type=content_type
                        )
                        if f:
                            f.close()

                        if verify_etag:
                            # Only reliable for single-part uploads (heuristic: size <= 64MB default MinIO part size).
                            stat = client.stat_object(bucket_name, object_name)
                            etag = getattr(stat, "etag", None)
                            if etag and "-" not in etag:
                                h = hashlib.md5()
                                with open(local_path, "rb") as vf:
                                    for block in iter(lambda: vf.read(1024 * 1024), b""):
                                        h.update(block)
                                if h.hexdigest() != etag.lower():
                                    raise ValueError(f"ETag mismatch for {object_name}: expected {etag}, got {h.hexdigest()}")
                        return object_name, {'status': 'ok', 'size': size}
                    except Exception as e:
                        attempt += 1
                        if f:
                            try:
                                f.close()
                            except Exception:
                                pass
                        if attempt < retries:
                            time.sleep(backoff_base * (2 ** (attempt - 1)))
                        else:
                            return object_name, {'status': 'error', 'error': e}

            results: dict[str, dict] = {}
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                future_map = {executor.submit(upload_one, p): p for p in file_paths}
                for fut in as_completed(future_map):
                    try:
                        obj_name, res = fut.result()
                        results[obj_name] = res
                    except Exception as e:
                        # Should rarely happen since we catch inside fetch, but just in case.
                        fallback_name = compute_object_name(future_map[fut])
                        results[fallback_name] = {'status': 'error', 'error': e}

            failures = [k for k, v in results.items() if v['status'] != 'ok']
            if failures:
                st.warning(f"{len(failures)} uploads failed.")
            else:
                st.success(f"Uploaded {len(results)} objects.")
            return results
        except Exception as e:
            st.error(f"Failed to upload objects to MinIO: {e}")
            return None

    elif framework_utils.platform() == "snowflake":
        try:
            session = snowflake_connections.get_snowpark_session()

            if db_schema is None:
                user_group = get_user_group(get_current_username())
                db_schema = f"{user_group}_group_db.curated_schema"
            stage_prefix = f"@{db_schema}.{bucket_name}_stage"

            def upload_one(local_path: str):
                object_name = compute_object_name(local_path)
                attempt = 0
                while attempt < retries:
                    f = None
                    try:
                        size = os.path.getsize(local_path)
                        f = open(local_path, "rb")
                        session.file.put_stream(
                            input_stream=f,
                            stage_location=f"{stage_prefix}/{object_name}",
                            auto_compress=False
                        )
                        if f:
                            f.close()

                        if verify_etag:
                            # List the directory to get md5.
                            dir_part = os.path.dirname(object_name)
                            list_location = stage_prefix + ("/" + dir_part if dir_part else "")
                            files = session.file.list(stage_location=list_location)
                            meta = next((m for m in files if m.get("name") == os.path.basename(object_name)), None)
                            remote_md5 = meta.get("md5") if meta else None
                            if remote_md5:
                                h = hashlib.md5()
                                with open(local_path, "rb") as vf:
                                    for block in iter(lambda: vf.read(1024 * 1024), b""):
                                        h.update(block)
                                if h.hexdigest().lower() != remote_md5.lower():
                                    raise ValueError(f"MD5 mismatch for {object_name}: expected {remote_md5}, got {h.hexdigest()}")
                        return object_name, {'status': 'ok', 'size': size}
                    except Exception as e:
                        attempt += 1
                        if f:
                            try:
                                f.close()
                            except Exception:
                                pass
                        if attempt < retries:
                            time.sleep(backoff_base * (2 ** (attempt - 1)))
                        else:
                            return object_name, {'status': 'error', 'error': e}

            results: dict[str, dict] = {}
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                future_map = {executor.submit(upload_one, p): p for p in file_paths}
                for fut in as_completed(future_map):
                    try:
                        obj_name, res = fut.result()
                        results[obj_name] = res
                    except Exception as e:
                        fallback_name = compute_object_name(future_map[fut])
                        results[fallback_name] = {'status': 'error', 'error': e}

            failures = [k for k, v in results.items() if v['status'] != 'ok']
            if failures:
                st.warning(f"{len(failures)} uploads failed.")
            else:
                st.success(f"Uploaded {len(results)} objects.")
            return results
        except Exception as e:
            st.error(f"Failed to upload objects to Snowflake stage: {e}")
            return None


#### 3. ORCHESTRATION FUNCTIONALITY ###############################################################


@st.cache_data()
def get_frontend_image_id():
    if framework_utils.platform() == "local":
        try:
            resp = requests.get("http://docker_orchestrator:8080/frontend_id", timeout=3)
            resp.raise_for_status()
            data = resp.json()
            return data.get("frontend_image_id")
        except Exception as e:
            st.error(f"Could not retrieve frontend image id: {e}")
            return None
    elif framework_utils.platform() == "snowflake":
        try:
            session = snowflake_connections.get_snowpark_session()
            return snowflake_orchestrator.frontend_id(username=get_current_username(), session=session)
        except Exception as e:
            st.error(f"Could not retrieve frontend image id: {e}")
            return None


def submit_job(job_id, blocking=True):
    if framework_utils.platform() == "local":
        try:
            update_job_status(job_id, "Submitted", "submission_time")
            if blocking:
                analysis_framework.run_local_analysis(job_id)  # This is the worker code.
            else:  # Asynchronous analysis that should execute analysis_framework.run_local_analysis(job_id) on a worker.
                # Submit to Docker orchestrator, which launches an ephemeral worker container
                resp = requests.post("http://docker_orchestrator:8080/jobs", json={"job_id": job_id}, timeout=10)
                resp.raise_for_status()
                data = resp.json()
                worker_image_id = data.get("worker_image_id")
                # Record the image used by the worker container for traceability
                if worker_image_id:
                    set_worker_image_id(job_id, worker_image_id)
            return True
        except Exception as e:
            st.error(f"Failed to submit job: {e}")
            return False
    elif framework_utils.platform() == "snowflake":
        try:
            update_job_status(job_id, "Submitted", "submission_time")
            if blocking:
                analysis_framework.run_local_analysis(job_id)  # This is the worker code.
            else:
                session = snowflake_connections.get_snowpark_session()
                worker_image_id = snowflake_orchestrator.submit_job(job_id=job_id, username=get_current_username(), session=session)
                if worker_image_id:
                    set_worker_image_id(job_id, worker_image_id)
            return True
        except Exception as e:
            st.error(f"Failed to submit job: {e}")
            return False


def shut_down_app():
    if framework_utils.platform() == "local":
        try:
            resp = requests.post("http://docker_orchestrator:8080/shutdown", timeout=10)
            resp.raise_for_status()
            st.success("Application is shutting down...")
            return True
        except Exception as e:
            st.error(f"Failed to shut down app: {e}")
            return False
    elif framework_utils.platform() == "snowflake":
        try:
            session = snowflake_connections.get_snowpark_session()
            snowflake_orchestrator.shutdown(username=get_current_username(), session=session)
            st.success("Application is shutting down...")
            return True
        except Exception as e:
            st.error(f"Failed to shut down app: {e}")
            return False


#### 4. OTHER FUNCTIONALITY #######################################################################


@st.cache_data()
def get_current_username():
    if framework_utils.platform() == "local":
        try:
            try:
                return os.getlogin()
            except (FileNotFoundError, OSError):
                pass
            username = os.getenv('APP_USER') or os.getenv('USER') or os.getenv('USERNAME') or os.getenv('LOGNAME')  # Delete "os.getenv('APP_USER') or " (and modify docker-compose.yml) in the future when we implement running the container as non-root.
            if username:
                return username
            return getpass.getuser()
        except Exception as e:
            st.error(f"Failed to get current username: {e}")
            return None
    elif framework_utils.platform() == "snowflake":
        try:
            return os.getenv("SNOWFLAKE_USER")
        except Exception as e:
            st.error(f"Failed to get current username: {e}")
            return None
