import streamlit as st
import os
import io
import psycopg2.pool  # Probably change to version 3 (psycopg) soon!
import minio
import getpass
import polars as pl
import requests
import atexit
import framework.utils as utils
import framework.analysis_framework as analysis_framework
import framework.snowflake_connections as snowflake_connections
import framework.snowflake_orchestrator as snowflake_orchestrator


ST_KEY_PREFIX_STARTUP = "startup.py__"


#### 1. DATABASE FUNCTIONALITY ####################################################################


DB_URL_GROUP = os.getenv('DB_URL_GROUP')
DB_URL_COMMON = os.getenv('DB_URL_COMMON')
APP_NAME = os.getenv('APP_NAME')


@st.cache_resource()
def get_database_pool(db_url):
    if utils.platform() == "local":
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
    elif utils.platform() == "snowflake":
        pass  # Like everything left as "pass", this is not needed on Snowflake.


# Note we could set this and the following up using @contextmanager as in 8/28/25 chat with GH Copilot, but keeping out for now for simplicity.
def get_database_connection(db_url):
    if utils.platform() == "local":
        try:
            db_pool = get_database_pool(db_url)
            if db_pool:
                return db_pool.getconn()
        except Exception as e:
            st.error(f"Failed to get database connection: {e}")
            return None
    elif utils.platform() == "snowflake":
        pass


def return_database_connection(conn, db_url):
    if utils.platform() == "local":
        try:
            db_pool = get_database_pool(db_url)
            if db_pool and conn:
                db_pool.putconn(conn)
            return True
        except Exception as e:
            st.error(f"Failed to return database connection: {e}")
            return False
    elif utils.platform() == "snowflake":
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
def set_up_database():
    if utils.platform() == "local":
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
    elif utils.platform() == "snowflake":
        pass


def write_archive_database_data(row_tuple):
    if utils.platform() == "local":
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
    elif utils.platform() == "snowflake":
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
    if utils.platform() == "local":
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
    elif utils.platform() == "snowflake":
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
    if utils.platform() == "local":
        try:
            conn = get_database_connection(DB_URL_GROUP)
            with conn.cursor() as cur:
                cur.execute(f"""
                    UPDATE {APP_NAME}_schema.app_sessions_table
                    SET explicit_shutdown_time = %s
                    WHERE app_session_id = %s
                """, (utils.get_timestamp(), app_session_id))
            conn.commit()
            return_database_connection(conn, DB_URL_GROUP)
            return True
        except Exception as e:
            st.error(f"Failed to set app session shutdown time: {e}")
            _rollback_and_return(conn, DB_URL_GROUP)
            return False
    elif utils.platform() == "snowflake":
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
    if utils.platform() == "local":
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
    elif utils.platform() == "snowflake":
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
    if utils.platform() == "local":
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
    elif utils.platform() == "snowflake":
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
    if utils.platform() == "local":
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
    elif utils.platform() == "snowflake":
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
    if utils.platform() == "local":
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
    elif utils.platform() == "snowflake":
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
    if utils.platform() == "local":
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
    elif utils.platform() == "snowflake":
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
    if utils.platform() == "local":
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
    elif utils.platform() == "snowflake":
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
    if utils.platform() == "local":
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
    elif utils.platform() == "snowflake":
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
    if utils.platform() == "local":
        try:
            conn = get_database_connection(DB_URL_GROUP)
            with conn.cursor() as cur:
                cur.execute(f"""
                    UPDATE {APP_NAME}_schema.jobs_table
                    SET job_status = %s, {time_column} = %s
                    WHERE job_id = %s
                """, (new_status, utils.get_timestamp(), job_id))
            conn.commit()
            return_database_connection(conn, DB_URL_GROUP)
            return True
        except Exception as e:
            st.error(f"Failed to update status of job {job_id} to {new_status} and update {time_column}: {e}")
            _rollback_and_return(conn, DB_URL_GROUP)
            return False
    elif utils.platform() == "snowflake":
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
    if utils.platform() == "local":
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
    elif utils.platform() == "snowflake":
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
    if utils.platform() == "local":
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
    elif utils.platform() == "snowflake":
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
    if utils.platform() == "local":
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
    elif utils.platform() == "snowflake":
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
    if utils.platform() == "local":
        try:
            conn = get_database_connection(DB_URL_GROUP)
            with conn.cursor() as cur:
                cur.execute(f"""
                    UPDATE {APP_NAME}_schema.app_sessions_table
                    SET explicit_shutdown_time = %s
                    WHERE app_session_id = %s
                """, (utils.get_timestamp(), app_session_id))
            conn.commit()
            return_database_connection(conn, DB_URL_GROUP)
            return True
        except Exception as e:
            st.error(f"Failed to record explicit shutdown time for app session {app_session_id}: {e}")
            _rollback_and_return(conn, DB_URL_GROUP)
            return False
    elif utils.platform() == "snowflake":
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
JOB_INPUTS_BUCKET_NAME = os.getenv('JOB_INPUTS_BUCKET_NAME')
JOB_OUTPUTS_BUCKET_NAME = os.getenv('JOB_OUTPUTS_BUCKET_NAME')
DATA_OBJECTS_BUCKET_NAME = os.getenv('DATA_OBJECTS_BUCKET_NAME')


@st.cache_resource()
def get_object_storage_client():
    if utils.platform() == "local":
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
    elif utils.platform() == "snowflake":
        pass


@st.cache_data()
def set_up_object_storage():
    if utils.platform() == "local":
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
            return True
        except Exception as e:
            st.error(f"Failed to set up object storage: {e}")
            return False
    elif utils.platform() == "snowflake":
        pass


def write_object_data(bucket_name, zip_id, zip_buffer):
    if utils.platform() == "local":
        try:
            client = get_object_storage_client()

            # Get the size of the zip buffer
            zip_buffer.seek(0, io.SEEK_END)
            zip_size = zip_buffer.tell()
            zip_buffer.seek(0)

            client.put_object(
                bucket_name=bucket_name,
                object_name=f"{zip_id}.zip",
                data=zip_buffer,
                length=zip_size,
                content_type='application/zip'
            )
            return True
        except Exception as e:
            st.error(f"Failed to write {zip_id}.zip to bucket {bucket_name}: {e}")
            return False
    elif utils.platform() == "snowflake":
        try:
            session = snowflake_connections.get_snowpark_session()
            zip_buffer.seek(0)
            results = session.file.put_stream(
                input_stream=zip_buffer,
                stage_location=f"@{get_user_group(get_current_username())}_group_db.{APP_NAME}_schema.{bucket_name}_stage/{zip_id}.zip",
                auto_compress=False,
            )
            return results
        except Exception as e:
            st.error(f"Failed to write {zip_id}.zip to bucket {bucket_name}: {e}")
            return None


# This could potentially be a lot of data, so we don't want to cache it using st.cache_data().
def download_object_data(bucket_name, zip_id):
    if utils.platform() == "local":
        response = None
        try:
            client = get_object_storage_client()
            object_name = f"{zip_id}.zip"
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
    elif utils.platform() == "snowflake":
        try:
            session = snowflake_connections.get_snowpark_session()
            stage_path = f"@{get_user_group(get_current_username())}_group_db.{APP_NAME}_schema.{bucket_name}_stage/{zip_id}.zip"
            bytes_io = session.file.get_stream(stage_location=stage_path)
            bytes_io.seek(0)
            return bytes_io
        except Exception as e:
            st.error(f"Failed to download {zip_id}.zip from bucket {bucket_name}: {e}")
            return None
        

def list_group_curated_object_data():
    if utils.platform() == "local":
        try:
            client = get_object_storage_client()
            objects = client.list_objects(DATA_OBJECTS_BUCKET_NAME)
            object_list = [obj.object_name for obj in objects]
            return object_list
        except Exception as e:
            st.error(f"Failed to list group curated data in {DATA_OBJECTS_BUCKET_NAME} bucket: {e}")
            return None
    elif utils.platform() == "snowflake":
        try:
            user_group = get_user_group(get_current_username())
            session = snowflake_connections.get_snowpark_session()
            stage_location = f"@{user_group}_group_db.curated_schema.{DATA_OBJECTS_BUCKET_NAME}_stage"
            files = session.file.list(stage_location=stage_location)
            object_list = [file['name'] for file in files]
            return object_list
        except Exception as e:
            st.error(f"Failed to list group curated data in {DATA_OBJECTS_BUCKET_NAME} stage: {e}")
            return None


#### 3. ORCHESTRATION FUNCTIONALITY ###############################################################


@st.cache_data()
def get_frontend_image_id():
    if utils.platform() == "local":
        try:
            resp = requests.get("http://docker_orchestrator:8080/frontend_id", timeout=3)
            resp.raise_for_status()
            data = resp.json()
            return data.get("frontend_image_id")
        except Exception as e:
            st.error(f"Could not retrieve frontend image id: {e}")
            return None
    elif utils.platform() == "snowflake":
        try:
            session = snowflake_connections.get_snowpark_session()
            return snowflake_orchestrator.frontend_id(username=get_current_username(), session=session)
        except Exception as e:
            st.error(f"Could not retrieve frontend image id: {e}")
            return None


def submit_job(job_id, blocking=True):
    if utils.platform() == "local":
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
    elif utils.platform() == "snowflake":
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
    if utils.platform() == "local":
        try:
            resp = requests.post("http://docker_orchestrator:8080/shutdown", timeout=10)
            resp.raise_for_status()
            st.success("Application is shutting down...")
            return True
        except Exception as e:
            st.error(f"Failed to shut down app: {e}")
            return False
    elif utils.platform() == "snowflake":
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
    if utils.platform() == "local":
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
    elif utils.platform() == "snowflake":
        try:
            return os.getenv("SNOWFLAKE_USER")
        except Exception as e:
            st.error(f"Failed to get current username: {e}")
            return None
