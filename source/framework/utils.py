import streamlit as st
import uuid
import datetime
import zoneinfo
import pathlib
import shutil
import pickle
import dill
import os
import zipfile
import io

ST_KEY_PREFIX_STARTUP = "startup.py__"
APP_TITLE = os.getenv("APP_TITLE")


@st.cache_data()
def _app_title_simple():
    return APP_TITLE.replace(" ", "_").lower()


# Not caching since the user could connect to the same Streamlit server by e.g. hitting refresh or opening a new tab at the same URL and would expect a new session directory.
# This function, session_dir(), and jobs_dir() below are the two places in the codebase that hardcode the local container directory to where any files are written in the app. Also, on Snowflake etc. we mount local storage at /tmp/multiplex_analysis_web_apps, so this is the isolated location where we can modify and understand these settings clearly, i.e., the only places where the app interacts with the local filesystem.
def session_dir():
    try:
        if ST_KEY_PREFIX_STARTUP + "app_session_id" not in st.session_state:
            # Worker environment fallback - use /tmp directory 
            worker_fallback_dir = f"/tmp/{_app_title_simple()}/worker_session_data"
            os.makedirs(worker_fallback_dir, exist_ok=True)
            # Also ensure output subdirectory exists for benchmark_collector
            output_dir = os.path.join(worker_fallback_dir, 'output')
            os.makedirs(output_dir, exist_ok=True)
            return worker_fallback_dir
        
        app_session_id = st.session_state[ST_KEY_PREFIX_STARTUP + "app_session_id"]

        session_dir = f"/tmp/{_app_title_simple()}/app_session_data/{app_session_id}"

        os.makedirs(session_dir, exist_ok=True)
        # Also ensure output subdirectory exists for benchmark_collector
        output_dir = os.path.join(session_dir, 'output')
        os.makedirs(output_dir, exist_ok=True)

        return session_dir
    except:
        # Ultimate fallback for worker environments where streamlit is not available
        worker_fallback_dir = f"/tmp/{_app_title_simple()}/worker_session_data"
        os.makedirs(worker_fallback_dir, exist_ok=True)
        # Also ensure output subdirectory exists for benchmark_collector
        output_dir = os.path.join(worker_fallback_dir, 'output')
        os.makedirs(output_dir, exist_ok=True)
        return worker_fallback_dir


@st.cache_data()
def jobs_dir():
    """Root directory for per-job temporary data.
    """
    jobs_dir = f"/tmp/{_app_title_simple()}/job_data"
    os.makedirs(jobs_dir, exist_ok=True)
    return jobs_dir


def get_unique_id():
    return uuid.uuid4().hex


def get_timestamp(timezone="America/New_York"):
    return datetime.datetime.now(zoneinfo.ZoneInfo(timezone))


def ensure_empty_directory(directory, create_if_missing=True):
    """Ensure a directory exists and is empty, creating it if necessary."""
    try:
        directory_path = pathlib.Path(directory)
        if directory_path.exists():
            shutil.rmtree(directory_path)

        if create_if_missing:
            directory_path.mkdir(parents=True, exist_ok=True)

        return True
    except Exception as e:
        st.error(f"Failed to ensure empty directory {directory}: {e}")
        return False


def serialize_dictionary_to_binary_files(dictionary, dict_name, directory, ignore_do_not_persist_flag=True):
    """Save the entire session state efficiently to the session directory"""
    try:
        serializable_dict = {}
        serializable_objects = {}
        unserializable_dict = {}
        unserializable_objects = {}

        for key, value in dictionary.items():
            if ignore_do_not_persist_flag or (not key.endswith("__do_not_persist")):
                try:
                    pickle.dumps(value)  # Test pickling
                    serializable_dict[key] = value
                    serializable_objects[key] = type(value).__name__
                except (TypeError, AttributeError, pickle.PicklingError):
                    unserializable_dict[key] = value
                    unserializable_objects[key] = type(value).__name__

        pkl_file = os.path.join(directory, f'{dict_name}.pkl')
        with open(pkl_file, 'wb') as f:
            f.write(pickle.dumps(serializable_dict))

        dill_file = os.path.join(directory, f'{dict_name}.dill')
        with open(dill_file, 'wb') as f:
            f.write(dill.dumps(unserializable_dict))

        return serializable_objects, unserializable_objects
    except Exception as e:
        st.error(f"Failed to serialize dictionary {dict_name} to directory {directory}: {e}")
        return None


def deserialize_binary_files_to_dictionary(dict_name, directory, dictionary=None):
    try:
        print(f"DEBUG: deserialize_binary_files_to_dictionary called with dict_name={dict_name}, directory={directory}", flush=True)
        
        if dictionary is None:
            dictionary = {}

        pkl_file = os.path.join(directory, f'{dict_name}.pkl')
        dill_file = os.path.join(directory, f'{dict_name}.dill')
        
        print(f"DEBUG: Looking for pkl_file: {pkl_file}, exists: {os.path.exists(pkl_file)}", flush=True)
        print(f"DEBUG: Looking for dill_file: {dill_file}, exists: {os.path.exists(dill_file)}", flush=True)
        
        if os.path.exists(pkl_file):
            print(f"DEBUG: Loading pickle file: {pkl_file}", flush=True)
            with open(pkl_file, 'rb') as f:
                pkl_data = pickle.loads(f.read())
                print(f"DEBUG: Loaded pickle data, type: {type(pkl_data)}, keys: {list(pkl_data.keys()) if isinstance(pkl_data, dict) else 'not a dict'}", flush=True)
                dictionary.update(pkl_data)

        if os.path.exists(dill_file):
            print(f"DEBUG: Loading dill file: {dill_file}", flush=True)
            with open(dill_file, 'rb') as f:
                dill_data = dill.loads(f.read())
                print(f"DEBUG: Loaded dill data, type: {type(dill_data)}, keys: {list(dill_data.keys()) if isinstance(dill_data, dict) else 'not a dict'}", flush=True)
                dictionary.update(dill_data)

        print(f"DEBUG: Final dictionary type: {type(dictionary)}, keys: {list(dictionary.keys()) if isinstance(dictionary, dict) else 'not a dict'}", flush=True)
        return dictionary
    except Exception as e:
        print(f"ERROR: Failed to deserialize binary files {dict_name}.pkl/.dill to dictionary in directory {directory}: {e}", flush=True)
        import traceback
        print(f"ERROR: Traceback: {traceback.format_exc()}", flush=True)
        st.error(f"Failed to deserialize binary files {dict_name}.pkl/.dill to dictionary in directory {directory}: {e}")
        return None


def zip_directory_to_buffer(directory, compresslevel=6):
    """Create a zip archive of a directory (includes hidden files and empty dirs; no symlink handling)."""
    try:
        main_path = pathlib.Path(directory)
        if not main_path.is_dir():
            st.error(f"Not a directory: {directory}")
            return None
        zip_buffer = io.BytesIO()
        with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED, compresslevel=compresslevel) as zip_file:
            all_paths = sorted(main_path.rglob('*'), key=lambda p: p.as_posix())
            for file_path in all_paths:
                rel = file_path.relative_to(main_path).as_posix()

                if file_path.is_dir():
                    if rel:
                        st_mode = file_path.stat().st_mode
                        mtime = file_path.stat().st_mtime
                        dt = datetime.datetime.fromtimestamp(mtime)
                        zinfo = zipfile.ZipInfo(
                            rel.rstrip('/') + '/',
                            date_time=(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)
                        )
                        zinfo.external_attr = st_mode << 16
                        zip_file.writestr(zinfo, '')
                    continue

                try:
                    zip_file.write(file_path, rel)  # retains original mtime automatically
                    zip_file.getinfo(rel).external_attr = file_path.stat().st_mode << 16
                except (OSError, IOError) as e:
                    st.warning(f"Skipping file {file_path}: {e}")

        zip_buffer.seek(0)
        return zip_buffer
    except Exception as e:
        st.error(f"Failed to zip directory {directory}: {e}")
        return None


def unzip_buffer_to_directory(zip_buffer, directory):
    try:
        if zip_buffer:
            with zipfile.ZipFile(zip_buffer, 'r') as zip_file:
                zip_file.extractall(directory)
        return True
    except Exception as e:
        st.error(f"Failed to unzip buffer to {directory}: {e}")
        return False


@st.cache_data()
def platform():
    app_platform = os.getenv("APP_PLATFORM")
    if app_platform in ("local", "snowflake"):
        return app_platform
    else:
        st.warning(f"APP_PLATFORM environment variable is not set to a valid value ('local' or 'snowflake'). Detected value: {app_platform}. Falling back to automatic detection.")
    if os.getenv("SNOWFLAKE_ACCOUNT") is not None:
        return "snowflake"
    else:
        return "local"
