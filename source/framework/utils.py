import streamlit as st
import uuid
import datetime
import zoneinfo
import pathlib
import shutil
import pickle
import os
import zipfile
import io
from fast_neighborhood_profiles import SpatialUMAP
import importlib
import operator


ST_KEY_PREFIX_STARTUP = "startup.py__"
APP_TITLE = os.getenv("APP_TITLE")


@st.cache_data()
def _app_title_simple():
    return APP_TITLE.replace(" ", "_").lower()


# Not caching since the user could connect to the same Streamlit server by e.g. hitting refresh or opening a new tab at the same URL and would expect a new session directory.
# This function, session_dir(), and jobs_dir() below are the two places in the codebase that hardcode the local container directory to where any files are written in the app. Also, on Snowflake etc. we mount local storage at /tmp/multiplex_analysis_web_apps, so this is the isolated location where we can modify and understand these settings clearly, i.e., the only places where the app interacts with the local filesystem.
def session_dir():

    if ST_KEY_PREFIX_STARTUP + "app_session_id" not in st.session_state:
        st.error("Session ID not found in session state; cannot return the session directory.")
        return None
    
    app_session_id = st.session_state[ST_KEY_PREFIX_STARTUP + "app_session_id"]

    session_dir = f"/tmp/{_app_title_simple()}/app_session_data/{app_session_id}"

    os.makedirs(session_dir, exist_ok=True)

    return session_dir


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
        
        save_dict = {}
        saved_objects_types = {}

        def print_dict_keys(d, indent=2):
            for k, v in d.items():
                print(" " * indent + f"{k}: {type(v).__name__}", flush=True)
                if isinstance(v, dict):
                    print_dict_keys(v, indent + 2)
        
        for key, value in dictionary.items():
            if ignore_do_not_persist_flag or (not key.endswith("__do_not_persist")):
                print(f"Serializing key (before deconstruction): {key}, type: {type(value).__name__}", flush=True)
                if isinstance(value, dict):
                    # Go through each key and if the corresponding value is a dict, print its keys. Do this recursively until there are no more dicts.
                    print_dict_keys(value)
                if key == "LAZYFRAMES":
                    val_copy = {lf_key: {k: v for k, v in lf_dict.items() if k != "lf"}
                                for lf_key, lf_dict in value.items()}
                elif key.endswith("__spatial_UMAP_results"):
                    spatial_umap = value["spatial_umap"]
                    # building_blocks_keys = ["um_per_px", "dist_bin_um", "dist_bin_px", "area_downsample", "arcs_radii", "arcs_masks", "counts", "areas", "cells", "x", "img_ellipse", "w", "h", "res", "cell_positions", "cell_labels", "region_ids", "species", "density", "umap_fit", "umap_test"]
                    building_blocks_keys = ["um_per_px", "dist_bin_um", "dist_bin_px", "area_downsample", "arcs_radii", "arcs_masks", "counts", "areas", "cells", "x", "img_ellipse", "w", "h", "res", "cell_positions", "cell_labels", "region_ids", "species", "density", "umap_test"]  # removed "umap_fit" to see if that's the pickling culprit
                    # save hyperparams instead?: params = umap_fit.get_params(deep=True). Note if we specifically set the seed then this should not return the random number generator which is likely the problem, i.e., as of now returning the hyperparams alone should still error out. Yes that was the problem of the write, without umap_fit we get a successful save of the job output data.
                    building_blocks = {attr: getattr(spatial_umap, attr) for attr in building_blocks_keys if hasattr(spatial_umap, attr)}
                    # Preserve other top-level keys without deepcopy
                    val_copy = {k: v for k, v in value.items() if k != "spatial_umap"}
                    val_copy["spatial_umap"] = building_blocks
                else:
                    val_copy = value
                print(f"Serializing key (after deconstruction): {key}, type: {type(val_copy).__name__}", flush=True)
                if isinstance(val_copy, dict):
                    print_dict_keys(val_copy)
                save_dict[key] = val_copy
                saved_objects_types[key] = type(val_copy).__name__

        pkl_file = os.path.join(directory, f'{dict_name}.pkl')
        with open(pkl_file, 'wb') as f:
            f.write(pickle.dumps(save_dict))

        return saved_objects_types
    except Exception as e:
        st.error(f"Failed to serialize dictionary {dict_name} to directory {directory}: {e}")
        return None


def deserialize_binary_files_to_dictionary(dict_name, directory, dictionary=None):
    try:
        if dictionary is None:
            dictionary = {}

        pkl_file = os.path.join(directory, f'{dict_name}.pkl')
        if os.path.exists(pkl_file):
            with open(pkl_file, 'rb') as f:
                dictionary.update(pickle.loads(f.read()))

        for key in dictionary:
            if key == "LAZYFRAMES":
                for lf_key in dictionary["LAZYFRAMES"]:

                    function_metadata = dictionary["LAZYFRAMES"][lf_key]["function_metadata"]
                    mod = importlib.import_module(function_metadata["module_name"])
                    function = operator.attrgetter(function_metadata["qualpath"])(mod)

                    input_dataset = dictionary["LAZYFRAMES"][lf_key]["input_dataset"]
                    params = dictionary["LAZYFRAMES"][lf_key]["params"]
                    if input_dataset is None:
                        result = function(**params)
                    elif input_dataset["type"] == "lf":
                        lf = dictionary["LAZYFRAMES"][input_dataset["keys"][0]]["lf"]
                        result = function(lf, **params)
                    elif input_dataset["type"] == "pandas_df":
                        pd_df = getattr(dictionary[input_dataset["keys"][0]][input_dataset["keys"][1]], input_dataset["keys"][2])  # Modify in the future; this is really specific to the format of sumap.cells on the run_spatial_umap.py page.
                        result = function(pd_df, **params)
                    if isinstance(result, tuple):
                        dictionary["LAZYFRAMES"][lf_key]["lf"] = result[0]
                        dictionary["LAZYFRAMES"][lf_key]["extras"] = result[1]
                    else:
                        dictionary["LAZYFRAMES"][lf_key]["lf"] = result
                        dictionary["LAZYFRAMES"][lf_key]["extras"] = None
            elif key.endswith("__spatial_UMAP_results"):
                print("AAAAAAAAAAAAAAAAAAAAAA")
                bb = dictionary[key]["spatial_umap"]
                print("before reconstruction:")
                print(dictionary[key]["spatial_umap"], flush=True)
                spatial_umap = SpatialUMAP.SpatialUMAP(dist_bin_um=bb["dist_bin_um"], um_per_px=bb["um_per_px"], area_downsample=bb["area_downsample"])
                for attr_key, attr_value in bb.items():
                    setattr(spatial_umap, attr_key, attr_value)
                dictionary[key]["spatial_umap"] = spatial_umap
                print("after reconstruction:")
                print(dictionary[key]["spatial_umap"], flush=True)
                print(dictionary[key]["spatial_umap"].area_downsample, flush=True)
                print("BBBBBBBBBBBBBBBBBBBB")

        return dictionary
    except Exception as e:
        st.error(f"Failed to deserialize binary file {dict_name}.pkl to dictionary in directory {directory}: {e}")
        return None


def zip_directory_to_buffer(directory, compresslevel=6, ignore_subdirs=[]):
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

                # Skip any ignored subdirectories (and their contents)
                if ignore_subdirs and any(rel == sub or rel.startswith(f"{sub}/") for sub in ignore_subdirs):
                    continue

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
