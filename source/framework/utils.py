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
import foundry_IO_lib
import benchmark_collector


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


def multiprint(message, functions):
    for function in functions:
        if function == print:
            function(message, flush=True)
        else:
            function(message)


def deconstruct_object(identifier, value, value_type):
    try:
        if value_type == "LazyFrame":
            lf_key = identifier
            return {"object_type": "Deconstructed LazyFrame", "lf_key": lf_key}  # Assuming the key in the dictionary is the last element in the path.
        elif value_type == "SpatialUMAP":
            spatial_umap = value
            component_keys = ["um_per_px", "dist_bin_um", "dist_bin_px", "area_downsample", "arcs_radii", "arcs_masks", "counts", "areas", "cells", "x", "img_ellipse", "w", "h", "res", "cell_positions", "cell_labels", "region_ids", "species", "density", "umap_test"]  # removed "umap_fit" to see if that's the pickling culprit
            # save hyperparams too?: params = umap_fit.get_params(deep=True). Note if we specifically set the seed then this should not return the random number generator which is likely the problem, i.e., as of now returning the hyperparams alone should still error out. Yes that was the problem of the write, without umap_fit we get a successful save of the job output data.
            components = {attr: getattr(spatial_umap, attr) for attr in component_keys if hasattr(spatial_umap, attr)}
            return {"object_type": "Deconstructed SpatialUMAP", "components": components}
        elif value_type == "foundry_IO_lib":
            fiol = value
            component_keys = ["onNIDAP"]
            components = {attr: getattr(fiol, attr) for attr in component_keys if hasattr(fiol, attr)}
            return {"object_type": "Deconstructed foundry_IO_lib", "components": components}
        elif value_type == "benchmark_collector":
            bc = value
            component_keys = ["benchmarkDF", "on_nidap", "benchmark_csv", "benchmark_project_path", "benchmark_dataset"]
            components = {attr: getattr(bc, attr) for attr in component_keys if hasattr(bc, attr)}
            return {"object_type": "Deconstructed benchmark_collector", "components": components}
        else:
            return value
    except Exception as e:
        multiprint(f"Failed to deconstruct object {identifier} of type {value_type}: {e}", (print, st.error))
        raise


def reconstruct_object(value, value_type, orig_dict):
    try:
        if value_type == "LazyFrame":
            lf_key = value["lf_key"]

            function_metadata = orig_dict["LAZYFRAMES"][lf_key]["function_metadata"]
            mod = importlib.import_module(function_metadata["module_name"])
            function = operator.attrgetter(function_metadata["qualpath"])(mod)

            input_dataset = orig_dict["LAZYFRAMES"][lf_key]["input_dataset"]
            params = orig_dict["LAZYFRAMES"][lf_key]["params"]
            if input_dataset is None:
                result = function(**params)
            elif input_dataset["type"] == "lf":
                lf = orig_dict["LAZYFRAMES"][input_dataset["keys"][0]]["lf"]
                result = function(lf, **params)
            elif input_dataset["type"] == "pandas_df":
                pd_df = getattr(orig_dict[input_dataset["keys"][0]][input_dataset["keys"][1]], input_dataset["keys"][2])  # Modify in the future; this is really specific to the format of sumap.cells on the run_spatial_umap.py page.
                result = function(pd_df, **params)
            return result  # Return the lazyframe.
        elif value_type == "SpatialUMAP":
            components = value["components"]
            spatial_umap = SpatialUMAP.SpatialUMAP(dist_bin_um=components["dist_bin_um"], um_per_px=components["um_per_px"], area_downsample=components["area_downsample"])
            for attr_key, attr_value in components.items():
                setattr(spatial_umap, attr_key, attr_value)
            return spatial_umap  # Return the spatial UMAP object.
        elif value_type == "foundry_IO_lib":
            components = value["components"]
            fiol = foundry_IO_lib.foundry_IO_lib()
            for attr_key, attr_value in components.items():
                setattr(fiol, attr_key, attr_value)
            return fiol  # Return the foundry_IO_lib object.
        elif value_type == "benchmark_collector":
            components = value["components"]
            if "fiol" in orig_dict:
                bc = benchmark_collector.benchmark_collector(orig_dict["fiol"])
            else:
                bc = benchmark_collector.benchmark_collector()
            for attr_key, attr_value in components.items():
                setattr(bc, attr_key, attr_value)
            return bc  # Return the benchmark_collector object.
        else:  # Functionality for this branch *should* be different than in deconstruct_object(). Overall, whether deconstructing or reconstructing, we should return a new object or the original one.
            raise ValueError(f"Unknown object type for reconstruction: {value_type}")
    except Exception as e:
        multiprint(f"Failed to reconstruct object of type {value_type}: {e}", (print, st.error))
        raise


def traverse_for_deconstruct(path, value, orig_dict, debug=False):
    """Recursively deconstruct objects, always building new containers."""
    try:
        key_str = ".".join(str(p) for p in path)
        value_type = type(value).__name__
        if debug:
            print(f"{key_str}: {value_type}", flush=True)

        try:
            if isinstance(value, dict):
                new_dict = {}
                for sub_key, sub_value in value.items():
                    new_dict[sub_key] = traverse_for_deconstruct(path + [sub_key], sub_value, orig_dict)
                return new_dict
            elif isinstance(value, list):
                return [traverse_for_deconstruct(path + [idx], sub_value, orig_dict) for idx, sub_value in enumerate(value)]
            elif isinstance(value, tuple):
                return tuple(traverse_for_deconstruct(path + [idx], sub_value, orig_dict) for idx, sub_value in enumerate(value))
            else:
                return deconstruct_object(identifier=path[-1], value=value, value_type=value_type)
        except RecursionError as e:
            print(f"{key_str}: RecursionError encountered (possible cyclic reference): {e}", flush=True)
            return value
    except Exception as e:
        multiprint(f"Failed to traverse for deconstruction at {key_str}: {e}", (print, st.error))
        raise


def traverse_for_reconstruct(path, value, orig_dict, debug=False):
    try:
        key_str = ".".join(str(p) for p in path)
        value_type = type(value).__name__
        if debug:
            print(f"{key_str}: {value_type}", flush=True)

        try:
            # Dict branch
            if isinstance(value, dict) and not ("object_type" in value and value["object_type"].startswith("Deconstructed ")):
                new_dict = None  # defer allocation
                changed = False
                for sub_key, sub_value in value.items():
                    result = traverse_for_reconstruct(path + [sub_key], sub_value, orig_dict)
                    if result is not sub_value and not changed:
                        # First change detected: materialize new_dict with prior unchanged items
                        new_dict = {k: (value[k] if k != sub_key else result)
                                    for k in value.keys()}
                        changed = True
                    elif changed:
                        new_dict[sub_key] = result
                return new_dict if changed else value

            # List branch
            elif isinstance(value, list):
                new_list = None
                changed = False
                for idx, sub_value in enumerate(value):
                    result = traverse_for_reconstruct(path + [idx], sub_value, orig_dict)
                    if result is not sub_value and not changed:
                        new_list = value[:]
                        new_list[idx] = result
                        changed = True
                    elif changed:
                        new_list[idx] = result
                return new_list if changed else value

            # Tuple branch
            elif isinstance(value, tuple):
                # We must build a list of results only if a change occurs
                temp = None
                changed = False
                for idx, sub_value in enumerate(value):
                    result = traverse_for_reconstruct(path + [idx], sub_value, orig_dict)
                    if result is not sub_value and not changed:
                        temp = list(value)
                        temp[idx] = result
                        changed = True
                    elif changed:
                        temp[idx] = result
                return tuple(temp) if changed else value

            # Leaf / deconstructed placeholder
            else:
                if isinstance(value, dict) and ("object_type" in value and value["object_type"].startswith("Deconstructed ")):
                    value_type = value["object_type"].removeprefix("Deconstructed ")
                    return reconstruct_object(value=value, value_type=value_type, orig_dict=orig_dict)
                return value

        except RecursionError as e:
            multiprint(f"{key_str}: RecursionError encountered (possible cyclic reference): {e}", (print, st.error))
            raise
    except Exception as e:
        multiprint(f"Failed to traverse for reconstruction at {key_str}: {e}", (print, st.error))
        raise


def serialize_dictionary_to_binary_files(dictionary, dict_name, directory, ignore_do_not_persist_flag=True):
    """Save the entire session state efficiently to the session directory"""
    try:
        transformed_dict = {}
        for key, value in dictionary.items():
            if ignore_do_not_persist_flag or (not key.endswith("__do_not_persist")):
                transformed_dict[key] = traverse_for_deconstruct([key], value, dictionary)

        pkl_file = os.path.join(directory, f'{dict_name}.pkl')
        with open(pkl_file, 'wb') as f:
            f.write(pickle.dumps(transformed_dict))
    except Exception as e:
        multiprint(f"Failed to serialize dictionary {dict_name} to directory {directory}: {e}", (print, st.error))
        raise


def deserialize_binary_files_to_dictionary(dict_name, directory, dictionary=None):
    try:
        if dictionary is None:
            dictionary = {}

        pkl_file = os.path.join(directory, f'{dict_name}.pkl')
        if os.path.exists(pkl_file):
            with open(pkl_file, 'rb') as f:
                dictionary.update(pickle.loads(f.read()))

        for key, value in dictionary.items():
            dictionary[key] = traverse_for_reconstruct([key], value, dictionary)

        return dictionary
    except Exception as e:
        multiprint(f"Failed to deserialize binary file {dict_name}.pkl to dictionary in directory {directory}: {e}", (print, st.error))
        raise


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
