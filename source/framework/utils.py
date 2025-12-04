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
        multiprint(f"Failed to ensure empty directory {directory}: {e}", (print,))
        raise


def multiprint(message, functions):
    for function in functions:
        if function == print:
            function(message, flush=True)
        else:
            function(message)


def deconstruct_object(path, value, value_type):
    try:
        if value_type == "LazyFrame":
            multiprint("Deconstructing LazyFrame object.", (print,))
            return {"object_type": "Deconstructed LazyFrame"}
        elif value_type == "SpatialUMAP":
            multiprint("Deconstructing SpatialUMAP object.", (print,))
            spatial_umap = value
            component_keys = ["um_per_px", "dist_bin_um", "dist_bin_px", "area_downsample", "arcs_radii", "arcs_masks", "counts", "areas", "cells", "x", "img_ellipse", "w", "h", "res", "cell_positions", "cell_labels", "region_ids", "species", "density", "umap_test"]  # removed "umap_fit" to see if that's the pickling culprit
            # save hyperparams too?: params = umap_fit.get_params(deep=True). Note if we specifically set the seed then this should not return the random number generator which is likely the problem, i.e., as of now returning the hyperparams alone should still error out. Yes that was the problem of the write, without umap_fit we get a successful save of the job output data.
            components = {attr: getattr(spatial_umap, attr) for attr in component_keys if hasattr(spatial_umap, attr)}
            return {"object_type": "Deconstructed SpatialUMAP", "components": components}
        elif value_type == "foundry_IO_lib":
            multiprint("Deconstructing foundry_IO_lib object.", (print,))
            fiol = value
            component_keys = ["onNIDAP"]
            components = {attr: getattr(fiol, attr) for attr in component_keys if hasattr(fiol, attr)}
            return {"object_type": "Deconstructed foundry_IO_lib", "components": components}
        elif value_type == "benchmark_collector":
            multiprint("Deconstructing benchmark_collector object.", (print,))
            bc = value
            component_keys = ["benchmarkDF", "on_nidap", "benchmark_csv", "benchmark_project_path", "benchmark_dataset"]
            components = {attr: getattr(bc, attr) for attr in component_keys if hasattr(bc, attr)}
            return {"object_type": "Deconstructed benchmark_collector", "components": components}
        elif value_type == "Platform":
            multiprint("Deconstructing Platform object. Nothing is actually being done since this will always get overwritten by the app that loads in the session state; see manage_sessions.load_session_state().", (print,))
            return {"object_type": "Deconstructed Platform"}
        elif value_type == "Standardized":
            multiprint("Deconstructing dataset_formats.Standardized object.", (print,))
            standardized = value
            component_keys = ["images_to_analyze", "phenotypes_to_analyze", "input_datafile", "sep", "data", "coord_units_in_microns", "min_coord_spacing_", "species_equivalents", "mapping_dict", "roi_width", "overlap", "phenotype_identification_tsv_file", "extra_cols_to_keep"]
            components = {attr: getattr(standardized, attr) for attr in component_keys if hasattr(standardized, attr)}
            return {"object_type": "Deconstructed Standardized", "components": components}
        else:
            return value
    except Exception as e:
        multiprint(f"Failed to deconstruct object in dictionary path {".".join(str(x) for x in path)} of type {value_type}: {e}", (print,))
        raise


def reconstruct_object(value, value_type):
    try:
        if value_type == "LazyFrame":
            multiprint("Skipping reconstruction of LazyFrame object in utils.reconstruct_object() since we do it afterward for all lazyframes at once in utils.deserialize_binary_files_to_dictionary().", (print,))
            return value
        elif value_type == "SpatialUMAP":
            multiprint("Reconstructing SpatialUMAP object.", (print,))
            components = value["components"]
            spatial_umap = SpatialUMAP.SpatialUMAP(dist_bin_um=components["dist_bin_um"], um_per_px=components["um_per_px"], area_downsample=components["area_downsample"])
            for attr_key, attr_value in components.items():
                setattr(spatial_umap, attr_key, attr_value)
            return spatial_umap  # Return the spatial UMAP object.
        elif value_type == "foundry_IO_lib":
            multiprint("Reconstructing foundry_IO_lib object.", (print,))
            components = value["components"]
            fiol = foundry_IO_lib.foundry_IO_lib()
            for attr_key, attr_value in components.items():
                setattr(fiol, attr_key, attr_value)
            return fiol  # Return the foundry_IO_lib object.
        elif value_type == "benchmark_collector":
            multiprint("Reconstructing benchmark_collector object.", (print,))
            components = value["components"]
            # Maybe it's good principle to not have orig_dict in here at all because there's no telling when various needed pieces may be updated. If it's needed it could instead be a sign of a poorly designed class.
            # if "fiol" in orig_dict and not isinstance(orig_dict["fiol"], dict):
            #     bc = benchmark_collector.benchmark_collector(orig_dict["fiol"])
            # else:
            bc = benchmark_collector.benchmark_collector()
            for attr_key, attr_value in components.items():
                setattr(bc, attr_key, attr_value)
            return bc  # Return the benchmark_collector object.
        elif value_type == "Platform":
            multiprint("Reconstructing Platform object. Nothing is actually being done since this will always get overwritten by the app that loads in the session state; see manage_sessions.load_session_state().", (print,))
            return None
        else:  # Functionality for this branch *should* be different than in deconstruct_object(). Overall, whether deconstructing or reconstructing, we should return a new object or the original one.
            raise ValueError(f"Unknown object type for reconstruction: {value_type}")
    except Exception as e:
        multiprint(f"Failed to reconstruct object of type {value_type}: {e}", (print,))
        raise


def traverse_for_deconstruct(path, value, debug=False):
    """Recursively deconstruct objects, always building new containers."""
    try:
        key_str = ".".join(str(p) for p in path)
        value_type = type(value).__name__
        if debug:
            multiprint(f"{key_str}: {value_type}", (print,))

        try:
            if isinstance(value, dict):
                new_dict = {}
                for sub_key, sub_value in value.items():
                    new_dict[sub_key] = traverse_for_deconstruct(path + [sub_key], sub_value, debug=debug)
                return new_dict
            elif isinstance(value, list):
                return [traverse_for_deconstruct(path + [idx], sub_value, debug=debug) for idx, sub_value in enumerate(value)]
            elif isinstance(value, tuple):
                return tuple(traverse_for_deconstruct(path + [idx], sub_value, debug=debug) for idx, sub_value in enumerate(value))
            else:
                return deconstruct_object(path=path, value=value, value_type=value_type)
        except RecursionError as e:
            multiprint(f"{key_str}: RecursionError encountered (possible cyclic reference): {e}", (print,))
            return value
    except Exception as e:
        multiprint(f"Failed to traverse for deconstruction at {key_str}: {e}", (print,))
        raise


def traverse_for_reconstruct(path, value, debug=False):
    try:
        key_str = ".".join(str(p) for p in path)
        value_type = type(value).__name__
        if debug:
            multiprint(f"{key_str}: {value_type}", (print,))

        try:
            # Dict branch
            if isinstance(value, dict) and not ("object_type" in value and value["object_type"].startswith("Deconstructed ")):
                new_dict = None  # defer allocation
                changed = False
                for sub_key, sub_value in value.items():
                    result = traverse_for_reconstruct(path + [sub_key], sub_value, debug=debug)
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
                    result = traverse_for_reconstruct(path + [idx], sub_value, debug=debug)
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
                    result = traverse_for_reconstruct(path + [idx], sub_value, debug=debug)
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
                    return reconstruct_object(value=value, value_type=value_type)
                return value

        except RecursionError as e:
            multiprint(f"{key_str}: RecursionError encountered (possible cyclic reference): {e}", (print,))
            return value
    except Exception as e:
        multiprint(f"Failed to traverse for reconstruction at {key_str}: {e}", (print,))
        raise


def serialize_dictionary_to_binary_files(dictionary, dict_name, directory, ignore_do_not_persist_flag=True):
    """Save the entire session state efficiently to the session directory"""
    try:
        transformed_dict = {}
        for key, value in dictionary.items():
            if ignore_do_not_persist_flag or (not key.endswith("__do_not_persist")):
                transformed_dict[key] = traverse_for_deconstruct([key], value)

        pkl_file = os.path.join(directory, f'{dict_name}.pkl')
        with open(pkl_file, 'wb') as f:
            f.write(pickle.dumps(transformed_dict))
    except Exception as e:
        multiprint(f"Failed to serialize dictionary {dict_name} to directory {directory}: {e}", (print,))
        raise


def deserialize_binary_files_to_dictionary(dict_name, directory, dictionary=None, extra_dict_to_load=None, topdir_for_lazyframe_data=None):
    try:
        if dictionary is None:
            dictionary = {}

        pkl_file = os.path.join(directory, f'{dict_name}.pkl')
        if os.path.exists(pkl_file):
            with open(pkl_file, 'rb') as f:
                dictionary.update(pickle.loads(f.read()))

        # Since old archives may store the old app session ID which were read in just above, but below the reconstruction may require the current app session ID, we update the dictionary so that reconstruction uses the current app session ID. Implementing this since benchmark collector uses session_dir() during its reconstruction which depends on the app session ID in the session state. Correspondingly commenting out this dictionary update in manage_sessions.load_session_state() since that's now done here.
        if extra_dict_to_load:
            dictionary.update(extra_dict_to_load)

        for key, value in dictionary.items():
            dictionary[key] = traverse_for_reconstruct([key], value)

        # We should always ensure that reconstructing lazyframes is a fast process; no long computations should be done here! E.g., the spatial UMAP lazyframe should be quickly computed from the spatial_umap.cells dataframe already present in the dictionary. The calculation of spatial_umap.cells could be a long operation.
        # Note also that lazyframe construction was assumed to depend on values in the session state (such as spatial_umap.cells) so rebuilding lazyframes after loading in the session state makes sense.
        if "LAZYFRAMES" in dictionary:

            # Ensure the location for storing the files from which the lazyframes could read exists.
            os.makedirs(os.path.join(topdir_for_lazyframe_data, "input"), exist_ok=True)

            # Build dependency order
            lf_meta = dictionary["LAZYFRAMES"]
            deps = {}
            for k, v in lf_meta.items():
                d = v.get("input_dataset")
                if d and d.get("type") == "lf":
                    deps[k] = {d["keys"][0]}
                else:
                    deps[k] = set()
            
            # Topological sort
            visited = set()
            order = []
            def visit(n):
                if n in visited:
                    return
                for prereq in deps.get(n, ()):
                    if prereq in lf_meta:
                        visit(prereq)
                visited.add(n)
                order.append(n)
            
            for lf_key in lf_meta.keys():
                visit(lf_key)
            
            # Reconstruct in order
            for lf_key in order:
                multiprint(f"Reconstructing LazyFrame: {lf_key}.", (print,))
                function_metadata = dictionary["LAZYFRAMES"][lf_key]["function_metadata"]
                mod = importlib.import_module(function_metadata["module_name"])
                function = operator.attrgetter(function_metadata["qualpath"])(mod)

                input_dataset = dictionary["LAZYFRAMES"][lf_key]["input_dataset"]
                params = dictionary["LAZYFRAMES"][lf_key]["params"]
                if input_dataset is None:
                    lf_out = function(**params, topdir=topdir_for_lazyframe_data)  # inject topdir_for_lazyframe_data here since there's no input so we must be downloading data (where we'll need to **store** hence the need) from which to create a lazyframe
                elif input_dataset["type"] == "lf":
                    lf = dictionary["LAZYFRAMES"][input_dataset["keys"][0]]["lf"]
                    lf_out = function(lf, **params)  # no need to inject topdir_for_lazyframe_data here, for now, since we're likely just transforming an already existing lazyframe that should already be set up to correctly reference an on-disk file
                elif input_dataset["type"] == "pandas_df":
                    pd_df = getattr(dictionary[input_dataset["keys"][0]][input_dataset["keys"][1]], input_dataset["keys"][2])  # Modify in the future; this is really specific to the format of sumap.cells on the run_spatial_umap.py page.
                    lf_out = function(pd_df, **params, topdir=topdir_for_lazyframe_data)  # inject topdir_for_lazyframe_data here since we're likely creating a file on disk (hence the need) from a pandas dataframe

                dictionary["LAZYFRAMES"][lf_key]["lf"] = lf_out

        return dictionary
    except Exception as e:
        multiprint(f"Failed to deserialize binary file {dict_name}.pkl to dictionary in directory {directory}: {e}", (print,))
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
                    multiprint(f"Skipping file {file_path}: {e}", (print,))
                    raise

        zip_buffer.seek(0)
        return zip_buffer
    except Exception as e:
        multiprint(f"Failed to zip directory {directory}: {e}", (print,))
        raise


def unzip_buffer_to_directory(zip_buffer, directory):
    try:
        if zip_buffer:
            with zipfile.ZipFile(zip_buffer, 'r') as zip_file:
                zip_file.extractall(directory)
        return True
    except Exception as e:
        multiprint(f"Failed to unzip buffer to {directory}: {e}", (print,))
        raise


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
