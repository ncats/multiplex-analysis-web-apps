# Import relevant libraries.
import streamlit as st
import subprocess
import os
import pathlib
import copy
import framework.utils as framework_utils
import framework.platform_abstraction as pa
import importlib

ST_KEY_PREFIX = "manage_sessions.py__"
ST_KEY_PREFIX_STARTUP = "startup.py__"
ST_KEY_PREFIX_APP = "app.py__"
ARCHIVES_BUCKET_NAME = os.getenv('ARCHIVES_BUCKET_NAME')


@st.cache_data()
def export_conda_environment():
    """Export the current conda environment and return the YAML content."""
    try:
        # Probably run this during startup and store the result (in a file) so no calls to subprocess are needed in the whole app.
        return subprocess.run(['micromamba', 'env', 'export'], capture_output=True, text=True, check=True).stdout
    except Exception as e:
        st.error(f"Failed to export conda environment: {e}")
        return None


def write_conda_environment(filename, directory):
    """Write the conda environment content to a YAML file in the session directory."""
    try:
        conda_env_content = export_conda_environment()
        yaml_filename = os.path.join(directory, filename)
        with open(yaml_filename, 'w') as f:
            f.write(conda_env_content)
        return True
    except Exception as e:
        st.error(f"Failed to write conda environment file {filename}: {e}")
        return False


def save_session_state():
    try:
        session_state_directory = framework_utils.session_dir()
        serializable_objects = framework_utils.serialize_dictionary_to_binary_files(st.session_state, "session_state", session_state_directory, ignore_do_not_persist_flag=False)

        # Write out what was serialized and how.
        info_file = os.path.join(session_state_directory, f'session_state_contents.txt')
        with open(info_file, 'w') as f:
            f.write("Serializable objects:\n")
            for key, value in serializable_objects.items():
                f.write(f"{key}: {value}\n")

        return True
    except Exception as e:
        st.error(f"Failed to save session state: {e}")
        return False


def load_session_state():
    """Load the session state from the session directory."""
    try:
        # Must save session_dir and not use utils.session_dir() in utils.deserialize... since the session_dir() depends on the session state, which would have just been cleared.
        session_dir = framework_utils.session_dir()

        # Back up app session-specific (i.e., startup.py-defined) variables we ultimately don't want to overwrite.
        keys_to_keep = [ST_KEY_PREFIX_STARTUP + "app_session_id", "previous_page_name", "current_page_name", ST_KEY_PREFIX_APP + "app_initialized", "platform"]
        startup_keys = {key: copy.deepcopy(st.session_state[key]) for key in keys_to_keep if key in st.session_state}

        # Delete everything in the session state.
        for key in list(st.session_state.keys()):
            del st.session_state[key]

        framework_utils.deserialize_binary_files_to_dictionary("session_state", session_dir, dictionary=st.session_state)

        # Restore the startup keys.
        st.session_state.update(startup_keys)

        return True
    except Exception as e:
        st.error(f"Failed to load session state: {e}")
        return False


def reset_session_state(extra_keys_to_keep=[], delete_input_dir=True):
    """Load the session state from the session directory."""
    try:
        # Back up app session-specific (i.e., startup.py-defined) variables we ultimately don't want to overwrite.
        keys_to_keep = [ST_KEY_PREFIX_STARTUP + "app_session_id", "previous_page_name", "current_page_name", ST_KEY_PREFIX_APP + "app_initialized", "platform"]

        # Delete everything in the session state but the keys to keep.
        for key in list(st.session_state.keys()):
            if (key not in keys_to_keep) and (key not in extra_keys_to_keep):
                del st.session_state[key]

        # Delete everything from the input and output directories.
        if delete_input_dir:
            framework_utils.ensure_empty_directory(os.path.join(framework_utils.session_dir(), "input"))
        framework_utils.ensure_empty_directory(os.path.join(framework_utils.session_dir(), "output"))

        return True
    except Exception as e:
        st.error(f"Failed to reset session state: {e}")
        return False


@st.cache_data()
def get_current_git_commit():
    try:
        # Probably run this during startup and store the result (in an environment variable) so no calls to subprocess are needed in the whole app.
        return subprocess.run(['git', '-c', 'safe.directory=/app', 'rev-parse', 'HEAD'], capture_output=True, text=True, check=True).stdout.strip()
        # return os.getenv('GIT_COMMIT')
    except Exception as e:
        st.error(f"Failed to get current git commit: {e}")
        return None


def write_dictionary_to_text_file(dictionary, dict_name, directory):
    """Write archive metadata to a text file with timestamp."""
    try:
        # Ensure we have a Path object
        session_path = pathlib.Path(directory)
        text_file = session_path / f"{dict_name}.txt"

        with open(text_file, "w") as f:
            for key, value in dictionary.items():
                f.write(f"{key}: {value}\n")

            # Add timestamp
            timestamp = framework_utils.get_timestamp()
            f.write(f"timestamp: {timestamp.strftime('%Y-%m-%d %H:%M:%S %Z')}\n")

        return True

    except Exception as e:
        st.error(f"Error writing dictionary {dict_name} to text file in directory {directory}: {e}")
        return False


def main():

    st.header("Save app session")

    with st.columns(2)[0]:
        st.write(f"Current session ID: {st.session_state[ST_KEY_PREFIX_STARTUP + 'app_session_id']}")

        # Allow user to describe the session state.
        key = ST_KEY_PREFIX + "session_description"
        if key not in st.session_state:
            st.session_state[key] = ""
        session_description = st.text_area("Describe the app session archive that will be saved:", key=key)

        # Allow the user to save the current session state.
        if st.button("Save app session"):
            username = pa.get_current_username()
            user_group = pa.get_user_group(username)
            current_git_commit = get_current_git_commit()
            container_image_id = pa.get_frontend_image_id()
            archive_id = framework_utils.get_unique_id()
            app_session_id = st.session_state[ST_KEY_PREFIX_STARTUP + "app_session_id"]
            archive_compatibility_id = pa.get_archive_compatibility_id(container_image_id)
            archive_metadata = {"username": username, "user_group": user_group, "session_description": session_description, "current_git_commit": current_git_commit, "container_image_id": container_image_id, "archive_id": archive_id, "app_session_id": app_session_id, "archive_compatibility_id": archive_compatibility_id}
            write_dictionary_to_text_file(archive_metadata, "archive_metadata", framework_utils.session_dir())  # Writes archive_metadata.txt to the session directory.
            write_conda_environment("environment.yml", framework_utils.session_dir())  # Writes environment.yml to the session directory.
            save_session_state()  # Writes session_state.pkl, session_state.dill, and session_state_contents.txt to the session directory.
            zip_buffer = framework_utils.zip_directory_to_buffer(framework_utils.session_dir(), ignore_subdirs=["input"])
            pa.write_archive_database_data(tuple(archive_metadata.values()))
            pa.upload_zip_object_data(ARCHIVES_BUCKET_NAME, archive_id, zip_buffer)
            pa.get_available_archives.clear()  # Do this to refresh the archive listing below.
            st.success("✅ App session saved successfully!")

    st.header("Load app session")

    if st.button("Refresh app session archive listing"):
        pa.get_available_archives.clear()

    # Display all archives that are compatible with the current frontend image.
    st.write(f"Archives compatible with the current app version (format: `{pa.get_archive_compatibility_id(pa.get_frontend_image_id())}`):")
    df = pa.get_available_archives(pa.get_archive_compatibility_id(pa.get_frontend_image_id()))
    key = ST_KEY_PREFIX + "archive_selection" + "__do_not_persist"
    if not df.is_empty():
        st.dataframe(df, on_select="rerun", selection_mode="single-row", key=key)
    else:
        st.info(f"No app session archives found.")

    if (key in st.session_state) and (st.session_state[key]["selection"]["rows"]):
        selected_row_index = st.session_state[key]["selection"]["rows"][0]
        selected_archive_id = df["Archive ID"][selected_row_index]
        st.write(f"Selected archive ID: {selected_archive_id}")

        if st.button("Load selected app session archive"):
            zip_buffer = pa.download_zip_object_data(ARCHIVES_BUCKET_NAME, selected_archive_id)
            framework_utils.ensure_empty_directory(framework_utils.session_dir())
            framework_utils.unzip_buffer_to_directory(zip_buffer, framework_utils.session_dir())
            load_session_state()
            os.makedirs(os.path.join(framework_utils.session_dir(), "input"), exist_ok=True)  # Ensure input directory exists since we deliberately exclude it when saving an archive.
            # st.rerun()  # Keeping this rerun because masking of errors here is less risky and it's really helpful to see the archive description just pop up when loading an archive.


if __name__ == "__main__":
    main()
