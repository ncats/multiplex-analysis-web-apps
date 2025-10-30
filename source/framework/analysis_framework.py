# This script exemplifies the general use case in which we have functions that potentially: (1) write results to memory (returning objects) and (2) write results to local disk.

import streamlit as st
import os
import framework.platform_abstraction as pa
import framework.utils as framework_utils
import framework.analysis_functions as analysis_functions

ST_KEY_PREFIX_STARTUP = "startup.py__"
JOB_INPUTS_BUCKET_NAME = os.getenv('JOB_INPUTS_BUCKET_NAME')
JOB_OUTPUTS_BUCKET_NAME = os.getenv('JOB_OUTPUTS_BUCKET_NAME')


@st.cache_data()
def get_available_async_resources():
    raw = os.getenv("AVAILABLE_COMPUTE_RESOURCES", "")
    if not raw.strip():
        return []
    return [r for r in raw.split() if r]


def initialize_job(function_name):
    try:
        job_id = framework_utils.get_unique_id()
        job_name = function_name
        submitter = pa.get_current_username()
        submitter_group = pa.get_user_group(submitter)
        app_session_id = st.session_state[ST_KEY_PREFIX_STARTUP + "app_session_id"]
        if not pa.log_job((job_id, job_name, submitter, submitter_group, app_session_id)):
            st.error(f"Failed to log job {job_id} to the jobs table, potentially due to duplicate job ID generation, which should never happen.")
            return None
        return job_id
    except Exception as e:
        st.error(f"Error occurred while initializing job: {e}")
        return None


def save_job_input_data(job_id, inputs):
    try:
        inputs_directory = os.path.join(framework_utils.session_dir(), "tmp_job_inputs")
        framework_utils.ensure_empty_directory(inputs_directory)  # Ensure the inputs directory is empty (create it if necessary).
        framework_utils.serialize_dictionary_to_binary_files(inputs, "inputs", inputs_directory)  # Writes inputs.pkl and inputs.dill to inputs_directory from the inputs dictionary.
        inputs_buffer = framework_utils.zip_directory_to_buffer(inputs_directory)  # Zip the entire inputs_directory to a buffer.
        framework_utils.ensure_empty_directory(inputs_directory, create_if_missing=False)  # Recursively delete the inputs_directory, including the directory itself.
        pa.upload_zip_object_data(JOB_INPUTS_BUCKET_NAME, job_id, inputs_buffer)  # Write the zip buffer to object storage.
        return True
    except Exception as e:
        st.error(f"Error occurred while saving job input data: {e}")
        return False


def load_job_input_data(job_id, job_dir):
    try:
        inputs_directory = os.path.join(job_dir, "inputs")
        inputs_buffer = pa.download_zip_object_data(JOB_INPUTS_BUCKET_NAME, job_id)
        framework_utils.ensure_empty_directory(inputs_directory)
        framework_utils.unzip_buffer_to_directory(inputs_buffer, inputs_directory)
        inputs = framework_utils.deserialize_binary_files_to_dictionary("inputs", inputs_directory)  # Loads inputs.pkl and inputs.dill from inputs_directory into an inputs dictionary.
        framework_utils.ensure_empty_directory(inputs_directory, create_if_missing=False)
        return inputs
    except Exception as e:
        st.error(f"Error occurred while loading job input data: {e}")
        return None


def save_job_output_data(job_id, outputs, job_dir):
    try:
        if outputs is not None:
            outputs_directory = os.path.join(job_dir, "outputs")
            os.makedirs(outputs_directory, exist_ok=True)  # We don't ensure this is an empty directory because it may potentially contain output files from a job.
            framework_utils.serialize_dictionary_to_binary_files(outputs, "outputs", outputs_directory)  # Writes outputs.pkl and outputs.dill to outputs_directory.
            outputs_buffer = framework_utils.zip_directory_to_buffer(outputs_directory)
            framework_utils.ensure_empty_directory(outputs_directory, create_if_missing=False)
            pa.upload_zip_object_data(JOB_OUTPUTS_BUCKET_NAME, job_id, outputs_buffer)
            return True
        else:
            return False
    except Exception as e:
        st.error(f"Error occurred while saving job output data: {e}")
        return None


def delete_serialized_files(dict_name, directory):
    try:
        pkl_file = os.path.join(directory, f'{dict_name}.pkl')
        os.remove(pkl_file)

        dill_file = os.path.join(directory, f'{dict_name}.dill')
        os.remove(dill_file)

        return True
    except Exception as e:
        st.error(f"Failed to delete serialized files {dict_name}.pkl/.dill in directory {directory}: {e}")
        return False


def load_job_output_data(job_id, outputs_directory):
    try:
        job_status = pa.get_job_status(job_id)

        outputs = None
        if job_status == "Submitted":
            pass
        elif job_status == "Running":
            pass
        elif job_status == "Completed":
            outputs_buffer = pa.download_zip_object_data(JOB_OUTPUTS_BUCKET_NAME, job_id)  # Creates a buffer of the job results. Buffer likely contains both .pkl/.dill files and any other output files the job may have created.
            framework_utils.unzip_buffer_to_directory(outputs_buffer, outputs_directory)  # Unzips the buffer to the session directory.
            outputs = framework_utils.deserialize_binary_files_to_dictionary("outputs", outputs_directory)  # Loads outputs.pkl and outputs.dill from the session directory into an "outputs" dictionary.
            delete_serialized_files("outputs", outputs_directory)  # Delete the .pkl/.dill files from the session directory.
        elif job_status == "Failed":
            pass
        else:
            pass

        return job_status, outputs
    except Exception as e:
        st.error(f"Error occurred while loading job output data: {e}")
        return None


def run_local_analysis(job_id):
    try:
        pa.update_job_status(job_id, "Running", "start_time")
        job_dir = os.path.join(framework_utils.jobs_dir(), job_id)
        function_name = pa.get_job_function_name(job_id)
        inputs = load_job_input_data(job_id, job_dir)
        outputs = analysis_functions.run_analysis_job(function_name, inputs, job_dir)
        success = save_job_output_data(job_id, outputs, job_dir)
        if success:
            pa.update_job_status(job_id, "Completed", "completion_time")
        else:
            pa.update_job_status(job_id, "Failed", "failure_time")
        return True
    except Exception as e:
        # Since this function can run on a worker (without Streamlit), echo the error both to the terminal and to the screen (if present).
        error_string = f"Error occurred while running local analysis for job {job_id}: {e}"
        st.error(error_string)
        print(error_string)
        return False


def run_analysis_job_wrapper(function_name, inputs, blocking=True, selected_compute_resource: str = None):
    """
    Run an analysis job by calling the specified function with the given inputs.
    """
    try:
        job_id = initialize_job(function_name)
        save_job_input_data(job_id, inputs)
        pa.submit_job(job_id, blocking=blocking, selected_compute_resource=selected_compute_resource)
        return job_id
    except Exception as e:
        st.error(f"Error occurred while running wrapper for analysis job {function_name}: {e}")
        return None


# TODO: Add ability to save archive (with description) at the end; see functionality in manage_sessions.py.
def job_submission(job_name, inputs, analysis_purpose, st_key_prefix):

    analysis_purpose_with_underscores = analysis_purpose.replace(" ", "_")

    with st.columns(1, border=True)[0]:

        st.subheader(analysis_purpose.capitalize())
    
        key = st_key_prefix + "async_analysis_" + analysis_purpose_with_underscores
        if key not in st.session_state:
            st.session_state[key] = True
        do_async_analysis = st.toggle(f"Run {analysis_purpose} asynchronously", key=key)

        available_async_resources = get_available_async_resources()
        if do_async_analysis and available_async_resources:
            selected_compute_resource = st.selectbox("Select compute resource:", available_async_resources, key=st_key_prefix + "selected_compute_resource_for_" + analysis_purpose_with_underscores)
        else:
            selected_compute_resource = None

        if st.button(f"Run {analysis_purpose}", type=("primary" if do_async_analysis else "secondary")):
            job_id = run_analysis_job_wrapper(job_name, inputs, blocking=not do_async_analysis, selected_compute_resource=selected_compute_resource)
            st.session_state["JOB_PENDING"] = {"job_id": job_id, "key_for_results": st_key_prefix + analysis_purpose_with_underscores + "_results"}
            # st.rerun()  # Remove to not mask any potential warnings/errors.
