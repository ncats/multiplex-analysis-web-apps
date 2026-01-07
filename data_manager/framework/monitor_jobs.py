import streamlit as st
import time
import framework.analysis_framework as analysis_framework
import framework.utils as framework_utils
import framework.platform_abstraction as pa
import os

ST_KEY_PREFIX_STARTUP = "startup.py__"
REFRESH_INTERVAL_SECONDS = int(os.getenv("MONITOR_JOBS_REFRESH_INTERVAL_SECONDS"))


def check_for_updates(c):
    if "JOB_PENDING" in st.session_state:
        job_info = st.session_state["JOB_PENDING"]
        job_id = job_info["job_id"]
        job_status, outputs = analysis_framework.load_job_output_data(job_id, framework_utils.session_dir())
        if job_status == "Completed":
            c.success(f"Job {job_id} completed.")  # Note completion time and duration.
            st.session_state[job_info["key_for_results"]] = outputs
            job_name = job_info["job_name"] if "job_name" in job_info else ""
            st.session_state["JOB_JUST_COMPLETED"] = job_name
            del st.session_state["JOB_PENDING"]
        else:
            c.info(f"Job {job_id} is still in progress. Status: {job_status}.")
    else:
        c.info("No job is currently pending.")


def main():
    c1 = st.empty()
    c2 = st.empty()
    while True:
        check_for_updates(c1)
        pa.get_jobs_table_data.clear()
        jobs_data = pa.get_jobs_table_data()
        c2.dataframe(jobs_data)
        time.sleep(REFRESH_INTERVAL_SECONDS)


if __name__ == "__main__":
    main()
