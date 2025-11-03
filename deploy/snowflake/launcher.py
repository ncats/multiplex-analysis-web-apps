# Import python packages
import streamlit as st
from snowflake.snowpark.context import get_active_session
import pandas as pd


def _get_delta(status):
    if status in ("RUNNING", "ACTIVE", "IDLE", "STARTED", "STARTING", "RESIZING", "STOPPING", "PENDING", "SUSPENDING"):
        delta = "1"
    elif status in ("SUSPENDED", "FAILED", "DONE"):
        delta = "-1"
    else:
        delta = None
    return delta


def start_app(session, db_str="data_app_db", username="andrewweisman"):
    session.sql(f"ALTER SERVICE {db_str}.app_runtime_schema.frontend_service_{username} RESUME;").collect()
    # do the following only if compute pool status isn't idle!
    compute_pool_status = session.sql(f"DESCRIBE COMPUTE POOL data_app_workers_compute_pool_{username};").collect()[0]["state"]
    if compute_pool_status not in ("IDLE", "ACTIVE", "RESIZING"):
        session.sql(f"ALTER COMPUTE POOL data_app_workers_compute_pool_{username} RESUME;").collect()
    st.info("App should be starting now...")


def show_endpoints(session, db_str="data_app_db", username="andrewweisman"):
    ingress_url = session.sql(f"SHOW ENDPOINTS IN SERVICE {db_str}.app_runtime_schema.frontend_service_{username};").collect()[0]["ingress_url"]
    full_url = f"https://{ingress_url}"
    st.write(f"**NOTE: You must open this link in a new tab by right-clicking and choosing something like \"Open link in new tab\" or doing a Ctrl+Click:**")
    st.write(f"{full_url}")


def show_objects(session, db_list=["data_app_db"], username="andrewweisman"):
    df_list = []
    for db_str in db_list:
        df_list.append(session.sql(f"show services in {db_str}.app_runtime_schema;").to_pandas().rename(columns={"\"status\"": "\"state\""}))
    df_list.append(session.sql(f"DESCRIBE COMPUTE POOL data_app_frontend_compute_pool_{username};").to_pandas())
    df_list.append(session.sql(f"DESCRIBE COMPUTE POOL data_app_workers_compute_pool_{username};").to_pandas())
    df_list.append(session.sql("show warehouses").to_pandas())
    df = pd.concat(df_list, ignore_index=True).sort_values("\"updated_on\"", ignore_index=True, ascending=False)
    # 🟢 for positive/active (delta "1"), 🔴 for negative (delta "-1"), none otherwise.
    if "\"state\"" in df.columns:
        def _add_state_emoji(val):
            if val is None:
                return val
            val_str = str(val)
            delta_val = _get_delta(val_str)
            if delta_val == "1":
                return f"🟢 {val_str}"
            elif delta_val == "-1":
                return f"🔴 {val_str}"
            return val_str

        df["\"state\""] = df["\"state\""].apply(_add_state_emoji)
    st.dataframe(df, hide_index=True, use_container_width=True)

def stop_frontend(session, stop_pool_too=True, db_str="data_app_db", username="andrewweisman"):
    session.sql(f"ALTER SERVICE {db_str}.app_runtime_schema.frontend_service_{username} SUSPEND;").collect()
    if stop_pool_too:
        session.sql(f"ALTER COMPUTE POOL data_app_frontend_compute_pool_{username} SUSPEND;").collect()
    st.info("Frontend should be stopping now...")

def stop_workers(session, username="andrewweisman"):
    # STUB: SELECT data_app_db.app_runtime_schema.job_service_andrewweisman_job_id_<JOB_ID>!SPCS_CANCEL_JOB();
    session.sql(f"ALTER COMPUTE POOL data_app_workers_compute_pool_{username} SUSPEND;").collect()
    st.info("Workers should be stopping now...")


def get_frontend_image_id(session, db_str="data_app_db", username="andrewweisman"):
    image_id = session.sql(f"show service containers in service {db_str}.app_runtime_schema.frontend_service_{username}").collect()[0]["image_digest"]
    st.write(image_id)


def main():

    st.title("App Launcher v2")
    st.set_page_config()

    # db_dict = {"Prototype data app": "data_app_db", "HALO Metadata Analyzer": "hma_db"}
    # db_dict = {"HALO Metadata Analyzer": "hma_db"}
    db_dict = {"Multiplex Analysis Web Apps": "mawa"}
    username = "scott_lawrence"  # Change this to your username

    app_names = list(db_dict.keys())
    db_list = list(db_dict.values())

    # Get the current credentials
    session = get_active_session()

    current_user = session.sql("select current_user();").collect()[0]["CURRENT_USER()"]
    current_role = session.sql("select current_role();").collect()[0]["CURRENT_ROLE()"]

    st.header("Information")
    st.text(f"Current user: {current_user}\nCurrent role: {current_role}\nStreamlit version: {st.__version__}")

    st.header("App control")
    chosen_key = st.selectbox("Select app to control:", app_names, key="chosen_key")
    db_str = db_dict[chosen_key]
    action_to_perform = st.selectbox("Select action to perform:", ["Start app", "Show URL", "Stop frontend", "Stop frontend (service only)", "Stop workers", "Retrieve frontend image ID"], key="action_to_perform")
    if st.button("Take action"):
        if action_to_perform == "Start app":
            start_app(session, db_str=db_str, username=username)
        elif action_to_perform == "Show URL":
            show_endpoints(session, db_str=db_str, username=username)
        elif action_to_perform == "Stop frontend":
            stop_frontend(session, db_str=db_str, username=username)
        elif action_to_perform == "Stop frontend (service only)":
            stop_frontend(session, stop_pool_too=False, db_str=db_str, username=username)
        elif action_to_perform == "Stop workers":
            stop_workers(session, username=username)
        elif action_to_perform == "Retrieve frontend image ID":
            get_frontend_image_id(session, db_str=db_str, username=username)

    st.header("General control")
    if st.button("Detect services and compute pools"):
        show_objects(session, db_list=db_list, username=username)


if __name__ == "__main__":
    main()
