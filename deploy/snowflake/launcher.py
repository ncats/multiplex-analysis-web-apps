# Import relevant libraries.
import streamlit as st
from snowflake.snowpark.context import get_active_session
import pandas as pd
import re


@st.cache_data()
def get_owner_role():
    session = get_active_session()
    return session.sql("select current_role();").collect()[0]["CURRENT_ROLE()"]  # Note current_user() doesn't work in Streamlit apps.


# Get services in the schema {chosen_app_shortname}_app_db.{user_group}_schema that have "_frontend_" in their name, returning a potentially empty list.
@st.cache_data()
def get_frontend_service_names(chosen_app_shortname, user_group, username):
    session = get_active_session()
    raw_services_df = session.sql(
        f"show services in schema {chosen_app_shortname}_app_db.{user_group}_schema;"
    ).to_pandas()
    frontend_services_df = raw_services_df[
        raw_services_df["name"].str.contains(f"_{username}_frontend_", case=False, na=False)
    ]
    return frontend_services_df["name"].tolist()  # potentially empty list


# This should be the same as in platform_abstraction.py in the Snowflake branch.
@st.cache_data()
def get_user_group(username):
    session = get_active_session()
    try:
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
def extract_user_segment(role: str) -> str | None:
    m = re.match(r"^data_apps_(.+?)_role$", role)
    return m.group(1) if m else None


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

    # Display the page title.
    st.title("App Launcher v2")

    # Get a list of apps subject to the new organization scheme.
    app_shortname_dict = {"Multiplex Analysis Web Apps": "mawa"}

    # Get the corresponding keys and values.
    app_titles = list(app_shortname_dict.keys())
    app_shortnames = list(app_shortname_dict.values())

    # Get app owner's role.
    current_role = get_owner_role()

    # Extract the username from the role.
    username = extract_user_segment(current_role)
    if username is None:
        st.error(f"Could not extract username from current role: {current_role}. Ensure you are using a role like data_apps_<username>_role.")
        return
    
    # Get the user group from the user_groups table in common_db.
    user_group = get_user_group(username)

    # Display some information.
    st.header("Information")
    st.text(f"Current user: {username}\nUser group: {user_group}\nCurrent role: {current_role}\nStreamlit version: {st.__version__}")

    # Start controling the app.
    st.header("App control")

    # Let the user choose the app to control.
    chosen_key = st.selectbox("Select app to control:", app_titles)
    chosen_app_shortname = app_shortname_dict[chosen_key]

    # Get the names of the frontend services available for the current user in their selected app.
    frontend_service_names = get_frontend_service_names(chosen_app_shortname, user_group, username)
    if not frontend_service_names:
        st.warning(f"No frontend services found in schema {chosen_app_shortname}_app_db.{user_group}_schema.")
        return
    
    # Allow the user to select which frontend service (i.e., version of the app) they'd like to control.
    chosen_frontend_service_name = st.selectbox("Select frontend service to control:", frontend_service_names)

    # Ask the user what they want to do.
    action_to_perform = st.selectbox("Select action to perform:", ["Start app", "Show URL", "Stop frontend", "Stop frontend (service only)", "Stop workers", "Retrieve frontend image ID"])

    #### PICK UP HERE WITH MODIFYING THE BELOW AND CONSIDERING WHETHER USERS WILL DETECT FRONTENDS FROM OTHER MEMBERS OF THEIR GROUP (probably taken care of now)!!!! ####
    if st.button("Take action"):
        if action_to_perform == "Start app":
            start_app(session, db_str=chosen_app_shortname, username=username)
        elif action_to_perform == "Show URL":
            show_endpoints(session, db_str=chosen_app_shortname, username=username)
        elif action_to_perform == "Stop frontend":
            stop_frontend(session, db_str=chosen_app_shortname, username=username)
        elif action_to_perform == "Stop frontend (service only)":
            stop_frontend(session, stop_pool_too=False, db_str=chosen_app_shortname, username=username)
        elif action_to_perform == "Stop workers":
            stop_workers(session, username=username)
        elif action_to_perform == "Retrieve frontend image ID":
            get_frontend_image_id(session, db_str=chosen_app_shortname, username=username)

    st.header("General control")
    if st.button("Detect services and compute pools"):
        show_objects(session, db_list=app_shortnames, username=username)


if __name__ == "__main__":
    main()
