# Import relevant libraries.
import streamlit as st
from snowflake.snowpark.context import get_active_session
import pandas as pd
import re


@st.cache_data()
def get_owner_role():
    session = get_active_session()
    return session.sql("select current_role();").collect()[0]["CURRENT_ROLE()"]  # Note current_user() doesn't work in Streamlit apps.


@st.cache_data()
def get_db_name(chosen_app_shortname):
    if chosen_app_shortname == "data_manager":
        return "data_manager_db"
    else:
        return f"{chosen_app_shortname}_app_db"


# Get services in the schema {chosen_app_shortname}_app_db.{user_group}_schema that have "_frontend_" in their name, returning a potentially empty list.
@st.cache_data()
def get_service_names(chosen_app_shortname, user_group, username, suffix=None, invert=False):
    if suffix is None:
        suffix = ""
    elif suffix == "frontend":
        suffix = "_frontend"
    elif suffix == "worker":
        suffix = "_worker"
    else:
        raise ValueError(f"Invalid suffix: {suffix}")
    session = get_active_session()
    raw_services_df = session.sql(
        f"show services in schema {get_db_name(chosen_app_shortname)}.{user_group}_schema;"
    ).to_pandas()
    mask = raw_services_df["\"name\""].str.lower().str.contains(f"_{username}{suffix}_", case=False, na=False)
    if invert:
        matches_df = raw_services_df[~mask]
    else:
        matches_df = raw_services_df[mask]
    return matches_df["\"name\""].tolist()  # potentially empty list


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
    m = re.match(r"^data_apps_(.+?)_role$", role.lower())
    return m.group(1) if m else None


def _get_delta(status):
    if status in ("RUNNING", "ACTIVE", "IDLE", "STARTED", "STARTING", "RESIZING", "STOPPING", "PENDING", "SUSPENDING"):
        delta = "1"
    elif status in ("SUSPENDED", "FAILED", "DONE"):
        delta = "-1"
    else:
        delta = None
    return delta


def start_app(chosen_app_shortname, user_group, service_name):
    session = get_active_session()
    session.sql(f"ALTER SERVICE {get_db_name(chosen_app_shortname)}.{user_group}_schema.{service_name} RESUME;").collect()
    st.info("App should be starting now...")


def show_endpoints(chosen_app_shortname, user_group, service_name):
    session = get_active_session()
    ingress_url = session.sql(f"SHOW ENDPOINTS IN SERVICE {get_db_name(chosen_app_shortname)}.{user_group}_schema.{service_name};").collect()[0]["ingress_url"]
    full_url = f"https://{ingress_url}"
    st.write(f"**NOTE: You must open this link in a new tab by right-clicking and choosing something like \"Open link in new tab\" or doing a Ctrl+Click:**")
    st.write(f"{full_url}")


def stop_frontend(chosen_app_shortname, user_group, service_name, stop_pool_too=True):
    session = get_active_session()
    session.sql(f"ALTER SERVICE {get_db_name(chosen_app_shortname)}.{user_group}_schema.{service_name} SUSPEND;").collect()
    if stop_pool_too:
        compute_pool_name = service_name.lower().removesuffix("_service") + "_compute_pool"
        session.sql(f"ALTER COMPUTE POOL {compute_pool_name} SUSPEND;").collect()
    st.info(f"Selected frontend {'(and compute pool)' if stop_pool_too else '(only)'} should be stopping now...")


def stop_workers(service_name):
    session = get_active_session()
    if "_frontend_" in service_name.lower():
        compute_pool_name = service_name.lower().replace("_frontend_", "_workers_").removesuffix("_service") + "_compute_pool"
        # STUB: SELECT data_app_db.app_runtime_schema.job_service_andrewweisman_job_id_<JOB_ID>!SPCS_CANCEL_JOB(); --> not needed with "STOP ALL" below since that stops all running jobs
        session.sql(f"ALTER COMPUTE POOL {compute_pool_name} STOP ALL;").collect()
        st.info(f"Worker pool {compute_pool_name} and corresponding jobs should be stopping now...")
    else:
        st.warning(f"Service name {service_name} is not associated with a worker pool.")


def write_app_image_id(chosen_app_shortname, user_group, service_name):
    session = get_active_session()
    image_id = session.sql(f"show service containers in service {get_db_name(chosen_app_shortname)}.{user_group}_schema.{service_name}").collect()[0]["image_digest"]
    st.write(image_id)


def show_objects(app_shortnames, user_group):
    session = get_active_session()
    df_list = []
    for app_shortname in app_shortnames:
        df_list.append(session.sql(f"show services in schema {get_db_name(app_shortname)}.{user_group}_schema;").to_pandas().rename(columns={"\"status\"": "\"state\""}))
    df_list.append(session.sql("show compute pools").to_pandas())
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
    st.dataframe(df, use_container_width=True)


def main():

    # Display the page title.
    st.title("App Launcher v2")

    # Get a list of apps subject to the new organization scheme.
    # app_shortname_dict = {"Data Manager": "data_manager", "Multiplex Analysis Web Apps": "mawa"}
    app_shortname_dict = {"Multiplex Analysis Web Apps": "mawa"}

    # Get the corresponding keys.
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

    # Get the names of the startable services available for the current user in their selected app.
    startable_service_names = get_service_names(chosen_app_shortname, user_group, username, suffix="worker", invert=True)
    if not startable_service_names:
        st.warning(f"No startable services found.")
        return
    
    # Allow the user to select which startable service (i.e., version of the app) they'd like to control. Make the default a minimal-compute-resource one, if available.
    key = "chosen_startable_service_name"
    if key not in st.session_state:
        minimal_services = [x for x in startable_service_names if "_1x_" in x.lower()]
        if minimal_services:
            st.session_state[key] = minimal_services[0]
        else:
            st.session_state[key] = startable_service_names[0]
    chosen_startable_service_name = st.selectbox("Select startable service to control:", startable_service_names, key=key)

    # Ask the user what they want to do.
    action_to_perform = st.selectbox("Select action to perform:", ["Start app", "Show URL", "Stop frontend (and compute pool)", "Stop frontend (service only)", "Stop corresponding workers", "Retrieve app image ID"])

    if st.button("Take action"):
        if action_to_perform == "Start app":
            start_app(chosen_app_shortname, user_group, chosen_startable_service_name)
        elif action_to_perform == "Show URL":
            show_endpoints(chosen_app_shortname, user_group, chosen_startable_service_name)
        elif action_to_perform == "Stop frontend (and compute pool)":
            stop_frontend(chosen_app_shortname, user_group, chosen_startable_service_name)
        elif action_to_perform == "Stop frontend (service only)":
            stop_frontend(chosen_app_shortname, user_group, chosen_startable_service_name, stop_pool_too=False)
        elif action_to_perform == "Stop corresponding workers":
            stop_workers(chosen_startable_service_name)
        elif action_to_perform == "Retrieve app image ID":
            write_app_image_id(chosen_app_shortname, user_group, chosen_startable_service_name)

    st.header("General control")
    if st.button("Detect all services, compute pools, and warehouses"):
        show_objects(app_shortnames, user_group)


if __name__ == "__main__":
    main()
