import streamlit as st
import framework.platform_abstraction as pa


def main():
    if st.button("Refresh data"):
        pa.get_user_groups_table_data.clear()
        pa.get_app_sessions_table_data.clear()
        pa.get_archives_table_data.clear()
        pa.get_jobs_table_data.clear()

    st.header("user_groups table contents")
    user_groups_data = pa.get_user_groups_table_data()
    st.dataframe(user_groups_data)

    st.header("app_sessions table contents")
    app_sessions_data = pa.get_app_sessions_table_data()
    st.dataframe(app_sessions_data)

    st.header("archives table contents")
    archives_data = pa.get_archives_table_data()
    st.dataframe(archives_data)

    st.header("jobs table contents")
    jobs_data = pa.get_jobs_table_data()
    st.dataframe(jobs_data)


if __name__ == "__main__":
    main()
