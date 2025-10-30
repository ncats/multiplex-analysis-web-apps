import re
import textwrap
from snowflake.snowpark import Session

USERNAME_RE = re.compile(r"^[A-Za-z0-9_]+$")
ID_RE = re.compile(r"^[A-Za-z0-9_]+$")


def _validate_identifier(value: str, pattern: re.Pattern, label: str):
    if not pattern.match(value):
        raise ValueError(f"Invalid {label}: {value}")
    return value


def submit_job(job_id: str, username: str, session: Session, selected_compute_resource: str):
    """
    Submit a Snowpark Container Services job (EXECUTE JOB SERVICE) and grant privileges.
    session: Snowpark Session
    """
    try:
        username_validated = _validate_identifier(username, USERNAME_RE, "username")
        job_id_validated = _validate_identifier(job_id, ID_RE, "job_id")

        compute_pool_name = f"data_app_workers_compute_pool_{username_validated}"
        job_service_name = f"data_app_db.app_runtime_schema.job_service_{username_validated}_job_id_{job_id_validated[:10]}"
        role_name = f"data_apps_{username_validated}_role"

        job_service_sql = textwrap.dedent(f"""
          EXECUTE JOB SERVICE
              IN COMPUTE POOL {compute_pool_name}
              NAME = {job_service_name}
              ASYNC = TRUE
              FROM SPECIFICATION $$
          spec:
            containers:                           # container list
              - name: worker
                image: /data_app_db/app_runtime_schema/image_repository/data_app:latest
                command:                            # optional list of strings
                  - /opt/conda/bin/python
                  - -c
                  - "import framework.analysis_framework as analysis_framework; analysis_framework.run_local_analysis('{job_id_validated}')"
                env:                                # optional
                  SNOWFLAKE_USER: {username_validated}             # must inject manually; Snowflake does not inject
                  SNOWFLAKE_WAREHOUSE: data_app_warehouse   # must inject manually; Snowflake does not inject
                  TZ: America/New_York
                  ARCHIVES_BUCKET_NAME: archives
                  OLD_ARCHIVES_BUCKET_NAME: oldarchives
                  JOB_INPUTS_BUCKET_NAME: inputs
                  JOB_OUTPUTS_BUCKET_NAME: outputs
                  DATA_OBJECTS_BUCKET_NAME: objects
                  APP_PLATFORM: snowflake
                  APP_NAME: mawa
                  APP_TITLE: "Multiplex Analysis Web Apps"
                  MONITOR_JOBS_REFRESH_INTERVAL_SECONDS: 5
                volumeMounts:                       # optional list
                  - name: tmp
                    mountPath: /tmp/multiplex_analysis_web_apps
                resources:                          # optional
                  requests:
                    memory: 28Gi
                    cpu: 6
                  limits:
                    memory: 28Gi
                    cpu: 6
            volumes:                               # optional volume list
              - name: tmp
                source: local
          $$
        """).strip()
        session.sql(job_service_sql).collect()
        session.sql(f"GRANT MONITOR, OPERATE ON SERVICE {job_service_name} TO ROLE {role_name}").collect()
        worker_image_id = session.sql(f"show service containers in service {job_service_name}").collect()[0]["image_digest"]
        return worker_image_id
    except Exception as e:
        print(f"Error submitting job {job_id} for user {username}: {e}")
        return None


# On Snowflake, since there are no database, object storage, or container orchestration resources to clean up like there are locally with Docker, this means we simply kill the frontend.
def shutdown(username: str, session: Session):
    try:
        username_validated = _validate_identifier(username, USERNAME_RE, "username")
        service_name = f"data_app_db.app_runtime_schema.frontend_service_{username_validated}"
        compute_pool_name = f"data_app_frontend_compute_pool_{username_validated}"
        session.sql(f"ALTER SERVICE {service_name} SUSPEND").collect()
        session.sql(f"ALTER COMPUTE POOL {compute_pool_name} SUSPEND").collect()
        return True
    except Exception as e:
        print(f"Error during shutdown: {e}")
        return False


def frontend_id(username: str, session: Session):
    try:
        username_validated = _validate_identifier(username, USERNAME_RE, "username")
        service_name = f"data_app_db.app_runtime_schema.frontend_service_{username_validated}"
        frontend_image_id = session.sql(f"show service containers in service {service_name}").collect()[0]["image_digest"]
        return frontend_image_id
    except Exception as e:
        print(f"Error retrieving frontend ID: {e}")
        return None
