import re
import textwrap
from snowflake.snowpark import Session

USERNAME_RE = re.compile(r"^[A-Za-z0-9_]+$")
ID_RE = re.compile(r"^[A-Za-z0-9_]+$")


def _parse_compute(resource: str):
    """
    Extract vcpu and gib integers from strings like:
    CPU_X64_XS_1vcpu_6gib_1x
    CPU_X64_L_28vcpu_116gib_14x
    HIGHMEM_X64_M_28vcpu_240gib_19x
    """
    m = re.search(r'_(\d+)vcpu_(\d+)gib_', resource)
    if not m:
        raise ValueError(f"Unrecognized compute resource format: {resource}")
    vcpu = int(m.group(1))
    gib = int(m.group(2))
    return vcpu, gib


def _validate_identifier(value: str, pattern: re.Pattern, label: str):
    if not pattern.match(value):
        raise ValueError(f"Invalid {label}: {value}")
    return value


def submit_job(job_id: str, username: str, session: Session, selected_compute_resource: str, group_name: str, app_shortname: str, image_name: str, image_tag: str, app_title: str):
    """
    Submit a Snowpark Container Services job (EXECUTE JOB SERVICE) and grant privileges.
    session: Snowpark Session
    """
    try:
        vcpu, gib = _parse_compute(selected_compute_resource)

        username_validated = _validate_identifier(username, USERNAME_RE, "username")
        job_id_validated = _validate_identifier(job_id, ID_RE, "job_id")

        compute_pool_name = f"{app_shortname}_{username_validated}_workers_{selected_compute_resource}_compute_pool"
        job_service_name = f"{app_shortname}_app_db.{group_name}_schema.{app_shortname}_{username_validated}_worker_{selected_compute_resource}_job_service_{job_id_validated[:10]}"
        role_name = f"data_apps_{username_validated}_role"

        job_service_sql = textwrap.dedent(f"""
          EXECUTE JOB SERVICE
            IN COMPUTE POOL {compute_pool_name}
            FROM @{app_shortname}_app_db.general_schema.general_stage SPECIFICATION_TEMPLATE_FILE='worker_service_spec.yaml'
            USING ( APP_SHORTNAME=>'{app_shortname}', IMAGE_NAME=>'{image_name}', IMAGE_TAG=>'{image_tag}', REQUESTS_MEMORY_GI=>{gib}, REQUESTS_CPU=>{vcpu}, LIMITS_MEMORY_GI=>{gib}, LIMITS_CPU=>{vcpu}, APP_TITLE=>'{app_title}', MONITOR_JOBS_REFRESH_INTERVAL_SECONDS=>5, SNOWFLAKE_USER=>'{username_validated}', COMPUTE_RESOURCE=>'{selected_compute_resource}', WAREHOUSE_SIZE=>'xs', JOB_ID_VALIDATED=>'{job_id_validated}' )
            NAME = {job_service_name}
            ASYNC = TRUE;
        """).strip()
        session.sql(job_service_sql).collect()
        session.sql(f"GRANT MONITOR, OPERATE ON SERVICE {job_service_name} TO ROLE {role_name}").collect()
        worker_image_id = session.sql(f"show service containers in service {job_service_name}").collect()[0]["image_digest"]
        return worker_image_id
    except Exception as e:
        print(f"Error submitting job {job_id} for user {username}: {e}")
        return None


# On Snowflake, since there are no database, object storage, or container orchestration resources to clean up like there are locally with Docker, this means we simply kill the frontend.
def shutdown(username: str, session: Session, app_shortname: str, group_name: str, compute_resource: str):
    try:

        username_validated = _validate_identifier(username, USERNAME_RE, "username")

        service_name = f"{app_shortname}_app_db.{group_name}_schema.{app_shortname}_{username_validated}_frontend_{compute_resource}_service"
        compute_pool_name = f"{app_shortname}_{username_validated}_frontend_{compute_resource}_compute_pool"

        session.sql(f"ALTER SERVICE {service_name} SUSPEND").collect()
        session.sql(f"ALTER COMPUTE POOL {compute_pool_name} SUSPEND").collect()
        return True
    except Exception as e:
        print(f"Error during shutdown: {e}")
        return False


def frontend_id(username: str, session: Session, app_shortname: str, group_name: str, compute_resource: str):
    try:
        username_validated = _validate_identifier(username, USERNAME_RE, "username")
        service_name = f"{app_shortname}_app_db.{group_name}_schema.{app_shortname}_{username_validated}_frontend_{compute_resource}_service"
        frontend_image_id = session.sql(f"show service containers in service {service_name}").collect()[0]["image_digest"]
        return frontend_image_id
    except Exception as e:
        print(f"Error retrieving frontend ID: {e}")
        return None
