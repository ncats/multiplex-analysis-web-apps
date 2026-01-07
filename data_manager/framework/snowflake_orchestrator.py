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
    m = re.search(r'^(\d+)vcpu_(\d+)gib_', resource)
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

        print("Submitting job with the following parameters:", flush=True)
        print(f"  job_id: {job_id}", flush=True)
        print(f"  username: {username}", flush=True)
        print(f"  session: {session}", flush=True)
        print(f"  selected_compute_resource: {selected_compute_resource}", flush=True)
        print(f"  group_name: {group_name}", flush=True)
        print(f"  app_shortname: {app_shortname}", flush=True)
        print(f"  image_name: {image_name}", flush=True)
        print(f"  image_tag: {image_tag}", flush=True)
        print(f"  app_title: {app_title}", flush=True)

        vcpu, gib = _parse_compute(selected_compute_resource)

        username_validated = _validate_identifier(username, USERNAME_RE, "username")
        job_id_validated = _validate_identifier(job_id, ID_RE, "job_id")

        compute_pool_name = f"{app_shortname}_{username_validated}_workers_{selected_compute_resource}_compute_pool"
        job_service_name = f"{app_shortname}_app_db.{group_name}_schema.{app_shortname}_{username_validated}_worker_{selected_compute_resource}_job_service_{job_id_validated[:10]}"
        role_name = f"data_apps_{username_validated}_role"

        print(f"Computed values:", flush=True)
        print(f"  vcpu: {vcpu}", flush=True)
        print(f"  gib: {gib}", flush=True)
        print(f"  username_validated: {username_validated}", flush=True)
        print(f"  job_id_validated: {job_id_validated}", flush=True)
        print(f"  compute_pool_name: {compute_pool_name}", flush=True)
        print(f"  job_service_name: {job_service_name}", flush=True)
        print(f"  role_name: {role_name}", flush=True)

        job_service_sql = textwrap.dedent(f"""
          EXECUTE JOB SERVICE
            IN COMPUTE POOL {compute_pool_name}
            FROM @{app_shortname}_app_db.general_schema.general_stage SPECIFICATION_TEMPLATE_FILE='worker_service_spec.yaml'
            USING ( APP_SHORTNAME=>'{app_shortname}', APP_TITLE=>' "{app_title}" ', MONITOR_JOBS_REFRESH_INTERVAL_SECONDS=>5, SNOWFLAKE_USER=>' "{username_validated}" ', COMPUTE_RESOURCE=>' "{selected_compute_resource}" ', IMAGE=>' "/{app_shortname}_app_db/general_schema/image_repository/{image_name}:{image_tag}" ', PYTHON_COMMAND=>' "import framework.analysis_framework as analysis_framework; analysis_framework.run_local_analysis(\\'{job_id_validated}\\')" ', SNOWFLAKE_WAREHOUSE=>' "{app_shortname}_{username_validated}_xs_warehouse" ', MOUNTPATH=>' "/tmp/{app_shortname}" ', MEMORY=>'{gib}Gi', CPU=>{vcpu} )
            NAME = {job_service_name}
            ASYNC = TRUE;
        """).strip()

        print(f"Executing job service SQL:\n{job_service_sql}", flush=True)

        session.sql(job_service_sql).collect()

        print(f"Job service {job_service_name} submitted successfully.", flush=True)

        session.sql(f"GRANT MONITOR, OPERATE ON SERVICE {job_service_name} TO ROLE {role_name}").collect()

        print(f"Granted privileges on job service {job_service_name} to role {role_name}.", flush=True)

        worker_image_id = session.sql(f"show service containers in service {job_service_name}").collect()[0]["image_digest"]

        print(f"Retrieved worker image ID: {worker_image_id}", flush=True)

        return worker_image_id
    except Exception as e:
        print(f"Error submitting job {job_id} for user {username}: {e}", flush=True)
        return None


# This is modified from the source/ version for the data manager.
def shutdown(username: str, session: Session, app_shortname: str, group_name: str, compute_resource: str):
    try:

        service_size_mapping = {"1vcpu_6gib_1x": "xs", "6vcpu_28gib_4x": "m"}

        service_size_str = service_size_mapping[compute_resource]

        username_validated = _validate_identifier(username, USERNAME_RE, "username")

        service_name = f"{app_shortname}_db.{group_name}_schema.{app_shortname}_{username_validated}_{service_size_str}_service"
        compute_pool_name = f"{app_shortname}_{username_validated}_{service_size_str}_compute_pool"

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
