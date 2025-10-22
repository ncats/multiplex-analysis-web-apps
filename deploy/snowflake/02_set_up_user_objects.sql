USE ROLE ACCOUNTADMIN;

CREATE COMPUTE POOL IF NOT EXISTS data_app_frontend_compute_pool_andrewweisman
    MIN_NODES = 1
    MAX_NODES = 2
    INSTANCE_FAMILY = CPU_X64_XS
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 1800;
GRANT USAGE, MONITOR, OPERATE ON COMPUTE POOL data_app_frontend_compute_pool_andrewweisman TO ROLE data_app_role;  -- at least initially, none of this may be needed, but doing anyway for completeness

CREATE COMPUTE POOL IF NOT EXISTS data_app_workers_compute_pool_andrewweisman
    MIN_NODES = 1
    MAX_NODES = 5
    INSTANCE_FAMILY = CPU_X64_M
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 600;
-- GRANT USAGE ON COMPUTE POOL data_app_workers_compute_pool_andrewweisman TO ROLE data_app_role;  -- at least initially, this is all that's actually needed
GRANT USAGE, MONITOR, OPERATE ON COMPUTE POOL data_app_workers_compute_pool_andrewweisman TO ROLE data_app_role;  -- doing this anyway for completeness

GRANT BIND SERVICE ENDPOINT ON ACCOUNT TO ROLE data_app_role; -- this seems to be necessary to create the service below

USE ROLE data_app_role;
USE WAREHOUSE data_app_warehouse;

-- **** ACTION: **** Push image to /data_app_db/app_runtime_schema/image_repository/data_app:latest using, e.g.:
-- docker compose build
-- docker compose up (to test locally)
-- docker tag frontend:latest nihnci-eval.registry.snowflakecomputing.com/data_app_db/app_runtime_schema/image_repository/data_app:latest (add a tag to the image we want to push to Snowflake)
-- snow spcs image-registry login --role data_app_role
-- docker push nihnci-eval.registry.snowflakecomputing.com/data_app_db/app_runtime_schema/image_repository/data_app:latest

-- If service already exists, drop it first.
DROP SERVICE IF EXISTS data_app_db.app_runtime_schema.frontend_service_andrewweisman;

CREATE SERVICE IF NOT EXISTS data_app_db.app_runtime_schema.frontend_service_andrewweisman
    IN COMPUTE POOL data_app_frontend_compute_pool_andrewweisman
    -- AUTO_SUSPEND_SECS = 1800
    AUTO_RESUME = FALSE  -- setting to false to hope that shutdown is actually first-time successful
    MIN_INSTANCES=1
    MAX_INSTANCES=1
    FROM SPECIFICATION $$
spec:
  containers:                           # container list
    - name: frontend
      image: /data_app_db/app_runtime_schema/image_repository/data_app:latest
      env:                                # optional
        SNOWFLAKE_USER: andrewweisman             # must inject manually; Snowflake does not inject
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
          memory: 6Gi
          cpu: 1
        limits:
          memory: 6Gi
          cpu: 1
  endpoints:                             # optional endpoint list
    - name: web
      port: 8501                     # specify this or portRange
      public: true
      protocol: HTTP
  volumes:                               # optional volume list
    - name: tmp
      source: local
serviceRoles:                   # Optional list of service roles
  - name: web_endpoint_service_role
    endpoints:
      - web
$$;
alter service data_app_db.app_runtime_schema.frontend_service_andrewweisman suspend;
alter compute pool data_app_frontend_compute_pool_andrewweisman suspend;
