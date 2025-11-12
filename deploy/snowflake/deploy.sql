-- Become admin.
use role accountadmin;

-- Create a warehouse for setup purposes.
CREATE WAREHOUSE IF NOT EXISTS setup_xs_warehouse
  WAREHOUSE_SIZE = 'XSMALL'
  AUTO_RESUME = TRUE
  INITIALLY_SUSPENDED = TRUE;
use warehouse setup_xs_warehouse;

-- Create a warehouse for general usage.
CREATE WAREHOUSE IF NOT EXISTS general_xs_warehouse
  WAREHOUSE_SIZE = 'XSMALL'
  AUTO_RESUME = TRUE
  INITIALLY_SUSPENDED = TRUE;

-- Do this because creating a warehouse such as above switches to that warehouse at least in the Snowflake VS Code extension.
use warehouse setup_xs_warehouse;


---------------- Database group_alpha_group_db. ---------------------------------------------------
-- Create the database and schemas.
create database if not exists group_alpha_group_db;
use database group_alpha_group_db;
create schema if not exists curated_schema;
create schema if not exists app_a_schema;

-- Create stages.
create stage if not exists curated_schema.objects_stage
  directory = ( enable = true );
create stage if not exists app_a_schema.archives_stage
  directory = ( enable = true );
create stage if not exists app_a_schema.inputs_stage
  directory = ( enable = true );
create stage if not exists app_a_schema.outputs_stage
  directory = ( enable = true );
create stage if not exists app_a_schema.oldarchives_stage
  directory = ( enable = true );

-- Create tables.
CREATE TABLE IF NOT EXISTS app_a_schema.app_sessions_table (
  id INTEGER IDENTITY PRIMARY KEY,
  app_session_id VARCHAR(255) UNIQUE NOT NULL,
  username VARCHAR(255),
  user_group VARCHAR(255),
  startup_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
  explicit_shutdown_time TIMESTAMP,
  container_image_id VARCHAR(255),
  compute_resource VARCHAR(255)
);
CREATE TABLE IF NOT EXISTS app_a_schema.archives_table (
  id INTEGER IDENTITY PRIMARY KEY,
  creator VARCHAR(255),
  user_group VARCHAR(255),
  archive_description VARCHAR,
  current_git_commit VARCHAR(255),
  container_image_id VARCHAR(255),
  archive_id VARCHAR(255) UNIQUE NOT NULL,
  app_session_id VARCHAR(255),
  archive_compatibility_id INTEGER,
  creation_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
CREATE TABLE IF NOT EXISTS app_a_schema.jobs_table (
  id INTEGER IDENTITY PRIMARY KEY,
  job_id VARCHAR(255) UNIQUE NOT NULL,
  job_name VARCHAR(255),
  job_status VARCHAR(255),
  submitter VARCHAR(255),
  submitter_group VARCHAR(255),
  app_session_id VARCHAR(255),
  worker_image_id VARCHAR(255),
  submission_time TIMESTAMP,
  start_time TIMESTAMP,
  completion_time TIMESTAMP,
  failure_time TIMESTAMP,
  compute_resource VARCHAR(255)
);

-- Create roles.
create database role if not exists curated_schema_rw_db_role;
create database role if not exists app_a_schema_rw_db_role;
create role if not exists data_apps_user_1_role; -- This is the account role that will be granted the service role for the app.

-- Grant curated_schema_rw_db_role privileges to see the database and schema.
GRANT USAGE ON DATABASE group_alpha_group_db
  TO DATABASE ROLE curated_schema_rw_db_role;
GRANT USAGE ON SCHEMA curated_schema
  TO DATABASE ROLE curated_schema_rw_db_role;

-- Grant curated_schema_rw_db_role privileges to read and write files from/to the stage.
GRANT READ ON STAGE curated_schema.objects_stage  -- Allow reading files from the stage (LIST/GET, COPY INTO <table> FROM @stage)
  TO DATABASE ROLE curated_schema_rw_db_role;
GRANT WRITE ON STAGE curated_schema.objects_stage  -- Allow writing files to the stage (PUT/REMOVE, COPY INTO @stage)
  TO DATABASE ROLE curated_schema_rw_db_role;

-- Grant app_a_schema_rw_db_role privileges to see the database and schema.
GRANT USAGE ON DATABASE group_alpha_group_db
  TO DATABASE ROLE app_a_schema_rw_db_role;
GRANT USAGE ON SCHEMA app_a_schema
  TO DATABASE ROLE app_a_schema_rw_db_role;

-- Grant app_a_schema_rw_db_role privileges to read and write files from/to the stages, except for the oldarchives_stage, which should be read-only.
GRANT READ ON STAGE app_a_schema.archives_stage  -- Allow reading files from the stage (LIST/GET, COPY INTO <table> FROM @stage)
  TO DATABASE ROLE app_a_schema_rw_db_role;
GRANT WRITE ON STAGE app_a_schema.archives_stage  -- Allow writing files to the stage (PUT/REMOVE, COPY INTO @stage)
  TO DATABASE ROLE app_a_schema_rw_db_role;
GRANT READ ON STAGE app_a_schema.inputs_stage  -- Allow reading files from the stage (LIST/GET, COPY INTO <table> FROM @stage)
  TO DATABASE ROLE app_a_schema_rw_db_role;
GRANT WRITE ON STAGE app_a_schema.inputs_stage  -- Allow writing files to the stage (PUT/REMOVE, COPY INTO @stage)
  TO DATABASE ROLE app_a_schema_rw_db_role;
GRANT READ ON STAGE app_a_schema.outputs_stage  -- Allow reading files from the stage (LIST/GET, COPY INTO <table> FROM @stage)
  TO DATABASE ROLE app_a_schema_rw_db_role;
GRANT WRITE ON STAGE app_a_schema.outputs_stage  -- Allow writing files to the stage (PUT/REMOVE, COPY INTO @stage)
  TO DATABASE ROLE app_a_schema_rw_db_role;
GRANT READ ON STAGE app_a_schema.oldarchives_stage  -- Allow reading files from the stage (LIST/GET, COPY INTO <table> FROM @stage)
  TO DATABASE ROLE app_a_schema_rw_db_role;

-- Grant app_a_schema_rw_db_role privileges to read and write data from/to the tables. Note that we probably don't need all of select, insert, update on all tables, but for simplicity we give all three here.
GRANT SELECT, INSERT, UPDATE  -- Table-level privileges (read + write)
  ON TABLE app_a_schema.app_sessions_table
  TO DATABASE ROLE app_a_schema_rw_db_role;

GRANT SELECT, INSERT, UPDATE
  ON TABLE app_a_schema.archives_table
  TO DATABASE ROLE app_a_schema_rw_db_role;

GRANT SELECT, INSERT, UPDATE
  ON TABLE app_a_schema.jobs_table
  TO DATABASE ROLE app_a_schema_rw_db_role;

-- Other required grants for the data_apps_user_1_role to robustly use the app/Snowflake.
GRANT USAGE ON WAREHOUSE general_xs_warehouse TO ROLE data_apps_user_1_role;
GRANT USAGE ON DATABASE group_alpha_group_db TO ROLE data_apps_user_1_role;
GRANT USAGE ON SCHEMA curated_schema TO ROLE data_apps_user_1_role;
GRANT USAGE ON SCHEMA app_a_schema TO ROLE data_apps_user_1_role;
GRANT DATABASE ROLE curated_schema_rw_db_role TO ROLE data_apps_user_1_role;
GRANT READ ON STAGE app_a_schema.oldarchives_stage TO ROLE data_apps_user_1_role;
GRANT WRITE ON STAGE app_a_schema.oldarchives_stage TO ROLE data_apps_user_1_role;
---------------- End database group_alpha_group_db. -----------------------------------------------


---------------- Database common_db. ---------------------------------------------------
-- Create the database and schemas.
create database if not exists common_db;
use database common_db;
create schema if not exists admin_schema;

-- Create tables.
-- I believe this table is now primarily/only really for indicating into which database the app should write given the user running the app.
CREATE TABLE IF NOT EXISTS admin_schema.user_groups_table (
  id INTEGER IDENTITY PRIMARY KEY,
  username VARCHAR(255) UNIQUE NOT NULL,
  user_group VARCHAR(255),
  user_added_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
  who_added VARCHAR(255),
  user_email VARCHAR(255)
);

-- SEE GITHUB README FOR WHAT DATA TO ADD TO THIS TABLE.

-- Create database roles.
create database role if not exists admin_schema_ro_db_role;

-- Grant admin_schema_ro_db_role privileges to see the database and schema.
GRANT USAGE ON DATABASE common_db
  TO DATABASE ROLE admin_schema_ro_db_role;
GRANT USAGE ON SCHEMA admin_schema
  TO DATABASE ROLE admin_schema_ro_db_role;

-- Grant admin_schema_ro_db_role privileges to read data from the tables.
GRANT SELECT  -- Table-level privileges (read only)
  ON TABLE admin_schema.user_groups_table
  TO DATABASE ROLE admin_schema_ro_db_role;
---------------- End database common_db. -----------------------------------------------


---------------- Database app_a_app_db. ---------------------------------------------------
-- Create the database and schemas.
create database if not exists app_a_app_db;
use database app_a_app_db;
create schema if not exists group_alpha_schema;
create schema if not exists general_schema;

-- Create a general stage for holding the code for the services.
create stage if not exists general_schema.general_stage
  directory = ( enable = true );

-- SEE GITHUB README FOR WHAT FILES TO UPLOAD TO THIS STAGE (the two service specification YAML files).

-- Create an image repository.
CREATE IMAGE REPOSITORY IF NOT EXISTS general_schema.image_repository;

-- SEE GITHUB README FOR HOW TO UPLOAD AN IMAGE TO THIS REPOSITORY.

-- Create table.
CREATE TABLE IF NOT EXISTS general_schema.image_metadata_table (
  id INTEGER IDENTITY PRIMARY KEY,
  image_id VARCHAR(255) UNIQUE NOT NULL,
  name VARCHAR(255),
  tag VARCHAR(255),
  git_commit VARCHAR(255),
  environment_yaml_file VARCHAR(255),
  archive_compatibility_id INTEGER,
  image_added_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
  who_added VARCHAR(255)
);

-- SEE GITHUB README FOR WHAT DATA TO ADD TO THIS TABLE.

-- Create roles.
create role if not exists app_a_group_alpha_role; -- This is the account role that owns and operates the app.
create database role if not exists general_schema_ro_db_role;
CREATE DATABASE ROLE IF NOT EXISTS group_alpha_schema_service_db_role;

-- Create a warehouse for the app.
CREATE WAREHOUSE IF NOT EXISTS app_a_user_1_xs_warehouse
  WAREHOUSE_SIZE = 'XSMALL'
  AUTO_RESUME = TRUE
  INITIALLY_SUSPENDED = TRUE;

-- Do this because creating a warehouse such as above switches to that warehouse at least in the Snowflake VS Code extension.
use warehouse setup_xs_warehouse;

-- Potential frontend compute pools.
CREATE COMPUTE POOL IF NOT EXISTS app_a_user_1_frontend_1vcpu_6gib_1x_compute_pool
    MIN_NODES = 1
    MAX_NODES = 1
    INSTANCE_FAMILY = CPU_X64_XS
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 600;
CREATE COMPUTE POOL IF NOT EXISTS app_a_user_1_frontend_3vcpu_13gib_2x_compute_pool
    MIN_NODES = 1
    MAX_NODES = 1
    INSTANCE_FAMILY = CPU_X64_S
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 600;
CREATE COMPUTE POOL IF NOT EXISTS app_a_user_1_frontend_6vcpu_28gib_4x_compute_pool
    MIN_NODES = 1
    MAX_NODES = 1
    INSTANCE_FAMILY = CPU_X64_M
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 600;
CREATE COMPUTE POOL IF NOT EXISTS app_a_user_1_frontend_6vcpu_58gib_5x_compute_pool
    MIN_NODES = 1
    MAX_NODES = 1
    INSTANCE_FAMILY = HIGHMEM_X64_S
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 600;
CREATE COMPUTE POOL IF NOT EXISTS app_a_user_1_frontend_14vcpu_58gib_7x_compute_pool
    MIN_NODES = 1
    MAX_NODES = 1
    INSTANCE_FAMILY = CPU_X64_SL
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 600;
CREATE COMPUTE POOL IF NOT EXISTS app_a_user_1_frontend_28vcpu_116gib_14x_compute_pool
    MIN_NODES = 1
    MAX_NODES = 1
    INSTANCE_FAMILY = CPU_X64_L
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 600;
CREATE COMPUTE POOL IF NOT EXISTS app_a_user_1_frontend_28vcpu_240gib_19x_compute_pool
    MIN_NODES = 1
    MAX_NODES = 1
    INSTANCE_FAMILY = HIGHMEM_X64_M
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 600;

-- Potential worker compute pools.
CREATE COMPUTE POOL IF NOT EXISTS app_a_user_1_workers_1vcpu_6gib_1x_compute_pool
    MIN_NODES = 1
    MAX_NODES = 5
    INSTANCE_FAMILY = CPU_X64_XS
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 600;
CREATE COMPUTE POOL IF NOT EXISTS app_a_user_1_workers_3vcpu_13gib_2x_compute_pool
    MIN_NODES = 1
    MAX_NODES = 5
    INSTANCE_FAMILY = CPU_X64_S
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 600;
CREATE COMPUTE POOL IF NOT EXISTS app_a_user_1_workers_6vcpu_28gib_4x_compute_pool
    MIN_NODES = 1
    MAX_NODES = 5
    INSTANCE_FAMILY = CPU_X64_M
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 600;
CREATE COMPUTE POOL IF NOT EXISTS app_a_user_1_workers_6vcpu_58gib_5x_compute_pool
    MIN_NODES = 1
    MAX_NODES = 5
    INSTANCE_FAMILY = HIGHMEM_X64_S
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 600;
CREATE COMPUTE POOL IF NOT EXISTS app_a_user_1_workers_14vcpu_58gib_7x_compute_pool
    MIN_NODES = 1
    MAX_NODES = 5
    INSTANCE_FAMILY = CPU_X64_SL
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 600;
CREATE COMPUTE POOL IF NOT EXISTS app_a_user_1_workers_28vcpu_116gib_14x_compute_pool
    MIN_NODES = 1
    MAX_NODES = 5
    INSTANCE_FAMILY = CPU_X64_L
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 600;
CREATE COMPUTE POOL IF NOT EXISTS app_a_user_1_workers_28vcpu_240gib_19x_compute_pool
    MIN_NODES = 1
    MAX_NODES = 5
    INSTANCE_FAMILY = HIGHMEM_X64_M
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 600;

-- Grant app_a_group_alpha_role privileges to see the database and schema.
GRANT USAGE ON DATABASE app_a_app_db
  TO ROLE app_a_group_alpha_role;
GRANT USAGE ON SCHEMA group_alpha_schema
  TO ROLE app_a_group_alpha_role;
GRANT USAGE ON SCHEMA general_schema TO ROLE app_a_group_alpha_role;

-- Grant general_schema_ro_db_role privileges to see the database and schema.
GRANT USAGE ON DATABASE app_a_app_db
  TO DATABASE ROLE general_schema_ro_db_role;
GRANT USAGE ON SCHEMA general_schema
  TO DATABASE ROLE general_schema_ro_db_role;

-- Grant group_alpha_schema_service_db_role privileges to see the database and schema.
GRANT USAGE ON DATABASE app_a_app_db
  TO DATABASE ROLE group_alpha_schema_service_db_role;
GRANT USAGE ON SCHEMA group_alpha_schema
  TO DATABASE ROLE group_alpha_schema_service_db_role;

-- Give the app user role the ability to even launch the app by granting access to the database and schema.
GRANT USAGE ON DATABASE app_a_app_db TO ROLE data_apps_user_1_role;
GRANT USAGE ON SCHEMA group_alpha_schema TO ROLE data_apps_user_1_role;

-- Give this account role the appropriate database roles.
grant database role group_alpha_group_db.app_a_schema_rw_db_role to role app_a_group_alpha_role;
grant database role group_alpha_group_db.curated_schema_rw_db_role to role app_a_group_alpha_role;
grant database role common_db.admin_schema_ro_db_role to role app_a_group_alpha_role;
grant database role app_a_app_db.general_schema_ro_db_role to role app_a_group_alpha_role;
GRANT DATABASE ROLE app_a_app_db.group_alpha_schema_service_db_role TO ROLE app_a_group_alpha_role;

-- Give the app user role the required database role (used in launcher.py).
grant database role common_db.admin_schema_ro_db_role to role data_apps_user_1_role;

-- Grant access to using the compute resources for the actual "service" role app_a_group_alpha_role.
GRANT USAGE ON WAREHOUSE app_a_user_1_xs_warehouse TO ROLE app_a_group_alpha_role;
GRANT USAGE, OPERATE ON COMPUTE POOL app_a_user_1_frontend_1vcpu_6gib_1x_compute_pool TO ROLE app_a_group_alpha_role;
GRANT USAGE, OPERATE ON COMPUTE POOL app_a_user_1_frontend_3vcpu_13gib_2x_compute_pool TO ROLE app_a_group_alpha_role;
GRANT USAGE, OPERATE ON COMPUTE POOL app_a_user_1_frontend_6vcpu_28gib_4x_compute_pool TO ROLE app_a_group_alpha_role;
GRANT USAGE, OPERATE ON COMPUTE POOL app_a_user_1_frontend_6vcpu_58gib_5x_compute_pool TO ROLE app_a_group_alpha_role;
GRANT USAGE, OPERATE ON COMPUTE POOL app_a_user_1_frontend_14vcpu_58gib_7x_compute_pool TO ROLE app_a_group_alpha_role;
GRANT USAGE, OPERATE ON COMPUTE POOL app_a_user_1_frontend_28vcpu_116gib_14x_compute_pool TO ROLE app_a_group_alpha_role;
GRANT USAGE, OPERATE ON COMPUTE POOL app_a_user_1_frontend_28vcpu_240gib_19x_compute_pool TO ROLE app_a_group_alpha_role;
GRANT USAGE ON COMPUTE POOL app_a_user_1_workers_1vcpu_6gib_1x_compute_pool TO ROLE app_a_group_alpha_role;
GRANT USAGE ON COMPUTE POOL app_a_user_1_workers_3vcpu_13gib_2x_compute_pool TO ROLE app_a_group_alpha_role;
GRANT USAGE ON COMPUTE POOL app_a_user_1_workers_6vcpu_28gib_4x_compute_pool TO ROLE app_a_group_alpha_role;
GRANT USAGE ON COMPUTE POOL app_a_user_1_workers_6vcpu_58gib_5x_compute_pool TO ROLE app_a_group_alpha_role;
GRANT USAGE ON COMPUTE POOL app_a_user_1_workers_14vcpu_58gib_7x_compute_pool TO ROLE app_a_group_alpha_role;
GRANT USAGE ON COMPUTE POOL app_a_user_1_workers_28vcpu_116gib_14x_compute_pool TO ROLE app_a_group_alpha_role;
GRANT USAGE ON COMPUTE POOL app_a_user_1_workers_28vcpu_240gib_19x_compute_pool TO ROLE app_a_group_alpha_role;

-- Grant monitor and operate on the compute resources to the user role data_apps_user_1_role so that users can monitor and operate the app.
GRANT MONITOR, OPERATE ON WAREHOUSE app_a_user_1_xs_warehouse TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON COMPUTE POOL app_a_user_1_frontend_1vcpu_6gib_1x_compute_pool TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON COMPUTE POOL app_a_user_1_frontend_3vcpu_13gib_2x_compute_pool TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON COMPUTE POOL app_a_user_1_frontend_6vcpu_28gib_4x_compute_pool TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON COMPUTE POOL app_a_user_1_frontend_6vcpu_58gib_5x_compute_pool TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON COMPUTE POOL app_a_user_1_frontend_14vcpu_58gib_7x_compute_pool TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON COMPUTE POOL app_a_user_1_frontend_28vcpu_116gib_14x_compute_pool TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON COMPUTE POOL app_a_user_1_frontend_28vcpu_240gib_19x_compute_pool TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON COMPUTE POOL app_a_user_1_workers_1vcpu_6gib_1x_compute_pool TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON COMPUTE POOL app_a_user_1_workers_3vcpu_13gib_2x_compute_pool TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON COMPUTE POOL app_a_user_1_workers_6vcpu_28gib_4x_compute_pool TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON COMPUTE POOL app_a_user_1_workers_6vcpu_58gib_5x_compute_pool TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON COMPUTE POOL app_a_user_1_workers_14vcpu_58gib_7x_compute_pool TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON COMPUTE POOL app_a_user_1_workers_28vcpu_116gib_14x_compute_pool TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON COMPUTE POOL app_a_user_1_workers_28vcpu_240gib_19x_compute_pool TO ROLE data_apps_user_1_role;

-- Grant read permissions on the table.
-- I'm not sure why I originally had this; commenting it out for the time being.
-- GRANT SELECT
--   ON TABLE general_schema.image_metadata_table
--   TO DATABASE ROLE general_schema_ro_db_role;

-- Grant permissions to allow the database role to create serivces (such as for submitting a job service) and to use the images in the image repository for doing so.
GRANT CREATE SERVICE ON SCHEMA group_alpha_schema
  TO DATABASE ROLE group_alpha_schema_service_db_role;
GRANT READ ON IMAGE REPOSITORY general_schema.image_repository
  TO DATABASE ROLE group_alpha_schema_service_db_role;

-- Allow the app to read the worker spec from this stage so it can launch a worker service.
GRANT READ ON STAGE general_schema.general_stage TO DATABASE ROLE general_schema_ro_db_role;

-- Allow the app role to bind the service endpoint to the account so that the app can be accessed from the web.
GRANT BIND SERVICE ENDPOINT ON ACCOUNT TO ROLE app_a_group_alpha_role;

-- Allow the app role to perform setup.
GRANT USAGE ON WAREHOUSE setup_xs_warehouse TO ROLE app_a_group_alpha_role;

-- We want the app owner to be the app role so switch to that prior to creating the app.
GRANT ROLE app_a_group_alpha_role TO ROLE accountadmin;
USE ROLE app_a_group_alpha_role;
USE WAREHOUSE setup_xs_warehouse;

-- Create the seven services, one with each set of compute resources.
-- CPU_X64_XS_1vcpu_6gib_1x
-- Note in the USING blocks if there are underscores, hyphens, spaces, etc., you need something like ' "bleh" ' instead of 'bleh'.
DROP SERVICE IF EXISTS group_alpha_schema.app_a_user_1_frontend_1vcpu_6gib_1x_service;
CREATE SERVICE group_alpha_schema.app_a_user_1_frontend_1vcpu_6gib_1x_service
  IN COMPUTE POOL app_a_user_1_frontend_1vcpu_6gib_1x_compute_pool
  FROM @app_a_app_db.general_schema.general_stage SPECIFICATION_TEMPLATE_FILE='frontend_service_spec.yaml'
  USING ( APP_SHORTNAME=>'app_a', APP_TITLE=>' "App A" ', MONITOR_JOBS_REFRESH_INTERVAL_SECONDS=>5, SNOWFLAKE_USER=>' "user_1" ', COMPUTE_RESOURCE=>' "1vcpu_6gib_1x" ', ALL_COMPUTE_RESOURCES=>' "1vcpu_6gib_1x 3vcpu_13gib_2x 6vcpu_28gib_4x 6vcpu_58gib_5x 14vcpu_58gib_7x 28vcpu_116gib_14x 28vcpu_240gib_19x" ', IMAGE=>' "/app_a_app_db/general_schema/image_repository/frontend:latest" ', SNOWFLAKE_WAREHOUSE=>' "app_a_user_1_xs_warehouse" ', MOUNTPATH=>' "/tmp/app_a" ', MEMORY=>'6Gi', CPU=>1, IMAGE_NAME=>' "frontend" ', IMAGE_TAG=>' "latest" ' )
  AUTO_RESUME = FALSE
  MIN_INSTANCES = 1
  MAX_INSTANCES = 1;
alter service group_alpha_schema.app_a_user_1_frontend_1vcpu_6gib_1x_service suspend;
alter compute pool app_a_user_1_frontend_1vcpu_6gib_1x_compute_pool suspend;

-- CPU_X64_S_3vcpu_13gib_2x
DROP SERVICE IF EXISTS group_alpha_schema.app_a_user_1_frontend_3vcpu_13gib_2x_service;
CREATE SERVICE group_alpha_schema.app_a_user_1_frontend_3vcpu_13gib_2x_service
  IN COMPUTE POOL app_a_user_1_frontend_3vcpu_13gib_2x_compute_pool
  FROM @app_a_app_db.general_schema.general_stage SPECIFICATION_TEMPLATE_FILE='frontend_service_spec.yaml'
  USING ( APP_SHORTNAME=>'app_a', APP_TITLE=>' "App A" ', MONITOR_JOBS_REFRESH_INTERVAL_SECONDS=>5, SNOWFLAKE_USER=>' "user_1" ', COMPUTE_RESOURCE=>' "3vcpu_13gib_2x" ', ALL_COMPUTE_RESOURCES=>' "1vcpu_6gib_1x 3vcpu_13gib_2x 6vcpu_28gib_4x 6vcpu_58gib_5x 14vcpu_58gib_7x 28vcpu_116gib_14x 28vcpu_240gib_19x" ', IMAGE=>' "/app_a_app_db/general_schema/image_repository/frontend:latest" ', SNOWFLAKE_WAREHOUSE=>' "app_a_user_1_xs_warehouse" ', MOUNTPATH=>' "/tmp/app_a" ', MEMORY=>'13Gi', CPU=>3, IMAGE_NAME=>' "frontend" ', IMAGE_TAG=>' "latest" ' )
  AUTO_RESUME = FALSE
  MIN_INSTANCES = 1
  MAX_INSTANCES = 1;
alter service group_alpha_schema.app_a_user_1_frontend_3vcpu_13gib_2x_service suspend;
alter compute pool app_a_user_1_frontend_3vcpu_13gib_2x_compute_pool suspend;

-- CPU_X64_M_6vcpu_28gib_4x
DROP SERVICE IF EXISTS group_alpha_schema.app_a_user_1_frontend_6vcpu_28gib_4x_service;
CREATE SERVICE group_alpha_schema.app_a_user_1_frontend_6vcpu_28gib_4x_service
  IN COMPUTE POOL app_a_user_1_frontend_6vcpu_28gib_4x_compute_pool
  FROM @app_a_app_db.general_schema.general_stage SPECIFICATION_TEMPLATE_FILE='frontend_service_spec.yaml'
  USING ( APP_SHORTNAME=>'app_a', APP_TITLE=>' "App A" ', MONITOR_JOBS_REFRESH_INTERVAL_SECONDS=>5, SNOWFLAKE_USER=>' "user_1" ', COMPUTE_RESOURCE=>' "6vcpu_28gib_4x" ', ALL_COMPUTE_RESOURCES=>' "1vcpu_6gib_1x 3vcpu_13gib_2x 6vcpu_28gib_4x 6vcpu_58gib_5x 14vcpu_58gib_7x 28vcpu_116gib_14x 28vcpu_240gib_19x" ', IMAGE=>' "/app_a_app_db/general_schema/image_repository/frontend:latest" ', SNOWFLAKE_WAREHOUSE=>' "app_a_user_1_xs_warehouse" ', MOUNTPATH=>' "/tmp/app_a" ', MEMORY=>'28Gi', CPU=>6, IMAGE_NAME=>' "frontend" ', IMAGE_TAG=>' "latest" ' )
  AUTO_RESUME = FALSE
  MIN_INSTANCES = 1
  MAX_INSTANCES = 1;
alter service group_alpha_schema.app_a_user_1_frontend_6vcpu_28gib_4x_service suspend;
alter compute pool app_a_user_1_frontend_6vcpu_28gib_4x_compute_pool suspend;

-- HIGHMEM_X64_S_6vcpu_58gib_5x
DROP SERVICE IF EXISTS group_alpha_schema.app_a_user_1_frontend_6vcpu_58gib_5x_service;
CREATE SERVICE group_alpha_schema.app_a_user_1_frontend_6vcpu_58gib_5x_service
  IN COMPUTE POOL app_a_user_1_frontend_6vcpu_58gib_5x_compute_pool
  FROM @app_a_app_db.general_schema.general_stage SPECIFICATION_TEMPLATE_FILE='frontend_service_spec.yaml'
  USING ( APP_SHORTNAME=>'app_a', APP_TITLE=>' "App A" ', MONITOR_JOBS_REFRESH_INTERVAL_SECONDS=>5, SNOWFLAKE_USER=>' "user_1" ', COMPUTE_RESOURCE=>' "6vcpu_58gib_5x" ', ALL_COMPUTE_RESOURCES=>' "1vcpu_6gib_1x 3vcpu_13gib_2x 6vcpu_28gib_4x 6vcpu_58gib_5x 14vcpu_58gib_7x 28vcpu_116gib_14x 28vcpu_240gib_19x" ', IMAGE=>' "/app_a_app_db/general_schema/image_repository/frontend:latest" ', SNOWFLAKE_WAREHOUSE=>' "app_a_user_1_xs_warehouse" ', MOUNTPATH=>' "/tmp/app_a" ', MEMORY=>'58Gi', CPU=>6, IMAGE_NAME=>' "frontend" ', IMAGE_TAG=>' "latest" ' )
  AUTO_RESUME = FALSE
  MIN_INSTANCES = 1
  MAX_INSTANCES = 1;
alter service group_alpha_schema.app_a_user_1_frontend_6vcpu_58gib_5x_service suspend;
alter compute pool app_a_user_1_frontend_6vcpu_58gib_5x_compute_pool suspend;

-- CPU_X64_SL_14vcpu_58gib_7x
DROP SERVICE IF EXISTS group_alpha_schema.app_a_user_1_frontend_14vcpu_58gib_7x_service;
CREATE SERVICE group_alpha_schema.app_a_user_1_frontend_14vcpu_58gib_7x_service
  IN COMPUTE POOL app_a_user_1_frontend_14vcpu_58gib_7x_compute_pool
  FROM @app_a_app_db.general_schema.general_stage SPECIFICATION_TEMPLATE_FILE='frontend_service_spec.yaml'
  USING ( APP_SHORTNAME=>'app_a', APP_TITLE=>' "App A" ', MONITOR_JOBS_REFRESH_INTERVAL_SECONDS=>5, SNOWFLAKE_USER=>' "user_1" ', COMPUTE_RESOURCE=>' "14vcpu_58gib_7x" ', ALL_COMPUTE_RESOURCES=>' "1vcpu_6gib_1x 3vcpu_13gib_2x 6vcpu_28gib_4x 6vcpu_58gib_5x 14vcpu_58gib_7x 28vcpu_116gib_14x 28vcpu_240gib_19x" ', IMAGE=>' "/app_a_app_db/general_schema/image_repository/frontend:latest" ', SNOWFLAKE_WAREHOUSE=>' "app_a_user_1_xs_warehouse" ', MOUNTPATH=>' "/tmp/app_a" ', MEMORY=>'58Gi', CPU=>14, IMAGE_NAME=>' "frontend" ', IMAGE_TAG=>' "latest" ' )
  AUTO_RESUME = FALSE
  MIN_INSTANCES = 1
  MAX_INSTANCES = 1;
alter service group_alpha_schema.app_a_user_1_frontend_14vcpu_58gib_7x_service suspend;
alter compute pool app_a_user_1_frontend_14vcpu_58gib_7x_compute_pool suspend;

-- CPU_X64_L_28vcpu_116gib_14x
DROP SERVICE IF EXISTS group_alpha_schema.app_a_user_1_frontend_28vcpu_116gib_14x_service;
CREATE SERVICE group_alpha_schema.app_a_user_1_frontend_28vcpu_116gib_14x_service
  IN COMPUTE POOL app_a_user_1_frontend_28vcpu_116gib_14x_compute_pool
  FROM @app_a_app_db.general_schema.general_stage SPECIFICATION_TEMPLATE_FILE='frontend_service_spec.yaml'
  USING ( APP_SHORTNAME=>'app_a', APP_TITLE=>' "App A" ', MONITOR_JOBS_REFRESH_INTERVAL_SECONDS=>5, SNOWFLAKE_USER=>' "user_1" ', COMPUTE_RESOURCE=>' "28vcpu_116gib_14x" ', ALL_COMPUTE_RESOURCES=>' "1vcpu_6gib_1x 3vcpu_13gib_2x 6vcpu_28gib_4x 6vcpu_58gib_5x 14vcpu_58gib_7x 28vcpu_116gib_14x 28vcpu_240gib_19x" ', IMAGE=>' "/app_a_app_db/general_schema/image_repository/frontend:latest" ', SNOWFLAKE_WAREHOUSE=>' "app_a_user_1_xs_warehouse" ', MOUNTPATH=>' "/tmp/app_a" ', MEMORY=>'116Gi', CPU=>28, IMAGE_NAME=>' "frontend" ', IMAGE_TAG=>' "latest" ' )
  AUTO_RESUME = FALSE
  MIN_INSTANCES = 1
  MAX_INSTANCES = 1;
alter service group_alpha_schema.app_a_user_1_frontend_28vcpu_116gib_14x_service suspend;
alter compute pool app_a_user_1_frontend_28vcpu_116gib_14x_compute_pool suspend;

-- HIGHMEM_X64_M_28vcpu_240gib_19x
DROP SERVICE IF EXISTS group_alpha_schema.app_a_user_1_frontend_28vcpu_240gib_19x_service;
CREATE SERVICE group_alpha_schema.app_a_user_1_frontend_28vcpu_240gib_19x_service
  IN COMPUTE POOL app_a_user_1_frontend_28vcpu_240gib_19x_compute_pool
  FROM @app_a_app_db.general_schema.general_stage SPECIFICATION_TEMPLATE_FILE='frontend_service_spec.yaml'
  USING ( APP_SHORTNAME=>'app_a', APP_TITLE=>' "App A" ', MONITOR_JOBS_REFRESH_INTERVAL_SECONDS=>5, SNOWFLAKE_USER=>' "user_1" ', COMPUTE_RESOURCE=>' "28vcpu_240gib_19x" ', ALL_COMPUTE_RESOURCES=>' "1vcpu_6gib_1x 3vcpu_13gib_2x 6vcpu_28gib_4x 6vcpu_58gib_5x 14vcpu_58gib_7x 28vcpu_116gib_14x 28vcpu_240gib_19x" ', IMAGE=>' "/app_a_app_db/general_schema/image_repository/frontend:latest" ', SNOWFLAKE_WAREHOUSE=>' "app_a_user_1_xs_warehouse" ', MOUNTPATH=>' "/tmp/app_a" ', MEMORY=>'240Gi', CPU=>28, IMAGE_NAME=>' "frontend" ', IMAGE_TAG=>' "latest" ' )
  AUTO_RESUME = FALSE
  MIN_INSTANCES = 1
  MAX_INSTANCES = 1;
alter service group_alpha_schema.app_a_user_1_frontend_28vcpu_240gib_19x_service suspend;
alter compute pool app_a_user_1_frontend_28vcpu_240gib_19x_compute_pool suspend;

-- Switch back to accountadmin role.
USE ROLE accountadmin;

-- Allow the app user to see and operate the frontend services.
GRANT MONITOR, OPERATE ON SERVICE group_alpha_schema.app_a_user_1_frontend_1vcpu_6gib_1x_service TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON SERVICE group_alpha_schema.app_a_user_1_frontend_3vcpu_13gib_2x_service TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON SERVICE group_alpha_schema.app_a_user_1_frontend_6vcpu_28gib_4x_service TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON SERVICE group_alpha_schema.app_a_user_1_frontend_6vcpu_58gib_5x_service TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON SERVICE group_alpha_schema.app_a_user_1_frontend_14vcpu_58gib_7x_service TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON SERVICE group_alpha_schema.app_a_user_1_frontend_28vcpu_116gib_14x_service TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON SERVICE group_alpha_schema.app_a_user_1_frontend_28vcpu_240gib_19x_service TO ROLE data_apps_user_1_role;

-- Allow the user to run the app from the web even though they have no access to the role that runs the app.
GRANT SERVICE ROLE group_alpha_schema.app_a_user_1_frontend_1vcpu_6gib_1x_service!web_endpoint_service_role TO ROLE data_apps_user_1_role;
GRANT SERVICE ROLE group_alpha_schema.app_a_user_1_frontend_3vcpu_13gib_2x_service!web_endpoint_service_role TO ROLE data_apps_user_1_role;
GRANT SERVICE ROLE group_alpha_schema.app_a_user_1_frontend_6vcpu_28gib_4x_service!web_endpoint_service_role TO ROLE data_apps_user_1_role;
GRANT SERVICE ROLE group_alpha_schema.app_a_user_1_frontend_6vcpu_58gib_5x_service!web_endpoint_service_role TO ROLE data_apps_user_1_role;
GRANT SERVICE ROLE group_alpha_schema.app_a_user_1_frontend_14vcpu_58gib_7x_service!web_endpoint_service_role TO ROLE data_apps_user_1_role;
GRANT SERVICE ROLE group_alpha_schema.app_a_user_1_frontend_28vcpu_116gib_14x_service!web_endpoint_service_role TO ROLE data_apps_user_1_role;
GRANT SERVICE ROLE group_alpha_schema.app_a_user_1_frontend_28vcpu_240gib_19x_service!web_endpoint_service_role TO ROLE data_apps_user_1_role;
---------------- End database app_a_app_db. -----------------------------------------------


---------------- Database app_launcher_db. ---------------------------------------------------
-- Create the database and schemas.
create database if not exists app_launcher_db;
use database app_launcher_db;
create schema if not exists group_alpha_schema;
create schema if not exists general_schema;

-- Create a general stage for holding the code for the launcher.
create stage if not exists general_schema.general_stage
  directory = ( enable = true );

-- SEE GITHUB README FOR WHAT FILES TO UPLOAD TO THIS STAGE (the streamlit app `launcher.py` and the environment `environment.yml`).

-- Create a warehouse for the app.
CREATE WAREHOUSE IF NOT EXISTS app_launcher_user_1_xs_warehouse
  WAREHOUSE_SIZE = 'XSMALL'
  AUTO_RESUME = TRUE
  INITIALLY_SUSPENDED = TRUE;

-- Do this because creating a warehouse such as above switches to that warehouse at least in the Snowflake VS Code extension.
use warehouse setup_xs_warehouse;

-- Grant data_apps_user_1_role privileges to see the database and schemas.
GRANT USAGE ON DATABASE app_launcher_db
  TO ROLE data_apps_user_1_role;
GRANT USAGE ON SCHEMA group_alpha_schema
  TO ROLE data_apps_user_1_role;
GRANT USAGE ON SCHEMA general_schema
  TO ROLE data_apps_user_1_role;

-- Grant access to the warehouse.
GRANT USAGE, OPERATE, MONITOR ON WAREHOUSE app_launcher_user_1_xs_warehouse TO ROLE data_apps_user_1_role;

-- Allow creation of Streamlit apps in the group_alpha_schema schema.
GRANT CREATE STREAMLIT ON SCHEMA group_alpha_schema TO ROLE data_apps_user_1_role;

-- Allow reading from the stage where the launcher code is stored.
GRANT READ ON STAGE general_schema.general_stage TO ROLE data_apps_user_1_role;

-- Allow the app user to use the setup warehouse to create the app.
GRANT USAGE ON WAREHOUSE setup_xs_warehouse TO ROLE data_apps_user_1_role;

-- We want the app owner to be the app user role so switch to that prior to creating the app.
GRANT ROLE data_apps_user_1_role TO ROLE accountadmin;
USE ROLE data_apps_user_1_role;
USE WAREHOUSE setup_xs_warehouse;

-- Create the Streamlit app.
CREATE OR REPLACE STREAMLIT group_alpha_schema.app_launcher_user_1_streamlit
  FROM @app_launcher_db.general_schema.general_stage
  MAIN_FILE = 'launcher.py'
  QUERY_WAREHOUSE = app_launcher_user_1_xs_warehouse
  TITLE = 'App Launcher (group_alpha)';

-- Switch back to accountadmin role.
USE ROLE accountadmin;

-- Revoke temporarily granted roles.
REVOKE CREATE STREAMLIT ON SCHEMA group_alpha_schema FROM ROLE data_apps_user_1_role;
REVOKE READ ON STAGE general_schema.general_stage FROM ROLE data_apps_user_1_role;
REVOKE USAGE ON WAREHOUSE setup_xs_warehouse FROM ROLE data_apps_user_1_role;
---------------- End database app_launcher_db. -----------------------------------------------


---------------- Database dmgr_db. ---------------------------------------------------

-- Create the database and schemas.
create database if not exists dmgr_db;
use database dmgr_db;
create schema if not exists group_alpha_schema;
create schema if not exists general_schema;

-- Create a general stage for holding the spec for the service.
create stage if not exists general_schema.general_stage
  directory = ( enable = true );

-- SEE GITHUB README FOR WHAT FILE TO UPLOAD TO THIS STAGE (not yet created nor actually present in the README).

-- Create an image repository.
CREATE IMAGE REPOSITORY IF NOT EXISTS general_schema.image_repository;

-- SEE GITHUB README FOR WHAT IMAGE TO UPLOAD TO THIS REPOSITORY (not yet created nor actually present in the README).

-- Create table.
CREATE TABLE IF NOT EXISTS general_schema.image_metadata_table (
  id INTEGER IDENTITY PRIMARY KEY,
  image_id VARCHAR(255) UNIQUE NOT NULL,
  name VARCHAR(255),
  tag VARCHAR(255),
  git_commit VARCHAR(255),
  environment_yaml_file VARCHAR(255),
  image_added_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
  who_added VARCHAR(255)
);

-- SEE GITHUB README FOR WHAT DATA TO ADD TO THIS TABLE (not yet actually present in the README).

-- Create roles.
create role if not exists dmgr_group_alpha_role; -- This is the account role that owns and operates the app.

-- Create a warehouse for the app.
CREATE WAREHOUSE IF NOT EXISTS dmgr_user_1_xs_warehouse
  WAREHOUSE_SIZE = 'XSMALL'
  AUTO_RESUME = TRUE
  INITIALLY_SUSPENDED = TRUE;

-- Do this because creating a warehouse such as above switches to that warehouse at least in the Snowflake VS Code extension.
use warehouse setup_xs_warehouse;

-- Create compute pool for the app.
CREATE COMPUTE POOL IF NOT EXISTS dmgr_user_1_xs_compute_pool
    MIN_NODES = 1
    MAX_NODES = 1
    INSTANCE_FAMILY = CPU_X64_XS
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 600;

-- Grant dmgr_group_alpha_role privileges to see the database and schema.
GRANT USAGE ON DATABASE dmgr_db
  TO ROLE dmgr_group_alpha_role;
GRANT USAGE ON SCHEMA group_alpha_schema
  TO ROLE dmgr_group_alpha_role;
GRANT USAGE ON SCHEMA general_schema
  TO ROLE dmgr_group_alpha_role;

-- Give the app user role the ability to even launch the app by granting access to the database and schema.
GRANT USAGE ON DATABASE dmgr_db TO ROLE data_apps_user_1_role;
GRANT USAGE ON SCHEMA group_alpha_schema TO ROLE data_apps_user_1_role;

-- Give this account role the appropriate database roles.
grant database role group_alpha_group_db.curated_schema_rw_db_role to role dmgr_group_alpha_role;
grant database role common_db.admin_schema_ro_db_role to role dmgr_group_alpha_role;

-- Grant access to the compute resources.
GRANT USAGE ON WAREHOUSE dmgr_user_1_xs_warehouse TO ROLE dmgr_group_alpha_role;
GRANT USAGE ON COMPUTE POOL dmgr_user_1_xs_compute_pool TO ROLE dmgr_group_alpha_role;

-- Grant resource management to the user.
GRANT MONITOR, OPERATE ON WAREHOUSE dmgr_user_1_xs_warehouse TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON COMPUTE POOL dmgr_user_1_xs_compute_pool TO ROLE data_apps_user_1_role;

-- Not sure why this wasn't here already.
GRANT MONITOR, OPERATE ON COMPUTE POOL dmgr_user_1_xs_compute_pool TO ROLE dmgr_group_alpha_role;

-- Allow the app role to create the service.
GRANT CREATE SERVICE ON SCHEMA group_alpha_schema TO ROLE dmgr_group_alpha_role;
GRANT READ ON IMAGE REPOSITORY general_schema.image_repository TO ROLE dmgr_group_alpha_role;
GRANT READ ON STAGE general_schema.general_stage TO ROLE dmgr_group_alpha_role;
GRANT BIND SERVICE ENDPOINT ON ACCOUNT TO ROLE dmgr_group_alpha_role;

-- Allow the app role to perform setup.
GRANT USAGE ON WAREHOUSE setup_xs_warehouse TO ROLE dmgr_group_alpha_role;

-- We want the app owner to be the app role so switch to that prior to creating the app.
GRANT ROLE dmgr_group_alpha_role TO ROLE accountadmin;
USE ROLE dmgr_group_alpha_role;
USE WAREHOUSE setup_xs_warehouse;

-- Create the app.
DROP SERVICE IF EXISTS group_alpha_schema.dmgr_user_1_xs_service;
CREATE SERVICE group_alpha_schema.dmgr_user_1_xs_service
  IN COMPUTE POOL dmgr_user_1_xs_compute_pool
  FROM @dmgr_db.general_schema.general_stage SPECIFICATION_TEMPLATE_FILE='snowflake_service_spec.yaml'
  USING ( APP_SHORTNAME=>'dmgr', APP_TITLE=>' "Data Manager" ', SNOWFLAKE_USER=>' "user_1" ', COMPUTE_RESOURCE=>' "1vcpu_6gib_1x" ', ALL_COMPUTE_RESOURCES=>' "1vcpu_6gib_1x 3vcpu_13gib_2x 6vcpu_28gib_4x 6vcpu_58gib_5x 14vcpu_58gib_7x 28vcpu_116gib_14x 28vcpu_240gib_19x" ', IMAGE=>' "/dmgr_db/general_schema/image_repository/data-manager:latest" ', SNOWFLAKE_WAREHOUSE=>' "dmgr_user_1_xs_warehouse" ', MOUNTPATH=>' "/tmp/dmgr" ', MEMORY=>'6Gi', CPU=>1, IMAGE_NAME=>' "data-manager" ', IMAGE_TAG=>' "latest" ' )
  AUTO_RESUME = FALSE
  MIN_INSTANCES = 1
  MAX_INSTANCES = 1;
alter service group_alpha_schema.dmgr_user_1_xs_service suspend;
alter compute pool dmgr_user_1_xs_compute_pool suspend;

-- Switch back to accountadmin role.
USE ROLE accountadmin;

-- Allow the app user to see and operate the service.
GRANT MONITOR, OPERATE ON SERVICE group_alpha_schema.dmgr_user_1_xs_service TO ROLE data_apps_user_1_role;

-- Allow the user to run the app from the web even though they have no access to the role that runs the app.
GRANT SERVICE ROLE group_alpha_schema.dmgr_user_1_xs_service!web_endpoint_service_role TO ROLE data_apps_user_1_role;

-- Allow app to read/write data to the required stages.
GRANT USAGE ON DATABASE group_alpha_group_db  -- Allow using the database
  TO ROLE dmgr_group_alpha_role;
GRANT USAGE ON SCHEMA group_alpha_group_db.curated_schema  -- Allow using the schema
  TO ROLE dmgr_group_alpha_role;
GRANT USAGE ON SCHEMA group_alpha_group_db.app_a_schema  -- Allow using the schema
  TO ROLE dmgr_group_alpha_role;
GRANT READ ON STAGE group_alpha_group_db.app_a_schema.oldarchives_stage  -- Allow reading files from the stage (LIST/GET, COPY INTO <table> FROM @stage)
  TO ROLE dmgr_group_alpha_role;
GRANT WRITE ON STAGE group_alpha_group_db.app_a_schema.oldarchives_stage  -- Allow writing files to the stage (PUT/REMOVE, COPY INTO @stage)
  TO ROLE dmgr_group_alpha_role;
---------------- End database dmgr_db. -----------------------------------------------


-- Assign the apps user role to the user. This is the one place (the argument of USER) that the real Snowflake username must be used. Other instances of "user_1" can be anything, as long as they have an entry in the user_groups table so we know which group they should be accessing. E.g., user_1_alpha should correspond to the group_alpha group and user_1_beta should correspond to the group_beta group in the user_groups table. Then this script will create e.g. (1) data_apps_user_1_alpha_role and assign it to user_1 and (2) data_apps_user_1_beta_role and assign it to user_1. Then, user_1 in Snowsight can select either role to access the app/data for either group.
GRANT ROLE data_apps_user_1_role TO USER user_1;
