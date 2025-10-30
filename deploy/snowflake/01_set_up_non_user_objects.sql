-- Become admin.
use role accountadmin;

-- Create a warehouse for setup purposes.
CREATE WAREHOUSE IF NOT EXISTS setup_xs_warehouse
  WAREHOUSE_SIZE = 'XSMALL'
  AUTO_RESUME = TRUE
  INITIALLY_SUSPENDED = TRUE;
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

-- Create database roles.
create database role if not exists curated_schema_rw_db_role;
create database role if not exists app_a_schema_rw_db_role;

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

-- Create an image repository.
CREATE IMAGE REPOSITORY IF NOT EXISTS general_schema.image_repository;

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

-- Create roles.
create role if not exists app_a_group_alpha_role; -- This is the account role that owns and operates the app.
create role if not exists data_apps_user_1_role; -- This is the account role that will be granted the service role for the app.
create database role if not exists general_schema_ro_db_role;
CREATE DATABASE ROLE IF NOT EXISTS group_alpha_schema_service_db_role;

-- Create a warehouse for the app.
CREATE WAREHOUSE IF NOT EXISTS app_a_user_1_xs_warehouse
  WAREHOUSE_SIZE = 'XSMALL'
  AUTO_RESUME = TRUE
  INITIALLY_SUSPENDED = TRUE;

-- Potential frontend compute pools.
CREATE COMPUTE POOL IF NOT EXISTS app_a_user_1_frontend_CPU_X64_XS_1vcpu_6gib_1x_compute_pool
    MIN_NODES = 1
    MAX_NODES = 1
    INSTANCE_FAMILY = CPU_X64_XS
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 600;
CREATE COMPUTE POOL IF NOT EXISTS app_a_user_1_frontend_CPU_X64_S_3vcpu_13gib_2x_compute_pool
    MIN_NODES = 1
    MAX_NODES = 1
    INSTANCE_FAMILY = CPU_X64_S
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 600;
CREATE COMPUTE POOL IF NOT EXISTS app_a_user_1_frontend_CPU_X64_M_6vcpu_28gib_4x_compute_pool
    MIN_NODES = 1
    MAX_NODES = 1
    INSTANCE_FAMILY = CPU_X64_M
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 600;
CREATE COMPUTE POOL IF NOT EXISTS app_a_user_1_frontend_HIGHMEM_X64_S_6vcpu_58gib_5x_compute_pool
    MIN_NODES = 1
    MAX_NODES = 1
    INSTANCE_FAMILY = HIGHMEM_X64_S
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 600;
CREATE COMPUTE POOL IF NOT EXISTS app_a_user_1_frontend_CPU_X64_SL_14vcpu_58gib_7x_compute_pool
    MIN_NODES = 1
    MAX_NODES = 1
    INSTANCE_FAMILY = CPU_X64_SL
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 600;
CREATE COMPUTE POOL IF NOT EXISTS app_a_user_1_frontend_CPU_X64_L_28vcpu_116gib_14x_compute_pool
    MIN_NODES = 1
    MAX_NODES = 1
    INSTANCE_FAMILY = CPU_X64_L
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 600;
CREATE COMPUTE POOL IF NOT EXISTS app_a_user_1_frontend_HIGHMEM_X64_M_28vcpu_240gib_19x_compute_pool
    MIN_NODES = 1
    MAX_NODES = 1
    INSTANCE_FAMILY = HIGHMEM_X64_M
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 600;

-- Potential worker compute pools.
CREATE COMPUTE POOL IF NOT EXISTS app_a_user_1_workers_CPU_X64_XS_1vcpu_6gib_1x_compute_pool
    MIN_NODES = 1
    MAX_NODES = 5
    INSTANCE_FAMILY = CPU_X64_XS
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 600;
CREATE COMPUTE POOL IF NOT EXISTS app_a_user_1_workers_CPU_X64_S_3vcpu_13gib_2x_compute_pool
    MIN_NODES = 1
    MAX_NODES = 5
    INSTANCE_FAMILY = CPU_X64_S
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 600;
CREATE COMPUTE POOL IF NOT EXISTS app_a_user_1_workers_CPU_X64_M_6vcpu_28gib_4x_compute_pool
    MIN_NODES = 1
    MAX_NODES = 5
    INSTANCE_FAMILY = CPU_X64_M
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 600;
CREATE COMPUTE POOL IF NOT EXISTS app_a_user_1_workers_HIGHMEM_X64_S_6vcpu_58gib_5x_compute_pool
    MIN_NODES = 1
    MAX_NODES = 5
    INSTANCE_FAMILY = HIGHMEM_X64_S
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 600;
CREATE COMPUTE POOL IF NOT EXISTS app_a_user_1_workers_CPU_X64_SL_14vcpu_58gib_7x_compute_pool
    MIN_NODES = 1
    MAX_NODES = 5
    INSTANCE_FAMILY = CPU_X64_SL
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 600;
CREATE COMPUTE POOL IF NOT EXISTS app_a_user_1_workers_CPU_X64_L_28vcpu_116gib_14x_compute_pool
    MIN_NODES = 1
    MAX_NODES = 5
    INSTANCE_FAMILY = CPU_X64_L
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 600;
CREATE COMPUTE POOL IF NOT EXISTS app_a_user_1_workers_HIGHMEM_X64_M_28vcpu_240gib_19x_compute_pool
    MIN_NODES = 1
    MAX_NODES = 5
    INSTANCE_FAMILY = HIGHMEM_X64_M
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 600;

-- **** CREATE THE 7 APPS like:
--   group_alpha_schema.app_a_user_1_frontend_CPU_X64_XS_1vcpu_6gib_1x_service
--   group_alpha_schema.app_a_user_1_frontend_CPU_X64_S_3vcpu_13gib_2x_service
--   group_alpha_schema.app_a_user_1_frontend_CPU_X64_M_6vcpu_28gib_4x_service
--   group_alpha_schema.app_a_user_1_frontend_HIGHMEM_X64_S_6vcpu_58gib_5x_service
--   group_alpha_schema.app_a_user_1_frontend_CPU_X64_SL_14vcpu_58gib_7x_service
--   group_alpha_schema.app_a_user_1_frontend_CPU_X64_L_28vcpu_116gib_14x_service
--   group_alpha_schema.app_a_user_1_frontend_HIGHMEM_X64_M_28vcpu_240gib_19x_service

-- Grant app_a_group_alpha_role privileges to see the database and schema.
GRANT USAGE ON DATABASE app_a_app_db
  TO ROLE app_a_group_alpha_role;
GRANT USAGE ON SCHEMA group_alpha_schema
  TO ROLE app_a_group_alpha_role;

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

-- This probably isn't needed since accountadmin is the one who handles the endpoints, but putting it here to keep it in mind.
-- GRANT BIND SERVICE ENDPOINT ON ACCOUNT TO ROLE app_a_group_alpha_role;

-- Change the owner of the apps to the service role app_a_group_alpha_role.
GRANT OWNERSHIP ON SERVICE group_alpha_schema.app_a_user_1_frontend_CPU_X64_XS_1vcpu_6gib_1x_service TO ROLE app_a_group_alpha_role COPY CURRENT GRANTS;
GRANT OWNERSHIP ON SERVICE group_alpha_schema.app_a_user_1_frontend_CPU_X64_S_3vcpu_13gib_2x_service TO ROLE app_a_group_alpha_role COPY CURRENT GRANTS;
GRANT OWNERSHIP ON SERVICE group_alpha_schema.app_a_user_1_frontend_CPU_X64_M_6vcpu_28gib_4x_service TO ROLE app_a_group_alpha_role COPY CURRENT GRANTS;
GRANT OWNERSHIP ON SERVICE group_alpha_schema.app_a_user_1_frontend_HIGHMEM_X64_S_6vcpu_58gib_5x_service TO ROLE app_a_group_alpha_role COPY CURRENT GRANTS;
GRANT OWNERSHIP ON SERVICE group_alpha_schema.app_a_user_1_frontend_CPU_X64_SL_14vcpu_58gib_7x_service TO ROLE app_a_group_alpha_role COPY CURRENT GRANTS;
GRANT OWNERSHIP ON SERVICE group_alpha_schema.app_a_user_1_frontend_CPU_X64_L_28vcpu_116gib_14x_service TO ROLE app_a_group_alpha_role COPY CURRENT GRANTS;
GRANT OWNERSHIP ON SERVICE group_alpha_schema.app_a_user_1_frontend_HIGHMEM_X64_M_28vcpu_240gib_19x_service TO ROLE app_a_group_alpha_role COPY CURRENT GRANTS;

-- Give this account role the appropriate database roles.
grant database role group_alpha_group_db.app_a_schema_rw_db_role to role app_a_group_alpha_role;
grant database role group_alpha_group_db.curated_schema_rw_db_role to role app_a_group_alpha_role;
grant database role common_db.admin_schema_ro_db_role to role app_a_group_alpha_role;
grant database role app_a_app_db.general_schema_ro_db_role to role app_a_group_alpha_role;
GRANT DATABASE ROLE app_a_app_db.group_alpha_schema_service_db_role TO ROLE app_a_group_alpha_role;

-- Grant access to using the compute resources for the actual "service" role app_a_group_alpha_role.
GRANT USAGE ON WAREHOUSE app_a_user_1_xs_warehouse TO ROLE app_a_group_alpha_role;
GRANT USAGE ON COMPUTE POOL app_a_user_1_frontend_CPU_X64_XS_1vcpu_6gib_1x_compute_pool TO ROLE app_a_group_alpha_role;
GRANT USAGE ON COMPUTE POOL app_a_user_1_frontend_CPU_X64_S_3vcpu_13gib_2x_compute_pool TO ROLE app_a_group_alpha_role;
GRANT USAGE ON COMPUTE POOL app_a_user_1_frontend_CPU_X64_M_6vcpu_28gib_4x_compute_pool TO ROLE app_a_group_alpha_role;
GRANT USAGE ON COMPUTE POOL app_a_user_1_frontend_HIGHMEM_X64_S_6vcpu_58gib_5x_compute_pool TO ROLE app_a_group_alpha_role;
GRANT USAGE ON COMPUTE POOL app_a_user_1_frontend_CPU_X64_SL_14vcpu_58gib_7x_compute_pool TO ROLE app_a_group_alpha_role;
GRANT USAGE ON COMPUTE POOL app_a_user_1_frontend_CPU_X64_L_28vcpu_116gib_14x_compute_pool TO ROLE app_a_group_alpha_role;
GRANT USAGE ON COMPUTE POOL app_a_user_1_frontend_HIGHMEM_X64_M_28vcpu_240gib_19x_compute_pool TO ROLE app_a_group_alpha_role;
GRANT USAGE ON COMPUTE POOL app_a_user_1_workers_CPU_X64_XS_1vcpu_6gib_1x_compute_pool TO ROLE app_a_group_alpha_role;
GRANT USAGE ON COMPUTE POOL app_a_user_1_workers_CPU_X64_S_3vcpu_13gib_2x_compute_pool TO ROLE app_a_group_alpha_role;
GRANT USAGE ON COMPUTE POOL app_a_user_1_workers_CPU_X64_M_6vcpu_28gib_4x_compute_pool TO ROLE app_a_group_alpha_role;
GRANT USAGE ON COMPUTE POOL app_a_user_1_workers_HIGHMEM_X64_S_6vcpu_58gib_5x_compute_pool TO ROLE app_a_group_alpha_role;
GRANT USAGE ON COMPUTE POOL app_a_user_1_workers_CPU_X64_SL_14vcpu_58gib_7x_compute_pool TO ROLE app_a_group_alpha_role;
GRANT USAGE ON COMPUTE POOL app_a_user_1_workers_CPU_X64_L_28vcpu_116gib_14x_compute_pool TO ROLE app_a_group_alpha_role;
GRANT USAGE ON COMPUTE POOL app_a_user_1_workers_HIGHMEM_X64_M_28vcpu_240gib_19x_compute_pool TO ROLE app_a_group_alpha_role;

-- Grant monitor and operate on the compute resources to the user role data_apps_user_1_role so that users can monitor and operate the app.
GRANT MONITOR, OPERATE ON WAREHOUSE app_a_user_1_xs_warehouse TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON COMPUTE POOL app_a_user_1_frontend_CPU_X64_XS_1vcpu_6gib_1x_compute_pool TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON COMPUTE POOL app_a_user_1_frontend_CPU_X64_S_3vcpu_13gib_2x_compute_pool TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON COMPUTE POOL app_a_user_1_frontend_CPU_X64_M_6vcpu_28gib_4x_compute_pool TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON COMPUTE POOL app_a_user_1_frontend_HIGHMEM_X64_S_6vcpu_58gib_5x_compute_pool TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON COMPUTE POOL app_a_user_1_frontend_CPU_X64_SL_14vcpu_58gib_7x_compute_pool TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON COMPUTE POOL app_a_user_1_frontend_CPU_X64_L_28vcpu_116gib_14x_compute_pool TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON COMPUTE POOL app_a_user_1_frontend_HIGHMEM_X64_M_28vcpu_240gib_19x_compute_pool TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON COMPUTE POOL app_a_user_1_workers_CPU_X64_XS_1vcpu_6gib_1x_compute_pool TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON COMPUTE POOL app_a_user_1_workers_CPU_X64_S_3vcpu_13gib_2x_compute_pool TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON COMPUTE POOL app_a_user_1_workers_CPU_X64_M_6vcpu_28gib_4x_compute_pool TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON COMPUTE POOL app_a_user_1_workers_HIGHMEM_X64_S_6vcpu_58gib_5x_compute_pool TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON COMPUTE POOL app_a_user_1_workers_CPU_X64_SL_14vcpu_58gib_7x_compute_pool TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON COMPUTE POOL app_a_user_1_workers_CPU_X64_L_28vcpu_116gib_14x_compute_pool TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON COMPUTE POOL app_a_user_1_workers_HIGHMEM_X64_M_28vcpu_240gib_19x_compute_pool TO ROLE data_apps_user_1_role;

-- Grant read permissions on the table.
GRANT SELECT
  ON TABLE general_schema.image_metadata_table
  TO DATABASE ROLE general_schema_ro_db_role;

-- Grant permissions to allow the database role to create serivces (such as for submitting a job service) and to use the images in the image repository for doing so.
GRANT CREATE SERVICE ON SCHEMA group_alpha_schema
  TO DATABASE ROLE group_alpha_schema_service_db_role;
GRANT USAGE, READ ON IMAGE REPOSITORY general_schema.image_repository
  TO DATABASE ROLE group_alpha_schema_service_db_role;

-- Allow the app to read the worker spec from this stage so it can launch a worker service.
GRANT READ ON STAGE general_schema.general_stage TO DATABASE ROLE general_schema_ro_db_role;
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

-- Create a warehouse for the app.
CREATE WAREHOUSE IF NOT EXISTS app_launcher_user_1_xs_warehouse
  WAREHOUSE_SIZE = 'XSMALL'
  AUTO_RESUME = TRUE
  INITIALLY_SUSPENDED = TRUE;

-- **** CREATE THE LAUNCHER STREAMLIT APP, group_alpha_schema.app_launcher_user_1_streamlit.

-- Grant app_a_group_alpha_role privileges to see the database and schema.
GRANT USAGE ON DATABASE app_launcher_db
  TO ROLE data_apps_user_1_role;
GRANT USAGE ON SCHEMA group_alpha_schema
  TO ROLE data_apps_user_1_role;

-- Change the owner of the launcher app to the user role data_apps_user_1_role.
GRANT OWNERSHIP ON STREAMLIT group_alpha_schema.app_launcher_user_1_streamlit TO ROLE data_apps_user_1_role COPY CURRENT GRANTS;

-- Grant access to the warehouse.
GRANT USAGE, OPERATE, MONITOR ON WAREHOUSE app_launcher_user_1_xs_warehouse TO ROLE data_apps_user_1_role;
---------------- End database app_launcher_db. -----------------------------------------------


---------------- Database data_manager_db. ---------------------------------------------------
-- Create the database and schemas.
create database if not exists data_manager_db;
use database data_manager_db;
create schema if not exists group_alpha_schema;
create schema if not exists general_schema;

-- Create a general stage for holding the code for the service.
create stage if not exists general_schema.general_stage
  directory = ( enable = true );

-- Create an image repository.
CREATE IMAGE REPOSITORY IF NOT EXISTS general_schema.image_repository;

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

-- Create roles.
create role if not exists data_manager_group_alpha_role; -- This is the account role that owns and operates the app.

-- Create a warehouse for the app.
CREATE WAREHOUSE IF NOT EXISTS data_manager_user_1_xs_warehouse
  WAREHOUSE_SIZE = 'XSMALL'
  AUTO_RESUME = TRUE
  INITIALLY_SUSPENDED = TRUE;

CREATE COMPUTE POOL IF NOT EXISTS data_manager_user_1_xs_compute_pool
    MIN_NODES = 1
    MAX_NODES = 1
    INSTANCE_FAMILY = CPU_X64_XS
    AUTO_RESUME = TRUE
    INITIALLY_SUSPENDED = TRUE
    AUTO_SUSPEND_SECS = 600;

-- **** CREATE THE APP, group_alpha_schema.data_manager_user_1_service.

-- Grant data_manager_group_alpha_role privileges to see the database and schema.
GRANT USAGE ON DATABASE data_manager_db
  TO ROLE data_manager_group_alpha_role;
GRANT USAGE ON SCHEMA group_alpha_schema
  TO ROLE data_manager_group_alpha_role;

-- Change the owner of the app to the service role data_manager_group_alpha_role.
GRANT OWNERSHIP ON SERVICE group_alpha_schema.data_manager_user_1_service TO ROLE data_manager_group_alpha_role COPY CURRENT GRANTS;

-- Give this account role the appropriate database roles.
grant database role group_alpha_group_db.curated_schema_rw_db_role to role data_manager_group_alpha_role;
grant database role common_db.admin_schema_ro_db_role to role data_manager_group_alpha_role;

-- Grant access to the compute resources.
GRANT USAGE ON WAREHOUSE data_manager_user_1_xs_warehouse TO ROLE data_manager_group_alpha_role;
GRANT USAGE ON COMPUTE POOL data_manager_user_1_xs_compute_pool TO ROLE data_manager_group_alpha_role;

-- Grant resource management to the user.
GRANT MONITOR, OPERATE ON WAREHOUSE data_manager_user_1_xs_warehouse TO ROLE data_apps_user_1_role;
GRANT MONITOR, OPERATE ON COMPUTE POOL data_manager_user_1_xs_compute_pool TO ROLE data_apps_user_1_role;
---------------- End database data_manager_db. -----------------------------------------------
