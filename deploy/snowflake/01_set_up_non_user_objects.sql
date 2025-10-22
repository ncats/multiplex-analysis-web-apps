USE ROLE ACCOUNTADMIN;

CREATE ROLE IF NOT EXISTS data_app_role;

CREATE WAREHOUSE IF NOT EXISTS data_app_warehouse
  WAREHOUSE_SIZE = 'XSMALL'
  AUTO_RESUME = TRUE
  INITIALLY_SUSPENDED = TRUE;
GRANT USAGE, OPERATE, MONITOR ON WAREHOUSE data_app_warehouse TO ROLE data_app_role;

CREATE DATABASE IF NOT EXISTS data_app_db;
GRANT OWNERSHIP ON DATABASE data_app_db TO ROLE data_app_role COPY CURRENT GRANTS;

-------------------------------------------------
-- Since I do not want myself (an admin) to have data_app_role, temporarily grant it to myself to do the setup, and revoke it in 03_set_up_user_role_permissions.sql.
-- Note that a better strategy long-term is to "Create a dedicated service / CI user (e.g. data_app_provisioner_user) that permanently has data_app_role (or a separate provisioning role that then grants ownership to the runtime role)."
-- As SECURITYADMIN (or ACCOUNTADMIN if early bootstrap)
GRANT ROLE data_app_role TO USER andrewweisman;
-------------------------------------------------

USE ROLE data_app_role;
USE WAREHOUSE data_app_warehouse;

create schema if not exists data_app_db.user_data_schema;
create schema if not exists data_app_db.app_data_schema;
create schema if not exists data_app_db.app_runtime_schema;

CREATE IMAGE REPOSITORY IF NOT EXISTS data_app_db.app_runtime_schema.image_repository;

CREATE STAGE IF NOT EXISTS data_app_db.app_data_schema.archives_stage
  DIRECTORY = ( ENABLE = TRUE );
CREATE STAGE IF NOT EXISTS data_app_db.app_data_schema.inputs_stage
  DIRECTORY = ( ENABLE = TRUE );
CREATE STAGE IF NOT EXISTS data_app_db.app_data_schema.outputs_stage
  DIRECTORY = ( ENABLE = TRUE );
CREATE STAGE IF NOT EXISTS data_app_db.app_data_schema.oldarchives_stage
  DIRECTORY = ( ENABLE = TRUE );

-- Tables
CREATE TABLE IF NOT EXISTS data_app_db.app_data_schema.user_groups_table (
  id INTEGER IDENTITY PRIMARY KEY,
  username VARCHAR(255) UNIQUE NOT NULL,
  user_group VARCHAR(255),
  user_added_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
  who_added VARCHAR(255),
  user_email VARCHAR(255)
);

MERGE INTO data_app_db.app_data_schema.user_groups_table t
USING (
  SELECT * FROM (
    VALUES
      ('andrew','dmap','andrew','andrew@example.com'),
      ('jessica','ABC Lab','andrew','jessica@example.com'),
      ('Tessa','dmap','andrew','tessa@example.com'),
      ('andrewweisman', 'dmap', 'andrewweisman', 'andrew.weisman@nih.gov')
  ) AS v(username, user_group, who_added, user_email)
) s
ON t.username = s.username
WHEN MATCHED THEN UPDATE SET
  user_group = s.user_group,
  who_added = s.who_added,
  user_email = s.user_email
WHEN NOT MATCHED THEN INSERT (username, user_group, who_added, user_email)
VALUES (s.username, s.user_group, s.who_added, s.user_email);

CREATE TABLE IF NOT EXISTS data_app_db.app_data_schema.app_sessions_table (
  id INTEGER IDENTITY PRIMARY KEY,
  app_session_id VARCHAR(255) UNIQUE NOT NULL,
  username VARCHAR(255),
  user_group VARCHAR(255),
  startup_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
  explicit_shutdown_time TIMESTAMP,
  container_image_id VARCHAR(255)
);

CREATE TABLE IF NOT EXISTS data_app_db.app_data_schema.archives_table (
  id INTEGER IDENTITY PRIMARY KEY,
  creator VARCHAR(255),
  user_group VARCHAR(255),
  archive_description VARCHAR,
  current_git_commit VARCHAR(255),
  container_image_id VARCHAR(255),
  archive_id VARCHAR(255) UNIQUE NOT NULL,
  app_session_id VARCHAR(255),
  creation_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS data_app_db.app_data_schema.jobs_table (
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
  failure_time TIMESTAMP
);
