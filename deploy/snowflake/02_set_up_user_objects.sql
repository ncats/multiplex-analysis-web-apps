-- **** ACTION: **** Push image to /data_app_db/app_runtime_schema/image_repository/data_app:latest using, e.g.:
-- docker compose build
-- docker compose up (to test locally)
-- docker tag frontend:latest nihnci-eval.registry.snowflakecomputing.com/data_app_db/app_runtime_schema/image_repository/data_app:latest (add a tag to the image we want to push to Snowflake)
-- snow spcs image-registry login --role data_app_role
-- docker push nihnci-eval.registry.snowflakecomputing.com/data_app_db/app_runtime_schema/image_repository/data_app:latest


----------------------------------------------
-- # PUT file://deploy/snowflake/frontend_service_spec.yaml @app_a_app_db.general_schema.general_stage;
-- If service already exists, drop it first.
DROP SERVICE IF EXISTS group_alpha_schema.app_a_user_1_frontend_CPU_X64_XS_1vcpu_6gib_1x_service;
CREATE SERVICE group_alpha_schema.app_a_user_1_frontend_CPU_X64_XS_1vcpu_6gib_1x_service
  IN COMPUTE POOL app_a_user_1_frontend_CPU_X64_XS_1vcpu_6gib_1x_compute_pool
  FROM @app_a_app_db.general_schema.general_stage SPECIFICATION_TEMPLATE_FILE='frontend_service_spec.yaml'
  USING ( APP_NAME=>'app_a', IMAGE_NAME=>'frontend', IMAGE_TAG=>'latest', REQUESTS_MEMORY_GI=>6, REQUESTS_CPU=>1, LIMITS_MEMORY_GI=>6, LIMITS_CPU=>1, APP_TITLE_ONE_WORD=>'app_a', APP_TITLE=>'App A', MONITOR_JOBS_REFRESH_INTERVAL_SECONDS=>5, SNOWFLAKE_USER=>'user_1', COMPUTE_RESOURCE=>'CPU_X64_XS_1vcpu_6gib_1x', WAREHOUSE_SIZE=>'xs', ALL_COMPUTE_RESOURCES=>'CPU_X64_XS_1vcpu_6gib_1x CPU_X64_S_3vcpu_13gib_2x CPU_X64_M_6vcpu_28gib_4x HIGHMEM_X64_S_6vcpu_58gib_5x CPU_X64_SL_14vcpu_58gib_7x CPU_X64_L_28vcpu_116gib_14x HIGHMEM_X64_M_28vcpu_240gib_19x' )
  AUTO_RESUME = FALSE
  MIN_INSTANCES = 1
  MAX_INSTANCES = 1;
alter service group_alpha_schema.app_a_user_1_frontend_CPU_X64_XS_1vcpu_6gib_1x_service suspend;
alter compute pool app_a_user_1_frontend_CPU_X64_XS_1vcpu_6gib_1x_compute_pool suspend;
