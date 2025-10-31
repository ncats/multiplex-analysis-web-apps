-- **** ACTION: **** Push image to /data_app_db/app_runtime_schema/image_repository/data_app:latest using, e.g.:
-- docker compose build
-- docker compose up (to test locally)
-- docker tag frontend:latest nihnci-eval.registry.snowflakecomputing.com/data_app_db/app_runtime_schema/image_repository/data_app:latest (add a tag to the image we want to push to Snowflake)
-- snow spcs image-registry login --role data_app_role
-- docker push nihnci-eval.registry.snowflakecomputing.com/data_app_db/app_runtime_schema/image_repository/data_app:latest


----------------------------------------------
-- # PUT file://deploy/snowflake/frontend_service_spec.yaml @app_a_app_db.general_schema.general_stage;
