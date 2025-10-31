# Deploying on Snowflake

## How to deploy on Snowflake

1. Make the following example substitutions in `deploy/snowflake/deploy.sql` and in the instructions that follow:
  * `group_alpha` --> `cil`
  * `app_a` --> `mawa`
  * `App A` --> `Multiplex Analysis Web Apps`
  * `user_1` --> `andrewweisman`
1. Step through the lines in that file.
1. upload files to stages
1. build image
1. push image to snowflake repo
1. update relevant tables
1. more?

## How to add a new deployment in general, e.g., Snowflake

1. Add setup `deploy.sql` script `deploy/snowflake`.
1. Add orchestration functionality (`snowflake_orchestrator.py`) to mimic that in `docker_orchestrator/main.py`.
1. Add "snowflake" branches in `platform_abstraction.py`.
1. Step through lines in the setup `deploy.sql` script in `deploy/snowflake`.

Note that the only existing code that is modified is `platform_abstraction.py`.
