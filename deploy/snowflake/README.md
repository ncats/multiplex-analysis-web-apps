# Deploying on Snowflake

## How to deploy on Snowflake

1. Step through the lines in that file.
1. upload files to stages
1. build image
1. push image to snowflake repo
1. potentially modify Streamlit launcher script
1. include diagram
1. combine with https://github.com/ncats/multiplex-analysis-web-apps/blob/full-stack/README.md and probably move this content to there
1. add instructions for data manager
1. see local instructions for missing snowflake instructions
1. more?

## How to add a new deployment in general, e.g., Snowflake

1. Add setup `deploy.sql` script `deploy/snowflake`.
1. Add orchestration functionality (`source/framework/snowflake_orchestrator.py`) to mimic that in `docker_orchestrator/main.py`.
    * If the orchestrator is not a separate container (like `snowflake_orchestrator.py`), it should be treated as such to preserve modularity. E.g., no usage of global variables such as via `streamlit` or `os.getenv()`.
1. Add "snowflake" branches in `platform_abstraction.py`.
1. Step through lines in the setup `deploy.sql` script in `deploy/snowflake`.

Note that the only existing code that is modified is `platform_abstraction.py`.
