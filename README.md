# Full Stack MAWA

## Build instructions

Build the `frontend` and `orchestrator` images using:

* `git clone git@github.com:ncats/multiplex-analysis-web-apps.git`.
* `cd multiplex-analysis-web-apps`
* `git checkout full-stack`.
* Ensure the clone contains the file `foundry_transforms_lib_python-0.881.0.tar.gz` in a `temp_vendor` subdirectory (not present in the repository by default).
* E.g., `IMAGE_TAG=2025-10-20-03 docker compose build`.
  * Unless you're already running on AMD64, include `--platform=linux/amd64` if you want to be able to use the same image on Snowflake (which we do). For testing locally on a Mac, this should work but should be a bit slower. Alternatively, you can leave off this extra argument for testing on a Mac, but know that the image will need to be rebuilt with the argument so that it works on Snowpark Container Services.

The other two images (`postgres` and `minio`) should be pulled when the multi-container app is launched, below.

## Run instructions

* E.g., `IMAGE_TAG=2025-10-20-03 docker compose up`.
* In a web browser go to http://localhost:8501.
