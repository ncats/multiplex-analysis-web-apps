# Full Stack MAWA

## Build instructions

Build the `frontend` and `orchestrator` images using:

* `git clone git@github.com:ncats/multiplex-analysis-web-apps.git`.
* `cd multiplex-analysis-web-apps`
* `git checkout full-stack`.
* Ensure the clone contains the file `foundry_transforms_lib_python-0.881.0.tar.gz` in a `temp_vendor` subdirectory (not present in the repository by default).
* `docker compose build`.
* Do something like: `tag=2025-10-20-01 && docker tag frontend:latest frontend:$tag && docker tag orchestrator:latest orchestrator:$tag`
  * Remember to adjust the date and increase the version each time this line is run.

The other two images (`postgres` and `minio`) should be pulled when the multi-container app is launched, below.

## Run instructions

* `docker compose up`.
* In a web browser go to http://localhost:8501.
