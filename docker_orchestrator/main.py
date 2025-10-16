import os
import threading
import time
import uuid
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
import docker
from contextlib import contextmanager  # added

app = FastAPI()

WORKER_IMAGE = os.getenv("WORKER_IMAGE")
WORKER_CPUS = os.getenv("WORKER_CPUS")  # e.g., "1.0"
WORKER_MEM = os.getenv("WORKER_MEM")    # e.g., "2g"

ENV_KEYS = [
    "DATABASE_URL",
    "DB_URL_GROUP",
    "DB_URL_COMMON",
    "TZ",
    "MINIO_ENDPOINT",
    "MINIO_ACCESS_KEY",
    "MINIO_SECRET_KEY",
    "ARCHIVES_BUCKET_NAME",
    "JOB_INPUTS_BUCKET_NAME",
    "JOB_OUTPUTS_BUCKET_NAME",
    "DATA_OBJECTS_BUCKET_NAME",
]

@contextmanager
def docker_client():
    client = docker.from_env()
    try:
        yield client
    finally:
        try:
            client.close()
        except Exception:
            pass

class SubmitRequest(BaseModel):
    job_id: str

@app.get("/health")
def health():
    return {"status": "ok"}

@app.post("/jobs")
def submit_job(req: SubmitRequest):
    if not WORKER_IMAGE:
        raise HTTPException(status_code=500, detail="WORKER_IMAGE not configured")

    with docker_client() as client:
        current_container_id = os.uname().nodename
        current_container = client.containers.get(current_container_id)
        networks = current_container.attrs.get("NetworkSettings", {}).get("Networks", {})
        network_name = next(iter(networks.keys())) if networks else None

        worker_env = {k: os.environ[k] for k in ENV_KEYS if k in os.environ}

        nano_cpus = int(float(WORKER_CPUS) * 1e9) if WORKER_CPUS else None
        mem_limit = WORKER_MEM if WORKER_MEM else None

        # Use the same conda env command as frontend image
        cmd = [
            "python", "-c",
            f"import framework.analysis_framework as analysis_framework; analysis_framework.run_local_analysis('{req.job_id}')"
        ]

        try:
            worker = client.containers.run(
                image=WORKER_IMAGE,
                command=cmd,
                name=f"job-worker-{req.job_id}-{uuid.uuid4().hex[:6]}",
                detach=True,
                auto_remove=True,
                network=network_name,
                environment=worker_env,
                mem_limit=mem_limit,
                nano_cpus=nano_cpus,
                working_dir="/app",
                labels={"app.role": "worker", "app.job_id": req.job_id},
            )
            # Resolve the image ID (content digest) for traceability
            worker.reload()
            image_id = worker.image.id if hasattr(worker, "image") else None
            return {"worker_container_id": worker.id, "worker_image_id": image_id}
        except Exception as e:
            raise HTTPException(status_code=500, detail=str(e))

@app.post("/shutdown")
def shutdown():
    with docker_client() as client:

        def has_running_workers() -> bool:
            return len(client.containers.list(filters={"name": "job-worker-"})) > 0

        def stop_if_exists(name: str):
            try:
                c = client.containers.get(name)
                if name == "streamlit":
                    c.kill(signal="SIGINT")  # Graceful shutdown for Streamlit
                else:
                    c.stop(timeout=5)
                return True
            except Exception:
                return False

        running_workers = has_running_workers()
        actions = []

        if running_workers:
            # Stop the frontend immediately; let workers finish; then stop rest asynchronously.
            if stop_if_exists("streamlit"):
                actions.append("sigint_sent:streamlit")
            else:
                actions.append("frontend_not_found")

            def wait_and_shutdown_all():
                with docker_client() as tclient:
                    # Poll until all workers complete
                    while True:
                        if not tclient.containers.list(filters={"name": "job-worker-"}):
                            break
                        time.sleep(5)
                    # Stop core services (adjust names if needed)
                    for svc in ["postgres", "minio"]:
                        try:
                            c = tclient.containers.get(svc)
                            c.stop(timeout=5)
                        except Exception:
                            pass
                    # Stop Docker orchestrator (self) last
                    try:
                        self_id = os.uname().nodename
                        self_container = tclient.containers.get(self_id)
                        self_container.stop(timeout=5)
                    except Exception:
                        pass

            threading.Thread(target=wait_and_shutdown_all, daemon=True).start()
            return {
                "mode": "workers_running",
                "actions": actions,
                "follow_up": "will_stop_all_after_workers_finish"
            }

        # No workers: stop everything now except self; then stop self asynchronously
        for svc in ["streamlit", "postgres", "minio"]:
            if stop_if_exists(svc):
                if svc == "streamlit":
                    actions.append(f"sigint_sent:{svc}")
                else:
                    actions.append(f"stopped:{svc}")
            else:
                actions.append(f"not_found:{svc}")

        def stop_self_later():
            with docker_client() as tclient:
                time.sleep(1.0)
                try:
                    self_id = os.uname().nodename
                    self_container = tclient.containers.get(self_id)
                    self_container.stop(timeout=5)
                except Exception:
                    pass

        threading.Thread(target=stop_self_later, daemon=True).start()
        return {
            "mode": "no_workers",
            "actions": actions,
            "docker_orchestrator": "stopping"
        }

@app.get("/frontend_id")
def frontend_id():
    with docker_client() as client:
        try:
            frontend = client.containers.get("streamlit")
            image_id = frontend.image.id if hasattr(frontend, "image") else None
            return {"frontend_id": "streamlit", "frontend_image_id": image_id}
        except Exception as e:
            raise HTTPException(status_code=500, detail=str(e))
