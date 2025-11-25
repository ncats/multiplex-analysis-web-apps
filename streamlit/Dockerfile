# Use Micromamba as base image.
FROM mambaorg/micromamba@sha256:800e7ade3ffe29c9a9ac2026163131495f8197c3852e572c5835beb4e8a33cd6

# Declare build arguments
ARG ENV_NAME
ARG ENV_PY_VER

# Set the working directory in the container.
WORKDIR /app

# Set up the environment.
COPY source/environment/${ENV_NAME}.lock .
COPY source/environment/requirements-${ENV_NAME}.txt .
RUN micromamba install -y --file ${ENV_NAME}.lock
RUN micromamba run python -m pip install -r requirements-${ENV_NAME}.txt
RUN micromamba clean --all --yes
COPY --chown=mambauser:mambauser temp_vendor/foundry_transforms_lib_python-0.881.0.tar.gz .
RUN tar -xzvf foundry_transforms_lib_python-0.881.0.tar.gz -C /opt/conda/lib/python${ENV_PY_VER}/site-packages/ && rm foundry_transforms_lib_python-0.881.0.tar.gz

# Copy source files to the container
COPY --chown=mambauser:mambauser source/ .

# Copy .git to get commit hash
COPY --chown=mambauser:mambauser .git .git

# optional (explicit; base image already sets this)
USER mambauser

# Expose the port that Streamlit runs on.
EXPOSE 8501

# Run the Streamlit app when the container starts.
CMD ["streamlit", "run", "app.py", "--server.address", "0.0.0.0", "--server.port", "8501"]
