# Use Micromamba as base image
# In the future, pin this version to ensure consistency.
FROM mambaorg/micromamba:latest

# Set the working directory in the container
WORKDIR /app

# Copy ONLY the environment file first (changes less frequently)
# COPY source/environment-leandro-compatible.yml .
# COPY source/environment-ana-20240814_to_20241219-compatible.yml .
# COPY source/environment-ana-20250605-compatible.yml .
# COPY source/environment-gmb-20240628_to_20240701-compatible.yml .
# COPY source/environment-gmb-20240917_to_20241003-compatible.yml .
COPY source/environment-dceg-compatible.yml .

# Create conda environment from environment.yml
# RUN micromamba install -y -f environment-leandro-compatible.yml && micromamba clean --all --yes
# RUN micromamba install -y -f environment-ana-20240814_to_20241219-compatible.yml && micromamba clean --all --yes
# RUN micromamba install -y -f environment-ana-20250605-compatible.yml && micromamba clean --all --yes
# RUN micromamba install -y -f environment-gmb-20240628_to_20240701-compatible.yml && micromamba clean --all --yes
# RUN micromamba install -y -f environment-gmb-20240917_to_20241003-compatible.yml && micromamba clean --all --yes
RUN micromamba install -y -f environment-dceg-compatible.yml && micromamba clean --all --yes

# "Install" foundry_transforms_lib_python by unpacking it into site-packages.
COPY temp_vendor/foundry_transforms_lib_python-0.881.0.tar.gz .
RUN tar -xzvf foundry_transforms_lib_python-0.881.0.tar.gz -C /opt/conda/lib/python3.12/site-packages/ && rm foundry_transforms_lib_python-0.881.0.tar.gz
# RUN tar -xzvf foundry_transforms_lib_python-0.881.0.tar.gz -C /opt/conda/lib/python3.9/site-packages/ && rm foundry_transforms_lib_python-0.881.0.tar.gz

# Copy source files to the container
COPY source/ .

# Copy .git to get commit hash
COPY .git .git

# Expose the port that Streamlit runs on
EXPOSE 8501

# Run the Streamlit app when the container starts
CMD ["streamlit", "run", "app.py", "--server.address", "0.0.0.0", "--server.port", "8501"]
