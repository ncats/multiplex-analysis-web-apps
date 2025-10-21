# Use Micromamba as base image
# In the future, pin this version to ensure consistency.
FROM mambaorg/micromamba:latest

# Set the working directory in the container
WORKDIR /app

# Copy ONLY the environment file first (changes less frequently)
COPY source/environment.yml .

# Create conda environment from environment.yml
RUN micromamba install -y -f environment.yml && micromamba clean --all --yes

# "Install" foundry_transforms_lib_python by unpacking it into site-packages.
COPY temp_vendor/foundry_transforms_lib_python-0.881.0.tar.gz .
RUN tar -xzvf foundry_transforms_lib_python-0.881.0.tar.gz -C /opt/conda/lib/python3.12/site-packages/ && rm foundry_transforms_lib_python-0.881.0.tar.gz

# Copy source files to the container
COPY source/ .

# Copy input directory to the container
COPY input/ ./input/

# Copy .git to get commit hash
COPY .git .git

# Expose the port that Streamlit runs on
EXPOSE 8501

# Run the Streamlit app when the container starts
CMD ["streamlit", "run", "app.py", "--server.address", "0.0.0.0", "--server.port", "8501"]
