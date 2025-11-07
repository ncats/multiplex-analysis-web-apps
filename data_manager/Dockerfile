# Use Micromamba as base image.
FROM mambaorg/micromamba:git-9a50bbf

# Set the working directory in the container.
WORKDIR /app

# Copy ONLY the environment file first (changes less frequently).
COPY data_manager/environment.yml .

# Create conda environment from environment.yml.
RUN micromamba install -y -f environment.yml && micromamba clean --all --yes

# Copy source files to the container.
COPY data_manager/ .

# Expose the port that Streamlit runs on.
EXPOSE 8501

# Run the Streamlit app when the container starts.
CMD ["streamlit", "run", "app.py", "--server.address", "0.0.0.0", "--server.port", "8501"]
