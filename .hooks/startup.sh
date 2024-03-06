#!/usr/bin/env bash

set -e

# This script is always executed after the Jupyter server starts.
# Its output will show up in startup logs.

# Create custom conda environment with ipykernel so it can be used from Jupyter
# mamba create -n my_env -y -q python==3.8 ipykernel palantir-sdk

pip install scikit-learn streamlit-extras squidpy split-file-reader st-pages
mamba install -y -q natsort
