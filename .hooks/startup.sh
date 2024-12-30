#!/usr/bin/env bash

set -e

# This script is always executed after the Jupyter server starts.
# Its output will show up in startup logs.

# Create custom conda environment with ipykernel so it can be used from Jupyter
# mamba create -n my_env -y -q python==3.8 ipykernel palantir-sdk
# mamba create -n my_env -y -q python==3.8 ipykernel foundry-transforms-lib-python

pip install scikit-learn streamlit-extras squidpy split-file-reader st-pages dill pympler objsize
mamba install -y -q natsort "foundry-transforms-lib-python>=0.578.0"
mamba install -y gcc_linux-64 gxx_linux-64 && mamba install -y python-annoy && mamba install -y hnswlib
