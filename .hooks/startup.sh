#!/usr/bin/env bash

set -e

# This script is always executed after the Jupyter server starts.
# Its output will show up in startup logs.

# Create custom conda environment with ipykernel so it can be used from Jupyter
# mamba create -n my_env -y -q python==3.8 ipykernel palantir-sdk

# # ---- Essentially complete the app "installation" in the Workspace ----------------------
# # During package installation
# DMAP_DASHBOARDS_DIR="$HOME/.dmap-dashboards"
# if [ ! -d $DMAP_DASHBOARDS_DIR ]; then
#     mkdir $DMAP_DASHBOARDS_DIR
# fi

# # During package installation on NIDAP
# ln -s /home/user/repo/config/platform/nidap/startup.py $DMAP_DASHBOARDS_DIR/sit_startup.py
# # ----------------------------------------------------------------------------------------

# pip install streamlit-javascript streamlit-extras squidpy split-file-reader asyncio
pip install streamlit-javascript streamlit-extras squidpy split-file-reader
