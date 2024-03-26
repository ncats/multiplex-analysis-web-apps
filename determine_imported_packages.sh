#/bin/bash

# See the end of this script for notes!

# Get a list of all imported packages in MAWA
imported_packages=$(grep -E -Hinr "^ *import |^ *from .+ import " | awk -v FS=: '{gsub(/^ */, "", $3); print $3}' | awk '{print $2}' | awk -v FS=. '{print $1}' | sort -u)

# Get a list of all local Python scripts in the current directory, which could potentially be modules
local_python_scripts=$(ls *.py | awk -v FS="\\\.py" '{print $1}')

# Output
echo "Imported packages in MAWA that are not local modules:"
echo "----"

# For every imported package...
for package in $imported_packages; do

    # Assume it's not a local script
    is_script=0

    # For every local script...
    for script in $local_python_scripts; do

        # If the package is a local script, then store that fact
        if [ "x$package" == "x$script" ]; then
            is_script=1
            break
        fi

    done

    # If the imported package is not a local script, then output it to the screen
    if [ $is_script -eq 0 ]; then
        echo $package
    fi

done

# Output
echo "----"
echo ""
echo "Now it's probably most efficient to go through this list manually to determine those that are actually used and that aren't part of the standard library."

# Script output on 2/8/24 at 11:50 AM EST:
# altair
# anndata
# ast
# contextlib
# dataset_format_conversion
# datetime
# foundry
# functools
# glob
# io
# json
# matplotlib
# multiprocessing
# natsort
# numpy
# operator
# os
# palantir
# pandas
# pathlib
# pickle
# plotly
# pprint
# random
# re
# requests
# scipy
# seaborn
# shutil
# skimage
# sklearn
# socket
# split_file_reader
# squidpy
# st_pages
# streamlit
# streamlit_extras
# streamlit_javascript
# string
# subprocess
# sys
# tarfile
# time
# tqdm
# umap
# urllib
# warnings
# yaml
# zipfile

# Andrew's manual list of requirements from the above list that needs to be installed to get MAWA working locally:
# anndata
# matplotlib
# natsort
# numpy
# palantir
# pandas
# plotly
# scipy
# seaborn
# skimage
# sklearn
# split_file_reader
# squidpy
# st_pages
# streamlit
# streamlit_extras
# streamlit_javascript
# tqdm
# umap
# yaml

# List of corresponding actual package names to install per README.md:
# anndata
# matplotlib
# natsort
# numpy
# palantir-sdk
# pandas
# plotly
# scipy
# seaborn
# scikit-image
# scikit-learn
# split-file-reader
# squidpy
# st-pages
# streamlit
# streamlit-extras
# streamlit-javascript
# tqdm
# umap-learn
# pyyaml
