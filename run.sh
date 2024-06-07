#!/bin/bash

input_dir="./input"
output_dir="./output"

if [ ! -d "$input_dir" ]; then
    mkdir "$input_dir"
fi
# if [ ! -d "$input_dir/phenotypes" ]; then
#     mkdir "$input_dir/phenotypes"
# fi
# if [ ! -d "$input_dir/annotations" ]; then
#     mkdir "$input_dir/annotations"
# fi
if [ ! -d "$output_dir" ]; then
    mkdir "$output_dir"
fi


streamlit run Multiplex_Analysis_Web_Apps.py
