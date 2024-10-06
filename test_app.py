import streamlit as st
import multiprocessing
import os
import pandas as pd
import squidpy as sq


import phenograph
import anndata
from annoy import AnnoyIndex
import dill
import hnswlib
import matplotlib.pyplot as plt
from natsort import natsorted
import numpy as np
import objsize
import parc
import parmap
import plotly.express as px
import plotnine as p9
import pynndescent
import yaml  # pyyaml
import scanpy as sc
import skimage
import sklearn
import scipy
import seaborn as sns
import setuptools_scm
import sklearn_ann
import split_file_reader
import squidpy
import streamlit_extras
from tqdm import tqdm
import umap



def main():

    st.title('Hello World!')
    st.write('This is a simple Streamlit app.')

    st.write(sq.__version__)

    st.selectbox('Select a number', [1, 2, 3], key='number')

    st.write(f'You selected {st.session_state["number"]}')

    # Write number of processors available
    st.write(f'Number of processors: {multiprocessing.cpu_count()}')

    st.write(pd.DataFrame(os.listdir('.')))

    uploaded_files = st.file_uploader(
        "Choose a CSV file", accept_multiple_files=True
    )

    st.write(len(uploaded_files))

    for uploaded_file in uploaded_files:
        st.write(uploaded_file.name)
        st.write(pd.read_csv(uploaded_file))


if __name__ == '__main__':
    main()
