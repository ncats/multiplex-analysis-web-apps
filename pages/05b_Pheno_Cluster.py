# Import relevant libraries
from ast import arg
from pyparsing import col
import streamlit as st
import app_top_of_page as top
import streamlit_dataframe_editor as sde
import streamlit as st 
import pandas as pd
import anndata as ad
import scanpy as sc
import seaborn as sns
import os
import matplotlib.pyplot as plt
import phenograph
import parc
from utag import utag
import numpy as np
import scanpy.external as sce
import plotly.express as px


# Functions 

# clusters differential expression
def phenocluster__diff_expr(adata, phenocluster__de_col, phenocluster__de_sel_groups, plot_column):
    sc.tl.rank_genes_groups(adata, groupby = phenocluster__de_col, method="wilcoxon", layer="counts")
    
    if "All" in phenocluster__de_sel_groups:
        phenocluster__de_results = sc.get.rank_genes_groups_df(adata, group=None)
    else:
        phenocluster__de_results = sc.get.rank_genes_groups_df(adata, group=phenocluster__de_sel_groups)
    with plot_column:
        st.dataframe(phenocluster__de_results, use_container_width=True)


def main():
    phenocluster__col1b, phenocluster__col2b  = st.columns([1, 6])
    with phenocluster__col1b:
        # differential expression
        phenocluster__de_col_options = list(st.session_state['phenocluster__clustering_adata'].obs.columns)
        st.selectbox('Select column for differential expression:', phenocluster__de_col_options, key='phenocluster__de_col')
        phenocluster__de_groups =  ["All"] +  list(pd.unique(st.session_state['phenocluster__clustering_adata'].obs[st.session_state['phenocluster__de_col']]))
        phenocluster__de_selected_groups = st.multiselect('Select group for differential expression table:', options = phenocluster__de_groups)
        st.session_state['phenocluster__de_sel_groups'] = phenocluster__de_selected_groups
        # Differential expression

        st.button('Run Differential Expression', on_click=phenocluster__diff_expr, args = [st.session_state['phenocluster__clustering_adata'], 
                                                                                            st.session_state['phenocluster__de_col'], 
                                                                                            st.session_state['phenocluster__de_sel_groups'],
                                                                                            phenocluster__col2b
                                                                                            ])

# Run the main function
if __name__ == '__main__':

    # Set page settings
    page_name = 'Your Page Name Here'
    st.set_page_config(layout='wide', page_title=page_name)
    st.title(page_name)
    
    # Run streamlit-dataframe-editor library initialization tasks at the top of the page
    st.session_state = sde.initialize_session_state(st.session_state)

    # Run Top of Page (TOP) functions
    st.session_state = top.top_of_page_reqs(st.session_state)
    
    # Call the main function
    main()

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)