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
        phenocluster__de_results[['pvals', 'pvals_adj']] = phenocluster__de_results[['pvals', 'pvals_adj']].applymap('{:.1e}'.format)
        st.dataframe(phenocluster__de_results, use_container_width=True)
        

# change cluster names
def phenocluster__edit_cluster_names(adata, edit_names_result):
    adata.obs['Edit_Cluster'] = adata.obs['Cluster'].map(edit_names_result.set_index('Cluster')['New_Name'])
    st.session_state['phenocluster__clustering_adata'] = adata

# make differential intensity plots    
def phenocluster__plot_diff_intensity(adata, groups, method, n_genes, plot_column):
    if "All" in groups:
        cur_groups = None
    else:
        cur_groups = groups
    
    if method == "Rank Plot":
        cur_fig = sc.pl.rank_genes_groups(adata, n_genes=n_genes, 
                                              groups=cur_groups, sharey=False)
    elif method == "Dot Plot":
        cur_fig = sc.pl.rank_genes_groups_dotplot(adata, n_genes=n_genes, 
                                                      groups=cur_groups)
    elif method == "Heat Map":
        cur_fig = sc.pl.rank_genes_groups_heatmap(adata, n_genes=n_genes, 
                                                      groups=cur_groups) 
    elif method == "Violin Plot":
        cur_fig = sc.pl.rank_genes_groups_stacked_violin(adata, n_genes=n_genes, 
                                                      groups=cur_groups, split = False)
        
    with plot_column:
         st.pyplot(fig = cur_fig, clear_figure=None, use_container_width=True)
     
    

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
        
    phenocluster__col3b, phenocluster__col4b  = st.columns([1, 6])
    phenocluster__col5b, phenocluster__col6b  = st.columns([1, 6])
    cur_clusters = list(pd.unique(st.session_state['phenocluster__clustering_adata'].obs["Cluster"]))
    edit_names_df = pd.DataFrame({"Cluster": cur_clusters, "New_Name": cur_clusters})
    
    with phenocluster__col3b:
        # Plot differential intensity
        phenocluster__dif_int_plot_methods = ["Rank Plot", "Dot Plot", "Heat Map", "Violin Plot"]
        st.selectbox('Select Plot Type:', phenocluster__dif_int_plot_methods, key='phenocluster__plot_diff_intensity_method')
        st.number_input(label = "Number of genes to plot", 
                                key = 'phenocluster__plot_diff_intensity_n_genes',
                                step = 1)
        
        
        st.button('Plot Markers', on_click=phenocluster__plot_diff_intensity, args = [st.session_state['phenocluster__clustering_adata'], 
                                                                                            st.session_state['phenocluster__de_sel_groups'],
                                                                                            st.session_state['phenocluster__plot_diff_intensity_method'],
                                                                                            st.session_state['phenocluster__plot_diff_intensity_n_genes'],
                                                                                            phenocluster__col4b
                                                                                            ])   
    
    with phenocluster__col6b:
        edit_clustering_names = st.data_editor(edit_names_df)
        st.session_state['phenocluster__edit_names_result'] = edit_clustering_names
    with phenocluster__col5b:
         #Edit cluster names
        st.button('Edit Clusters Names', on_click=phenocluster__edit_cluster_names, args = [st.session_state['phenocluster__clustering_adata'], 
                                                                                            st.session_state['phenocluster__edit_names_result']
                                                                                            ])
    
        

# Run the main function
if __name__ == '__main__':

    # Set page settings
    page_name = 'Differential Intensity'
    st.set_page_config(layout='wide', page_title=page_name)
    st.title(page_name)
    st.set_option('deprecation.showPyplotGlobalUse', False)
    
    # Run streamlit-dataframe-editor library initialization tasks at the top of the page
    st.session_state = sde.initialize_session_state(st.session_state)

    # Run Top of Page (TOP) functions
    st.session_state = top.top_of_page_reqs(st.session_state)
    
    # Call the main function
    main()

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)