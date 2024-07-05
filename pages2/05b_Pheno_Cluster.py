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
import streamlit_dataframe_editor as sde
import basic_phenotyper_lib as bpl
import nidap_dashboard_lib as ndl 



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

def phenocluster__edit_cluster_names_2(adata, edit_names_result):
    edit_names_result_2 = edit_names_result.reconstruct_edited_dataframe()
    adata.obs['Edit_Cluster'] = adata.obs['Cluster'].map(dict(zip(edit_names_result_2['Cluster'].to_list(), edit_names_result_2['New_Name'].to_list())))
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

def data_editor_change_callback():
    '''
    data_editor_change_callback is a callback function for the streamlit data_editor widget
    which updates the saved value of the user-created changes after every instance of the 
    data_editor on_change method. This ensures the dashboard can remake the edited data_editor
    when the user navigates to a different page.
    '''

    st.session_state.df = bpl.assign_phenotype_custom(st.session_state.df, st.session_state['phenocluster__edit_names_result_2a'].reconstruct_edited_dataframe())

    # Create Phenotypes Summary Table based on 'phenotype' column in df
    st.session_state.pheno_summ = bpl.init_pheno_summ(st.session_state.df)

    # Perform filtering
    st.session_state.df_filt = ndl.perform_filtering(st.session_state)

    # Set Figure Objects based on updated df
    st.session_state = ndl.setFigureObjs(st.session_state, st.session_state.pointstSliderVal_Sel)  
          
def phenocluster__plotly_umaps_b(adata, umap_cur_col, umap_cur_groups, umap_color_col, plot_column):
    with plot_column:
        subcol1, subcol2 = st.columns(2)
        for i, umap_cur_group in enumerate(umap_cur_groups):
            if umap_cur_group == "All":
                subDat = adata
            else:
                subDat = adata[adata.obs[umap_cur_col] == umap_cur_group]
            umap_coords = subDat.obsm['X_umap']
            df = pd.DataFrame(umap_coords, columns=['UMAP_1', 'UMAP_2'])
            clustersList = list(subDat.obs[umap_color_col] )
            df[umap_color_col] = clustersList
            df[umap_color_col] = df[umap_color_col].astype(str)
            # Create the seaborn plot
            fig = px.scatter(df, 
             x="UMAP_1", 
             y="UMAP_2", 
             color=umap_color_col, 
             title="UMAP " + umap_cur_group
             #color_discrete_sequence=px.colors.sequential.Plasma
             )
            fig.update_traces(marker=dict(size=3)) # Adjust the size of the dots
            fig.update_layout(
                title=dict(
                    text="UMAP " + umap_cur_group,
                    x=0.5, # Center the title
                    xanchor='center',
                    yanchor='top'
                ),
                legend=dict(
                    orientation="h",
                    yanchor="top",
                    y=-0.2,
                    xanchor="right",
                    x=1
                )
            )
            if i % 2 == 0:
                subcol1.plotly_chart(fig, use_container_width=True)
            else:
                subcol2.plotly_chart(fig, use_container_width=True)

def spatial_plots_cust_2b(adata, umap_cur_col, umap_cur_groups, umap_color_col, plot_column):
    with plot_column:
        subcol3, subcol4 = st.columns(2)
        for i, umap_cur_group in enumerate(umap_cur_groups):
            if umap_cur_group == "All":
                subDat = adata
            else:
                subDat = adata[adata.obs[umap_cur_col] == umap_cur_group]
            umap_coords = subDat.obs[['Centroid X (µm)_(standardized)', 'Centroid Y (µm)_(standardized)']]
            df = pd.DataFrame(umap_coords).reset_index().drop('index', axis = 1)
            clustersList = list(subDat.obs[umap_color_col] )
            df[umap_color_col] = clustersList
            df[umap_color_col] = df[umap_color_col].astype(str)
            fig = px.scatter(df, 
             x="Centroid X (µm)_(standardized)", 
             y="Centroid Y (µm)_(standardized)", 
             color=umap_color_col, 
             title="Spatial " + umap_cur_group
             #color_discrete_sequence=px.colors.sequential.Plasma
             )
            fig.update_traces(marker=dict(size=3)) # Adjust the size of the dots
            fig.update_layout(
                title=dict(
                    text="Spatial " + umap_cur_group,
                    x=0.5, # Center the title
                    xanchor='center',
                    yanchor='top'
                ),
                legend=dict(
                    orientation="h",
                    yanchor="top",
                    y=-0.2,
                    xanchor="right",
                    x=1
                )
            )
            if i % 2 == 0:
                subcol3.plotly_chart(fig, use_container_width=True)
            else:
                subcol4.plotly_chart(fig, use_container_width=True)


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
    st.session_state['phenocluster__edit_names_df'] = edit_names_df
    
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
        #st.table(st.session_state['phenocluster__edit_names_df'])
        #edit_clustering_names = st.data_editor(edit_names_df)
        #st.session_state['phenocluster__edit_names_result'] = edit_clustering_names
        if 'phenocluster__edit_names_result_2' not in st.session_state:
            st.session_state['phenocluster__edit_names_result_2'] = sde.DataframeEditor(df_name='phenocluster__edit_names_result_2a', default_df_contents=st.session_state['phenocluster__edit_names_df'])
            
        #st.session_state['phenocluster__edit_names_result_2'].dataframe_editor(on_change=data_editor_change_callback, reset_data_editor_button_text='Reset New Clusters Names')
        st.session_state['phenocluster__edit_names_result_2'].dataframe_editor(reset_data_editor_button_text='Reset New Clusters Names')
        
    with phenocluster__col5b:
         #Edit cluster names
        st.button('Edit Clusters Names', on_click=phenocluster__edit_cluster_names_2, args = [st.session_state['phenocluster__clustering_adata'], 
                                                                                            st.session_state['phenocluster__edit_names_result_2']
                                                                                            ])
    phenocluster__col7b, phenocluster__col8b, phenocluster__col9b = st.columns([1, 3, 3])
    def make_all_plots_2():
        spatial_plots_cust_2b(st.session_state['phenocluster__clustering_adata'], 
                              st.session_state['phenocluster__umap_cur_col'], 
                              st.session_state['phenocluster__umap_cur_groups'],
                              st.session_state['phenocluster__umap_color_col_2'],
                              phenocluster__col9b
                              )
    # make umaps plots
        if 'X_umap' in st.session_state['phenocluster__clustering_adata'].obsm.keys():
            phenocluster__plotly_umaps_b(st.session_state['phenocluster__clustering_adata'], 
                                            st.session_state['phenocluster__umap_cur_col'], 
                                            st.session_state['phenocluster__umap_cur_groups'],
                                            st.session_state['phenocluster__umap_color_col_2'],
                                            phenocluster__col8b
                                            )
    with phenocluster__col7b:
            # select column for umap coloring
            st.session_state['phenocluster__umeta_columns'] = list(st.session_state['phenocluster__clustering_adata'].obs.columns)
            if 'Edit_Cluster' in st.session_state['phenocluster__umeta_columns']:
                st.session_state['phenocluster__umap_color_col_index'] = st.session_state['phenocluster__umeta_columns'].index('Edit_Cluster')
            else:
                st.session_state['phenocluster__umap_color_col_index'] = st.session_state['phenocluster__umeta_columns'].index(st.session_state['phenocluster__umap_color_col'])
                
            st.session_state['phenocluster__umap_color_col_2'] = st.selectbox('Select column to color groups:', 
                                                                    st.session_state['phenocluster__umeta_columns'],
                                                                    index=st.session_state['phenocluster__umap_color_col_index']
                                                                    )
            # select column for umap subsetting
            st.session_state['phenocluster__umap_cur_col'] = st.selectbox('Select column to subset plots:', 
                                                                    st.session_state['phenocluster__umeta_columns'], key='phenocluster__umap_col_dropdown_subset'
                                                                    )
            # list of available subsetting options
            umap_cur_groups=  ["All"] + list(pd.unique(st.session_state['phenocluster__clustering_adata'].obs[st.session_state['phenocluster__umap_cur_col']]))
            umap_sel_groups = st.multiselect('Select groups to be plotted',
                                                                            options = umap_cur_groups)
            st.session_state['phenocluster__umap_cur_groups'] = umap_sel_groups
            
            st.button('Make Plots' , on_click=make_all_plots_2)
        

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