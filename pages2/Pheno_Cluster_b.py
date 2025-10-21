# Import relevant libraries
from ast import arg
from pyparsing import col
import streamlit as st
import streamlit_dataframe_editor as sde
import streamlit as st 
import pandas as pd
import anndata as ad
import scanpy as sc
import seaborn as sns
import os
import matplotlib.pyplot as plt
import phenograph
#import parc
import numpy as np
import scanpy.external as sce
import plotly.express as px
import streamlit_dataframe_editor as sde
import basic_phenotyper_lib as bpl
import nidap_dashboard_lib as ndl 
import plotnine
from plotnine import *
import math


# Functions 

# clusters differential expression
def phenocluster__diff_expr(adata, phenocluster__de_col, phenocluster__de_sel_groups, only_positive):
    sc.tl.rank_genes_groups(adata, groupby = phenocluster__de_col, method="wilcoxon", layer="counts")
    
    if "All" in phenocluster__de_sel_groups:
        phenocluster__de_results = sc.get.rank_genes_groups_df(adata, group=None)
    else:
        phenocluster__de_results = sc.get.rank_genes_groups_df(adata, group=phenocluster__de_sel_groups)
        
    phenocluster__de_results[['pvals', 'pvals_adj']] = phenocluster__de_results[['pvals', 'pvals_adj']].applymap('{:.1e}'.format)
    #st.dataframe(phenocluster__de_results, use_container_width=True)
    if only_positive:
        phenocluster__de_results_filt = phenocluster__de_results[phenocluster__de_results['logfoldchanges'] > 0].reset_index(drop=True)
        st.session_state['phenocluster__de_results'] = phenocluster__de_results_filt
    else:
        st.session_state['phenocluster__de_results'] = phenocluster__de_results
    
    #st.session_state['phenocluster__de_markers'] = pd.unique(st.session_state['phenocluster__de_results']["names"])
        

# change cluster names
def phenocluster__edit_cluster_names(adata, edit_names_result):
    adata.obs['Edit_Cluster'] = adata.obs['Cluster'].map(edit_names_result.set_index('Cluster')['New_Name'])
    st.session_state['phenocluster__clustering_adata'] = adata   

def phenocluster__edit_cluster_names_2(adata, edit_names_result):
    edit_names_result_2 = edit_names_result.reconstruct_edited_dataframe()
    adata.obs['Edit_Cluster'] = adata.obs['Cluster'].map(dict(zip(edit_names_result_2['Cluster'].to_list(), edit_names_result_2['New_Name'].to_list())))
    st.session_state['phenocluster__clustering_adata'] = adata
    
# make differential intensity plots    
def phenocluster__plot_diff_intensity(adata, groups, method, n_genes, cur_col):
    if "All" in groups:
        cur_groups = None
    else:
        cur_groups = groups
    
    if method == "Rank Plot":
        cur_fig = sc.pl.rank_genes_groups(adata, n_genes=n_genes, 
                                              groups=cur_groups, sharey=False, fontsize=20, show=False)
    elif method == "Dot Plot":
        cur_fig = sc.pl.rank_genes_groups_dotplot(adata, n_genes=n_genes, 
                                                      groups=cur_groups)
    elif method == "Heat Map":
        # cur_fig = sc.pl.rank_genes_groups_heatmap(adata, n_genes=n_genes, 
        #                                               groups=cur_groups) 
        #sc.pp.normalize_total(adata)
        #sc.pp.log1p(adata)
        #sc.pp.scale(adata)

        if "Edit_Cluster" in adata.obs.columns:
            adata_sub  = adata[adata.obs['Edit_Cluster'].isin(cur_groups)]
            top_names = pd.unique(st.session_state['phenocluster__de_results'].groupby('group')['names'].apply(lambda x: x.head(n_genes)))
            cluster_group = "Edit_Cluster"
            cur_fig = sc.pl.heatmap(adata_sub, top_names, groupby=cluster_group, swap_axes=False)
        else:
            adata_sub  = adata[adata.obs['Cluster'].isin(cur_groups)]
            top_names = pd.unique(st.session_state['phenocluster__de_results'].groupby('group')['names'].apply(lambda x: x.head(n_genes)))
            cluster_group = "Cluster"
            cur_fig = sc.pl.heatmap(adata_sub, top_names, groupby=cluster_group, swap_axes=False)
    elif method == "Violin Plot":
        cur_fig = sc.pl.rank_genes_groups_stacked_violin(adata, n_genes=n_genes, 
                                                      groups=cur_groups, split = False)
    # elif method == "UMAP" and "X_umap" in adata.obsm.keys():
    #     adata_sub  = adata[adata.obs['Cluster'].isin(cur_groups)]
    #     top_names = pd.unique(st.session_state['phenocluster__de_results'].groupby('group')['names'].apply(lambda x: x.head(n_genes)))
    #     with cur_col:
    #        cur_fig  = st.pyplot(sc.pl.umap(adata, color=[*top_names], legend_loc="on data",frameon=True,
    #                    ncols=3, show=True, 
    #                    wspace = 0.2 ,save  = False), use_container_width = True , clear_figure = True)
        
def phenocluster__plot_diff_intensity_2(adata, groups, method, n_genes, 
                                        plot_column, normalize_total,
                                        log_normalize, z_normalize, 
                                        cluster_group):
    if "All" in groups:
        cur_groups = list(pd.unique(st.session_state['phenocluster__de_results']["group"]))
    else:
        cur_groups = groups
    adata_norm = adata.copy()
    if normalize_total:
        sc.pp.normalize_total(adata_norm)
    if log_normalize:
        sc.pp.log1p(adata_norm)
    if z_normalize:
        sc.pp.scale(adata_norm)
    
    adata_sub  = adata_norm[adata_norm.obs[cluster_group].isin(cur_groups)]
    top_names = pd.unique(st.session_state['phenocluster__de_results'].groupby('group')['names'].apply(lambda x: x.head(n_genes)))

    with plot_column:
        if method == "Heat Map":
            obs_df = adata_sub.to_df().fillna(0)
            obs_df[cluster_group] = adata_sub.obs[cluster_group]
            mean_per_group = obs_df.groupby(cluster_group).mean()
            matrix_avg = mean_per_group.stack().reset_index()
            matrix_avg.columns = [cluster_group, "Marker", "Intensity"]
            plot_mat = matrix_avg[matrix_avg["Marker"].isin(top_names)].reset_index(drop=True)
            plot_mat["Marker"] = pd.Categorical(plot_mat["Marker"], categories=top_names[::-1], ordered=True)
            plotnine.options.figure_size = (5, 4)
            plot = (
                ggplot(plot_mat, aes(cluster_group, "Marker"))
                + geom_tile(mapping = aes(fill = "Intensity"))
                + scale_fill_distiller(type = 'div', palette = 'RdYlBu')
                + theme(axis_text_x=element_text(rotation=0, hjust=0.5, size=8))
                + theme(axis_text_y=element_text(rotation=0, hjust=1, size=8))
                + theme(axis_title_x = element_blank(), axis_title_y = element_text(angle=90))
                + theme(text=element_text(size=8))
                )
            
            st.pyplot(ggplot.draw(plot), use_container_width=False)
        
        elif method == "UMAP":
            obs_df = adata_norm[:, top_names].to_df().reset_index(drop=True)
            umap_coords = adata_norm.obsm['X_umap']
            umap_df = pd.DataFrame(umap_coords, columns=['UMAP_1', 'UMAP_2']).reset_index(drop=True)
            umap_df["Cells_Id"] = umap_df.index
            obs_df = pd.concat([obs_df, umap_df], axis=1)
            # Get the column names excluding 'UMAP_1' and 'UMAP_2'
            columns = top_names
            # Calculate the number of rows needed for the subplots
            n_rows = int(np.ceil(len(columns) / 4))
            # Create a figure with subplots
            fig, axs = plt.subplots(n_rows, 4, figsize=(20, 5*n_rows))
            # Flatten the axes array to make it easier to iterate over
            axs = axs.flatten()
            for ax, col in zip(axs, columns):
                plot = sns.scatterplot(data=obs_df, x='UMAP_1', y='UMAP_2', hue=col, 
                                       palette='viridis', ax=ax, s=3)
                ax.set_title(col)
                plot.legend(loc='upper left', bbox_to_anchor=(1, 1))
                ax.set_aspect('equal', 'box')
            # Remove any unused subplots
            for ax in axs[len(columns):]:
                ax.remove()
            plt.tight_layout()
            st.pyplot(fig, use_container_width=True)
        elif method == "Rank Plot":
            #obs_df = st.session_state['phenocluster__de_results'][st.session_state['phenocluster__de_results']['group'].isin(cur_groups)].reset_index(drop=True)
            top_n_df = st.session_state['phenocluster__de_results'].groupby('group').head(n_genes).reset_index(drop=True)
            top_n_df["ranking"]= top_n_df.groupby('group')['names'].cumcount()
            
            # Get the unique groups
            groups = top_n_df['group'].unique()

            # Calculate the number of rows needed for the subplots
            rows = math.ceil(len(groups) / 2)
            cols = 2

            # Create a figure with subplots
            fig, axes = plt.subplots(rows, cols, figsize=(15, 5*rows))

            # Flatten the axes array for easier indexing
            axes = axes.flatten()

            # For each group, create a scatter plot
            for i, group in enumerate(groups):
                df_group = top_n_df[top_n_df['group'] == group]
                
                # Create scatter plot
                sns.scatterplot(x='ranking', y='scores', data=df_group, ax=axes[i], alpha=0)
                
                # Add text labels
                for _, row in df_group.iterrows():
                    axes[i].text(row['ranking'], row['scores'], row['names'], 
                                ha='center', va='bottom', rotation='vertical', fontsize=12)  # Increase fontsize here
                
                # Set the x-axis range and ticks
                axes[i].set_xlim(df_group['ranking'].min()-2, df_group['ranking'].max()+2)
                axes[i].set_xticks(range(int(df_group['ranking'].min()-2), int(df_group['ranking'].max()+2), 2))
                axes[i].set_ylim(df_group['scores'].min()-2, df_group['scores'].max()+10)
                
                # Set titles and labels
                axes[i].set_title(group, fontsize=16)  # Increase fontsize here
                axes[i].set_xlabel('ranking', fontsize=14)  # Increase fontsize here
                axes[i].set_ylabel('scores', fontsize=14)  # Increase fontsize here
                axes[i].tick_params(axis='x', labelsize=14)  # Increase fontsize here
                axes[i].tick_params(axis='y', labelsize=14) 

            # Remove any unused subplots
            for j in range(i+1, len(axes)):
                fig.delaxes(axes[j])

            # Adjust the layout and show the plot
            plt.tight_layout()
            # Show the plot
            st.pyplot(fig, use_container_width=True)
        
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
                ),
                xaxis=dict(
                    scaleanchor="y",
                    scaleratio=1)
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
                ),
                xaxis=dict(
                    scaleanchor="y",
                    scaleratio=1)
            )
            if i % 2 == 0:
                subcol3.plotly_chart(fig, use_container_width=True)
            else:
                subcol4.plotly_chart(fig, use_container_width=True)

def phenocluster__add_edit_clusters_to_input_df():
    if "phenocluster__phenotype_cluster_cols" in st.session_state:
        cur_df = st.session_state['input_dataset'].data
        cur_df = cur_df.drop(columns=st.session_state["phenocluster__phenotype_cluster_cols"])
        st.session_state['input_dataset'].data = cur_df
    
    st.session_state['input_dataset'].data["Phenotype_Cluster"] = 'Phenotype ' + st.session_state['phenocluster__clustering_adata'].obs['Edit_Cluster'].astype(str)
    dummies = pd.get_dummies(st.session_state['phenocluster__clustering_adata'].obs['Edit_Cluster'], prefix='Phenotype Cluster').astype(int)
    #dummies = dummies.replace({1: '+', 0: '-'})
    cur_df = pd.concat([st.session_state['input_dataset'].data, dummies], axis=1)
    st.session_state['input_dataset'].data = cur_df
    new_cluster_cols = list(dummies.columns)
    st.session_state["phenocluster__phenotype_cluster_cols"] = new_cluster_cols

def main():
    if 'phenocluster__clustering_adata' not in st.session_state:
        st.error("Please run the clustering step first",
                     icon="🚨")
        return
    
    if "Cluster" not in st.session_state['phenocluster__clustering_adata'].obs.columns:
        st.error("Please run the clustering step first",
                     icon="🚨")
        return
              
    phenocluster__col1b, phenocluster__col2b  = st.columns([2, 6])
    phenocluster__col3b, phenocluster__col4b  = st.columns([2, 6])
    sc.set_figure_params(figsize=(10, 10), fontsize = 16)
    if 'phenocluster__dif_int_plot_methods' not in st.session_state:
        st.session_state['phenocluster__dif_int_plot_methods'] = ["Rank Plot", "Heat Map", "UMAP"]
        
    with phenocluster__col1b:
        # differential expression
        phenocluster__de_col_options = list(st.session_state['phenocluster__clustering_adata'].obs.columns)
        st.selectbox('Select column for differential intensity:', phenocluster__de_col_options, key='phenocluster__de_col')
        phenocluster__de_groups =  ["All"] +  list(pd.unique(st.session_state['phenocluster__clustering_adata'].obs[st.session_state['phenocluster__de_col']]))
        phenocluster__de_selected_groups = st.multiselect('Select group for differential intensity table:', options = phenocluster__de_groups)
        st.session_state['phenocluster__de_sel_groups'] = phenocluster__de_selected_groups
        st.checkbox('Only Positive Markers', key='phenocluster__de_only_positive')
        # Differential expression

        st.button('Run Differential Intensity', on_click=phenocluster__diff_expr, args = [st.session_state['phenocluster__clustering_adata'], 
                                                                                            st.session_state['phenocluster__de_col'], 
                                                                                            st.session_state['phenocluster__de_sel_groups'],
                                                                                            st.session_state['phenocluster__de_only_positive']
                                                                                            ])
    if 'phenocluster__de_results' in st.session_state:
        with phenocluster__col2b:
            st.dataframe(st.session_state['phenocluster__de_results'], use_container_width=True)
            
        with phenocluster__col3b:
            # Plot differential intensity
            st.selectbox('Select Plot Type:', st.session_state['phenocluster__dif_int_plot_methods'], key='phenocluster__plot_diff_intensity_method')
            st.number_input(label = "Number of markers to plot", 
                                key = 'phenocluster__plot_diff_intensity_n_genes',
                                step = 1)
            st.toggle("Normalize total intensity", key='phenocluster__plot_normalize_total_intensity')
            st.toggle("Log normalize", key='phenocluster__plot_log_normalize')
            st.toggle("Z-score normalize columns", key='phenocluster__plot_zscore_normalize')
            
            # phenocluster__plot_diff_intensity(st.session_state['phenocluster__clustering_adata'], 
            #                                   st.session_state['phenocluster__de_sel_groups'],
            #                                   st.session_state['phenocluster__plot_diff_intensity_method'],
            #                                   st.session_state['phenocluster__plot_diff_intensity_n_genes'],
            #                                   phenocluster__col4b)
            st.button('Plot Markers', on_click=phenocluster__plot_diff_intensity_2, args = [st.session_state['phenocluster__clustering_adata'], 
                                                                                          st.session_state['phenocluster__de_sel_groups'],
                                                                                          st.session_state['phenocluster__plot_diff_intensity_method'],
                                                                                          st.session_state['phenocluster__plot_diff_intensity_n_genes'],
                                                                                            phenocluster__col4b,
                                                                                            st.session_state['phenocluster__plot_normalize_total_intensity'],
                                                                                            st.session_state['phenocluster__plot_log_normalize'],
                                                                                            st.session_state['phenocluster__plot_zscore_normalize'],
                                                                                            st.session_state['phenocluster__de_col']
                                                                                          ])   
            
    # make plots for differential intensity markers 
    # if 'phenocluster__diff_intensity_plot' in st.session_state:
    #     with phenocluster__col4b:
    #         cur_fig = st.session_state['phenocluster__diff_intensity_plot']
            #st.pyplot(fig = cur_fig, use_container_width = True, clear_figure = False)
    
    phenocluster__col5b, phenocluster__col6b  = st.columns([2, 6])
    
    cur_clusters = list(pd.unique(st.session_state['phenocluster__clustering_adata'].obs["Cluster"]))
    edit_names_df = pd.DataFrame({"Cluster": cur_clusters, "New_Name": cur_clusters})
    st.session_state['phenocluster__edit_names_df'] = edit_names_df
    
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
    phenocluster__col7b, phenocluster__col8b= st.columns([2, 6])
    def make_all_plots_2():
        spatial_plots_cust_2b(st.session_state['phenocluster__clustering_adata'], 
                              st.session_state['phenocluster__umap_cur_col'], 
                              st.session_state['phenocluster__umap_cur_groups'],
                              st.session_state['phenocluster__umap_color_col_2'],
                              phenocluster__col8b
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
            
            st.button('Make Spatial Plots' , on_click=spatial_plots_cust_2b, 
                      args = [st.session_state['phenocluster__clustering_adata'], 
                              st.session_state['phenocluster__umap_cur_col'], 
                              st.session_state['phenocluster__umap_cur_groups'],
                              st.session_state['phenocluster__umap_color_col_2'],
                              phenocluster__col8b]
                      )
            
            if 'X_umap' in st.session_state['phenocluster__clustering_adata'].obsm.keys():
                st.button("Make UMAPS", on_click=phenocluster__plotly_umaps_b,
                          args= [st.session_state['phenocluster__clustering_adata'], 
                                            st.session_state['phenocluster__umap_cur_col'], 
                                            st.session_state['phenocluster__umap_cur_groups'],
                                            st.session_state['phenocluster__umap_color_col_2'],
                                            phenocluster__col8b]
                                            )
            
            st.button('Add Edited Clusters to Input Data', on_click=phenocluster__add_edit_clusters_to_input_df)
    

# Run the main function
if __name__ == '__main__':
    main()
