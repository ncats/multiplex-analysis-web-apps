'''
NIDAP Dashboard Library (NDL) is a series of functions
that supports the types of analysis activities needed by the 
NCI researchers with their data on NIDAP.
'''

import os
import time
from copy import copy
import numpy as np
import pandas as pd
import altair as alt
from natsort import natsorted
from pathlib import Path
from datetime import datetime
alt.data_transformers.disable_max_rows()

# Import relevant libraries
import basic_phenotyper_lib as bpl                  # Useful functions for cell phenotyping
from foundry_IO_lib import foundry_IO_lib           # Foundry Input/Output Class
from benchmark_collector import benchmark_collector # Benchmark Collector Class
from neighborhood_profiles import NeighborhoodProfiles, UMAPDensityProcessing
import PlottingTools as umPT

def identify_col_type(col):
    '''
    Quick and dirty column identifaction function
    '''
    dtypes = col.dtypes
    n_uni = col.nunique()

    if n_uni == 2:
        return 'bool'
    if n_uni == 1:
        return 'constant'
    elif dtypes == 'category':
        return 'category'
    elif dtypes == 'object':
        return 'object'
    else:
        return 'not_bool'

def init_session_state(session_state):
    """
    Initialize session_state values for streamlit processing
    """

    session_state.init             = True
    session_state.init_phenotyping = True

    # Create an instance of the foundry IO Library
    session_state.fiol = foundry_IO_lib()
    session_state.bc   = benchmark_collector(session_state.fiol)

    # Analysis Settings
    session_state.marker_pre  = 'Phenotype '

    # Create a default dataframe (df_default)
    # This is a placeholder for figures/tables before
    # real data is loaded
    df_dict = {'Slide ID': ['imagenum_Demo'],
               'Cell X Position': [0],
               'Cell Y Position': [0],
               'Index': [1],
               'Phenotype a': [1]}

    df_default                = pd.DataFrame(data = df_dict)
    session_state.reqFeatures = df_default.columns[:-1]

    # Features for filtering
    session_state.file_format = 'Native'
    session_state.SEL_feat = []
    session_state.CHK_feat = []

    # session_state['phenotyping_micron_coordinate_units'] = 0.25

    # Dataset Dictionaries of files in each unstructure dataset
    # session_state.files_dict = {}
    # for dataset in session_state.usDatasetPaths:
    #     session_state.files_dict[dataset] = load_listofFiles(session_state.fiol, dataset)

    # List of DataSets to save CSV to
    session_state.OutputCSVPaths_U = './output'

    # List of DataSets to save PNGS to
    session_state.OutputPNGPaths = session_state.OutputCSVPaths_U

    session_state.files_to_export = pd.DataFrame(columns = ['Item Name', 'File Name', 'Date Time Added'])
    session_state.files_to_export['Date Time Added'] = pd.to_datetime(session_state.files_to_export['Date Time Added'])

    # Selected Dataset Meta Information
    session_state.selectProj = ''
    session_state.datafile   = ''

    # Default Upload File value
    session_state.uploaded_file = None

    # Output file information
    session_state.df_update_filename_U = ''
    session_state.pheno_assign_filename_U = 'phenotype_summary'
    session_state.NeighborProPng = 'NeighborhoodProfile'

    # Image settings
    session_state.figsize = (12, 12)
    session_state.imgFileSuffixText = ''
    session_state.colorPalette = np.array([[131, 201, 255], [125, 239, 161], [109, 63, 192], [255, 171, 171],
                                           [0, 104, 201], [41, 176, 157], [213, 218, 229], [255, 43, 43],
                                           [232, 197, 77], [255, 209, 106], [131, 201, 255], [160, 120, 148], 
                                           [178, 70, 181], [255, 135, 0], [0, 104, 0], [82, 166, 56],
                                           [166, 63, 152], [141, 166, 42], [130, 86, 133], [133, 54, 23],
                                           [9, 89, 133], [240, 135, 228], [240, 188, 176], [113, 32, 240],
                                           [57, 240, 223], [95, 166, 94], [94, 94, 65], [94, 51, 51],
                                           [50, 94, 32], [252, 226, 17]])/256

    ### LOAD DATA BUTTON ###
    # 'Load Data' Button for our default Dataset
    session_state = loadDataButton(session_state, df_default, 'None', 'None')

    ## Some error handling text:
    # Place holder for possible error message
    session_state.errmsg_def2row = '''
                                       
                                   '''
    # Error message when selected csv does not include the requisite features
    session_state.errmsg_wrongCol = ''':green[Please only select .csv files that include columns listed in the
                                       Required DATASET Features (About Page). See Dante or Andrew for help]'''

    # Inital Phenotyping state
    session_state.phenotyping_completed = False

    # Reset all of the neighborhood Profiles settings
    session_state = reset_neigh_profile_settings(session_state)

    # General Neighborhood Profile Page Settings
    session_state.cpu_pool_size = 7
    session_state.umap_subset_toggle = False
    session_state.umap_subset_per = 20
    session_state.area_filter_per = 0.001
    session_state.clust_minmax = [1, 10]
    session_state.toggle_clust_diff = False
    session_state.appro_feat = False

    # Set data_loaded = False.
    # This needs to happen at the end to counteract the 'loadDataButton' action
    session_state.data_loaded = False
    return session_state

def reset_neigh_profile_settings(session_state):
    '''
    Resets all the variables required for neighborhood
    profiles analysis 
    '''

    print('Resetting Neighborhood Profiles Analysis Settings')

    # Define the checkpoint directory
    session_state.checkpoint_dir = './output/checkpoints/neighborhood_profiles'
    if not os.path.exists(session_state.checkpoint_dir):
        os.makedirs(session_state.checkpoint_dir)

    # Has the UMAP been completed yet?
    session_state.density_completed = False
    session_state.umap_completed    = False
    session_state.cluster_completed = False
    session_state.UMAPFigType       = 'Density'

    # UMAP Lineage Display
    session_state.lineageDisplayToggle = 'Phenotypes'
    session_state.lineageDisplayToggle_clus = 'Phenotypes'

    # Unfiltered dropdown default options
    session_state.defLineageOpt    = 'All Phenotypes'
    session_state.defumapOutcomes  = 'phenotype'
    session_state.definciOutcomes  = 'Cell Counts'

    # Default UMAP dropdown options
    session_state.umapPheno = [session_state.defLineageOpt]
    session_state.umapMarks = [session_state.defLineageOpt]
    session_state.umaplineages = [session_state.defLineageOpt]
    session_state.umapOutcomes = [session_state.defumapOutcomes]

    # Default Incidence dropdown options
    session_state.outcomes     = [session_state.defumapOutcomes]
    session_state.inciOutcomes = [session_state.definciOutcomes]

    # Default UMAPInspect settings
    session_state.umapInspect_Ver = session_state.defLineageOpt
    session_state.umapInspect_Feat = session_state.defumapOutcomes

    # Default UMAP differences settings
    session_state.diffUMAPSel_Ver  = session_state.defLineageOpt
    session_state.diffUMAPSel_Feat = session_state.defumapOutcomes

    # Default Incidence settings
    session_state.inciPhenoSel   = session_state.defLineageOpt
    session_state.inciOutcomeSel = session_state.definciOutcomes
    session_state.Inci_Value_display = 'Count Differences'

    # Default Cluster_Dict()
    session_state.cluster_dict = {0: 'No Cluster'}

    # Neighborhood Profiles Line Plot Settings
    session_state.compare_clusters_as = 'Difference'
    session_state.palette_dict = 'bwr'

    return session_state

def load_listofFiles(fiol, projectPath):
    """
    Identify datasets available within Unstructured Dataset. 
    Calling function from Foundry IO Library (FIOL)
    """
    return fiol.load_listofFiles(projectPath)

def loadDataButton(session_state, df_import, projectName, fileName):
    """
    All the required data processing steps when 'Load Data' Button is pressed
    """
    print('Loading Data')

    # Meta Data
    session_state.selectProj = projectName # Project Name
    session_state.datafile   = fileName    # File Name
    session_state.df_update_filename_U = session_state.datafile + '_updated'

    # Identify Markers in the dataset
    session_state.bc.startTimer()
    session_state.marker_names = bpl.identify_marker_columns(df_import, session_state.marker_pre)
    session_state.bc.printElapsedTime(msg = 'Identifying Marker Names')

    # Set Phenotyping Elements
    session_state.bc.startTimer()
    session_state = set_phenotyping_elements(session_state, df_import)
    session_state.bc.printElapsedTime(msg = 'Setting Phenotying Elements')

    # Data has now undergone enough transformation to be called 'LOADED'
    session_state.data_loaded = True

    # Analysis Setting Init
    session_state.loaded_marker_names = session_state.marker_names
    session_state.marker_multi_sel = session_state.marker_names
    session_state.pointstSliderVal_Sel = 100
    session_state.calcSliderVal  = 100
    session_state.selected_nClus = 1         # Clustering (If applicable)
    session_state.NormHeatRadio  = 'No Norm' # Heatmap Radio

    # Initalize Filtering Settings
    session_state.SEL_feat_widg = []
    session_state.CHK_feat_widg = []
    session_state.SEL_feat = session_state.SEL_feat_widg + ['Slide ID']
    session_state.CHK_feat = session_state.CHK_feat_widg + ['has_pos_mark']

    # if session_state.file_format == 'REEC':
    #     session_state.SEL_feat.extend(['tNt'])
    #     session_state.CHK_feat.extend(['GOODNUC'])

    # All filter categories
    features4filter = session_state.SEL_feat + session_state.CHK_feat
    # Create variables in session state
    for feature in features4filter:
        session_state[eval('"uni" + feature')] = natsorted(session_state.df_raw[feature].unique())    # Unique Values
        if feature in session_state.CHK_feat:
            session_state[eval('"sel" + feature')] = 0
        else:
            session_state[eval('"sel" + feature')] = session_state[eval('"uni" + feature')][0] # Selected Value (default)

    # Slide ID Progression Initializeion
    session_state['idxSlide ID'] = 0
    session_state['numSlide ID'] = len(session_state['uniSlide ID'])
    session_state['uniSlide ID_short'] = session_state['uniSlide ID']
    session_state['selSlide ID_short'] = session_state['uniSlide ID_short'][0]

    session_state.prog_left_disabeled = True
    session_state.prog_right_disabeled = False
    if session_state['numSlide ID'] == 1:
        session_state.prog_right_disabeled = True

    # Perform Filtering
    session_state.bc.startTimer()
    df_plot = perform_filtering(session_state)
    # session_state.bc.printElapsedTime(msg = 'Performing Filtering')

    # Set Figure Objects
    session_state.bc.startTimer()
    session_state = setFigureObjs(session_state, df_plot)
    session_state.pointstSliderVal_Sel = session_state.calcSliderVal
    # session_state.bc.printElapsedTime(msg = 'Setting Figure Objects')

    session_state.bc.set_value_df('file', fileName)
    session_state.bc.set_value_df('nSlides', session_state['numSlide ID'])
    session_state.bc.set_value_df('nCells', df_import.shape[0])
    session_state.bc.set_value_df('CellsxSlide', [[session_state.df.loc[session_state.df['Slide ID'] == x, :].shape[0] for x in session_state['uniSlide ID']]])

    # Identify the size of the smallest image
    dataset_img_sizes = [group.shape[0] for ind, group in df_import.groupby('Slide ID')]
    session_state.datafile_min_img_size = min(dataset_img_sizes)

    return session_state

def set_phenotyping_elements(session_state, df_orig):
    """
    To be run each time new data is loaded using the 'Load Data' method
    """

    # Perform pre-processing (phenotying columns, pheno_assign table, pheno_summ table)
    session_state.df_raw, \
    session_state.df, \
    session_state.spec_summ, \
    session_state.pheno_summ = bpl.preprocess_df(df_orig, session_state.marker_names, session_state.marker_pre, session_state.bc)

    # Initalize Custom Phenotyping Variables
    session_state.spec_summ_load       = session_state.spec_summ # Default version that is loaded
    session_state.spec_summ_dataeditor = session_state.spec_summ # Default version that is used for custom phenotyping table

    if 'dataeditor__do_not_persist' in session_state:
        del session_state.dataeditor__do_not_persist

    # Initalize Phenotyping Settings (Radio BUttons)
    session_state.noPhenoOpt = 'Not Selected'
    session_state.phenoMeth  = 'Species'                         # Default when first loaded
    session_state.selected_phenoMeth = session_state.noPhenoOpt  # Default when first loaded

    return session_state

def load_dataset(fiol, dataset_path, files_dict, file_path, loadCompass=False):
    """
    Load selected Dataset (Either from NIDAP or locally). Calling functions from 
    the Foundry IO Library (FIOL).
    
    returns a PANDAS dataframe (df)
    """
    return fiol.load_dataset(dataset_path, files_dict, file_path, loadCompass)

def updatePhenotyping(session_state):
    '''
    Function that is run when changes are made to the phenotyping settings
    of the apps
    Args:
        session_state: Streamlit data structure

    Returns:
        session_state: Streamlit data structure
    '''

    # Create session_state.df
    session_state.df = assign_phenotype_col(session_state.df_raw,
                                            session_state.spec_summ_load,
                                            session_state.selected_phenoMeth,
                                            session_state.marker_names)

    # Initalize Species Summary Table
    session_state.spec_summ    = bpl.init_pheno_assign(session_state.df)

    if 'dataeditor__do_not_persist' in session_state:
        del session_state.dataeditor__do_not_persist

    # session_state.spec_summ_load       = session_state.spec_summ
    session_state.spec_summ_dataeditor = session_state.spec_summ

    # Create Phenotypes Summary Table based on 'phenotype' column in df
    session_state.pheno_summ = bpl.init_pheno_summ(session_state.df)

    # Filtered dataset
    df_plot = perform_filtering(session_state)

    # Update and reset Figure Objects
    session_state = setFigureObjs(session_state, df_plot)

    return session_state

def assign_phenotype_col(df_raw, spec_summ_load, phenoMeth, marker_names):
    """
    Assign a new column to the raw dataset (df_raw) called 'phenotype' based on the 
    phenotyping method selected.

    This function is called within the updatePhenotyping function. It is also responsible
    for creating the 'df' dataframe that is used in the rest of the analysis.

    Args:
        df_raw: Raw dataset
        spec_summ_load: Species Summary Table
        phenoMeth: Phenotyping Method
        marker_names: Marker Names
    
    Returns:
        df: Updated dataset
    """

    df = df_raw.copy()

    if phenoMeth != 'Custom':
        if phenoMeth == 'Species':
            allow_compound_species=True
        elif phenoMeth == 'Marker':
            allow_compound_species=False

        # If compound species are allowed (multiple positive markers are allowed),
        # then we have Will's "exclusive" case; otherwise, it's possible cells are
        # overlapping and we must duplicate the coordinates for the rows having
        # multiple positive markers
        df = bpl.remove_compound_species(df, marker_names, allow_compound_species=allow_compound_species)

        # Assign phenotype column to dataframe based on species name
        df = bpl.assign_phenotype_species(df)
    else:
        # Assign phenotype column to dataframe based on species summary
        df = bpl.assign_phenotype_custom(df, spec_summ_load)

    return df

def perform_filtering(session_state):
    """
    Sets up the filter dictionaries to be used in the 
    filter_dataset step, and then returns the filter_dataset step.
    I suppose I could have combined these two functions, but I like
    keeping the inputs to filter_dataset simple, and using the high level
    session_state input in the top-level function.
    """

    # Create dictionaries of filter types
    session_state = init_filter_struct(session_state,
                                       session_state.SEL_feat,
                                       session_state.CHK_feat)

    # Filter the dataset
    return filter_dataset(session_state.df, session_state.SELdict, session_state.CHKdict)

def init_filter_struct(session_state, SEL_feat, CHK_feat):
    """
    Initalize filtering data structures
    """

    SELdict = dict()
    CHKdict = dict()

    for key in SEL_feat:
        SELdict[f'{key}'] = session_state[eval('"sel" + key')]

    for key in CHK_feat:
        CHKdict[f'{key}'] = session_state[eval('"sel" + key')]

    session_state.SELdict = SELdict
    session_state.CHKdict = CHKdict

    return session_state

def filter_dataset(df, SELdict, CHKdict):
    """
    filter_dataset creates a filtered dataframe based on the an input dataframe (df)
    and dictionaries of filter values. This function is agnostic to the number of unique
    filter/feature combinations in the dictionaries.
    
    SELdict is a dictionary of feature keys and values for filters that are always only 
    1 of many selected from a selectbox
    
    CHKdict is a dictionary of feature keys and values for filters that either check for 
    TRUE or ignored based on a checkbox
    """

    df_filt = df

    # Select box filters
    for selfilt in SELdict:
        df_filt = df_filt[(df_filt[selfilt] == SELdict[selfilt])]

    # Check box filters
    for chkfilt in CHKdict:
        if CHKdict[chkfilt] is True:
            df_filt = df_filt[(df_filt[chkfilt] == CHKdict[chkfilt])]

    return df_filt

def date_time_adjust(df, field):
    """
    Make datetime adjustments to a dataframe (to prevent errors)
    """
    df[field] = pd.to_datetime(df[field])
    return df

def check_upload_df(df, reqFeatures, marker_pre):
    """
    Check if the file to be loaded has apporpriate column names
    """

    hasReqCol  = all(item in df.columns for item in reqFeatures)
    hasMarkers = df.filter(regex='^{}'.format(marker_pre))

    up_file_rdy = False
    if (hasReqCol) & (not hasMarkers.empty):
        up_file_rdy = True
    elif ~(hasReqCol):
        print('Does not have required columns')
    elif hasMarkers.empty:
        print('Marker style is not setup')
    return up_file_rdy

def export_results_dataset(fiol, df, path, filename, saveCompass=False, type = 'S'):
    """
    Export Results/Updated Dataset. Calling functions from the Foundry IO Library (FIOL)
    """
    fiol.export_results_dataset(df, path, filename, saveCompass, type)

def setFigureObjs(session_state, df_plot, InSliderVal = None):
    """
    Organize Figure Objects to be used in plotting
    """

    title = [f'DATASET: {session_state.datafile}',
             f'PHENO METHOD: {session_state.selected_phenoMeth}',
             f'SLIDE ID: {session_state["selSlide ID_short"]}']

    session_state.phenoOrder = list(session_state.pheno_summ.loc[session_state.pheno_summ['phenotype_count'].index, 'phenotype'])

    # num_points
    targ_cell_count = 150000

    num_points = df_plot.shape[0]
    if (num_points > targ_cell_count) & (InSliderVal is None):
        n = targ_cell_count

        calc_slider_val = int(np.ceil(100*n/num_points))
        df_plot = df_plot.sample(n)
        session_state.plotPointsCustom = False
    elif InSliderVal is not None:

        calc_slider_val = InSliderVal
        df_plot = df_plot.sample(frac = calc_slider_val/100)
        session_state.plotPointsCustom = True
    else:
        n = num_points
        calc_slider_val = 100
        session_state.plotPointsCustom = False

    session_state.calcSliderVal = calc_slider_val
    session_state.drawnPoints = df_plot.shape[0]

    # Seaborn
    session_state.phenoFig, session_state.ax = bpl.draw_scatter_fig(figsize=session_state.figsize)
    session_state.phenoFig = bpl.scatter_plot(df_plot, session_state.phenoFig, session_state.ax, title,
                                              xVar = 'Cell X Position', yVar = 'Cell Y Position', hueVar='phenotype',
                                              hueOrder=session_state.phenoOrder)

    # Altair
    # session_state.chart = drawAltairObj(df_plot, title, session_state.phenoOrder, session_state.phenoFig, session_state.ax)

    return session_state

def setFigureObjs_UMAP(session_state, palette = 'tab20'):
    """
    Organize Figure Objects to be used in plotting but for clustering

    Args:
        session_state: Streamlit data structure

    Returns:
        session_state: Streamlit data structure
    """

    title = [f'DATASET: {session_state.datafile}',
             f'PHENO METHOD: {session_state.selected_phenoMeth}',
             f'SLIDE ID: {session_state["selSlide ID_short"]}']

    clust_order = sorted(session_state.spatial_umap.df_umap['clust_label'].unique())
    # Seaborn
    session_state.seabornFig_clust, session_state.ax = bpl.draw_scatter_fig(figsize=session_state.figsize)
    session_state.seabornFig_clust = bpl.scatter_plot(session_state.spatial_umap.df_umap_filt, session_state.seabornFig_clust, session_state.ax, title,
                                                      xVar = 'Cell X Position', yVar = 'Cell Y Position', hueVar = 'clust_label',
                                                      hueOrder=clust_order, palette = palette)

    # Altair
    session_state.altairFig_clust = drawAltairObj(session_state.spatial_umap.df_umap_filt, title, clust_order, session_state.seabornFig_clust, session_state.ax, legendCol = 'clust_label')

    return session_state

def setFigureObjs_UMAPDifferences(session_state):
    '''
    Organize parts of dataframes to be used in the figures
    that are the results of Neighborhood Profile analyses

    Args:
        session_state: Streamlit data structure

    Returns:
        session_state: Streamlit data structure
    '''

    title = [f'DATASET: {session_state.datafile}',
             f'PHENO METHOD: {session_state.selected_phenoMeth}']

    # Full UMAP
    udp_full = session_state.udp_full

    # Inspection UMAP properties
    session_state.umap_ins_msg = None

    udp_ins_raw = copy(udp_full)
    udp_ins_raw.filter_by_lineage(session_state.lineageDisplayToggle, session_state.umapInspect_Ver, session_state.defLineageOpt)

    # Filter by Feature for Inspection
    if session_state.umapInspect_Feat != session_state.defumapOutcomes:

        split_dict_full_ins = udp_ins_raw.split_df_by_feature(session_state.umapInspect_Feat)
        if split_dict_full_ins['appro_feat']:
            # Perform Density Calculations for each Condition
            udp_fals = UMAPDensityProcessing(session_state.npf, split_dict_full_ins['df_umap_fals'], xx=udp_ins_raw.xx, yy=udp_ins_raw.yy)
            udp_true = UMAPDensityProcessing(session_state.npf, split_dict_full_ins['df_umap_true'], xx=udp_ins_raw.xx, yy=udp_ins_raw.yy)

            ## Set Feature Labels
            udp_fals.set_feature_label(session_state.umapInspect_Feat, split_dict_full_ins['fals_msg'])
            udp_true.set_feature_label(session_state.umapInspect_Feat, split_dict_full_ins['true_msg'])

            udp_fals.cluster_dict = udp_ins_raw.cluster_dict
            udp_true.cluster_dict = udp_ins_raw.cluster_dict
            udp_fals.palette_dict = udp_ins_raw.palette_dict
            udp_true.palette_dict = udp_ins_raw.palette_dict

            udp_ins = udp_true
        else:
            udp_ins = udp_ins_raw
            session_state.umap_ins_msg = 'Please choose a boolean or numerical feature'
    else:
        udp_ins = udp_ins_raw

    # Full UMAP figures colored by Density
    if session_state.UMAPFigType == 'Density':

        # All UMAP Figure
        session_state.UMAPFig = udp_full.UMAPdraw_density()

        # UMAP for Lineage/Outcome Inspection
        session_state.UMAPFigInsp = udp_ins.UMAPdraw_density()

    # Full UMAP figures colored by clust_label
    elif session_state.UMAPFigType == 'Clusters':

        # All UMAP Figure
        session_state.UMAPFig = udp_full.umap_draw_clusters()

        # UMAP for Lineage/Outcome Inspection
        session_state.UMAPFigInsp = udp_ins.umap_draw_clusters()

    # Difference UMAP properties
    draw_diff = False
    session_state.umap_diff_msg = None

    udp_diff_raw = copy(udp_full)
    udp_diff_raw.filter_by_lineage(session_state.lineageDisplayToggle, session_state.diffUMAPSel_Ver, session_state.defLineageOpt)

    # Filter by Feature for Inspection
    if session_state.diffUMAPSel_Feat != session_state.defumapOutcomes:
        split_dict_full_diff = udp_diff_raw.split_df_by_feature(session_state.diffUMAPSel_Feat)

        if split_dict_full_diff['appro_feat']:

            # Perform Density Calculations for each Condition
            udp_fals = UMAPDensityProcessing(session_state.npf, split_dict_full_diff['df_umap_fals'], xx=udp_diff_raw.xx, yy=udp_diff_raw.yy)
            udp_true = UMAPDensityProcessing(session_state.npf, split_dict_full_diff['df_umap_true'], xx=udp_diff_raw.xx, yy=udp_diff_raw.yy)

            udp_fals.cluster_dict = udp_diff_raw.cluster_dict
            udp_true.cluster_dict = udp_diff_raw.cluster_dict
            udp_fals.palette_dict = udp_diff_raw.palette_dict
            udp_true.palette_dict = udp_diff_raw.palette_dict

            ## Copy over
            udp_diff = copy(udp_fals)
            ## Perform difference calculation
            udp_diff.dens_mat = np.log10(udp_fals.dens_mat) - np.log10(udp_true.dens_mat)
            ## Rerun the min/max calcs
            udp_diff.umap_summary_stats()

            ## Set Feature Labels
            udp_fals.set_feature_label(session_state.diffUMAPSel_Feat, split_dict_full_diff['fals_msg'])
            udp_true.set_feature_label(session_state.diffUMAPSel_Feat, split_dict_full_diff['true_msg'])
            udp_diff.set_feature_label(session_state.diffUMAPSel_Feat, 'Difference')

            draw_diff = True
        else:
            udp_fals = udp_diff_raw
            udp_true = udp_diff_raw
            udp_diff = udp_diff_raw

            session_state.umap_diff_msg = 'Please choose a boolean or numerical feature'

    else:
        udp_fals = udp_diff_raw
        udp_true = udp_diff_raw
        udp_diff = udp_diff_raw

    session_state.UMAPFigDiff0_Dens = udp_fals.UMAPdraw_density()
    session_state.UMAPFigDiff1_Dens = udp_true.UMAPdraw_density()
    session_state.UMAPFigDiff2_Dens = udp_diff.UMAPdraw_density(diff=draw_diff)

    session_state.UMAPFigDiff0_Clus = udp_fals.umap_draw_clusters()
    session_state.UMAPFigDiff1_Clus = udp_true.umap_draw_clusters()
    session_state.UMAPFigDiff2_Clus = udp_diff.umap_draw_clusters()

    return session_state

def set_figure_objs_clusters_analyzer(session_state):
    '''
    Sets the figures for the Cluster Analyzer Page

    Args:
        session_state: Streamlit data structure
    
    Returns:
        session_state: Streamlit data structure
    '''

    title = [f'DATASET: {session_state.datafile}',
             f'PHENO METHOD: {session_state.selected_phenoMeth}']

    ######## Heatmap/Incidence #########
    df_umap = session_state.spatial_umap.df_umap
    list_clusters = list(session_state.cluster_dict.values())
    list_clusters.remove('No Cluster')

    ### Cluster/Phenotype Heatmap ###
    if session_state.NormHeatRadio == 'Norm within Clusters':
        norm_axis = 0
    elif session_state.NormHeatRadio == 'Norm within Phenotypes':
        norm_axis = 1
    else:
        norm_axis = None

    session_state.heatmapfig = bpl.createHeatMap(df_umap, session_state.pheno_summ['phenotype'], title, norm_axis)

    ### Incidence Line Graph ###
    # Filter by the lineage
    df_umap = filterLineage4UMAP(df_umap, session_state.lineageDisplayToggle_clus, session_state.defLineageOpt, session_state.inciPhenoSel)

    # Set up incidence dataframe
    comp_thresh = None
    inci_df = pd.DataFrame()
    inci_df.index = list_clusters
    inci_df['counts'] = 0
    inci_df['featureCount1'] = 0 # True Condition
    inci_df['featureCount0'] = 0 # False Condition

    # Not Cell Counts
    if session_state.inciOutcomeSel != session_state.definciOutcomes:
        col = df_umap[session_state.inciOutcomeSel]
        if identify_col_type(col) == 'not_bool':
            comp_thresh = 0
            df_umap['chosen_feature'] = df_umap.apply(lambda row: 1 if row[session_state.inciOutcomeSel] >= comp_thresh else 0, axis = 1)
        elif identify_col_type(col) == 'bool':
            df_umap['chosen_feature'] = df_umap[session_state.inciOutcomeSel]
        else:
            df_umap['chosen_feature'] = df_umap[session_state.inciOutcomeSel]

        # Compute the Difference
        for clust_label, group in df_umap.groupby('clust_label'):
            if clust_label != 'No Cluster':
                inci_df.loc[clust_label, 'counts'] = group['chosen_feature'].count()
                inci_df.loc[clust_label, 'featureCount1'] = sum(group['chosen_feature'] == 1)
                inci_df.loc[clust_label, 'featureCount0'] = sum(group['chosen_feature'] == 0)

        inci_df['Count Differences'] = inci_df['featureCount1'] - inci_df['featureCount0']

        sumf1 = sum(inci_df['featureCount1'])
        sumf0 = sum(inci_df['featureCount0'])

        inci_df['Percentages']  = 100*inci_df['featureCount1']/sumf1
        inci_df['Percentages0'] = 100*inci_df['featureCount0']/sumf0

        inci_df['Percentages1_adj'] = 100*(inci_df['featureCount1'] + 1)/(sumf1 + 1*session_state.selected_nClus)
        inci_df['Percentages0_adj'] = 100*(inci_df['featureCount0'] + 1)/(sumf0 + 1*session_state.selected_nClus)

        inci_df['Ratios'] = np.log10(inci_df['Percentages1_adj']/inci_df['Percentages0_adj'])
    # Cell Counts
    else:
        for clust_label, group in df_umap.groupby('clust_label'):
            if clust_label != 'No Cluster':
                inci_df.loc[clust_label, 'counts'] = group['Slide ID'].count()

    # Title
    inci_title = ['Incidence by Cluster']

    # Draw Incidence Figure
    session_state.inciFig = bpl.drawIncidenceFigure(inci_df, inci_title,
                                                    phenotype  = session_state.inciPhenoSel,
                                                    feature    = session_state.inciOutcomeSel,
                                                    displayas  = session_state.Inci_Value_display,
                                                    comp_thresh = comp_thresh)

    return session_state

def filterLineage4UMAP(df, display_toggle, defVal, dropVal):
    '''
    Function for filtering UMAP function based on Phenotypes or Markers
    '''
    if dropVal != defVal:
        if display_toggle == 'Phenotypes':
            df = df.loc[df['Lineage'] == dropVal, :]
        elif display_toggle == 'Markers':
            df = df.loc[df['species_name_short'].str.contains(dropVal), :]

    return df

def drawAltairObj(df, title, sortOrder, fig, ax = None, legendCol='phenotype'):
    """
    Draw Altair Objects
    """
    ## Draw the Scatter Plot
    # Wrap the Title
    wrapTitle = wrapTitleText(title)

    if ax is not None:
        minXLim, maxXLim = ax.get_xlim()
        minYLim, maxYLim = ax.get_ylim()
        bbox = fig.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        width, height = bbox.width*fig.dpi*0.75, bbox.height*fig.dpi*0.75
    else:
        minXLim = df['CentroidX'].min() - 50
        minYLim = df['CentroidY'].min() - 50
        maxXLim = df['CentroidX'].max() + 50
        maxYLim = df['CentroidY'].max() + 50
        width, height = 750, 750

    numLgdCol = 4
    # if len(sortOrder) >= numLgdCol:
    #     sortOrderTran = np.array(sortOrder).reshape(-1, numLgdCol).T.flatten().reshape(-1, numLgdCol).T.flatten()
    # else:
    sortOrderTran = sortOrder

    # Altair Visualization
    selection = alt.selection_point(fields=[legendCol], bind='legend')
    chart = alt.Chart(df).mark_circle(size=3).encode(
            alt.X('CentroidX:Q',
                    scale=alt.Scale(domain=(minXLim, maxXLim)),
                    title='CentroidX (\u03BCm)'),
            alt.Y('CentroidY:Q',
                    scale=alt.Scale(domain=(minYLim, maxYLim)),
                    title='CentroidY (\u03BCm)'),
            color= alt.Color(legendCol, scale=alt.Scale(domain = sortOrder, scheme = 'category20'),
                                          sort=sortOrderTran,
                                          legend=alt.Legend(
                                                            orient='bottom',
                                                            columns = numLgdCol)),
            order=alt.Order('color_phenotype_sort_index:Q'),
            opacity=alt.condition(selection, alt.value(1), alt.value(0.2)),
            tooltip=[legendCol]
            ).properties(width=width,height=height, title=wrapTitle
            ).interactive().add_params(selection)

    # Histogram
    # chart = alt.Chart(st.session_state.df).mark_bar().encode(
    #     alt.X("phenotype:N", bin=True),
    #     y='count()',
    # )

    return chart

def wrapTitleText(title):
    """
    Helps with Wrapping text
    """
    char_lim = 70
    wrap_title = []
    for x in title:
        while len(x) > char_lim:
            x1 = x[:char_lim]
            x2 = x[char_lim:]
            wrap_title.append(x1)
            x = x2
        wrap_title.append(x)

    return wrap_title

def add_item_export_list(session_state, item_name, file_name):
    tempdf = pd.DataFrame(data = {'Item Name' : [item_name],
                                  'File Name' : [file_name],
                                  'Date Time Added': [datetime.now()]})
    session_state.files_to_export = pd.concat([session_state.files_to_export, tempdf]).reset_index(drop=True)

def read_markdown_file(markdown_file):
    '''
    Simple markdown reading function
    '''
    return Path(markdown_file).read_text()

def save_csv(df, df_name):
    '''
    Simple method for saving csv to the output folder
    '''

    output_folder = 'output'
    df.to_csv(f'{output_folder}/{df_name}_{time.strftime("%Y%m%d-%H%M%S")}.csv')

def save_png(img_obj, fig_type, suffix = None):
    '''
    Simple method for saving png to the output folder
    '''

    output_folder = 'output'
    if suffix is not None:
        suffix = '_' + suffix
    file_name_full = f'{output_folder}/{fig_type}_{time.strftime("%Y%m%d-%H%M%S")}{suffix}.png'
    # Save as a png in the local directory using the Matplotlib 'savefig' method
    img_obj.savefig(file_name_full)

def save_png_dataset(fiol, datafile, pngFileName, pltFig):
    """
    Save png image to dataset. Calling functions from the Foundry IO Library (FIOL)

    Args:
        fiol (obj): Foundry IO Library object for handling Palantir SDK calls.
        datafile (str): Path to the dataset that the image will be saved to
        pngFileName (str): Filename for the image, not included the suffix (added later)
        pltfig (obj): Matplotlib figure object to be save as png
    """
    fiol.save_png_dataset(datafile, pngFileName, pltFig)
