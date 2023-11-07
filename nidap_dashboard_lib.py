'''
NIDAP Dashboard Library (NDL) is a series of functions
that supports the types of analysis activities needed by the 
NCI researchers with their data on NIDAP.
'''

import os
import time
import yaml
import numpy as np
import pandas as pd
import altair as alt
from pathlib import Path
alt.data_transformers.disable_max_rows()

# Import relevant libraries
import basic_phenotyper_lib as bpl                  # Useful functions for cell phenotyping
from foundry_IO_lib import foundry_IO_lib           # Foundry Input/Output Class
from benchmark_collector import benchmark_collector # Benchmark Collector Class

def init_session_state(session_state, settings_yaml_file):
    """
    Initialize session_state values for streamlit processing
    """

    session_state.init          = True
    session_state.cpu_pool_size = 8

    # Create an instance of the foundry IO Library
    session_state.fiol = foundry_IO_lib()
    session_state.bc   = benchmark_collector(session_state.fiol)

    # Set benchmarking values
    session_state.bc.set_value_df('counts_multiprocess', True)
    session_state.bc.set_value_df('cpu_pool_size', session_state.cpu_pool_size)

    # Set the directory configurations
    d = os.path.dirname(os.path.abspath(__file__))
    settings_yaml_path = os.path.join(d, settings_yaml_file)
    with open(settings_yaml_path, mode='rt') as file:
        settings = yaml.load(file, yaml.UnsafeLoader)

    # Analysis Settings
    session_state.marker_pre  = settings['analysis']['marker_pre']

    # df Default
    df_dict = {}
    for feature in settings['def_df']:
        df_dict[feature] = settings['def_df'][feature]

    df_default                = pd.DataFrame(data = df_dict)
    session_state.reqFeatures = df_default.columns[:-1]

    # Features for filtering
    session_state.SEL_feat = settings['filt_feat']['SEL_feat']
    session_state.CHK_feat = settings['filt_feat']['CHK_feat']

    # Features for Outcomes Analysis
    session_state.outcomes_BOOL = settings['analysis']['outcomes_BOOL']
    session_state.outcomes_nBOOL = settings['analysis']['outcomes_nBOOL']
    session_state.outcomes_nBOOL_thresh = settings['analysis']['outcomes_nBOOL_thresh']
    session_state.outcomes = session_state.outcomes_BOOL + session_state.outcomes_nBOOL

    # Predefined Project Paths
    session_state.projectPaths   = settings['dir']['projectPaths']
    session_state.usDatasetPaths = settings['dir']['usDatasetPaths']

    # Dataset Dictionaries of files in each unstructure dataset
    session_state.files_dict = {}
    for dataset in session_state.usDatasetPaths:
        session_state.files_dict[dataset] = load_listofFiles(session_state.fiol, dataset)

    # List of DataSets to save CSV to
    session_state.OutputCSVPaths_S = settings['dir']['OutputCSVPaths_S']
    session_state.OutputCSVPaths_U = settings['dir']['OutputCSVPaths_U']

    # List of DataSets to save PNGS to
    session_state.OutputPNGPaths = session_state.OutputCSVPaths_U

    # Selected Dataset Meta Information
    session_state.selectProj = ''
    session_state.datafileS  = ''
    session_state.datafileU  = ''
    session_state.datafile   = ''

    # Default Upload File value
    session_state.uploaded_file = None

    # Output file information
    session_state.df_update_filename_S = ''
    session_state.df_update_filename_U = ''
    session_state.pheno_assign_filename_S = 'phenotype_summary'
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

    # Completing UMAP analysis
    # Has the UMAP been completed yet?
    session_state.umapCompleted = False
    session_state.UMAPFigType = 'Density'

    # UMAP Lineage Display
    session_state.lineageDisplayToggle = 'Phenotypes'

    # Unfiltered dropdown default options
    session_state.defLineageOpt   = 'All Phenotypes'
    session_state.defumapOutcomes = 'No Outcome'
    session_state.definciOutcomes  = 'Cell Counts'

    # Default UMAP dropdown options
    session_state.umapPheno = [session_state.defLineageOpt]
    session_state.umapMarks = [session_state.defLineageOpt]
    session_state.umaplineages = [session_state.defLineageOpt]
    session_state.umapOutcomes = [session_state.defumapOutcomes]

    # Default Incidence dropdown options
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

    # Set data_loaded = False. 
    # This needs to happen at the end to counteract the 'loadDataButton' action
    session_state.data_loaded = False

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

    loadDataSt = time.time()
    session_state.data_loaded = True
    # Prepare dataframe(s) for use in the dashboard and note the time
    session_state.df_raw, \
    session_state.df, \
    session_state.marker_names, \
    session_state.spec_summ, \
    session_state.assign_pheno = prepare_data(df_import)
    prepDataSp = time.time()

    # Meta Data
    session_state.selectProj = projectName # Project Name
    session_state.datafile   = fileName    # File Name
    session_state.spec_summ_load = session_state.spec_summ.copy() # Default version that is loaded
    session_state.spec_summ_dataeditor = session_state.spec_summ.copy() # Default version that is used for custom phenotyping table

    # Default Phenotyping Method (Radio Button)
    session_state.noPhenoOpt = 'Not Selected'
    session_state.phenoMeth = 'Species'                          # Default when first loaded
    session_state.selected_phenoMeth = session_state.noPhenoOpt  # Default when first loaded

    # Clustering (If applicable)
    session_state.selected_nClus = 1

    # Note the time after completing above processing steps
    resetVarsSp = time.time()

    # Filtering
    # All filter categories
    features4filter = session_state.SEL_feat + session_state.CHK_feat
    # Create variables in session state
    for feature in features4filter:
        session_state[eval('"uni" + feature')] = sorted(session_state.df_raw[feature].unique())    # Unique Values
        if feature in session_state.CHK_feat:
            session_state[eval('"sel" + feature')] = 0
        else:
            session_state[eval('"sel" + feature')] = session_state[eval('"uni" + feature')][0] # Selected Value (default)

    # Filtered dataset based on default filter settings and note the time
    session_state.df_filt = perform_filtering(session_state) 
    setfiltSp = time.time()

    # Draw Points Slider
    session_state.pointstSliderVal_Sel = 100
    session_state.calcSliderVal = 100

    # Output File names
    session_state.df_update_filename_S = session_state.datafile + '_updated'
    session_state.df_update_filename_U = session_state.datafile + '_updated'

    # Heatmap Radio
    session_state.NormHeatRadio = 'No Norm'

    # Set Figure Objects and note the time
    session_state = setFigureObjs(session_state) # First view of the filtered datasets (Seaborn/Altair)
    session_state.pointstSliderVal_Sel = session_state.calcSliderVal
    setFigSp = time.time()

    loadDataTD = {'prepData': np.round(prepDataSp - loadDataSt, 3),
                  'resetVar': np.round(resetVarsSp - prepDataSp, 3),
                  'setFilt': np.round(setfiltSp - resetVarsSp, 3),
                  'setFig': np.round(setFigSp - setfiltSp, 3)}
    # print(loadDataTD)

    return session_state

def prepare_data(df_orig):
    """
    To be run each time the 'Load Data' button is hit
    """

    # Time the components of the prepare_data step
    prepDataSt = time.time()

    # Fix column names if needed and note the time
    df_orig = fix_df_cols(df_orig)
    fx_colsSp = time.time()

    # Set df_raw as the baseline dataframe
    df_raw = df_orig.copy()

    # Perform pre-processing (based on app-specific needs)
    df_raw, marker_names, spec_summ, assign_pheno = bpl.preprocess_df(df_raw)
    procDFSP = time.time()

    # Make a copy of df_raw as df
    df = df_raw.copy()

    prepDataTD = {'fx_cols': np.round(fx_colsSp - prepDataSt, 3),
                  'procDF': np.round(procDFSP - fx_colsSp, 3)}
    # print(prepDataTD)

    return df_raw, df, marker_names, spec_summ, assign_pheno

def fix_df_cols(df):
    """
    Fixing column names with goofy (illegal) characters
    """

    char_to_replace = {' (': '_', '(': '_', '): ': '_',  ')': '_', ': ': '_', ':': '_', ' ': '_'}

    colList = list(df.columns)
    newColDict = {}
    for col in colList:
        colNew = col
        for key, value in char_to_replace.items():
            # Replace key character with value character in string
            colNew = colNew.replace(key, value)
        newColDict[col] = colNew
    df = df.rename(newColDict, axis=1)
    return df

def load_dataset(fiol, dataset_path, files_dict, file_path, loadCompass=False):
    """
    Load selected Dataset (Either from NIDAP or locally). Calling functions from 
    the Foundry IO Library (FIOL).
    
    returns a PANDAS dataframe (df)
    """
    return fiol.load_dataset(dataset_path, files_dict, file_path, loadCompass)

def updatePhenotyping(session_state):
    '''
    
    '''
    session_state.selected_phenoMeth = session_state.phenoMeth
    session_state.df = changePhenoMeth(session_state.df_raw,
                                       session_state.spec_summ_load,
                                       session_state.selected_phenoMeth,
                                       session_state.marker_names)

    session_state.spec_summ = bpl.init_species_summary(session_state.df)
    session_state.assign_pheno = bpl.init_assign_pheno(session_state.df)
    session_state.spec_summ_cuTb = session_state.spec_summ.copy()

    # Perform filtering
    session_state.df_filt = perform_filtering(session_state)

    # Set Figure Objects
    session_state = setFigureObjs(session_state)

    return session_state

def changePhenoMeth(df_raw, spec_summ_load, phenoMeth, marker_names):
    """
    Toggle for changing the phenotyping method
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

        # Add a pretty phenotype column to the dataframe
        df = bpl.species_as_phenotype(df)
    else:
        # Update df based on previously saved species_summ
        df = bpl.update_df_phenotype(df, spec_summ_load)

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
    session_state = init_filter_struct(session_state)

    # Filter the dataset
    return filter_dataset(session_state.df, session_state.SELdict, session_state.CHKdict)

def init_filter_struct(session_state):
    """
    Initalize filtering data structures
    """

    SELdict = dict()
    CHKdict = dict()

    for key in session_state.SEL_feat:
        SELdict['{}'.format(key)] = session_state[eval('"sel" + key')]

    for key in session_state.CHK_feat:
        CHKdict['{}'.format(key)] = session_state[eval('"sel" + key')]

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

    df_filt = df.copy()

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
    return up_file_rdy

def export_results_dataset(fiol, df, path, filename, saveCompass=False, type = 'S'):
    """
    Export Results/Updated Dataset. Calling functions from the Foundry IO Library (FIOL)
    """
    fiol.export_results_dataset(df, path, filename, saveCompass, type)

def setFigureObjs(session_state, InSliderVal = None):
    """
    Organize Figure Objects to be used in plotting
    """

    title = [f'PROJECT Path: {session_state.selectProj}', 
             f'DATASET: {session_state.datafile}',
             f'PHENO METHOD: {session_state.selected_phenoMeth}']

    session_state.phenoOrder = list(session_state.assign_pheno.loc[session_state.assign_pheno['phenotype_count'].index, 'phenotype'])

    # NumPoints
    targCellCount = 150000 
    df_plot = session_state.df_filt.copy()

    numPoints = session_state.df_filt.shape[0]
    if (numPoints > targCellCount) & (InSliderVal is None):
        n = targCellCount

        calcSliderVal = int(np.ceil(100*n/numPoints))
        df_plot = df_plot.sample(n)
        session_state.plotPointsCustom = False
    elif InSliderVal is not None:

        calcSliderVal = InSliderVal
        df_plot = df_plot.sample(frac = calcSliderVal/100)
        session_state.plotPointsCustom = True
    else:
        n = numPoints
        calcSliderVal = 100
        session_state.plotPointsCustom = False

    session_state.calcSliderVal = calcSliderVal
    session_state.drawnPoints = df_plot.shape[0]

    # Seaborn
    session_state.phenoFig, session_state.ax = bpl.draw_scatter_fig(figsize=session_state.figsize)
    session_state.phenoFig = bpl.scatter_plot(df_plot, session_state.phenoFig, session_state.ax, title, 
                                              xVar = 'CentroidX', yVar = 'CentroidY', hueVar='phenotype', 
                                              hueOrder=session_state.phenoOrder)

    # Altair
    session_state.chart = drawAltairObj(df_plot, title, session_state.phenoOrder, session_state.phenoFig, session_state.ax)

    return session_state

def setFigureObjs_UMAP(session_state):
    """
    Organize Figure Objects to be used in plotting but for clustering
    """

    title = [f'PROJECT Path: {session_state.selectProj}', 
             f'DATASET: {session_state.datafile}',
             f'PHENO METHOD: {session_state.selected_phenoMeth}']

    clustered_cells = session_state.spatial_umap.cells.loc[session_state.spatial_umap.cells.loc[:, 'umap_test'] == True, :]
    clustOrder = sorted(clustered_cells['clust_label'].unique())
    # Seaborn
    session_state.seabornFig_clust, session_state.ax = bpl.draw_scatter_fig(figsize=session_state.figsize)
    session_state.seabornFig_clust = bpl.scatter_plot(clustered_cells, session_state.seabornFig_clust, session_state.ax, title,
                                                      xVar = 'CentroidX', yVar = 'CentroidY', hueVar = 'clust_label',
                                                      hueOrder=clustOrder)

    # Altair
    session_state.altairFig_clust = drawAltairObj(clustered_cells, title, clustOrder, session_state.seabornFig_clust, session_state.ax, legendCol = 'clust_label')

    return session_state

def setFigureObjs_UMAPDifferences(session_state):
    import matplotlib.pyplot as plt

    title = [f'PROJECT Path: {session_state.selectProj}',
             f'DATASET: {session_state.datafile}',
             f'PHENO METHOD: {session_state.selected_phenoMeth}']

    dfUMAP = pd.DataFrame(data = session_state.spatial_umap.umap_test, columns = ['X', 'Y'])
    dfUMAP['Cluster'] = session_state.spatial_umap.cells['clust_label'].values[session_state.spatial_umap.cells['umap_test']]
    dfUMAP['Lineage'] = session_state.spatial_umap.cells['Lineage'].values[session_state.spatial_umap.cells['umap_test']]
    for outcome in session_state.outcomes:
        dfUMAP[outcome] = session_state.spatial_umap.cells[outcome].values[session_state.spatial_umap.cells['umap_test']]
    clustOrder = sorted(dfUMAP['Cluster'].unique())

    n_bins = 200
    xx = np.linspace(np.min(dfUMAP['X']), np.max(dfUMAP['X']), n_bins + 1)
    yy = np.linspace(np.min(dfUMAP['Y']), np.max(dfUMAP['Y']), n_bins + 1)
    n_pad = 30

    minXY = dfUMAP[['X', 'Y']].min()-1
    maxXY = dfUMAP[['X', 'Y']].max()+1

    dfUMAPI = dfUMAP.copy()
    dfUMAPD = dfUMAP.copy()

    # Lineage filtering
    dfUMAPI = filterLineage4UMAP(session_state, dfUMAPI, session_state.defLineageOpt, session_state.umapInspect_Ver)
    dfUMAPD = filterLineage4UMAP(session_state, dfUMAPD, session_state.defLineageOpt, session_state.diffUMAPSel_Ver)

    vlim = .97

    # Inspection UMAP properties
    if session_state.umapInspect_Feat != session_state.defumapOutcomes:
        w_Ins = dfUMAPI[session_state.umapInspect_Feat]
        w_Ins, dfUMAPI = bpl.preprocess_weighted_umap(w_Ins, dfUMAPI)
    else:
        w_Ins = None

    # Difference UMAP properties
    if session_state.diffUMAPSel_Feat != session_state.defumapOutcomes:
        w = dfUMAPD[session_state.diffUMAPSel_Feat]
        if session_state.diffUMAPSel_Feat in session_state.outcomes_nBOOL:
            compThresh = 0
            w = np.array(w > compThresh).astype('int')

            featComp1 = f'>= {compThresh}'
            featComp2 = f'< {compThresh}'

            dfUMAPD_A = dfUMAPD.loc[dfUMAPD[session_state.diffUMAPSel_Feat] >= compThresh, :]
            dfUMAPD_B = dfUMAPD.loc[dfUMAPD[session_state.diffUMAPSel_Feat] < compThresh, :]
            dfUMAPD_AB = dfUMAPD_B.copy()
        else:
            featComp1 = '= 1'
            featComp2 = '= 0'

            dfUMAPD_A = dfUMAPD.loc[dfUMAPD[session_state.diffUMAPSel_Feat] == 1, :]
            dfUMAPD_B = dfUMAPD.loc[dfUMAPD[session_state.diffUMAPSel_Feat] == 0, :]
            dfUMAPD_AB = dfUMAPD_B.copy()

        w, dfUMAPD = bpl.preprocess_weighted_umap(w, dfUMAPD)

        w_DiffA = w
        w_DiffB = max(w) - w
        w_Diff  = w_DiffA - w_DiffB

        feat_label0 = f'{session_state.diffUMAPSel_Feat} {featComp1} '
        feat_label1 = f'{session_state.diffUMAPSel_Feat} {featComp2} '
        feat_label2 = None

    else:
        w_DiffA = None
        w_DiffB = None
        w_Diff  = None

        feat_label0 = None
        feat_label1 = None
        feat_label2 = None

        dfUMAPD_A = dfUMAPD.copy()
        dfUMAPD_B = dfUMAPD.copy()
        dfUMAPD_AB = dfUMAPD.copy()

    feat_labels = [feat_label0, feat_label1, feat_label2]
    dfUMAPDs = [dfUMAPD_A, dfUMAPD_B, dfUMAPD_AB]

    # UMAP colored by Density
    if session_state.UMAPFigType == 'Density':
        w = None

        # All UMAP Figure
        session_state.UMAPFig     = bpl.UMAPdraw_density(dfUMAP, bins = [xx, yy], w=None, n_pad=n_pad, vlim=vlim)

        # UMAP for Lineage/Outcome Inspection
        session_state.UMAPFigInsp = bpl.UMAPdraw_density(dfUMAPI, bins = [xx, yy], w=w_Ins, n_pad=n_pad, vlim=vlim)

    # UMAP colored by Clusters
    elif session_state.UMAPFigType == 'Clusters':
        # Make a new dataframe to send to the phenotyping library scatterplot function

        # All UMAP Figure
        session_state.UMAPFig, session_state.UMAPax = bpl.draw_scatter_fig(figsize=session_state.figsize)
        session_state.UMAPFig = bpl.scatter_plot(dfUMAP, session_state.UMAPFig, session_state.UMAPax, title,
                                                 xVar = 'X', yVar = 'Y', hueVar = 'Cluster',
                                                 hueOrder = clustOrder,
                                                 xLim = [minXY[0], maxXY[0]], yLim = [minXY[1], maxXY[1]], boxoff=True, clusters_label = True)
        
        # UMAP for Lineage/Outcome Inspection
        session_state.UMAPFigInsp, session_state.UMAPInspax = bpl.draw_scatter_fig(figsize=session_state.figsize)
        session_state.UMAPFigInsp = bpl.scatter_plot(dfUMAPI, session_state.UMAPFigInsp, session_state.UMAPInspax, title,
                                                 xVar = 'X', yVar = 'Y', hueVar = 'Cluster',
                                                 hueOrder = clustOrder, 
                                                 xLim = [minXY[0], maxXY[0]], yLim = [minXY[1], maxXY[1]], boxoff=True, clusters_label = True)

    # UMAP Difference Figures
    session_state.UMAPFigDiff0_Dens = bpl.UMAPdraw_density(dfUMAPD, bins = [xx, yy], w=w_DiffA, n_pad=n_pad, vlim=vlim, feat = feat_label0)
    session_state.UMAPFigDiff1_Dens = bpl.UMAPdraw_density(dfUMAPD, bins = [xx, yy], w=w_DiffB, n_pad=n_pad, vlim=vlim, feat = feat_label1)
    session_state.UMAPFigDiff2_Dens = bpl.UMAPdraw_density(dfUMAPD, bins = [xx, yy], w=w_Diff, n_pad=n_pad, vlim=vlim, diff = True)

    # UMAP Difference Figures
    for i in range(3):
        fig, ax = bpl.draw_scatter_fig(figsize=session_state.figsize)

        fig = bpl.scatter_plot(dfUMAPDs[i], fig, ax, title,
                                xVar = 'X', yVar = 'Y', hueVar = 'Cluster',
                                hueOrder = clustOrder, boxoff=True, 
                                xLim = [minXY[0], maxXY[0]], yLim = [minXY[1], maxXY[1]],
                                feat = feat_labels[i], small_ver = True, clusters_label = True)
        
        session_state[eval('"UMAPFigDiff" + str(i) + "_Clus"')] = fig
        session_state[eval('"UMAPax" + str(i)')] = ax

    ### Cluster/Phenotype Heatmap ###
    if session_state.NormHeatRadio == 'Norm within Clusters':
        normAxis = 0
    elif session_state.NormHeatRadio == 'Norm within Phenotypes':
        normAxis = 1
    else:
        normAxis = None

    session_state.heatmapfig = bpl.createHeatMap(session_state.spatial_umap, title, normAxis)

    ### Incidence Line Graph ###
    # Create a cell df based on the cells that performed in the UMAP_test
    cellsUMAP = session_state.spatial_umap.cells.loc[session_state.spatial_umap.cells['umap_test'] == True, :]
    # Filter by the lineage
    if session_state.inciPhenoSel != session_state.defLineageOpt:
        cellsUMAP = cellsUMAP.loc[cellsUMAP['Lineage'] == session_state.inciPhenoSel, :]
    
    # Not Cell Counts
    if session_state.inciOutcomeSel != session_state.definciOutcomes:
        # Remake nonboolean variable into a boolean.
        if session_state.inciOutcomeSel in session_state.outcomes_nBOOL:
            compThresh = 0
            cellsUMAP[session_state.inciOutcomeSel] = cellsUMAP.apply(lambda row: 1 if row[session_state.inciOutcomeSel] >= compThresh else 0, axis = 1)
        else:
            compThresh = None

        # Computer the Difference
        if session_state.Inci_Value_display == 'Count Differences':
            inciDF = cellsUMAP.groupby('clust_label')[session_state.inciOutcomeSel].agg(lambda x: sum(x) - (len(x) -sum(x)))
        elif session_state.Inci_Value_display == 'Ratios':
            inciDF = cellsUMAP.groupby('clust_label')[session_state.inciOutcomeSel].agg(lambda x: np.log10((sum(x==1))/(sum(x==0) + 1)))
        elif session_state.Inci_Value_display == 'Percentages':
            inciDF = cellsUMAP.groupby('clust_label')[session_state.inciOutcomeSel].agg(lambda x: 100*sum(x)/len(x))
    
    # Cell Counts
    else:
        compThresh = None
        inciDF = cellsUMAP.groupby('clust_label')['ID'].count()

    # Make a Common x range
    commonIdx = np.arange(0, session_state.selected_nClus, 1.0)
    inciTitle = [f'Incidence by Cluster']


    # Draw Incidence Figure
    session_state.inciFig = bpl.drawIncidenceFigure(commonIdx, inciDF, inciTitle, phenotype=session_state.inciPhenoSel,
                                                    outcome = session_state.inciOutcomeSel, compThresh=compThresh, 
                                                    displayas=session_state.Inci_Value_display)

    return session_state

def filterLineage4UMAP(session_state, df, defVal, dropVal):
    '''
    Function for filtering UMAP function based on Phenotypes or Markers
    '''
    if dropVal != defVal:
        if session_state.lineageDisplayToggle == 'Phenotypes':
            df = df.loc[df['Lineage'] == dropVal, :]
        elif session_state.lineageDisplayToggle == 'Markers':
            df = df.loc[df['Lineage'].str.contains(f'{dropVal}'), :]

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

def read_markdown_file(markdown_file):
    '''
    Simple markdown reading function
    '''
    return Path(markdown_file).read_text()

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
