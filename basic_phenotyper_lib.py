'''
The Basic Phenotyping Library (BPL) handles the basic functions
required for phenotyping 

'''

import time
import math
import warnings
import numpy as np
import pandas as pd
import umap  # slow
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans # K-Means
from SpatialUMAP import SpatialUMAP
import PlottingTools as umPT
import utils

warnings.simplefilter(action='ignore', category= FutureWarning)
warnings.filterwarnings("ignore", message=".*The 'nopython' keyword.*")
pd.options.mode.chained_assignment = None  # default='warn'

def preprocess_df(df_orig, marker_names, marker_col_prefix, bc):
    '''Perform some preprocessing on our dataset to apply tranforms
    and collect meta-data
    '''

    # Step 1: Initalize the columns needed for phenotyping
    bc.startTimer()
    df_raw = init_pheno_cols(df_orig, marker_names, marker_col_prefix)
    bc.printElapsedTime(msg = '     Initializing Phenotyping Columns')

    # Step 2: Intialize Phenotying Assignments dataframe
    bc.startTimer()
    pheno_assign = init_pheno_assign(df_raw)
    bc.printElapsedTime(msg = '     Initializing Phenotyping Assignments Table')

    # Step 3: Intialize Phenotype Summary Dataframe (based on Phenotying Assigments)
    bc.startTimer()
    pheno_summ = init_pheno_summ(df_raw)
    bc.printElapsedTime(msg = '     Initializing Phenotype Summary Table')

    return df_raw, df_raw, pheno_assign, pheno_summ

def identify_marker_columns(df, marker_col_prefix):
    '''
    Identify the marker columns that are present in the base dataset
    '''
    # Create a dataframe of just the marker columns
    df_markers = df.filter(regex='^{}'.format(marker_col_prefix))

    # Get a list of the markers in the datafile
    return [x.removeprefix(marker_col_prefix) for x in df_markers.columns]

def init_pheno_cols(df, marker_names, marker_col_prefix):
    """Add a column to the dataframe containing a string of the marker bits in the same order as the 
    also-returned marker_names list.

    Args:
        df (Pandas dataframe): Dataframe containing data from the input dataset
        marker_col_prefix (str): Prefix to the marker columns after which is only the marker name

    Returns:
        Pandas dataframe: Dataframe with the marker bits column appended
        list: List of marker names in the dataset
    """

    # Create a dataframe of just the marker columns
    marker_cols = [marker_col_prefix + x for x in marker_names]
    df_markers = df[marker_cols]

    # Add a column to the original dataframe containing a concatenation of the bits in the marker columns to a single string, e.g., '0110'
    # Previously called Species String
    # This was previously really slow. Code basically taken from new_phenotyping_lib.py
    marker_cols_first_row = df_markers.iloc[0, :].to_list()  # get just the first row of marker values
    if (0 not in marker_cols_first_row) and (1 not in marker_cols_first_row):

        # Null values in df_markers will break the .map() step so check for and remove them here
        ser_num_of_null_rows_in_each_column = df_markers.isnull().sum()
        if ser_num_of_null_rows_in_each_column.sum() != 0:

            # For the time being, import Streamlit so warnings can be rendered. Otherwise, this file does not import streamlit and it should remain that way but this is a minimal fix for the time being
            import streamlit as st

            st.warning('Null values have been detected in the phenotype columns. Next time, please check for and remove null rows in the datafile unification step (File Handling > Datafile Unification). We are removing them for you now but it would be *much* better to do this in the Datafile Unifier now! Otherwise, downstream functionality may not work. Here are the numbers of null rows found in each column containing them:')
            ser_num_of_null_rows_in_each_column.name = 'Number of null rows'
            st.write(ser_num_of_null_rows_in_each_column[ser_num_of_null_rows_in_each_column != 0])

            # Perform the operation
            row_count_before = len(df)
            df = df.dropna(subset=marker_cols)
            row_count_after = len(df)

            # Display a success message
            st.write(f'{row_count_before - row_count_after} rows deleted')

            # Update df_markers
            df_markers = df[marker_cols]

        df_markers = df_markers.map(lambda x: {'+': '1', '-': '0'}[x[-1]])
    df['mark_bits'] = df_markers.astype(str).apply(''.join, axis='columns')  # efficiently create a series of strings that are the columns (in string format) concatenated together

    # Add a column of prettier names for the species, e.g., 'VIM- ECAD+ COX2+ NOS2-'
    df['species_name_long'] = df['mark_bits'].apply(lambda mark_bits: ' '.join([marker_name + ('+' if marker_bit == '1' else '-') for marker_name, marker_bit in zip(marker_names, mark_bits)]))

    # Add a column dropping the negative markers from these pretty names, e.g., 'ECAD+ COX2+'
    def species_name_long_to_short(species_name_long):
        x = '+ '.join([marker_names[iy] for iy, y in enumerate([x for x in species_name_long if x in ('+', '-')]) if y == '+']) + '+'
        species_name_short = x if len(x) != 1 else 'Other'
        return species_name_short
    # This can possibly be made faster (if it's correct) via but I haven't tested it:
        # marker_indices = [i for i, x in enumerate(species_name_long) if x == '+']
        # if not marker_indices:
        #     return 'Other'
        # return ' + '.join(marker_names[i] for i in marker_indices) + '+'
    df['species_name_short'] = df['species_name_long'].apply(species_name_long_to_short)

    # Create a new column called 'has pos mark' identifying which species_name_shorts are not Other
    df['has_pos_mark'] = True
    df.loc[df['species_name_short'] == 'Other', 'has_pos_mark'] = False

    # Create phenotype column and assign a value of 'unassigned'
    df['phenotype'] = 'unassigned'

    # Create cluster column and assign a value of '-1'
    df['cluster'] = -1

    # Return the dataframe with the marker bits column appended as well as the list of marker names
    return df

def init_pheno_assign(df):
    '''For each unique species (elsewhere called "exclusive" phenotyping), 
    generate information concerning their prevalence in a new dataframe.

    Args:
        df (Pandas dataframe): Dataframe containing data from the input dataset, 
                               including a "mark_bits" column

    Returns:
        spec_summ (Pandas dataframe): Dataframe containing the value counts 
                                      of each "exclusive" species
    '''

    st_init_species = time.time()
    spec_summ = df[['species_name_short', 'phenotype', 'species_name_long']]
    sp_init_species = time.time()
    elapsed = round(sp_init_species - st_init_species, 3)
    print(f'        Initalizing Phenotying Assignments: {elapsed}s')

    # This line seems to throw a TypeError: unhashable type: 'numpy.ndarray' error
    spec_summ['species_count'] = spec_summ['species_name_short'].groupby(spec_summ['species_name_short']).transform('count')
    spec_summ = spec_summ.drop_duplicates().reset_index(drop=True)

    # The above seems a bit inefficient and should probably be replaced with something like this:
    # spec_summ = spec_summ['species_name_short'].value_counts().reset_index()
    # spec_summ.columns = ['species_name_short', 'species_count']

    sp_species_count = time.time()
    elapsed_counts = round(sp_species_count - sp_init_species, 3)
    print(f'        Phenotying Assignments Counts Calculations: {elapsed_counts}s')

    spec_summ['species_percent'] = [round(100*x/sum(spec_summ['species_count']), 2) for x in spec_summ['species_count']]
    sp_species_per = time.time()
    elapsed_per = round(sp_species_per - sp_species_count, 3)
    print(f'        Phenotying Assignments Percents Calculations: {elapsed_per}s')

    spec_summ = spec_summ.sort_values(by='species_count', ascending= False).reset_index(drop=True)
    sp_species_sort = time.time()
    elapsed_sort = round(sp_species_sort - sp_species_per, 3)
    print(f'        Phenotying Assignments sorting: {elapsed_sort}s')

    # Return the created dataframe
    return spec_summ

def init_pheno_summ(df):
    '''For each unique species (elsewhere called "exclusive" phenotyping),
    generate information concerning their prevalence in a new dataframe.

    Args:
        df (Pandas dataframe): Dataframe containing data from the input dataset,
                                including a "mark_bits" column

    Returns:
        assign_pheno (Pandas dataframe): Dataframe containing the value counts of
                                        each "exclusive" species
    '''

    assign_pheno = df[['phenotype', 'species_name_short', 'species_name_long']].groupby(by='phenotype', as_index = False).agg(lambda x: np.unique(list(x)))

    assign_pheno['phenotype_count'] = [sum(df['phenotype'] == x) for x in assign_pheno.phenotype]
    assign_pheno['phenotype_percent'] = [round(100*x/sum(assign_pheno['phenotype_count']), 2) for x in assign_pheno['phenotype_count']]
    assign_pheno = assign_pheno.sort_values(by='phenotype_count', ascending=False)

    return assign_pheno

def remove_compound_species(df, marker_names, allow_compound_species=True):
    '''
    For each compound species ('Species int' not just a plain power of two), add each 
    individual phenotype to the end of the dataframe individually 
    and then delete the original compound entry
    '''

    # Remove compound species if requested
    if not allow_compound_species:

        marker_names = np.array(marker_names)

        df['mark_bits'] = [[bool(int(y)) for y in x] for x in df['mark_bits']]
        df['species_name_short'] = [marker_names[x] for x in df['mark_bits']]
        df = df.explode('species_name_short', ignore_index = True)

        df = df.dropna(subset=['species_name_short']).reset_index(drop=True)
        df['species_name_long'] = df['species_name_short']

    return df

def incorporate_phenotype_identifications(df_objects_orig, tsv_file=None):
    """Take hyphen-separated phenotypes ("reduced_marker_set" below) specified by biologists 
    in Excel for sets of surface markers exported from the plotting map and incorporate these 
    identifications into the main data object
    Workflow:
        (1) Run tci.TIMECellInteraction() as usual
        (2) Export the plotting map using utils.export_plotting_map()
        (3) Ask biologists to fill out the resulting CSV file in Excel in order to identify the phenotypes
        (4) Delete the .pkl files created so far
        (5) Run the entire pipeline again from scratch, this time adding, e.g., "phenotype_identification_tsv_file='gmb_phenotype_ids_as_of_2022-06-06.tsv'" in the arguments to tci.preprocess_dataset()

    Originally taken from utilities/utils.py.

    Args:
        df_objects_orig (Pandas dataframe): Input dataframe
        tsv_file (str): Name of the TSV file generated by pasting into a text editor columns 
                        from Excel in the order specified in the read_csv() method below

    NOTE: We're not yet using this function though we should!
    """

    if tsv_file is not None:

        # Import relevant libraries
        import pandas as pd

        # Define a hardcoded constant that in the future I can put in the same place as
        # phenotype_identification_tsv_file='gmb_phenotype_ids_as_of_2022-06-06.tsv'
        species_id_col = 'Species int'

        # Read in the data from the filled-out-by-collaborators TSV file specifying the phenotype identifications 
        # from the surface markers, setting the 'species_id' as both the index and keeping the values in a column, 
        # dropping the unnecessary "species_perc" column
        df_phenotypes = pd.read_csv(tsv_file, sep='\t', names=['species_id', 'species_count', 'species_percent', 'positive_markers', 'reduced_marker_set']).set_index('species_id', drop=False).drop('species_percent', axis='columns')

        # Initialize a dataframe containing just the final set of species, using as initial data the "first" rows of the phenotype data grouped by species
        df_reduced_species = df_phenotypes[df_phenotypes['reduced_marker_set'].apply(lambda x: '-' not in x)].groupby(by='reduced_marker_set').first()

        # Add a column to the phenotypes dataframe with the species IDs to which we want to copy the object data
        df_phenotypes['ids_to_map'] = df_phenotypes['reduced_marker_set'].apply(lambda x: [df_reduced_species['species_id'].loc[y] for y in x.split('-')])

        # Clean up the dataframe of the final phenotypes, dropping the "positive_markers" column (since it only contains 
        # one of multiple possible sets), initializing the species counts to zero, and setting the index to the species_id,
        # which is unique though corresponds to the single set of positive markers, which again is only one of multiple possible sets
        df_reduced_species = df_reduced_species.drop('positive_markers', axis='columns')
        df_reduced_species['species_count'] = 0
        df_reduced_species = df_reduced_species.reset_index().set_index('species_id')

        # Store the value counts of the original set of species in the original dataset
        orig_species_value_counts = df_objects_orig[species_id_col].value_counts()

        # Check that the total number of objects in the original dataset is consistent using multiple measures. 
        # Note that this can get flagged if the dataset is inconsistent with the TSV file used for phenotype identification,
        # which hardcodes some properties of the original dataset, like the total number of cells of each naiive phenotype.
        tot_number_of_objects_orig = [orig_species_value_counts.sum(), df_phenotypes['species_count'].sum(), len(df_objects_orig)]
        if len(set(tot_number_of_objects_orig)) != 1:
            print('WARNING: Inconsistent total number of original objects! ({})'.format(tot_number_of_objects_orig))

        # Initialize some variables
        size_of_orig_data = 0  # stores total number of objects in original dataset
        size_of_new_data  = 0  # stores total number of objects in new dataset
        list_objects_new  = [] # stores dataframes to be combined into the new dataset

        # For each row in the original set of phenotypes...
        for row in df_phenotypes.iterrows():

            # Get the current species ID and phenotype data for that ID, in particular, 
            # the species IDs to which to copy the current set of object data
            species_id, row_data = row
            ids_to_map = row_data['ids_to_map']

            # Print out what we're doing
            print('Now copying data for objects with positive surface markers {} to the phenotype(s) {}'.format(row_data['positive_markers'], row_data['reduced_marker_set']))

            # Get a boolean series of the original object data that correspond to the current species ID
            object_is_species = df_objects_orig[species_id_col] == species_id

            # Determine the total number of matches for the current species ID
            curr_num_objects = object_is_species.sum()

            # Check that the current number of object matches is consistent with what was previously calculated using the value_counts() method
            assert curr_num_objects == orig_species_value_counts[species_id], 'ERROR: Current number of species matches is inconsistent! ({}, {})'.format(curr_num_objects, orig_species_value_counts[species_id])

            # Store a copy of the object data for the current species ID
            curr_object_data = df_objects_orig[object_is_species].copy()

            # Update the total number of original and new objects
            size_of_orig_data = size_of_orig_data + curr_num_objects
            size_of_new_data = size_of_new_data + curr_num_objects * len(ids_to_map)

            # For each species ID to which to map the current object data...
            for id_to_map in ids_to_map:

                # Print out what we're doing
                print('  surface markers {} (ID {}) --> phenotype {} (ID {}) '.format(row_data['positive_markers'], species_id, df_phenotypes.loc[id_to_map, 'reduced_marker_set'], id_to_map))

                # Update the total numbers of the new, mapped species
                df_reduced_species.loc[id_to_map, 'species_count'] = df_reduced_species.loc[id_to_map, 'species_count'] + curr_num_objects

                # Overwrite the species IDs for the current data dataframe with the new species ID
                curr_object_data[species_id_col] = id_to_map

                # Save the current object data with updated species ID to the dataframe
                list_objects_new.append(curr_object_data.copy())

        # Create new dataframe with the new phenotypes
        df_objects_new = pd.concat(list_objects_new)

        # Sort the new species holder by decreasing frequency in the entire dataset
        df_reduced_species = df_reduced_species.sort_values(by='species_count', ascending=False)

        # Run checks on the sizes of the old and new datasets
        assert size_of_orig_data == len(df_objects_orig), 'ERROR: Inconsistent total number of original objects (second check)! ({}, {})'.format(size_of_orig_data, len(df_objects_orig))
        assert size_of_new_data == len(df_objects_new), 'ERROR: Inconsistent total number of new objects! ({}, {})'.format(size_of_new_data, len(df_objects_new))

        # Visually compare the species counts calculated two different ways
        print('NOTE: Visually ensure the following two structures are the same (automate this in the future)!')
        print(df_reduced_species)
        print(df_objects_new[species_id_col].value_counts())

    else:  # if tsv_file is None
        df_objects_new = df_objects_orig

    # Save the calculated data as properties of the class object;
    # actually as this doesn't seem to be working,
    # just return these four variables and overwrite them in the class instantiation;
    # I may have run into this issue before
    return df_objects_new

def assign_phenotype_species(df):
    """Add a "phenotype" column to df based on the from the marker columns 
    via the calculated species integer.

    Args:
        df (Pandas dataframe): Dataframe containing the "species_name_short" column

    Returns:
        df (Pandas dataframe): Input dataframe with an added "phenotype" column 
                               appended or overwritten
    """

    # Create a "phenotype" column based on the from the "species_name_short"
    df['phenotype'] = df['species_name_short']

    # Return dataframe with phenotype column
    return df

def assign_phenotype_custom(df, spec_summ):
    '''Add a "phenotype" column to df based on a species summary dataframe 
    which identfies custom phenotype assignments from custom phenotyping

    Args:
        df (Pandas dataframe): Dataframe containing the "species_name_short" column
        spec_summ (Pandas dataframe): Dataframe containing the custom assignments 
                                      of phenotypes based on "species_name_short"

    Returns:
        df (Pandas dataframe): Input dataframe with an added "phenotype" column 
                               appended or overwritten
    '''
    # Create a "phenotype" column mapped from a species summary dataset
    df['phenotype'] = df['species_name_short'].map(dict(zip(spec_summ['species_name_short'].to_list(), spec_summ['phenotype'].to_list())))

    df.loc[df['phenotype'].isna(), 'phenotype'] = 'unassigned'

    # Return dataframe with phenotype column
    return df

def load_previous_species_summary(filename):
    '''
    Load previous species summary file
    '''

    spec_summ_load = pd.read_csv(filename)

    spec_summ = spec_summ_load

    return spec_summ

def draw_pheno_summ_bar_fig(pheno_summ, omit_other):
    
    pheno_order = pheno_summ['phenotype'].tolist()

    slc_bg   = '#0E1117'  # Streamlit Background Color
    slc_text = '#FAFAFA'  # Streamlit Text Color

    # Read in the tab20 palette from seaborn and convert it to a format suitable for plotly
    palette = sns.color_palette('tab20')[0:len(pheno_summ)]
    palette = [f'rgba({int(r*255)}, {int(g*255)}, {int(b*255)}, 1)' for r, g, b in palette]

    if omit_other:
        
        pheno_summ = pheno_summ[pheno_summ['phenotype'] != 'Other']
        idx = pheno_order.index('Other')
        pheno_order = pheno_order[:idx] + pheno_order[idx+1:]
        palette = palette[:idx] + palette[idx+1:]

    fig = go.Figure(
        data=[
            go.Bar(
                x=pheno_summ['phenotype'],
                y=pheno_summ['phenotype_count'],
                marker_color=palette,
                hovertemplate='<b>Phenotype:</b> %{x}<br><b>Count:</b> %{y}<extra></extra>'
            )
        ]
    )
    fig.update_layout(
        title='Phenotype Counts',
        xaxis_title='Phenotype',
        yaxis_title='Count',
        plot_bgcolor=slc_bg,
        paper_bgcolor=slc_bg,
        font=dict(color=slc_text)
    )

    return fig

def draw_scatter_fig(figsize=(12, 12)):
    '''
    Setup Scatter plot figure and axes
    '''

    slc_bg   = '#0E1117'  # Streamlit Background Color
    slc_text = '#FAFAFA'  # Streamlit Text Color
    slc_bg2  = '#262730'  # Streamlit Secondary Background Color

    fig = plt.figure(figsize=figsize, facecolor = slc_bg)
    ax = fig.add_subplot(1, 1, 1, facecolor = slc_bg)

    ax.spines['left'].set_color(slc_text)
    ax.spines['bottom'].set_color(slc_text)
    ax.tick_params(axis='x', colors=slc_text, which='both')
    ax.tick_params(axis='y', colors=slc_text, which='both')

    return fig, ax

def scatter_plot(df, fig, ax, figTitle, xVar, yVar, hueVar, hueOrder, xLim = None, yLim = None, boxoff = False, small_ver = False, feat = None, clusters_label = None, figname='scatter_plot.png', dpi=200, saveFlag=0, palette = 'tab20'):
    """Create a 2D scatter plot and color the points by a specific variable in a dataframe

    Args:
        df (PANDAS dataframe): Dataframe including the coordinates and phenotypes.
        fig (MATPLOTLib Figure Obj): Pregenerated matplotlib figure object
        ax (MATPLOTLib Axes Obj): Pregenerated matplotlib axes object
        figTitle (string): Title to be displayed at the top of the figure
        xVar (string): Name of the dataframe variable representing the X coordinates
        yVar (string): Name of the dataframe variable representing the Y coordinates
        hueVar (string): Name of the dataframe variable representing the distinguishing (color) characteristic
        hueOrder (list of strings): Order of the legend entries
        
        figname (str, optional): Name of the figure to create. Defaults to 'scatter_plot.png'.
        dpi (int, optional): Dots per inch of the saved image. Defaults to 200.
        saveFlag (Boolean, optional): Boolean Flag to save the figure to disk using matplotlib methods. Defaults to FALSE
    """

    figTitle = wrap_title_text(figTitle)
    plot_title = ''
    for i in figTitle:
        plot_title = plot_title + i + '\n'

    slc_bg   = '#0E1117'  # Streamlit Background Color
    slc_text = '#FAFAFA'  # Streamlit Text Color
    slc_bg2 = '#262730'  # Streamlit Secondary Background Color

    # Create the scatter plot
    sns.scatterplot(df,
                    x = xVar,
                    y = yVar,
                    s = 15,
                    hue = hueVar,
                    hue_order = hueOrder,
                    linewidth = 0,
                    palette = palette,
                    ax = ax)

    bbox = ax.get_yticklabels()[-1].get_window_extent()
    x,_ = ax.transAxes.inverted().transform([bbox.x0, bbox.y0])

    # Set the Title
    ax.set_frame_on(False) # Turn off the Frame

    if xVar == 'Cell X Position':
        ax.set_title(plot_title, fontsize = 14, color = slc_text, ha='left', x=x, wrap=True)
        ax.set_xlabel(r'Centroid X ($\mu m)$', fontsize = 14, color = slc_text)
        ax.set_ylabel(r'Centroid Y ($\mu m)$', fontsize = 14, color = slc_text)
        ax.set_aspect(1)       # Set the Aspect Ratio
    else:
        ax.set_xlabel('')
        ax.set_ylabel('')

    if boxoff:
        for sp in ax.spines:
            ax.spines[sp].set_visible(False)
        ax.set(xticks=[], yticks=[])

    if xLim is not None:
        ax.set_xlim(xLim[0], xLim[1])
    else:
        xLim = ax.get_xlim()

    if yLim is not None:
        ax.set_ylim(yLim[0], yLim[1])
    else:
        yLim = ax.get_ylim()

    if small_ver is True:
        lgd_fontsize = 20
        lgd_markscale = 6
    else:
        lgd_fontsize = 10
        lgd_markscale = 6

    # Put the legend outside of the plot
    ax.legend(bbox_to_anchor = (-0.05, -0.1),
              loc = 'upper left',
              fontsize = lgd_fontsize,
              markerscale = lgd_markscale,
              borderaxespad = 0,
              ncols = 4,
              facecolor = slc_bg2,
              edgecolor = slc_bg2,
              labelcolor = slc_text)

    if clusters_label:
        ax.text(0.80*xLim[1], yLim[0], 'Clusters', c = slc_text, fontsize = 25)

    if feat is not None:
        ax.text(xLim[0], 0.93*yLim[1], feat, c = slc_text, fontsize = 30)

    # Save the figure to disk
    if saveFlag:
        fig.savefig(figname, dpi=dpi, bbox_inches='tight')
    return fig

def wrap_title_text(title):
    '''
    Helps with wrapping title text around a 75 character limit
    '''
    char_lim =75
    wrap_title = []
    for x in title:
        while len(x) > char_lim:
            x1 = x[:char_lim]
            x2 = x[char_lim:]
            wrap_title.append(x1)
            x = x2
        wrap_title.append(x)

    return wrap_title

def setup_Spatial_UMAP(df, marker_names, pheno_order, smallest_image_size):
    '''
    Setup the requirements for running spatial UMAP

    Args:
        df (Pandas dataframe): Dataframe containing the data
        marker_names (list): List of marker names
        pheno_order (list): List of phenotype order
        smallest_image_size (int): The size of the smallest image in the dataset
    
    Returns:
        SpatialUMAP: SpatialUMAP object
    '''

    # Initialize the SpatialUMAP object
    # dist_bin_um (np.array): Array of distances in microns
    # um_per_px (float): Microns per pixel
    # area_downsample (float): Area downsample
    spatial_umap = SpatialUMAP(dist_bin_um=np.array([25, 50, 100, 150, 200]), um_per_px=1.0, area_downsample=1.0)
    spatial_umap.cells = df
    spatial_umap.patients = spatial_umap.makeDummyClinic(10)

    # Set Lineage and sort
    spatial_umap.cells['Lineage'] = spatial_umap.cells['phenotype']
    spatial_umap.cells['Lineage'] = spatial_umap.cells['Lineage'].astype("category")
    spatial_umap.cells['Lineage'] = spatial_umap.cells['Lineage'].cat.set_categories(pheno_order)
    # spatial_umap.cells = spatial_umap.cells.sort_values(["Lineage"])

    # Assign pheno_order
    spatial_umap.phenoLabel = pheno_order

    # Set regions
    spatial_umap.cells['TMA_core_id'] = spatial_umap.cells['Slide ID']
    # Set sample number
    if 'Sample_number' not in spatial_umap.cells:
        spatial_umap.cells['Sample_number'] = np.ones(spatial_umap.cells.shape[0])
    print(f'There are {spatial_umap.cells["TMA_core_id"].unique().size} images in this dataset ')

    # Define the number of species we will be working with (how many different get_dummies)
    spatial_umap.species = pheno_order
    spatial_umap.markers = sorted(marker_names)
    spatial_umap.markers = [x + '+' for x in spatial_umap.markers]
    spatial_umap.num_species = len(spatial_umap.species)
    spatial_umap.num_markers = len(spatial_umap.markers)

    # set explicitly as numpy array the cell coordinates (x, y)
    # Notice here that I needed to change the script to CentroidX, CentroidY
    spatial_umap.cell_positions = spatial_umap.cells[['Cell X Position', 'Cell Y Position']].values
    # set explicitly as one hot data frame the cell labels
    spatial_umap.cell_labels = pd.get_dummies(spatial_umap.cells['Lineage'])
    spatial_umap.cell_labels = spatial_umap.cell_labels[spatial_umap.species]
    # set the region is to be analyzed (a TMA core is treated similar to a region of a interest)
    spatial_umap.region_ids = spatial_umap.cells.TMA_core_id.unique()
    # default cluster values
    spatial_umap.cells['clust_label'] = 'No Cluster'

    spatial_umap.elbow_fig = None
    spatial_umap.smallest_image_size = smallest_image_size

    # sets flags for analysis processing
    spatial_umap.phenotyping_completed = True
    spatial_umap.density_completed     = False
    spatial_umap.umap_completed        = False
    spatial_umap.cluster_completed     = False

    return spatial_umap

def perform_density_calc(spatial_umap, bc, calc_areas, cpu_pool_size = 1, area_threshold = 0.001):
    '''
    Calculate the cell counts, cell areas,
    perform the cell densities and cell proportions analyses.

    This is using Andrew's code to calculate the cell counts

    Args:
        spatial_umap (SpatialUMAP): SpatialUMAP object
        bc (benchmark_collector): Benchmark Collector object
        cpu_pool_size (int): Number of CPUs to use for parallel processing
        area_threshold (float): Area threshold to use for cell areas

    Returns:
        SpatialUMAP: SpatialUMAP object with the cell counts, cell areas, 
                    cell densities and cell proportions analyses performed
    '''

    # clear metrics
    spatial_umap.clear_counts()
    spatial_umap.clear_areas()

    print(f'cpu_pool_size set to {cpu_pool_size}\n')

    # get the counts per cell and save to pickle file
    print('Starting Cell Counts process')
    bc.startTimer()
    spatial_umap.get_counts_And(cpu_pool_size=cpu_pool_size)
    bc.printElapsedTime(f'Calculating Counts for {len(spatial_umap.cells)} cells')

    # get the areas of cells and save to pickle file
    print(f'\nStarting Cell Areas process with area threshold of {area_threshold}')
    bc.startTimer()
    spatial_umap.get_areas(calc_areas, area_threshold, pool_size=cpu_pool_size)
    bc.printElapsedTime(f'Calculating Areas for {len(spatial_umap.cells)} cells')

    # calculate density based on counts of cells / area of each arc examine
    spatial_umap.calc_densities(area_threshold)
    # calculate proportions based on species counts/# cells within an arc
    spatial_umap.calc_proportions(area_threshold)

    spatial_umap.density_completed = True

    return spatial_umap

def perform_spatialUMAP(spatial_umap, bc, umap_subset_per_fit, umap_subset_toggle, umap_subset_per):
    '''
    Perform the spatial UMAP analysis

    Args:
        spatial_umap (spatial_umap): spatial_umap object
        bc (benchmark_collector): Benchmark Collector object
        UMAPStyle (str): Style of UMAP to use
    
    Returns:
        spatial_umap: spatial_umap object with the UMAP analysis performed
    '''

    min_image_size = spatial_umap.smallest_image_size
    n_fit = int(min_image_size*umap_subset_per_fit/100)
    n_tra = n_fit + int(min_image_size*umap_subset_per/100)

    # set training and "test" cells for umap training and embedding, respectively
    print('Setting Train/Test Split')
    spatial_umap.set_train_test(n_fit=n_fit, n_tra = n_tra, groupby_label = 'TMA_core_id', seed=54321, umap_subset_toggle = umap_subset_toggle)

    # fit umap on training cells
    print('Fitting Model')
    spatial_umap.umap_fit = umap.UMAP().fit(spatial_umap.density[spatial_umap.cells['umap_train'].values].reshape((spatial_umap.cells['umap_train'].sum(), -1)))
    bc.printElapsedTime(f'      Fitting {np.sum(spatial_umap.cells["umap_train"] == 1)} points to a model', split = True)

    # Transform test cells based on fitted model
    print('Transforming Data')
    spatial_umap.umap_test = spatial_umap.umap_fit.transform(spatial_umap.density[spatial_umap.cells['umap_test'].values].reshape((spatial_umap.cells['umap_test'].sum(), -1)))
    bc.printElapsedTime(f'      Transforming {np.sum(spatial_umap.cells["umap_test"] == 1)} points with the model', split = True)

    spatial_umap.umap_completed = True

    return spatial_umap

def kmeans_calc(clust_data, n_clusters = 5, random_state = None):
    '''
    Perform KMeans clustering on sets of 2D data

    Args:
        clust_data (numpy array): Data to be clustered
        nClus (int): Number of clusters to use
        random_state (int): Random state to use
    
    Returns:
        kmeans_obj: KMeans obj created from KMeans
    '''

    print(f'Starting KMeans Calculation for {n_clusters} clusters')
    # Create KMeans object for a chosen cluster
    kmeans_obj = KMeans(n_clusters = n_clusters,
                        init ='k-means++',
                        max_iter = 300,
                        n_init = 50,
                        random_state = random_state)

    # Fit the data to the KMeans object
    kmeans_obj.fit(clust_data)

    print(f'...Completed KMeans Calculation for {n_clusters} clusters')

    return kmeans_obj

def umap_clustering(spatial_umap, n_clusters, clust_minmax, cpu_pool_size = 8):
    '''
    perform clustering for the UMAP data using KMeans

    Args:
        spatial_umap (spatial_umap): spatial_umap object
        n_clusters (int): Number of clusters to use
        clust_minmax (tuple): Tuple of min and max clusters to use

    Returns:
        spatial_umap: spatial_umap object with the clustering performed
    '''
    # Reset the cluster labels just in case
    spatial_umap.df_umap.loc[:, 'clust_label'] = -1

    clust_range = range(clust_minmax[0], clust_minmax[1]+1)

    kwargs_list = []
    for clust in clust_range:
        kwargs_list.append(
            (
                spatial_umap.umap_test,
                clust
            )
        )

    results = utils.execute_data_parallelism_potentially(kmeans_calc,
                                                         kwargs_list,
                                                         nworkers = cpu_pool_size,
                                                         task_description='KMeans Clustering',
                                                         use_starmap=True)
    # mp_start_method = mp.get_start_method()
    # # Create a pool of worker processes
    # with mp.get_context(mp_start_method).Pool(processes=cpu_pool_size) as pool:
    #     results = pool.starmap(kmeans_calc, kwargs_list)

    wcss = [x.inertia_ for x in results]

    # Create WCSS Elbow Plot
    spatial_umap.elbow_fig = draw_wcss_elbow_plot(clust_range, wcss, n_clusters)

    # Identify the kmeans obj that matches the selected cluster number
    kmeans_obj_targ = results[n_clusters-1]

    spatial_umap.cluster_dict = dict()
    for i in range(n_clusters):
        spatial_umap.cluster_dict[i+1] = f'Cluster {i+1}'
    spatial_umap.cluster_dict[0] = 'No Cluster'

    spatial_umap.palette_dict = dict()
    for i in range(n_clusters):
        spatial_umap.palette_dict[f'Cluster {i+1}'] = sns.color_palette('tab20')[i]
    spatial_umap.palette_dict['No Cluster'] = 'white'

    # Assign values to cluster_label column in df_umap
    spatial_umap.df_umap.loc[:, 'clust_label'] = [spatial_umap.cluster_dict[key] for key in (kmeans_obj_targ.labels_+1)]

    return spatial_umap

def draw_wcss_elbow_plot(clust_range, wcss, sel_clus):
    '''
    Calculate possible clusters and plot the elbow plot

    Args:
        clust_range (list): List of cluster values
        wcss (list): List of within-cluster sum of squares
        sel_clus (int): Selected cluster value

    Returns:
        fig: Matplotlib figure object
    '''

    # Streamlit Theming
    slc_bg   = '#0E1117'  # Streamlit Background Color
    slc_text = '#FAFAFA'  # Streamlit Text Color
    slc_bg2  = '#262730'  # Streamlit Secondary Background Color

    fig = plt.figure(figsize = (5,5), facecolor = slc_bg)
    ax = fig.add_subplot(1,1,1, facecolor = slc_bg)
    ax.set_xlabel('Number of Clusters', fontsize = 10, color = slc_text)
    ax.set_ylabel('WCSS', fontsize = 10, color = slc_text)
    ax.set_xlim(0, clust_range[-1])
    ax.set_xticks(np.linspace(0, clust_range[-1], clust_range[-1]+1))

    plt.plot(clust_range, wcss)
    # plt.axvline(sel_clus, linestyle='--', color='r')

    ax.spines['left'].set_color(slc_text)
    ax.spines['bottom'].set_color(slc_text)
    ax.spines['top'].set_color(slc_bg)
    ax.spines['right'].set_color(slc_bg)
    ax.tick_params(axis='x', colors=slc_text, which='both')
    ax.tick_params(axis='y', colors=slc_text, which='both')
    return fig

def draw_heatmap_fig(df, pheno_list, title, norm_axis = None):
    '''
    Create a heatmap of the phenotypes and clusters

    Args:
        df:
        pheno_list:
        title:
        norm_axis:

    Returns:
        fig: Matplotlib figure
    '''
    # Create heatmap df
    heatmap_df = pd.DataFrame()

    for clust_label, group in df.groupby('clust_label'):
        clust_value_counts = group['Lineage'].value_counts()
        clust_value_counts.name = f'{clust_label}'

        heatmap_df = pd.concat([heatmap_df, pd.DataFrame([clust_value_counts])])

    # Fix the NA
    heatmap_df[heatmap_df.isna()] = 0
    heatmap_df = heatmap_df.astype('int')

    # Rearrange Columsn in order of prevalences
    heatmap_df = heatmap_df.loc[:, pheno_list]

    if norm_axis == 0:
        heatmap_title = 'Phenotype/Cluster Heatmap: Normalized within Clusters'
        heatmap_df = round(heatmap_df.div(heatmap_df.sum(axis=1), axis=0), 3)
        vmin = 0
        vmax = 1
    elif norm_axis == 1:
        heatmap_title = 'Phenotype/Cluster Heatmap: Normalized within Phenotypes'
        heatmap_df = round(heatmap_df.div(heatmap_df.sum(axis=0), axis=1), 3)
        vmin = 0
        vmax = 1
    else:
        heatmap_title = 'Phenotype/Cluster Heatmap: '
        vmin = heatmap_df.min().min()
        vmax = heatmap_df.max().max()
    title.append(heatmap_title)

    fig_title = wrap_title_text(title)
    plot_title = ''
    for i in fig_title:
        plot_title = plot_title + i + '\n'

    # Define Output Variables
    phenotypes = heatmap_df.columns
    clusters = heatmap_df.index

    # Theme Styles
    slc_bg   = '#0E1117'  # Streamlit Background Color
    slc_text = '#FAFAFA'  # Streamlit Text Color
    slc_bg2  = '#262730'  # Streamlit Secondary Background Color

    fig = plt.figure(figsize = (12,12), facecolor = slc_bg)
    ax = fig.add_subplot(1,1,1, facecolor = slc_bg)
    im = ax.imshow(heatmap_df, cmap = 'inferno', vmin = vmin, vmax = vmax)

    bbox = ax.get_yticklabels()[-1].get_window_extent()
    x, _ = ax.transAxes.inverted().transform([bbox.x0, bbox.y0])

    # Show all ticks and label them with the respective list entries
    ax.set_title(plot_title, fontsize = 20, loc = 'left', color = slc_text, x=3*x, wrap=True)
    ax.set_xticks(np.arange(len(phenotypes)), labels = phenotypes, fontsize = 14, color = slc_text)
    ax.set_yticks(np.arange(len(clusters)), labels = clusters, fontsize = 14, color = slc_text)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
            rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    for i, cluster in enumerate(range(len(clusters))):
        for j, phenotype in enumerate(range(len(phenotypes))):
            value = heatmap_df.iloc[cluster, phenotype]
            if value >= 0.85*vmax:
                text_color = 'b'
            else:
                text_color = 'w'
            ax.text(j, i, value, ha="center", va="center", color=text_color)

    fig.tight_layout()

    return fig

def draw_neigh_profile_fig(spatial_umap, ax, sel_clus, cmp_clus = None, cmp_style = None, hide_other = False, hide_no_cluster = False, legend_flag = True):
    '''
    draw_neigh_profile_fig is the method that draws the neighborhood profile
    line plots

    Args:
        spatial_umap: The spatial UMAP object containing the data to plot
        ax: The matplotlib axis to draw the plot on
        sel_clus: The selected cluster to highlight
        cmp_clus: The cluster to compare against (optional)
        cmp_style: The comparison style to use (optional)
        hide_other: Whether to hide other phenotypes (optional)
        hide_no_cluster: Whether to hide cells with no cluster (optional)
        legend_flag: Whether to show the legend (optional)

    Returns:
        None
    '''

    dens_df_mean_base = spatial_umap.dens_df_mean
    if hide_other:
        dens_df_mean_base = dens_df_mean_base.loc[dens_df_mean_base['phenotype'] != 'Other', :]
    if hide_no_cluster:
        dens_df_mean_base = dens_df_mean_base.loc[dens_df_mean_base['clust_label'] != 'No Cluster', :]

    maxdens_df   = 1.05*max(dens_df_mean_base['density_mean'] + dens_df_mean_base['density_sem'])
    dens_df_mean_sel = dens_df_mean_base.loc[dens_df_mean_base['clust_label'] == sel_clus, :].reset_index(drop=True)
    ylim = [1, maxdens_df]
    dens_df_mean = dens_df_mean_sel.copy()
    cluster_title = f'{sel_clus}'

    if cmp_clus is not None:
        dens_df_mean_cmp = dens_df_mean_base.loc[dens_df_mean_base['clust_label'] == cmp_clus, :].reset_index(drop=True)

        dens_df_mean = dens_df_mean_cmp.copy()
        dens_df_mean['density_sem'] = 0
        # Subtraction Compare Style
        if cmp_style == 'Difference':
            dens_df_mean['density_mean'] = dens_df_mean_sel['density_mean'] - dens_df_mean_cmp['density_mean']

            range_values = [min(dens_df_mean['density_mean']), max(dens_df_mean['density_mean'])]
            top_range = 1.05*max(abs(range_values[0]), abs(range_values[1]))
            ylim = [-top_range, top_range]
            cluster_title = f'{sel_clus} - {cmp_clus}'
        # Ratio Compare Style
        elif cmp_style == 'Ratio':
            dens_df_mean['density_mean'] = dens_df_mean_sel['density_mean'] / dens_df_mean_cmp['density_mean']

            range_values = [dens_df_mean.loc[np.isfinite(dens_df_mean['density_mean']), 'density_mean'].min(),
                            dens_df_mean.loc[np.isfinite(dens_df_mean['density_mean']), 'density_mean'].max()]
            if range_values[0] < 0:
                ymin = 1.05*range_values[0]
                ymax = 1.05*range_values[1]
            else:
                ymin = 0.95*range_values[0]
                ymax = 1.05*range_values[1]
            ylim = np.array([ymin, ymax])
            cluster_title = f'{sel_clus} / {cmp_clus}'
    else:
        cmp_style = None

    if not np.all([math.isfinite(x) for x in ylim]):
        ylim = [1, 10]

    umPT.plot_mean_neighborhood_profile(ax = ax,
                                        dist_bin = spatial_umap.dist_bin_um,
                                        pheno_order = spatial_umap.phenoLabel,
                                        npf_dens_mean = dens_df_mean,
                                        cluster_title = cluster_title,
                                        cmp_style = cmp_style,
                                        max_dens = ylim,
                                        leg_flag = legend_flag)

def preprocess_weighted_umap(w, df_umap):
    '''
    Perform perprocessing on UMAP data and weights

    w will be the values from a specific feature, 
    not any and all features

    Weights are essentially any chosen feature of the data beyond
    the x/y coordinates, lineage, and cluster number

    Args:
        w (numpy array): Weights for the UMAP data
        df_umap (Pandas dataframe): UMAP data

    Returns:
        w (numpy array): Preprocessed weights
        df_umap (Pandas dataframe): Preprocessed UMAP data2003
    '''

    # Check for NaN in the w and remove them
    not_nan = ~np.isnan(w)
    w = w[not_nan]                    # Remove NaNs from w
    df_umap = df_umap.loc[not_nan, :] # Remove NaNs from df_umap

    # Raise all values of w about 0
    if np.any(w < 0):
        w = w + min(w)*-1.2

    # Apply the log to the weights
    w = np.log(0.1 * w + 0.1)
    w -= np.min(w)

    return w, df_umap

def UMAPdraw_density(d, bins, w, n_pad, vlim, feat = None, diff = False, legendtype = 'colorbar', figsize=(12, 12)):
    '''
    UMAPdraw_density is the method that draws the UMAP density plot

    Args:
        d (numpy array): UMAP data
        bins (int): Number of bins to use
        w (numpy array): Weights for the UMAP data
        n_pad (int): Padding for the bins
        vlim (list): Limits for the colorbar
        feat (str): Feature to display
        diff (bool): Flag to display the difference
        legendtype (str): Type of legend to display
        figsize (tuple): Size of the figure

    Returns:
        umap_fig (MATPLOTLIB Figure Obj): UMAP density plot figure
    '''

    # Streamlit Theming
    slc_bg   = '#0E1117'  # Streamlit Background Color
    slc_text = '#FAFAFA'  # Streamlit Text Color
    slc_bg2  = '#262730'  # Streamlit Secondary Background Color

    # color maps
    cmap_viridis = plt.get_cmap('viridis').copy()
    cmap_viridis.set_under('white')
    cmap_magma = plt.get_cmap('magma').copy()
    cmap_magma.set_under('white')
    cmap_bwr = plt.get_cmap('bwr').copy()

    # Set up Figure
    umap_fig = plt.figure(figsize=figsize, facecolor = slc_bg)
    ax = umap_fig.add_subplot(1, 1, 1, facecolor = slc_bg)

    if w is None and diff is False:
        cmap = cmap_viridis
        circle_type = None
    elif diff is False:
        cmap = cmap_magma
        circle_type = None
    else:
        cmap = cmap_bwr
        circle_type = 'arch'

    umPT.plot_2d_density(d, bins=bins, w=w, n_pad=n_pad, ax=ax, cmap=cmap,
                         vlim = vlim, circle_type = circle_type, legendtype = legendtype)

    x_lim = ax.get_xlim()
    y_lim = ax.get_ylim()

    ax.text(0.82*x_lim[1], 0.03*y_lim[1], 'Density', c = slc_text, fontsize = 25)

    if feat is not None:
        if cmap == cmap_bwr:
            ax.text(x_lim[0], 0.90*y_lim[1], feat, c = 'black', fontsize = 30)
        else:
            ax.text(x_lim[0], 0.90*y_lim[1], feat, c = slc_text, fontsize = 30)

    return umap_fig

def draw_incidence_fig(inci_df, fig_title, phenotype = 'All Phenotypes', feature = 'Cell Counts', displayas = 'Counts Difference', msg_tags = ['=0', '=1'], show_raw_counts= False):
    '''
    Draws the line plot figure which describes the incideces of a 
    selected features

    Args:
        inci_df (pd.DataFrame): DataFrame containing the incidence data
        fig_title (str): Title of the figure
        phenotype (str): Phenotype to display
        feature (str): Feature to display
        displayas (str): How to display the data (Counts Difference, Ratios, Percentages)
        comp_thresh (float): Comparison threshold
        show_raw_counts (bool): Whether to show raw counts

    Returns:
        inci_fig (MATPLOTLIB Figure Obj): Incidence figure
    '''

    inci_fig = go.Figure()

    slc_bg   = '#0E1117'  # Streamlit Background Color
    slc_text = '#FAFAFA'  # Streamlit Text Color
    slc_bg2  = '#262730'  # Streamlit Secondary Background Color
    slc_ylw  = '#F6EB61'  # Streamlit Yellow Color
    slc_red  = '#FF4B4B'  # Streamlit Red Color

    up_tag = msg_tags[1]
    dn_tag = msg_tags[0]

    inci_df = inci_df.round(2)

    anno2 = False
    if feature != 'Cell Counts':

        df = inci_df[displayas]

        dfmin = df.loc[(df != np.nan)].min()
        dfmax = df.loc[(df != np.nan)].max()
        up_limit = max(-1*dfmin, dfmax)
        if displayas == 'Count Differences':
            anno2 = True

            df_up = inci_df['featureCount1']
            df_dn = inci_df['featureCount0']

            dfmin = df_dn.loc[(df_dn != np.nan)].max()
            dfmax = df_up.loc[(df_up != np.nan)].max()
            up_limit = max(dfmin, dfmax)

            if up_limit < 2:
                up_limit = 2
            ylim = [-1.05*up_limit, 1.05*up_limit]
            feature_pos = [1, up_limit*.95]
            feature_text = f'{feature}{up_tag}'
            feature_pos2 = [1, -up_limit*.95]
            feature_text2 =f'{feature}{dn_tag}'
            hover_template = '<b>Cluster:</b> %{x}<br><b>Count Difference:</b> %{y}<extra></extra>'
            outcome_suff = ' (Counts)'
        elif displayas == 'Ratios':
            anno2 = True

            df_up = inci_df['Percentages1_adj_log']
            df_dn = inci_df['Percentages0_adj_log']

            dfmin = df_dn.loc[(df_dn != np.nan)].max()
            dfmax = df_up.loc[(df_up != np.nan)].max()
            up_limit = max(dfmin, dfmax)

            ylim = [-1.05*up_limit, 1.05*up_limit]
            feature_pos = [1, up_limit*.95]
            feature_text = f'{feature}{up_tag}'
            feature_pos2 = [1, -up_limit*.95]
            feature_text2 =f'{feature}{dn_tag}'
            hover_template = '<b>Cluster:</b> %{x}<br><b>Ratio:</b> %{y}<extra></extra>'
            outcome_suff = ' Ratio (log10)'
        elif displayas == 'Percentages':

            df_up = inci_df['Percentages']
            df_dn = inci_df['Percentages0']

            dfmin = df_dn.loc[(df_dn != np.nan)].max()
            dfmax = df_up.loc[(df_up != np.nan)].max()
            up_limit = dfmax

            ylim = [-1.05, 1.05*up_limit]
            feature_pos = [1, up_limit*.95]
            feature_text = f'{feature}{up_tag}'
            hover_template = '<b>Cluster:</b> %{x}<br><b>Percentage:</b> %{y}<extra></extra>'
            outcome_suff = ' (%)'
    else:
        df = inci_df['counts']

        dfmin = df.min()
        dfmax = df.max()
        up_limit = max(-1*dfmin, dfmax)
        ylim = [0, dfmax*1.05]

        feature_pos = [1, up_limit*.95]
        feature_text = f'{feature}'
        hover_template = '<b>Cluster:</b> %{x}<br><b>Count:</b> %{y}<extra></extra>'
        outcome_suff = ' (Counts)'

    if feature != 'Cell Counts':
        if show_raw_counts:
            inci_fig.add_trace(go.Bar(
                x=df_up.index,
                y=df_up.values,
                name=f"{phenotype}{up_tag}",
                marker=dict(color=slc_ylw),
                hovertemplate=hover_template,
                hoverlabel=dict(
                bgcolor=slc_ylw,
                bordercolor=slc_ylw,
                font=dict(color=slc_bg)
                ),
                opacity=0.65,
                offsetgroup='1',
                showlegend=False,
                text=[f"<b>{y:,}</b>" for y in df_up.values],
                textposition='inside',
                textfont=dict(color=slc_bg, size=14)
            ))
            inci_fig.add_trace(go.Bar(
                x=df_dn.index,
                y=-df_dn.values,
                name=f"{phenotype}{dn_tag}",
                marker=dict(color=slc_red),
                hovertemplate=hover_template,
                hoverlabel=dict(
                bgcolor=slc_red,
                bordercolor=slc_red,
                font=dict(color=slc_text)
                ),
                opacity=0.65,
                offsetgroup='1',
                showlegend=False,
                text=[f"<b>{y:,}</b>" for y in df_dn.values],
                textposition='inside',
                textfont=dict(color=slc_bg, size=14)
            ))

        inci_fig.add_trace(go.Scatter(
            x=df.index,
            y=df.values,
            mode='lines+markers',
            name=phenotype,
            line=dict(color='#0E86D4', width=2.5),
            marker=dict(color='#0E86D4', size=14),
            hovertemplate=hover_template,
            hoverlabel=dict(
                bgcolor='#0E86D4',
                bordercolor='#0E86D4',
                font=dict(color=slc_text)
            )
        ))
    else:
        inci_fig.add_trace(go.Bar(
            x=df.index,
            y=df.values,
            name=phenotype,
            marker=dict(color='#0E86D4'),
            hovertemplate=hover_template,
            hoverlabel=dict(
                bgcolor='#0E86D4',
                bordercolor='#0E86D4',
                font=dict(color=slc_text)
            ),
            text=[f"<b>{int(y):,}</b>" for y in df.values],
            textposition='inside',
            textfont=dict(color=slc_bg, size=14)
        ))

    annotations = [
        dict(
            x=feature_pos[0],
            y=feature_pos[1],
            text=feature_text,
            showarrow=False,
            font=dict(size=40, color=slc_text),
            xanchor='center',
            yanchor='bottom',
            opacity=0.3,
        )
    ]

    # Add a second annotation if anno2 is True
    if anno2:
        annotations.append(
            dict(
                x=feature_pos2[0],
                y=feature_pos2[1],
                text=feature_text2,
                showarrow=False,
                font=dict(size=40, color=slc_text),
                xanchor='center',
                yanchor='bottom',
                opacity=0.3,
            )
        )

    inci_fig.update_layout(
        title=dict(
            text=fig_title,
            font=dict(size=25),
            x=0.08,  # Align title to the left (vertical axis)
            xanchor='left'
        ),
        xaxis_title="Cluster #",
        yaxis_title=f'{feature}{outcome_suff}',
        plot_bgcolor=slc_bg,
        paper_bgcolor=slc_bg,
        font=dict(color=slc_text, size=14),
        xaxis=dict(showgrid=True, gridcolor=slc_bg2, showline=False, zeroline=False),
        yaxis=dict(
            showgrid=True,
            gridcolor=slc_bg2,
            showline=False,
            zeroline=False,
            range=ylim
        ),
        width=2000,
        height=1000,
        legend=dict(
            title=None,
            bgcolor=slc_bg2,
            bordercolor=slc_bg2,
            borderwidth=1,
            orientation='v',
            x=0.85,
            y=1,
            xanchor='left',
            yanchor='top'
        ),
        showlegend=True,
        annotations=annotations,
        shapes=[
            dict(
                type='line',
                xref='paper',
                x0=0,
                x1=1,
                yref='y',
                y0=0,
                y1=0,
                line=dict(
                    color=slc_text,
                    width=2,
                    dash='dash'
                ),
                opacity=0.7,
                layer='below'
            )
        ]
    )

    return inci_fig
