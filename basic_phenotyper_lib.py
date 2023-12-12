import time
import numpy as np
import pandas as pd
import warnings
warnings.simplefilter(action='ignore', category= FutureWarning)
warnings.filterwarnings("ignore", message=".*The 'nopython' keyword.*")
pd.options.mode.chained_assignment = None  # default='warn'
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans # K-Means

from benchmark_collector import benchmark_collector # Benchmark Collector Class

def preprocess_df(df_orig, marker_names, marker_col_prefix):
    '''Perform some preprocessing on our dataset to apply tranforms
    and collect meta-data
    '''

    # Create the bench mark collector obj
    bc = benchmark_collector()

    # Set df_raw as the baseline dataframe
    df_raw = df_orig.copy()
    
    # Step 1: Initalize the columns needed for phenotyping
    bc.startTimer()
    df_raw = init_pheno_cols(df_raw, marker_names, marker_col_prefix)
    bc.printElapsedTime(msg = '     Initializing Phenotyping Columns')

    # Step 2: Intialize Phenotying Assignments dataframe
    bc.startTimer()
    pheno_assign = init_pheno_assign(df_raw)
    bc.printElapsedTime(msg = '     Initializing Phenotyping Assignments Table')

    # Step 3: Intialize Phenotype Summary Dataframe (based on Phenotying Assigments)
    bc.startTimer()
    pheno_summ = init_pheno_summ(df_raw)
    bc.printElapsedTime(msg = '     Initializing Phenotype Summary Table')

    # Make a copy of df_raw as df
    df = df_raw.copy()

    return df_raw, df, pheno_assign, pheno_summ

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
    df['mark_bits'] = df_markers.apply(lambda row: ''.join((str(1) if str(x)[-1] == '+' else str(0) for x in row)), axis='columns')

    # Add a column of prettier names for the species, e.g., 'VIM- ECAD+ COX2+ NOS2-'
    df['species_name_long'] = df['mark_bits'].apply(lambda mark_bits: ' '.join([marker_name + ('+' if marker_bit == '1' else '-') for marker_name, marker_bit in zip(marker_names, mark_bits)]))

    # Add a column dropping the negative markers from these pretty names, e.g., 'ECAD+ COX2+'
    df['species_name_short'] = df['species_name_long'].apply(lambda species_name_long: '+ '.join([marker_names[iy] for iy, y in enumerate([x for x in species_name_long if x in ('+', '-')]) if y == '+']) + '+').apply(lambda x: (x if len(x) != 1 else 'Other'))

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
    """For each unique species (elsewhere called "exclusive" phenotyping), generate information concerning their prevalence in a new dataframe.

    Args:
        df (Pandas dataframe): Dataframe containing data from the input dataset, including a "mark_bits" column

    Returns:
        Pandas dataframe: Dataframe containing the value counts of each "exclusive" species
    """

    st_init_species = time.time()
    spec_summ = df[['species_name_short', 'phenotype', 'species_name_long']]
    sp_init_species = time.time()
    elapsed = round(sp_init_species - st_init_species, 3)
    print(f'        Initalizing Phenotying Assignments: {elapsed}s')

    spec_summ['species_count'] = spec_summ['species_name_short'].groupby(spec_summ['species_name_short']).transform('count')
    spec_summ = spec_summ.drop_duplicates().reset_index(drop=True)
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
    """For each unique species (elsewhere called "exclusive" phenotyping), generate information concerning their prevalence in a new dataframe.

    Args:
        df (Pandas dataframe): Dataframe containing data from the input dataset, including a "mark_bits" column

    Returns:
        assign_pheno (Pandas dataframe): Dataframe containing the value counts of each "exclusive" species
    """

    assign_pheno = df[['phenotype', 'species_name_short', 'species_name_long']].groupby(by='phenotype', as_index = False).agg(lambda x: np.unique(list(x)))

    assign_pheno['phenotype_count'] = [sum(df['phenotype'] == x) for x in assign_pheno.phenotype]
    assign_pheno['phenotype_percent'] = [round(100*x/sum(assign_pheno['phenotype_count']), 2) for x in assign_pheno['phenotype_count']]
    assign_pheno = assign_pheno.sort_values(by='phenotype_percent', ascending=False)

    return assign_pheno

def remove_compound_species(df, marker_names, allow_compound_species=True):
    '''
    For each compound species ('Species int' not just a plain power of two), add each 
    individual phenotype to the end of the dataframe individually and then delete the original compound entry
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
        tsv_file (str): Name of the TSV file generated by pasting into a text editor columns from Excel 
                        in the order specified in the read_csv() method below

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

            # Store a copy of the object data for the current species ID... copy is indeed needed or else I'd get a typical Pandas warning!
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

                # Save the current object data with updated species ID to the dataframe holder that we will later combine
                # into the final dataset... note the .copy() is absolutely needed or else when there are multiple phenotypes 
                # per set of surface markers, copies of the dataframes are added to list_objects_new instead of ones with different species IDs!!
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
    """Add a "phenotype" column to df based on the from the marker columns via the calculated species integer.

    Args:
        df (Pandas dataframe): Dataframe containing the "species_name_short" column

    Returns:
        df (Pandas dataframe): Input dataframe with an added "phenotype" column appended or overwritten
    """

    # Create a "phenotype" column based on the from the "species_name_short"
    df['phenotype'] = df['species_name_short']

    # Return dataframe with phenotype column
    return df

def assign_phenotype_custom(df, spec_summ):
    """Add a "phenotype" column to df based on a species summary dataframe which identfies custom phenotype assignments from custom phenotyping

    Args:
        df (Pandas dataframe): Dataframe containing the "species_name_short" column
        spec_summ (Pandas dataframe): Dataframe containing the custom assignments of phenotypes based on "species_name_short"

    Returns:
        df (Pandas dataframe): Input dataframe with an added "phenotype" column appended or overwritten
    """
    # Create a "phenotype" column mapped from a species summary dataset
    df['phenotype'] = df['species_name_short'].map(dict(zip(spec_summ['species_name_short'].to_list(), spec_summ['phenotype'].to_list())))
    
    df.loc[df['phenotype'].isna(), 'phenotype'] = 'unassigned'

    # Return dataframe with phenotype column
    return df

def load_previous_species_summary(filename):

    spec_summ_load = pd.read_csv(filename)

    spec_summ = spec_summ_load

    return spec_summ

def draw_scatter_fig(figsize=(12, 12)):

    SlBgC  = '#0E1117'  # Streamlit Background Color
    SlTC   = '#FAFAFA'  # Streamlit Text Color
    Sl2BgC = '#262730'  # Streamlit Secondary Background Color

    fig = plt.figure(figsize=figsize, facecolor = SlBgC)
    ax = fig.add_subplot(1, 1, 1)

    ax.spines['left'].set_color(SlTC)
    ax.spines['bottom'].set_color(SlTC)
    ax.tick_params(axis='x', colors=SlTC, which='both')
    ax.tick_params(axis='y', colors=SlTC, which='both')

    return fig, ax
    
def scatter_plot(df, fig, ax, figTitle, xVar, yVar, hueVar, hueOrder, xLim = None, yLim = None, boxoff = False, small_ver = False, feat = None, clusters_label = None, figname='scatter_plot.png', dpi=200, saveFlag=0):
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

    figTitle = wrapTitleText(figTitle)
    pltTitle = ''
    for i in figTitle:
        pltTitle = pltTitle + i + '\n'

    SlBgC  = '#0E1117'  # Streamlit Background Color
    SlTC   = '#FAFAFA'  # Streamlit Text Color
    Sl2BgC = '#262730'  # Streamlit Secondary Background Color

    # Create the scatter plot
    sns.scatterplot(df,
                    x = xVar,
                    y = yVar,
                    s = 4,
                    hue = hueVar,
                    hue_order = hueOrder,
                    palette = 'tab20',
                    ax = ax)
    
    bbox = ax.get_yticklabels()[-1].get_window_extent()
    x,_ = ax.transAxes.inverted().transform([bbox.x0, bbox.y0])

    # Set the Title
    ax.set_frame_on(False) # Turn off the Frame

    if xVar == 'Cell X Position':
        ax.set_title(pltTitle, fontsize = 14, color = SlTC, ha='left', x=x, wrap=True)
        ax.set_xlabel('Centroid X ('r'$\mu m)$', fontsize = 14, color = SlTC)
        ax.set_ylabel('Centroid Y ('r'$\mu m)$', fontsize = 14, color = SlTC)
        ax.set_aspect(1)       # Set the Aspect Ratio
    else:
        ax.set_xlabel('')
        ax.set_ylabel('')

    if boxoff:
        [ax.spines[sp].set_visible(False) for sp in ax.spines]
        ax.set(xticks=[], yticks=[])

    if xLim is not None:
        ax.set_xlim(xLim[0], xLim[1])
    else:
        xLim = ax.get_xlim()
    
    if yLim is not None:
        ax.set_ylim(yLim[0], yLim[1])
    else:
        yLim = ax.get_ylim()

    if small_ver == True:
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
              facecolor = Sl2BgC,
              edgecolor = Sl2BgC,
              labelcolor = SlTC)

    if clusters_label:
        ax.text(0.80*xLim[1], yLim[0], 'Clusters', c = SlTC, fontsize = 25)

    if feat is not None:
        ax.text(xLim[0], 0.93*yLim[1], feat, c = SlTC, fontsize = 30)


    # Save the figure to disk
    if saveFlag:
        fig.savefig(figname, dpi=dpi, bbox_inches='tight')
    return fig

# Helps with Wrapping text
def wrapTitleText(title):
    charLim =70
    wrapTitle = []
    for x in title:
        while len(x) > charLim:
            x1 = x[:charLim]
            x2 = x[charLim:]
            wrapTitle.append(x1)
            x = x2
        wrapTitle.append(x)

    return wrapTitle

def setup_Spatial_UMAP(df, marker_names, phenoOrder, cpu_pool_size = 1):
    # Import Libraries REMOVE THIS APPEND LATER
    import pandas as pd
    import numpy as np
    import matplotlib as mpl
    from SpatialUMAP import SpatialUMAP
        
    spatial_umap = SpatialUMAP(dist_bin_um=np.array([25, 50, 100, 150, 200]), um_per_px=0.5, area_downsample=.2)
    spatial_umap.cells = df
    spatial_umap.patients = spatial_umap.makeDummyClinic(10)

    # Set Lineage and sort
    spatial_umap.cells['Lineage'] = spatial_umap.cells['phenotype']
    spatial_umap.cells['Lineage'] = spatial_umap.cells['Lineage'].astype("category")
    spatial_umap.cells['Lineage'] = spatial_umap.cells['Lineage'].cat.set_categories(phenoOrder)
    # spatial_umap.cells = spatial_umap.cells.sort_values(["Lineage"])

    # Set regions
    spatial_umap.cells['TMA_core_id'] = spatial_umap.cells['Slide ID']
    # Set sample number
    if 'Sample_number' not in spatial_umap.cells:
        spatial_umap.cells['Sample_number'] = np.ones(spatial_umap.cells.shape[0])
    print(f'There are {spatial_umap.cells["Sample_number"].unique().size} samples in this dataset ')

    # Define the number of species we will be working with (how many different get_dummies)
    spatial_umap.species = sorted(spatial_umap.cells['Lineage'].unique())
    spatial_umap.markers = sorted(marker_names)
    spatial_umap.markers = [x + '+' for x in spatial_umap.markers]
    spatial_umap.num_species = len(spatial_umap.species)
    spatial_umap.num_markers = len(spatial_umap.markers)

    # set explicitly as numpy array the cell coordinates (x, y)
    # Notice here that I needed to change the script to CentroidX, CentroidY
    spatial_umap.cell_positions = spatial_umap.cells[['Cell X Position', 'Cell Y Position']].values
    # set explicitly as one hot data frame the cell labels
    spatial_umap.cell_labels = pd.get_dummies(spatial_umap.cells['Lineage'])
    # set the region is to be analyzed (a TMA core is treated similar to a region of a interest)
    spatial_umap.region_ids = spatial_umap.cells.TMA_core_id.unique()
    # clear metrics
    spatial_umap.clear_counts()
    spatial_umap.clear_areas()

    print(f'cpu_pool_size set to {cpu_pool_size}\n')

    # get the counts per cell and save to pickle file
    print('Starting Cell Counts process')
    spatial_umap.get_counts(pool_size=cpu_pool_size)

    # get the areas of cells and save to pickle file
    area_threshold = 0.001
    print('Starting Cell Areas process')
    spatial_umap.get_areas(area_threshold, pool_size=cpu_pool_size)

    # calculate density based on counts of cells / area of each arc examine
    spatial_umap.calc_densities(area_threshold)
    # calculate proportions based on species counts/# cells within an arc
    spatial_umap.calc_proportions(area_threshold)

    spatial_umap.phenoLabel = phenoOrder

    spatial_umap.cells['clust_label'] = -1

    return spatial_umap

def perform_spatialUMAP(spatial_umap, UMAPStyle):
    import umap
    from benchmark_collector import benchmark_collector # Benchmark Collector Class

    bc = benchmark_collector()
    # set training and "test" cells for umap training and embedding, respectively
    print('Setting Train/Test Split')
    spatial_umap.set_train_test(n=2500, groupByLabel = 'Sample_number', seed=54321)
    
    # fit umap on training cells
    bc.startTimer()
    print('Fitting Model')
    spatial_umap.umap_fit = umap.UMAP().fit(spatial_umap.density[spatial_umap.cells['umap_train'].values].reshape((spatial_umap.cells['umap_train'].sum(), -1)))
    bc.printElapsedTime(f'      Fitting {np.sum(spatial_umap.cells["umap_train"] == 1)} points to a model')
    
    # apply umap embedding on test cells
    bc.startTimer()
    print('Transforming Data')
    spatial_umap.umap_test = spatial_umap.umap_fit.transform(spatial_umap.density[spatial_umap.cells['umap_test'].values].reshape((spatial_umap.cells['umap_test'].sum(), -1)))
    bc.printElapsedTime(f'      Transforming {np.sum(spatial_umap.cells["umap_test"] == 1)} points with the model')
    return spatial_umap

def KMeans_calc(spatial_umap, nClus):
    # Create KMeans object for a chosen cluster
    kmeansObj = KMeans(n_clusters = nClus, init ='k-means++', max_iter=300, n_init=10, random_state=42)
    # Fit the data to the KMeans object
    kmeansObj.fit(spatial_umap.umap_test)
    # Extract out the labels of trained data from the fit model
    return kmeansObj

def measure_possible_clust(spatial_umap, clust_minmax):
    clust_range = range(clust_minmax[0], clust_minmax[1])
    wcss = [] # Within-Cluster Sum of Squares
    for nClus in clust_range: 
        # Perform clustering for chosen 
        kmeansObj = KMeans_calc(spatial_umap, nClus)
        # Append Within-Cluster Sum of Squares measurement
        wcss.append(kmeansObj.inertia_)
    return list(clust_range), wcss

def perform_clusteringUMAP(spatial_umap, nClus):

    # Reperform clustering for chosen 
    kmeansObj = KMeans_calc(spatial_umap, nClus)
    # Add cluster label column to cells dataframe
    spatial_umap.cells.loc[spatial_umap.cells.loc[:, 'umap_test'] == True, 'clust_label'] = kmeansObj.labels_
    # With the cluster labels assigned, perform mean calculations
    spatial_umap.mean_measures()

    return spatial_umap

def draw_wcss_elbow_plot(clust_range, wcss, selClus):

    SlBgC  = '#0E1117'  # Streamlit Background Color
    SlTC   = '#FAFAFA'  # Streamlit Text Color
    Sl2BgC = '#262730'  # Streamlit Secondary Background Color

    fig = plt.figure(figsize = (5,5), facecolor = SlBgC)
    ax = fig.add_subplot(1,1,1, facecolor = SlBgC) 
    ax.set_xlabel('Number of Clusters', fontsize = 10, color = SlTC)
    ax.set_ylabel('WCSS', fontsize = 10, color = SlTC)
    # ax.set_xticks(np.linspace(0, 21, 22))
    plt.plot(clust_range, wcss)
    plt.axvline(selClus, linestyle='--', color='r')

    ax.spines['left'].set_color(SlTC)
    ax.spines['bottom'].set_color(SlTC)
    ax.spines['top'].set_color(SlBgC)
    ax.spines['right'].set_color(SlBgC)
    ax.tick_params(axis='x', colors=SlTC, which='both')
    ax.tick_params(axis='y', colors=SlTC, which='both')
    return fig

def createHeatMap(df, phenoList, title, normAxis = None):
    
    # Create heatmap df
    heatMapDf = pd.DataFrame()

    for clust_label, group in df.groupby('clust_label'):
        clust_value_counts = group['Lineage'].value_counts()
        clust_value_counts.name = f'Cluster {clust_label}'

        heatMapDf = pd.concat([heatMapDf, pd.DataFrame([clust_value_counts])])

    # Fix the NA
    heatMapDf[heatMapDf.isna()] = 0
    heatMapDf = heatMapDf.astype('int')

    # Rearrange Columsn in order of prevalences
    heatMapDf = heatMapDf.loc[:, phenoList]

    if normAxis == 0:
        heatMapTitle = 'Phenotype/Cluster Heatmap: Normalized within Clusters'
        heatMapDf = round(heatMapDf.div(heatMapDf.sum(axis=1), axis=0), 3)
        vmin = 0
        vmax = 1
    elif normAxis == 1:
        heatMapTitle = 'Phenotype/Cluster Heatmap: Normalized within Phenotypes'
        heatMapDf = round(heatMapDf.div(heatMapDf.sum(axis=0), axis=1), 3)
        vmin = 0
        vmax = 1
    else:
        heatMapTitle = 'Phenotype/Cluster Heatmap: '
        vmin = heatMapDf.min().min()
        vmax = heatMapDf.max().max()
    title.append(heatMapTitle)

    figTitle = wrapTitleText(title)
    pltTitle = ''
    for i in figTitle:
        pltTitle = pltTitle + i + '\n'

    # Define Output Variables
    phenotypes = heatMapDf.columns
    clusters = heatMapDf.index

    # Theme Styles
    SlBgC  = '#0E1117'  # Streamlit Background Color
    SlTC   = '#FAFAFA'  # Streamlit Text Color
    Sl2BgC = '#262730'  # Streamlit Secondary Background Color

    fig = plt.figure(figsize = (12,12), facecolor = SlBgC)
    ax = fig.add_subplot(1,1,1, facecolor = SlBgC) 
    im = ax.imshow(heatMapDf, cmap = 'inferno', vmin = vmin, vmax = vmax)

    bbox = ax.get_yticklabels()[-1].get_window_extent()
    x, _ = ax.transAxes.inverted().transform([bbox.x0, bbox.y0])

    # Show all ticks and label them with the respective list entries
    ax.set_title(pltTitle, fontsize = 20, loc = 'left', color = SlTC, x=3*x, wrap=True)
    ax.set_xticks(np.arange(len(phenotypes)), labels = phenotypes, fontsize = 14, color = SlTC)
    ax.set_yticks(np.arange(len(clusters)), labels = clusters, fontsize = 14, color = SlTC)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
            rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    for i, cluster in enumerate(range(len(clusters))):
        for j, phenotype in enumerate(range(len(phenotypes))):
            value = heatMapDf.iloc[cluster, phenotype]
            if value >= 0.85*vmax:
                text_color = 'b'
            else:
                text_color = 'w'
            ax.text(j, i, value, ha="center", va="center", color=text_color)

    fig.tight_layout()

    return fig

def neighProfileDraw(spatial_umap, selClus, figsize=(14, 16)):
    import PlottingTools as umPT

    SlBgC  = '#0E1117'  # Streamlit Background Color
    SlTC   = '#FAFAFA'  # Streamlit Text Color
    Sl2BgC = '#262730'  # Streamlit Secondary Background Color

    NeiProFig = plt.figure(figsize=figsize, facecolor = SlBgC)
    ax = NeiProFig.add_subplot(1, 1, 1, facecolor = SlBgC)

    umPT.plot_mean_neighborhood_profile(ax, spatial_umap.dist_bin_um, spatial_umap.dens_df, selClus, maxDens=spatial_umap.maxdens_df, legF=1)

    return NeiProFig

def preprocess_weighted_umap(w, dfUMAP):

    # Raise everything about 0
    if np.any(w < 0):
        w = w + min(w)*-1.2

    # Check for NaN in the w and remove them
    notNAN = ~np.isnan(w)
    w = w[notNAN]
    dfUMAP = dfUMAP.loc[notNAN, :]

    # Apply the log
    w = np.log(0.1 * w + 0.1)
    w -= np.min(w)

    return w, dfUMAP

def UMAPdraw_density(df, bins, w, n_pad, vlim, feat = None, diff = False, figsize=(12, 12)):
    import PlottingTools as umPT

    # Streamlit Theming
    SlBgC  = '#0E1117'  # Streamlit Background Color
    SlTC   = '#FAFAFA'  # Streamlit Text Color
    Sl2BgC = '#262730'  # Streamlit Secondary Background Color

    # color maps
    cmap_viridis = plt.get_cmap('viridis').copy()
    cmap_viridis.set_under('white')
    cmap_magma = plt.get_cmap('magma').copy()
    cmap_magma.set_under('white')
    cmap_bwr = plt.get_cmap('bwr').copy()

    # Set up Figure
    UMAPFig = plt.figure(figsize=figsize, facecolor = SlBgC)
    ax = UMAPFig.add_subplot(1, 1, 1, facecolor = SlBgC)

    if w is None:
        cmap = cmap_viridis
        circle_type = None
    elif diff == False:
        cmap = cmap_magma
        circle_type = None
    else:
        cmap = cmap_bwr
        circle_type = 'arch'

    umPT.plot_2d_density(df['X'], df['Y'], bins=bins, w=w, n_pad=n_pad, 
                         ax=ax, cmap=cmap, vlim = vlim, circle_type = circle_type)
    
    xLim = ax.get_xlim()
    yLim = ax.get_ylim()
    
    ax.text(0.82*xLim[1], 0.03*yLim[1], 'Density', c = SlTC, fontsize = 25)

    if feat is not None:
        ax.text(xLim[0], 0.93*yLim[1], feat, c = SlTC, fontsize = 30)

    return UMAPFig

def drawIncidenceFigure(df, figTitle, phenotype = 'All Phenotypes', feature = 'Cell Counts', displayas = 'Counts Difference', compThresh = None, figsize=(12,12)):
    import PlottingTools as umPT

    SlBgC  = '#0E1117'  # Streamlit Background Color
    SlTC   = '#FAFAFA'  # Streamlit Text Color
    Sl2BgC = '#262730'  # Streamlit Secondary Background Color

    figTitle = wrapTitleText(figTitle)
    pltTitle = ''
    for i in figTitle:
        pltTitle = pltTitle + i + '\n'

    inciFig = plt.figure(figsize=figsize, facecolor = SlBgC)
    ax = inciFig.add_subplot(1, 1, 1, facecolor = SlBgC)

    if compThresh is not None:
        upTag = f' >= {compThresh}'
        dnTag = f' < {compThresh}'
    else:
        upTag = f' = 1'
        dnTag = f' = 0'

    if feature != 'Cell Counts':

        df = df[displayas]

        dfmin = df.loc[(df != np.nan)].min()
        dfmax = df.loc[(df != np.nan)].max()
        upLimit = max(-1*dfmin, dfmax)
        if displayas == 'Count Differences':
            if upLimit < 2:
                upLimit = 2
            ax.set_ylim([-1.05*upLimit, 1.05*upLimit])
            plt.axhline(y = 0, color = SlTC, linestyle = 'dashed', alpha = 0.7)
            ax.text(0.5, upLimit*.95, f'{feature}{upTag}', c = SlTC, fontsize = 30, alpha = 0.3)
            ax.text(0.5, -upLimit*.95, f'{feature}{dnTag}', c = SlTC, fontsize = 30, alpha = 0.3)
            outcome_suff = ' (Counts)'
        elif displayas == 'Ratios':
            ax.set_ylim([-1.05*upLimit, 1.05*upLimit])
            plt.axhline(y = 0, color = SlTC, linestyle = 'dashed', alpha = 0.7)
            ax.text(0.5, upLimit*.95, f'{feature}{upTag}', c = SlTC, fontsize = 30, alpha = 0.3)
            ax.text(0.5, -upLimit*.95, f'{feature}{dnTag}', c = SlTC, fontsize = 30, alpha = 0.3)
            outcome_suff = ' Ratio (log10)'
        elif displayas == 'Percentages':
            ax.set_ylim([-1.05, 1.05*upLimit])
            plt.axhline(y = 0, color = SlTC, linestyle = 'dashed', alpha = 0.7)
            ax.text(0.5, upLimit*.95, f'{feature}{upTag}', c = SlTC, fontsize = 30, alpha = 0.3)
            outcome_suff = ' (%)'

    else:
        df = df['counts']

        dfmin = df.min()
        dfmax = df.max()
        upLimit = max(-1*dfmin, dfmax)
        limrange = dfmax-dfmin
        liminc = limrange/8
        ax.set_ylim([dfmin-(liminc*0.1), dfmax + (liminc*0.1)])
        ax.text(0.5, upLimit*.95, f'{feature}', c = SlTC, fontsize = 30, alpha = 0.3)
        plt.axhline(y = 0, color = SlTC, linestyle = 'dashed', alpha = 0.7)
        outcome_suff = ' (Counts)'

    umPT.plot_incidence_line(ax, df, phenotype)
    
    # Reset xticks after 
    ax.set_xticks(df.index)
    ax.set_title(pltTitle, fontsize = 20, loc = 'left', color = SlTC)
    ax.set_xlabel('Cluster #', fontsize = 14, color = SlTC)
    ax.set_ylabel(f'{feature}{outcome_suff}', fontsize = 14, color = SlTC)

    return inciFig