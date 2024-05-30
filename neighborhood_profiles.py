'''
Set of scripts responsible for running the neighborhood profiles analysis

Class NeighborhoodProfiles:
    Organization of the methods and attributes that are required to run
    the neighborhood profiles analysis

Class UMAPDensityProcessing:
    Individual processing of UMAP density matrices
'''

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans # K-Means
import umap
from SpatialUMAP import SpatialUMAP
from scipy import ndimage as ndi

import basic_phenotyper_lib as bpl  # Useful functions for cell phenotyping
import nidap_dashboard_lib as ndl   # Useful functions for dashboards connected to NIDAP
from benchmark_collector import benchmark_collector # Benchmark Collector Class
import PlottingTools as umPT

class NeighborhoodProfiles:
    '''
    Organization of the methods and attributes that are required to run
    the neighborhood profiles analysis
    '''
    def __init__(self, bc = None):


        self.clust_minmax = [0, 40]
        self.reset_neigh_profile_settings()

        if bc is None:
            bc = benchmark_collector()
        self.bc = bc

        # Spectrogram Plotting Settings
        self.n_bins = 100
        self.n_pad  = 0
        self.vlim = .97
        self.xx = 0
        self.yy = 0

        self.w = None
        self.w_DiffA = None
        self.w_DiffB = None
        self.w_Diff  = None
        self.df_umapD = None

        self.umap_completed = False

        self.spatial_umap = None
        self.outcomes     = None
        self.df_umap      = pd.DataFrame()
        self.df_umap_filt = pd.DataFrame()

        # Streamlit Coloring (slc)
        self.slc_bg   = '#0E1117' # Streamlit Color -Background
        self.slc_text = '#FAFAFA' # Streamlit Color -Text
        self.slc_bg2  = '#262730' # Streamlit Color -Secondary Background

    def reset_neigh_profile_settings(self):
        '''
        Resets all the variables required for neighborhood
        profiles analysis 
        '''

        print('Resetting Neighborhood Profiles Analysis Settings')

        # Has the UMAP been completed yet?
        self.cell_counts_completed = False
        self.umap_completed        = False
        self.clustering_completed  = False
        self.UMAPFigType           = 'Density'

        # UMAP Lineage Display
        self.lineageDisplayToggle      = 'Phenotypes'
        self.lineageDisplayToggle_clus = 'Phenotypes'

        # Unfiltered dropdown default options
        self.defLineageOpt    = 'All Phenotypes'
        self.defumapOutcomes  = 'No Outcome'
        self.definciOutcomes  = 'Cell Counts'

        # Default UMAP dropdown options
        self.umapPheno    = [self.defLineageOpt]
        self.umapMarks    = [self.defLineageOpt]
        self.umaplineages = [self.defLineageOpt]
        self.umapOutcomes = [self.defumapOutcomes]

        # Default Incidence dropdown options
        self.inciOutcomes = [self.definciOutcomes]

        # Default UMAPInspect settings
        self.umapInspect_Ver = self.defLineageOpt
        self.umapInspect_Feat = self.defumapOutcomes

        # Default UMAP differences settings
        self.diffUMAPSel_Ver  = self.defLineageOpt
        self.diffUMAPSel_Feat = self.defumapOutcomes

        # Default Incidence settings
        self.inciPhenoSel   = self.defLineageOpt
        self.inciOutcomeSel = self.definciOutcomes
        self.Inci_Value_display = 'Count Differences'

    def setup_spatial_umap(self, df, marker_names, pheno_order):
        '''
        Silly I know. I will fix it later
        '''

        self.spatial_umap = bpl.setup_Spatial_UMAP(df, marker_names, pheno_order)

    def perform_density_calc(self, cpu_pool_size = 1):
        '''
        Calculate the cell counts, cell areas,
        perform the cell densities and cell proportions analyses.

        This is using Andrew's code to calculate the cell counts

        Args:
            spatial_umap (SpatialUMAP): SpatialUMAP object
            bc (benchmark_collector): Benchmark Collector object
            cpu_pool_size (int): Number of CPUs to use for parallel processing

        Returns:
            SpatialUMAP: SpatialUMAP object with the cell counts, cell areas, 
                        cell densities and cell proportions analyses performed
        '''

        # clear metrics
        self.spatial_umap.clear_counts()
        self.spatial_umap.clear_areas()

        print(f'cpu_pool_size set to {cpu_pool_size}\n')

        # get the counts per cell and save to pickle file
        print('Starting Cell Counts process')
        self.bc.startTimer()
        self.spatial_umap.get_counts_And()
        self.bc.printElapsedTime(f'Calculating Counts for {len(self.spatial_umap.cells)} cells')

        # get the areas of cells and save to pickle file
        area_threshold = 0.001
        print('\nStarting Cell Areas process')
        self.spatial_umap.get_areas(area_threshold, pool_size=cpu_pool_size)

        # calculate density based on counts of cells / area of each arc examine
        self.spatial_umap.calc_densities(area_threshold)
        # calculate proportions based on species counts/# cells within an arc
        self.spatial_umap.calc_proportions(area_threshold)

    def perform_spatial_umap(self, session_state, umap_style = 'density'):
        '''
        Perform the spatial UMAP analysis

        Args:
            spatial_umap (spatial_umap): spatial_umap object
            bc (benchmark_collector): Benchmark Collector object
            UMAPStyle (str): Style of UMAP to use
        
        Returns:
            spatial_umap: spatial_umap object with the UMAP analysis performed
        '''

        # set training and "test" cells for umap training and embedding, respectively
        print('Setting Train/Test Split')
        self.spatial_umap.set_train_test(n=2500, groupby_label = 'TMA_core_id', seed=54321)

        # fit umap on training cells
        self.bc.startTimer()
        print('Fitting Model')
        self.spatial_umap.umap_fit = umap.UMAP().fit(self.spatial_umap.density[self.spatial_umap.cells['umap_train'].values].reshape((self.spatial_umap.cells['umap_train'].sum(), -1)))
        self.bc.printElapsedTime(f'      Fitting {np.sum(self.spatial_umap.cells["umap_train"] == 1)} points to a model')

        # Transform test cells based on fitted model
        self.bc.startTimer()
        print('Transforming Data')
        self.spatial_umap.umap_test = self.spatial_umap.umap_fit.transform(self.spatial_umap.density[self.spatial_umap.cells['umap_test'].values].reshape((self.spatial_umap.cells['umap_test'].sum(), -1)))
        self.bc.printElapsedTime(f'      Transforming {np.sum(self.spatial_umap.cells["umap_test"] == 1)} points with the model')

        # Identify all of the features in the dataframe
        self.outcomes = self.spatial_umap.cells.columns

        self.spatial_umap.prepare_df_umap_plotting(self.outcomes)

        # Setup the session_state default parameters

        # List of possible UMAP Lineages as defined by the completed UMAP
        session_state.umapPheno = [session_state.defLineageOpt]
        session_state.umapPheno.extend(session_state.pheno_summ['phenotype'])
        session_state.umapMarks = [session_state.defLineageOpt]
        session_state.umapMarks.extend(self.spatial_umap.markers)
        session_state.umapMarks.extend(['Other'])

        # List of possible outcome variables as defined by the config yaml files
        session_state.umapOutcomes = [session_state.defumapOutcomes]
        session_state.umapOutcomes.extend(self.outcomes)
        session_state.inciOutcomes = [session_state.definciOutcomes]
        session_state.inciOutcomes.extend(self.outcomes)

        # Perform possible cluster variations with the completed UMAP
        # session_state.bc.startTimer()
        # session_state.clust_range, session_state.wcss = bpl.measure_possible_clust(session_state.spatial_umap, clust_minmax)
        # session_state.bc.printElapsedTime(msg = 'Calculating possible clusters')

        session_state.wcss_calc_completed = True
        session_state.umap_completed = True
        session_state.umapCompleted = True

        return session_state

    @staticmethod
    def kmeans_calc(spatial_umap, n_clusters = 5):
        '''
        Perform KMeans clustering on the spatial UMAP data

        Args:
            spatial_umap (spatial_umap): spatial_umap object
            nClus (int): Number of clusters to use
        
        Returns:
            kmeans_obj: KMeans obj created from KMeans
        '''
        # Create KMeans object for a chosen cluster
        kmeans_obj = KMeans(n_clusters = n_clusters,
                            init ='k-means++',
                            max_iter = 300,
                            n_init = 10,
                            random_state = 42)
        # Fit the data to the KMeans object
        kmeans_obj.fit(spatial_umap.umap_test)

        return kmeans_obj

    def measure_possible_clust(self, clust_minmax = [0, 40]):
        '''
        Method for measuring the within-cluster sum of squares for
        a range of cluster values

        Args:
            clust_minmax (list): List of min and max cluster values to use

        Returns:
            clust_range (list): List of cluster values
            wcss (list): List of within-cluster sum of squares
        '''
        clust_range = range(clust_minmax[0], clust_minmax[1])
        wcss = [] # Within-Cluster Sum of Squares
        for n_clusters in clust_range:
            # Perform clustering for chosen
            kmeans_obj = self.kmeans_calc(self.spatial_umap, n_clusters)
            # Append Within-Cluster Sum of Squares measurement
            wcss.append(kmeans_obj.inertia_)
        return list(clust_range), wcss

    def perform_clustering(self, n_clusters):
        '''
        perform clustering for the UMAP data using KMeans

        Args:
            spatial_umap (spatial_umap): spatial_umap object
            n_clusters (int): Number of clusters to use

        Returns:
            spatial_umap: spatial_umap object with the clustering performed
        '''

        # Reperform clustering for chosen cluster value
        kmeans_obj = self.kmeans_calc(self.spatial_umap, n_clusters)
        # Add cluster label column to cells dataframe
        self.spatial_umap.df_umap.loc[:, 'clust_label'] = kmeans_obj.labels_
        self.spatial_umap.df_umap.loc[:, 'cluster'] = kmeans_obj.labels_
        self.spatial_umap.df_umap.loc[:, 'Cluster'] = kmeans_obj.labels_

        # After assigning the cluster labels, perform mean measure calculations
        self.spatial_umap.mean_measures()

    def filter_and_plot(self, session_state):
        '''
        function to update the filtering and the figure plotting
        '''
        if self.umap_completed:
            self.df_umap_filt = self.df_umap.loc[self.df_umap['Slide ID'] == session_state['selSlide ID'], :]
            session_state = ndl.setFigureObjs_UMAP(session_state)

        return session_state

    def preprocess_weighted_umap(self, w, df_umap):
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

    def setup_spectrogram_settings(self):
        self.xx = np.linspace(np.min(self.df_umap['X']), np.max(self.df_umap['X']), self.n_bins + 1)
        self.yy = np.linspace(np.min(self.df_umap['Y']), np.max(self.df_umap['Y']), self.n_bins + 1)

        self.minXY = self.df_umap[['X', 'Y']].min()-1
        self.maxXY = self.df_umap[['X', 'Y']].max()+1

    def perform_df_umap_diff(self, sel_feat):
        '''
        Create the weights which contribute to the 
        '''

        self.setup_spectrogram_settings()

        feat_comp1 = '= 1'
        feat_comp2 = '= 0'

        self.w = self.df_umap[sel_feat]
        self.w, self.df_umapD = self.preprocess_weighted_umap(self.w, self.df_umap)

        self.w_DiffA = self.w
        self.w_DiffB = max(self.w) - self.w
        self.w_Diff  = self.w_DiffA - self.w_DiffB

        feat_label0 = f'{sel_feat} {feat_comp1} '
        feat_label1 = f'{sel_feat} {feat_comp2} '

        self.UMAPFigDiff0_Dens = self.UMAPdraw_density(self.df_umapD, bins = [self.xx, self.yy], w=self.w_DiffA, n_pad=self.n_pad, vlim=self.vlim, feat = feat_label0)
        self.UMAPFigDiff1_Dens = self.UMAPdraw_density(self.df_umapD, bins = [self.xx, self.yy], w=self.w_DiffB, n_pad=self.n_pad, vlim=self.vlim, feat = feat_label1)
        self.UMAPFigDiff2_Dens = self.UMAPdraw_density(self.df_umapD, bins = [self.xx, self.yy], w=self.w_Diff,  n_pad=self.n_pad, vlim=self.vlim, diff = True)

    def prepare_umap_density(self, X, Y=None, w=None, gaussian_sigma=0.5):
        '''
        create the density matrix for the UMAP data
        '''

        if Y is not None:
            if w is not None:
                b, _, _ = np.histogram2d(X, Y, bins=self.n_bins)
                b = ndi.gaussian_filter(b.T, sigma=gaussian_sigma)

                s, _, _ = np.histogram2d(X, Y, bins=self.n_bins, weights=w)
                s = ndi.gaussian_filter(s.T, sigma=gaussian_sigma)

                d = np.zeros_like(b)
                # d[b > 0] = s[b > 0] / b[b > 0]
                d = s
                d = ndi.gaussian_filter(d, sigma=gaussian_sigma)
            else:
                d, _, _ = np.histogram2d(X, Y, bins=self.n_bins)
                d /= np.sum(d)
                d = ndi.gaussian_filter(d.T, sigma=gaussian_sigma)
        else:
            d = X

        self.d = d

    @staticmethod
    def UMAPdraw_density(d, bins, w, n_pad, vlim, feat = None, diff = False, figsize=(12, 12)):
        '''
        Draw the UMAP density data
        '''

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
        umap_fig = plt.figure(figsize=figsize, facecolor = SlBgC)
        ax = umap_fig.add_subplot(1, 1, 1, facecolor = SlBgC)

        if w is None:
            cmap = cmap_viridis
            circle_type = None
        elif diff is False:
            cmap = cmap_magma
            circle_type = None
        else:
            cmap = cmap_bwr
            circle_type = 'arch'

        umPT.plot_2d_density(d, bins=bins, w=w, n_pad=n_pad, 
                            ax=ax, cmap=cmap, vlim = vlim, circle_type = circle_type)
        
        xLim = ax.get_xlim()
        yLim = ax.get_ylim()
        
        ax.text(0.82*xLim[1], 0.03*yLim[1], 'Density', c = SlTC, fontsize = 25)

        if feat is not None:
            ax.text(xLim[0], 0.93*yLim[1], feat, c = SlTC, fontsize = 30)

        return umap_fig
class UMAPDensityProcessing():
    '''
    Individual processing of UMAP density matrices
    '''
    def __init__(self, npf, df, xx=None, yy=None):
        self.df = df

        self.n_bins = npf.n_bins
        self.n_pad  = npf.n_pad
        self.vlim   = npf.vlim

        # Streamlit Coloring (slc)
        self.slc_bg   = npf.slc_bg   # Streamlit Color -Background
        self.slc_text = npf.slc_text # Streamlit Color -Text
        self.slc_bg2  = npf.slc_bg2  # Streamlit Color -Secondary Background

        # Preset Summary Stats
        self.dfmin    = 0
        self.dfmax    = 0
        self.dens_min = 0
        self.dens_max = 0
        self.minabs   = 0

        # Feature Label
        self.feat_label = None
        self.cluster_dict = None
        self.palette_dict = None

        if xx is not None:
            self.xx = xx
            self.yy = yy
        else:
            self.xx = np.linspace(np.min(self.df['X']), np.max(self.df['X']), self.n_bins + 1)
            self.yy = np.linspace(np.min(self.df['Y']), np.max(self.df['Y']), self.n_bins + 1)

        self.prepare_umap_density(self.df['X'], self.df['Y'])

    def prepare_umap_density(self, x, y, w=None):
        '''
        create the density matrix from the UMAP data
        '''

        self.dens_mat, \
        self.bin_indices_df_group,\
        self.empty_bin_ind = umPT.plot_2d_density(x, y, bins = [self.xx, self.yy],
                                                  w = w, return_matrix = True)

        self.umap_summary_stats()

    def umap_summary_stats(self):
        '''
        Identify the minimum and maximum values of the density matrix
        '''

        self.dfmin  = self.df[['X', 'Y']].min()
        self.dfmax  = self.df[['X', 'Y']].max()

        self.dens_min = np.min(self.dens_mat)
        self.dens_max = np.max(self.dens_mat)
        self.minabs   = np.min([np.abs(self.dens_min), np.abs(self.dens_max)])

    def filter_by_lineage(self, display_toggle, drop_val, default_val):
        '''
        Function for filtering UMAP function based on Phenotypes or Markers

        Args:
            display_toggle (str): Toggle to display as Phenotypes or Markers
            drop_val (str): Value selected from the drop value
            default_val (str): Default Value of phenotyping or markers

        Returns:
            None
        '''
        if drop_val != default_val:
            if display_toggle == 'Phenotypes':
                self.df = self.df.loc[self.df['Lineage'] == drop_val, :]
            elif display_toggle == 'Markers':
                self.df = self.df.loc[self.df['species_name_short'].str.contains(drop_val), :]

    def split_df_by_feature(self, feature):
        '''
        split_df_by_feature takes in a feature from a dataframe
        and first identifies if the feature is boolean, if it contains 
        float values, or neither. If its a boolean, it will split the
        dataframe between values of 0 and 1 for the selected feature.
        If the feature is a float, it will split the dataframe based on
        the median value of the feature. If the feature is neither boolean
        nor float, it will not split the dataframe. 

        In all cases this function will return a dictionary of the outcome
        of the split with the most importannt value being, appro_feat, 
        which will be True if the feature is appropriate for splitting, and
        False if not.

        Args:
            feature (str): Feature to split the dataframe by

        Returns:
            split_dict (dict): Dictionary of the outcomes of splitting
             the dataframe
        '''

        split_dict = dict()
        # Idenfify the column type that is splitting the UMAP
        col_type = ndl.identify_col_type(self.df[feature])

        if col_type == 'not_bool':
            # Identify UMAP by Condition
            median = np.round(self.df[feature].median(), 2)
            split_dict['df_umap_fals'] = self.df.loc[self.df[feature] <= median, :]
            split_dict['df_umap_true'] = self.df.loc[self.df[feature] > median, :]
            split_dict['fals_msg']   = f'<= {median}'
            split_dict['true_msg']   = f'> {median}'
            split_dict['appro_feat'] = True
        elif col_type == 'bool':
            # Identify UMAP by Condition
            values = self.df[feature].unique()
            split_dict['df_umap_fals'] = self.df.loc[self.df[feature] == values[0], :]
            split_dict['df_umap_true'] = self.df.loc[self.df[feature] == values[1], :]
            split_dict['fals_msg']   = f'= {values[0]}'
            split_dict['true_msg']   = f'= {values[1]}'
            split_dict['appro_feat'] = True
        else:
            split_dict['appro_feat'] = False

        return split_dict

    def set_feature_label(self, feature, feat_label):
        '''
        Setting feature label
        '''
        self.feat_label = f'{feature} {feat_label}'

    def UMAPdraw_density(self, w=None, diff = False, figsize=(12, 12), legendtype = 'colorbar'):
        '''
        Calls the UMAPdraw_density function from PlottingTools.py
        '''

        return bpl.UMAPdraw_density(d = self.dens_mat,
                                    bins = [self.xx, self.yy],
                                    w = w,
                                    n_pad = self.n_pad,
                                    vlim = self.vlim,
                                    feat = self.feat_label,
                                    diff = diff,
                                    figsize = figsize,
                                    legendtype = legendtype)

    def umap_draw_clusters(self, figsize = (12, 12)):
        '''
        Draw the UMAP colored by clusters
        '''

        umap_clust_fig, ax = bpl.draw_scatter_fig(figsize = figsize)
        umap_clust_fig = bpl.scatter_plot(self.df, umap_clust_fig, ax, 'Clusters',
                                          xVar = 'X', yVar = 'Y', hueVar='clust_label',
                                          xLim = [self.dfmin[0], self.dfmax[0]],
                                          yLim = [self.dfmin[1], self.dfmax[1]],
                                          hueOrder = self.cluster_dict.values(),
                                          palette  = self.palette_dict)

        return umap_clust_fig

    def filter_density_matrix(self, cutoff= 0.01, empty_bin_ind = None):
        '''
        Filter the density matrix based on the cutoff value
        to create a binary mask of values that are above or below
        the cutoff value. 

        This takes a list of empty bin indices to filter out any
        bins that are meant to be empty (0) no matter what the actual
        value is in the bin.

        Args:
            cutoff (float): Cutoff value to use for filtering
            empty_bin_ind (list): List of empty bin indices
        
        Returns:
            None
        '''

        dens_mat_shape = self.dens_mat.shape

        # Filtering and Masking
        for x_bin in range(dens_mat_shape[0]):
            for y_bin in range(dens_mat_shape[1]):
                if tuple([x_bin, y_bin]) in empty_bin_ind:
                    self.dens_mat[x_bin, y_bin] = 0
                else:
                    if self.dens_mat[x_bin, y_bin] > cutoff:
                        self.dens_mat[x_bin, y_bin] = 1
                    elif self.dens_mat[x_bin, y_bin] < -cutoff:
                        self.dens_mat[x_bin, y_bin] = -1

    def perform_clustering(self, dens_mat_cmp, num_clus_0, num_clus_1):
        '''
        Sets up clustering
        '''

        print(f'Performing Clustering with {num_clus_0} clusters for Negative Condition and {num_clus_1} clusters for Positive Condition')
        # Perform k-menas clustering for the Negative Condition
        kmeans_obj_cond0 = KMeans(n_clusters = num_clus_0,
                                  init ='k-means++',
                                  max_iter = 300,
                                  n_init = 50)
        
        # Perform k-menas clustering for the Positive Condition
        kmeans_obj_cond1 = KMeans(n_clusters = num_clus_1,
                                  init ='k-means++',
                                  max_iter = 300,
                                  n_init = 50)

        # Identify the indices of the negative condition
        cond0_ind = np.nonzero(dens_mat_cmp == 1)
        cells_cond0 = np.vstack(cond0_ind).T
        # Fit the negative condition to the negative kmeans object
        kmeans_obj_cond0.fit(cells_cond0)

        # Identify the indices of the positive condition
        cond1_ind = np.nonzero(dens_mat_cmp == -1)
        cells_cond1 = np.vstack(cond1_ind).T
        # Fit the positive condition to the positive kmeans object
        kmeans_obj_cond1.fit(cells_cond1)

        # Replace the labels in the density matrix with the cluster labels
        self.dens_mat[cond0_ind] = kmeans_obj_cond0.labels_ + 1
        self.dens_mat[cond1_ind] = -kmeans_obj_cond1.labels_ - 1

        unique_fals, counts_fals = np.unique(self.dens_mat[cond0_ind], return_counts=True)
        unique_true, counts_true = np.unique(self.dens_mat[cond1_ind], return_counts=True)

        unique_set_fals = pd.DataFrame(data = {'vals': unique_fals, 'counts': counts_fals}).sort_values('counts', ascending = False, ignore_index = True).astype(int)
        unique_set_true = pd.DataFrame(data = {'vals': unique_true, 'counts': counts_true}).sort_values('counts', ascending = False, ignore_index = True).astype(int)

        # Set up the cluster dictionary
        self.cluster_dict = dict()
        self.cluster_dict[0] = 'No Cluster'
        for i in unique_set_fals.index:
            self.cluster_dict[unique_set_fals.vals[i]] = f'False Cluster {i+1}'
        for i in unique_set_true.index:
            self.cluster_dict[unique_set_true.vals[i]] = f'True Cluster {i+1}'

        set_blues = sns.color_palette('Blues_r', 10)
        set_reds = sns.color_palette('Reds_r', 10)

        # Set Palette Dictionary
        self.palette_dict = dict()
        self.palette_dict['No Cluster'] = 'white'
        for i in unique_set_fals.index:
            self.palette_dict[f'False Cluster {i+1}'] = set_reds[i]
        for i in unique_set_true.index:
            self.palette_dict[f'True Cluster {i+1}'] = set_blues[i]

    def kmean_calc(self, dens_mat, n_clusters, cond):
        '''
        Perform clustering on the density matrix

        Args:
            dens_mat (numpy array): Density matrix
            n_clusters (int): Number of clusters to use
            cond (int): Condition to use for clustering
        
        Returns:
            kmeans_obj: KMeans object created from KMeans
        '''

        kmeans_obj = KMeans(n_clusters = n_clusters,
                            init ='k-means++',
                            max_iter = 300,
                            n_init = 50)

        # Identify the indices of the target condition
        cond_ind = np.nonzero(dens_mat == cond)
        cells_cond = np.vstack(cond_ind).T
        # Fit the condition to the kmeans object
        kmeans_obj.fit(cells_cond)

        return kmeans_obj
