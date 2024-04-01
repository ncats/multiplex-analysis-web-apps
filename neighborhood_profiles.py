'''
Set of scripts responsible for running the neighborhood profiles analysis
'''

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans # K-Means
import umap
from SpatialUMAP import SpatialUMAP
from scipy import ndimage as ndi

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
        self.n_bins = 200
        self.n_pad = 40
        self.vlim = .97

        self.w = None
        self.w_DiffA = None
        self.w_DiffB = None
        self.w_Diff  = None
        self.df_umapD = None

        self.umap_completed = False

        self.spatial_umap = None
        self.outcomes = None
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
        self.lineageDisplayToggle = 'Phenotypes'
        self.lineageDisplayToggle_clus = 'Phenotypes'

        # Unfiltered dropdown default options
        self.defLineageOpt    = 'All Phenotypes'
        self.defumapOutcomes  = 'No Outcome'
        self.definciOutcomes  = 'Cell Counts'

        # Default UMAP dropdown options
        self.umapPheno = [self.defLineageOpt]
        self.umapMarks = [self.defLineageOpt]
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
        Setup the requirements for running spatial UMAP
        '''

        spatial_umap = SpatialUMAP(dist_bin_um=np.array([25, 50, 100, 150, 200]), um_per_px=0.5, area_downsample=.2)
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
        spatial_umap.cell_labels = spatial_umap.cell_labels[spatial_umap.species]
        # set the region is to be analyzed (a TMA core is treated similar to a region of a interest)
        spatial_umap.region_ids = spatial_umap.cells.TMA_core_id.unique()
        # default cluster values
        spatial_umap.cells['clust_label'] = -1

        self.spatial_umap = spatial_umap

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
        self.df_umap = self.spatial_umap.df_umap

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

        session_state = self.filter_and_plot(session_state)

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

    @staticmethod
    def plot_2d_density(d, bins=200, n_pad=40, ax=None, cmap=plt.get_cmap('viridis'), vlim=np.array([0.001, 0.98]), circle_type='bg', box_off=True, return_matrix=False):
        '''
        plot_2d_density(X, Y, bins, n_pad, w, ax, gaussian_sigma, cmap, vlim, circle_type, box_off, return_matrix)
        is a method for drawing 2D histograms figures. In this particular instance, we are plotting the outputs of the UMAP.
        
        Parameters:
        X: X-values
        Y: Y-values
        bins: Bins. Vector tickMin and tickMax for X and Y direction.
        n_pad: number of elements to pad around a bin to make the spectrogram a little prettier. (I think)
        w: weights used to augment the density maps based on secondary phenotypes
        ax: Axes handle for the current figure
        gaussian_sigma: Sigma value for a gaussian filter with weight values
        cmap: color-map to use for the 2D Histogram
        vlim: Value range that the colormap should cover
        circle_type: Type of colormesh to return. ('bg', 'arch'). Will need to investigate further
        box_off: Flag to turn box around the figure off
        return_matrix: Flag for returning matrix. Default = 0 to not return density matrix, but intstead draws the figure.

        Returns:
        
        '''

        if return_matrix:
            return d
        else:
            if d[d > 0].shape == (0,):
                vmin = 0
                vmax = 1
                vlim = [vmin, vmax]
            else:
                if ~all((d < 0)[0]):
                    if np.isscalar(vlim):
                        vlim = np.array([0, np.quantile(d[d > 0].flatten(), vlim)])
                    else:
                        if np.all((vlim < 1) & (vlim > 0)):
                            vlim = np.quantile(d[d > 0].flatten(), vlim)
                else:
                    vlim = np.zeros(200)

            if ax is None:
                _, ax = plt.subplots()

            if np.isscalar(bins):
                n_bins = bins
            else:
                n_bins = len(bins[0]) - 1

            if circle_type == 'bg':
                extend = 'max'
                c = np.meshgrid(np.arange(2 * n_pad + n_bins), np.arange(2 * n_pad + n_bins))
                c = np.sqrt(((c[0] - ((2 * n_pad + n_bins) / 2)) ** 2) + ((c[1] - ((2 * n_pad + n_bins) / 2)) ** 2)) < (0.95 * ((2 * n_pad + n_bins) / 2))
                ax.pcolormesh(np.pad(d, [n_pad, n_pad]) + c, vmin=1, vmax=1 + vlim[1], cmap=cmap, shading='gouraud', alpha=1)
                # ax.pcolormesh(np.log10(np.pad(d, [n_pad, n_pad]) + c + 1), vmin=np.log10(2), vmax=np.log10(2 + vlim[1]), cmap=cmap, shading='gouraud', alpha=1)
            elif circle_type == 'arch':
                extend = 'both'
                if any((d<0)[0]):
                    vmin = np.quantile(d[d < 0].flatten(), 0.03)
                else:
                    vmin = 0.03
                if any((d>0)[0]):
                    vmax = np.quantile(d[d > 0].flatten(), 0.97)
                else:
                    vmax = 0.97
                c = (n_bins / 2)
                ax.add_artist(plt.Circle((c + n_pad, c + n_pad), 0.95 * (c + n_pad), color='black', fill=False))
                ax.pcolormesh(np.pad(d, [n_pad, n_pad]), vmin=vmin, vmax=vmax, cmap=cmap, shading='gouraud', alpha=1)
            else:
                extend = 'max'
                ax.pcolormesh(np.pad(d, [n_pad, n_pad]), vmin=0, vmax=vlim[1], cmap=cmap, shading='gouraud', alpha=1)

            cax = ax.inset_axes([0.95, 0.1, 0.01, 0.85])
            plt_cmap(ax=cax, cmap=cmap, extend=extend, width=0.01)

            if box_off is True:
                [ax.spines[sp].set_visible(False) for sp in ax.spines]
                ax.set(xticks=[], yticks=[])


@staticmethod
def plt_cmap(ax, cmap, extend, width, ylabel = None):
    '''plt_cmap(ax, cmap, extend, width, ylabel) draws a colorbar for the current colormap at the correct
    axes location, and with the correct label.

    Parameters:
    ax: Matplotlib axes handle
    cmap: Matplotlib colormap
    extend: {'neither', 'both', 'min', 'max'}  
            Make pointed end(s) for out-of-range values (unless 'neither'). 
            These are set for a given colormap using the colormap set_under and set_over methods.
    width: Width of the colorbar in Figure coordinates. '0.01' suggested value
    ylabel: String of the colomap label

    Returns:

    '''
    cb = mpl.colorbar.Colorbar(ax=ax, cmap=cmap, extend=extend)
    cb.set_ticks([])
    pos = ax.get_position().bounds
    ax.set_position([pos[0], pos[1], width, pos[3]])
    if ylabel is not None:
        ax.set(ylabel=ylabel)

class UMAPDensityProcessing():
    '''
    Individual processing of UMAP density values
    '''
    def __init__(self, npf, df):
        self.df = df

        self.n_bins = npf.n_bins
        self.n_pad = npf.n_pad
        self.vlim = npf.vlim

        # Streamlit Coloring (slc)
        self.slc_bg   = npf.slc_bg   # Streamlit Color -Background
        self.slc_text = npf.slc_text # Streamlit Color -Text
        self.slc_bg2  = npf.slc_bg2  # Streamlit Color -Secondary Background

        self.prepare_umap_density(self.df['X'], self.df['Y'])


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

    def UMAPdraw_density(self, d, bins, w=None, n_pad, vlim, feat = None, diff = False, figsize=(12, 12)):
        '''
        Draw the UMAP density data
        '''

        # color maps
        cmap_viridis = plt.get_cmap('viridis').copy()
        cmap_viridis.set_under('white')
        cmap_magma = plt.get_cmap('magma').copy()
        cmap_magma.set_under('white')
        cmap_bwr = plt.get_cmap('bwr').copy()

        # Set up Figure
        umap_fig = plt.figure(figsize=figsize, facecolor = self.slc_bg)
        ax = umap_fig.add_subplot(1, 1, 1, facecolor = self.slc_bg)

        if w is None:
            cmap = cmap_viridis
            circle_type = None
        elif diff is False:
            cmap = cmap_magma
            circle_type = None
        else:
            cmap = cmap_bwr
            circle_type = 'arch'

        self.plot_2d_density(d, bins=bins, w=w, n_pad=n_pad,
                                ax=ax, cmap=cmap, vlim = vlim, circle_type = circle_type)

        x_lim = ax.get_xlim()
        y_lim = ax.get_ylim()

        ax.text(0.82*x_lim[1], 0.03*y_lim[1], 'Density', c = self.slc_text, fontsize = 25)

        if feat is not None:
            ax.text(x_lim[0], 0.93*y_lim[1], feat, c = self.slc_text, fontsize = 30)

        return umap_fig

    @staticmethod
    def plot_2d_density(d, bins=200, n_pad=40, ax=None, cmap=plt.get_cmap('viridis'), vlim=np.array([0.001, 0.98]), circle_type='bg', box_off=True, return_matrix=False):
        '''
        plot_2d_density(X, Y, bins, n_pad, w, ax, gaussian_sigma, cmap, vlim, circle_type, box_off, return_matrix)
        is a method for drawing 2D histograms figures. In this particular instance, we are plotting the outputs of the UMAP.
        
        Parameters:
        X: X-values
        Y: Y-values
        bins: Bins. Vector tickMin and tickMax for X and Y direction.
        n_pad: number of elements to pad around a bin to make the spectrogram a little prettier. (I think)
        w: weights used to augment the density maps based on secondary phenotypes
        ax: Axes handle for the current figure
        gaussian_sigma: Sigma value for a gaussian filter with weight values
        cmap: color-map to use for the 2D Histogram
        vlim: Value range that the colormap should cover
        circle_type: Type of colormesh to return. ('bg', 'arch'). Will need to investigate further
        box_off: Flag to turn box around the figure off
        return_matrix: Flag for returning matrix. Default = 0 to not return density matrix, but intstead draws the figure.

        Returns:
        
        '''

        if return_matrix:
            return d
        else:
            if d[d > 0].shape == (0,):
                vmin = 0
                vmax = 1
                vlim = [vmin, vmax]
            else:
                if ~all((d < 0)[0]):
                    if np.isscalar(vlim):
                        vlim = np.array([0, np.quantile(d[d > 0].flatten(), vlim)])
                    else:
                        if np.all((vlim < 1) & (vlim > 0)):
                            vlim = np.quantile(d[d > 0].flatten(), vlim)
                else:
                    vlim = np.zeros(200)

            if ax is None:
                _, ax = plt.subplots()

            if np.isscalar(bins):
                n_bins = bins
            else:
                n_bins = len(bins[0]) - 1

            if circle_type == 'bg':
                extend = 'max'
                c = np.meshgrid(np.arange(2 * n_pad + n_bins), np.arange(2 * n_pad + n_bins))
                c = np.sqrt(((c[0] - ((2 * n_pad + n_bins) / 2)) ** 2) + ((c[1] - ((2 * n_pad + n_bins) / 2)) ** 2)) < (0.95 * ((2 * n_pad + n_bins) / 2))
                ax.pcolormesh(np.pad(d, [n_pad, n_pad]) + c, vmin=1, vmax=1 + vlim[1], cmap=cmap, shading='gouraud', alpha=1)
                # ax.pcolormesh(np.log10(np.pad(d, [n_pad, n_pad]) + c + 1), vmin=np.log10(2), vmax=np.log10(2 + vlim[1]), cmap=cmap, shading='gouraud', alpha=1)
            elif circle_type == 'arch':
                extend = 'both'
                if any((d<0)[0]):
                    vmin = np.quantile(d[d < 0].flatten(), 0.03)
                else:
                    vmin = 0.03
                if any((d>0)[0]):
                    vmax = np.quantile(d[d > 0].flatten(), 0.97)
                else:
                    vmax = 0.97
                c = (n_bins / 2)
                ax.add_artist(plt.Circle((c + n_pad, c + n_pad), 0.95 * (c + n_pad), color='black', fill=False))
                ax.pcolormesh(np.pad(d, [n_pad, n_pad]), vmin=vmin, vmax=vmax, cmap=cmap, shading='gouraud', alpha=1)
            else:
                extend = 'max'
                ax.pcolormesh(np.pad(d, [n_pad, n_pad]), vmin=0, vmax=vlim[1], cmap=cmap, shading='gouraud', alpha=1)

            cax = ax.inset_axes([0.95, 0.1, 0.01, 0.85])
            plt_cmap(ax=cax, cmap=cmap, extend=extend, width=0.01)

            if box_off is True:
                [ax.spines[sp].set_visible(False) for sp in ax.spines]
                ax.set(xticks=[], yticks=[])