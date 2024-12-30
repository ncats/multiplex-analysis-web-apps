'''
Author: Alex Baras, MD, PhD (https://github.com/alexbaras)
NCATS Maintainer: Dante J Smith, PhD (https://github.com/djsmith17)
'''
import time
import multiprocessing as mp
from functools import partial
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy import optimize
from scipy import ndimage as ndi
from scipy.spatial import ConvexHull
from sklearn.metrics.pairwise import euclidean_distances
from skimage import draw as skdraw, transform as sktran
np.seterr(divide='ignore', invalid='ignore')

import utils

class SpatialUMAP:
    '''
    SpatialUMAP() is a class that handles the cell count 
    processing and concentric circle area measurements for 
    cell density analysis in SpatialUMAP studies.

    Methods:
    * __init__
    * construct_arcs
    * process_cell_areas
    * process_cell_counts
    * clear_counts
    * clear_areas
    * start_pool
    * close_pool
    * process_region_counts
    * process_region_areas
    * get_counts
    * get_areas
    * set_train_test

    Parameters:
    dist_bin_um: Concentric radial distance bins for assessing areas of 
        cell densities
    um_per_px: Conversion of micrometers to pixels
    area_downsample: Downsample rate of area? Not exactly sure.

    '''

    @staticmethod
    def get_dataframes(results):
        for df in results:
            yield df

    @staticmethod
    def construct_arcs(dist_bin_px):
        '''
        construct_arcs() creates a boolean mask of concentric circles
        based on the distance bins in pixels.

        Args:
            dist_bin_px (np.array): distance bins in pixels

        Returns:
            np.array: boolean mask of concentric circles
        '''
        # set bool mask of the arcs
        arcs = np.zeros([int(2 * dist_bin_px[-1]) + 1] * 2 + [len(dist_bin_px), ], dtype=bool)
        for i in range(len(dist_bin_px)):
            # circle based on radius
            rr, cc = skdraw.disk(center=(np.array(arcs.shape[:2]) - 1) / 2, radius=dist_bin_px[i] + 1, shape=arcs.shape[:2])
            arcs[rr, cc, i] = True
        # difference logic to produce arcs
        return np.stack([arcs[:, :, 0]] + [arcs[:, :, i] != arcs[:, :, i - 1] for i in range(1, arcs.shape[2])], axis=2)

    @staticmethod
    def process_cell_areas(idx, region_id, cell_positions, dist_bin_px, img_mask, arcs):
        '''Processing the cell_area information'''

        # Print the image name
        # print(f'Calculating neighborhood area for image {region_id} (cells)...')

        # Record the start time
        start_time = time.time()

        region_area_df = pd.DataFrame()
        for i, idx_i in enumerate(idx):
            # true bounds to match arcs
            bounds = np.array([cell_positions[i].astype(int) - dist_bin_px[-1].astype(int), dist_bin_px[-1].astype(int) + 1 + cell_positions[i].astype(int)]).T
            # actual coordinate slices given tissue image
            coords = np.stack([np.maximum(0, bounds[:, 0]), np.array([np.minimum(a, b) for a, b in zip(np.array(img_mask.shape) - 1, bounds[:, 1])])], axis=1)
            # padded extract
            areas = np.pad(img_mask[tuple(map(lambda x: slice(*x), coords))], (bounds - coords) * np.array([-1, 1])[np.newaxis, :], mode='constant', constant_values=0)
            # area in square pixels
            areas = (areas[:, :, np.newaxis] & arcs).sum(axis=(0, 1))
            # append to df
            region_area_df = pd.concat([region_area_df, pd.DataFrame(data = {'idx': idx_i, 'area': areas})], axis=0)

        # Print the time taken to calculate the neighbor counts for the current image
        # print(f'  ...finished calculating neighborhood areas for image {region_id} (cells) in {time.time() - start_time:.2f} seconds')

        # return i and areas
        return region_area_df

    @staticmethod
    def process_cell_counts(i, cell_positions, cell_labels, dist_bin_px):
        # squared distance
        counts = np.sum(np.square(cell_positions[i][np.newaxis, :] - cell_positions), axis=1)
        # inequalities around arcs
        counts = counts[np.newaxis, :] <= np.square(np.concatenate([[0], dist_bin_px]))[:, np.newaxis]
        # matmul to counts
        counts = np.diff(np.matmul(counts.astype(int), cell_labels.astype(int)), axis=0)
        # return index and counts
        return counts

    def per_image_cell_counts_euc(self, image, cell_positions, cell_labels, targ_labels, dist_bin_px):
        '''
        per_image_cell_counts_euc() returns the number of cells within a given image

        Parameters:
            cell_positions (pd.DataFrame): DataFrame containing the cell positions
            cell_labels (np.array): labels of the cells
            targ_labels (np.array): labels of the cells to be counted
            dist_bin_px (np.array): distance bins in pixels
        '''

        start_time = time.time()
        print(f'Starting analysis for image {image}')
        # calculate pairwise distances between all cells in the image
        dist_st_time = time.time()
        distances = euclidean_distances(cell_positions)
        dist_end_time = (time.time() - dist_st_time) / 60
        print(f'Finished distance calculation for image {image} ({len(cell_positions)} cells) in {dist_end_time:.2f} minutes')

        image_counts = None
        for i in range(len(distances)):
            counts = self.euclidian_counts(i, distances, cell_labels, targ_labels, dist_bin_px)
            if image_counts is not None:
                image_counts = np.vstack((image_counts, counts))
            else:
                image_counts = counts

        comp_time = (time.time() - start_time) / 60
        print(f'Finished analysis for image {image} in {comp_time:.2f} minutes')
        return image_counts

    @staticmethod
    def euclidian_counts(idx, distances, cell_labels, targ_labels, dist_bin_px):
        '''
        euclidian_counts() returns the number of cells within a given 
        distance of a given cell.

        Parameters:
            idx (int): index of the cell to be counted
            distances (np.array): pairwise distances between cells
            cell_labels (np.array): labels of the cells
            targ_labels (np.array): labels of the cells to be counted
            dist_bin_px (np.array): distance bins in pixels
        '''

        idx_counts = None
        dist_bin_px = np.concatenate([[0], dist_bin_px])
        for i in range(len(dist_bin_px)-1):
            present_cells = cell_labels[(distances[idx] > dist_bin_px[i]) & (distances[idx] <= dist_bin_px[i+1])]
            these_counts = [sum(present_cells == label) for label in targ_labels]

            if idx_counts is not None:
                idx_counts = np.vstack((idx_counts, these_counts))
            else:
                idx_counts = np.array(these_counts)

        return idx_counts[np.newaxis, :]

    def calculate_density_matrix_for_all_images(self, cpu_pool_size = 8):
        """
        Calculate the density matrix for all images.

        Args:
            image_names (numpy.ndarray): The array of image names.
            df (pandas.DataFrame): The dataframe containing the data for all images.
            phenotypes (numpy.ndarray): The array of phenotypes.
            phenotype_column_name (str): The name of the column containing the phenotype information.
            image_column_name (str): The name of the column containing the image information.
            coord_column_names (list): The list of column names containing the coordinate information.
            radii (numpy.ndarray): The array of radii.
            range_strings (list): The list of range strings.
            debug_output (bool, optional): Whether to print debug output.
            num_cpus_to_use (int, optional): The number of CPUs to use. Defaults to 1.

        Returns:
            pandas.DataFrame: The dataframe containing the density matrix for all images.
        """

        df          = self.cells
        phenotypes  = self.species
        radii       = np.concatenate([[0], self.dist_bin_px])

        coord_column_names = ['Cell X Position', 'Cell Y Position']
        phenotype_column_name = 'Lineage'
        image_column_name     = 'Slide ID'
        image_names = df[image_column_name].unique()
        num_ranges = len(radii) - 1
        range_strings = [f'({radii[iradius]}, {radii[iradius + 1]}]' for iradius in range(num_ranges)]

        # Initialize keyword arguments
        kwargs_list = []

        # Loop through the images
        for image in image_names:

            # Create a dictionary for the variables
            kwargs_list.append(
                (
                    df[df[image_column_name] == image][[phenotype_column_name] + coord_column_names].copy(),
                    image,
                    coord_column_names,
                    phenotypes,
                    radii,
                    phenotype_column_name
                )
            )

        # Create a pool of worker processes
        with mp.Pool(processes=cpu_pool_size) as pool:
            results = pool.starmap(utils.fast_neighbors_counts_for_block2, kwargs_list)

        df_density_matrix = pd.concat(self.get_dataframes(results))
        full_array = None
        for ii, phenotype in enumerate(phenotypes):
            cols2Use = [f'{phenotype} in {x}' for x in range_strings]
            array_set = df_density_matrix.loc[:, cols2Use].to_numpy()
            if full_array is None:
                full_array = array_set
            else:
                full_array = np.dstack((full_array, array_set))

        full_array_nan = np.isnan(full_array)
        full_array[full_array_nan] = 0

        # Concatenate the results into a single dataframe
        return full_array

    def __init__(self, dist_bin_um, um_per_px, area_downsample):
        # microns per pixel
        self.um_per_px = um_per_px
        # distance arcs
        self.dist_bin_um = dist_bin_um
        # in pixels
        self.dist_bin_px = self.dist_bin_um / self.um_per_px
        # num of expected species. Reset this value in a top level script for your use-case
        self.num_species = 5
        # downsampling factor for area calculations
        self.area_downsample = area_downsample
        self.arcs_radii = (self.dist_bin_px * self.area_downsample).astype(int)
        # Create a boolean mask of concentric circles based on the distance bins in pixels
        self.arcs_masks = SpatialUMAP.construct_arcs(self.arcs_radii)

        # Attributes to be created in higher level script
        self.cells          = pd.DataFrame()
        self.cell_positions = pd.DataFrame()
        self.cell_labels    = pd.DataFrame()
        self.region_ids     = np.array([])
        self.pool = None
        self.species = None
        self.counts = None
        self.areas = None

        self.phenoLabel = None
        self.umap_test  = np.array([])
        self.patients   = np.array([])

        self.density    = None
        self.proportion = None

        # Mean Densities
        self.dens_df = pd.DataFrame()
        self.prop_df = pd.DataFrame()
        self.dens_df_mean = pd.DataFrame(data = {'clust_label': ['No Cluster'],
                                                 'phenotype': ['Other'], 
                                                 'dist_bin': [25], 
                                                 'density_mean': [0],
                                                 'density_sem': [0]})
        self.dens_df_se = pd.DataFrame()
        self.maxdens_df = pd.DataFrame()

        # UMAP Data for plotting
        self.df_umap = None

    def clear_counts(self):
        self.counts = np.empty((self.cell_positions.shape[0], len(self.dist_bin_um), self.num_species))

    def clear_areas(self):
        self.areas = np.empty((self.cell_positions.shape[0], len(self.dist_bin_um)))

    def start_pool(self, processes):
        start_method = mp.get_start_method()
        if start_method == 'fork':
            start_method = 'forkserver'  # to prevent crashing resulting in "Stopping..."
        self.pool = mp.get_context(start_method).Pool(processes)

    def close_pool(self):
        self.pool.close()
        self.pool.join()
        del self.pool

    def process_region_counts(self, region_id, pool_size):
        '''
        Process_region_counts
        '''
        # get indices of cells from this region
        idx = np.where(region_id == self.cells['TMA_core_id'])[0]
        # get counts if there are cells in region
        if len(idx) > 0:
            # partial for picklable fn for pool for process with data from this region
            args = dict(cell_positions=self.cell_positions[idx],
                        cell_labels=self.cell_labels.values[idx],
                        dist_bin_px=self.dist_bin_px)
            pool_map_fn = partial(SpatialUMAP.process_cell_counts, **args)
            chunk_size = 10000
            idxchunk = [idx[i:i + chunk_size] for i in range(0, len(idx), chunk_size)]
            # process
            with mp.Pool(pool_size) as pool:
                for i, chunk in enumerate(idxchunk):
                    chunk_range = range(0, len(chunk), 1)
                    chunk_range2 = [x + chunk_size*i for x in chunk_range]
                    results = list(pool.map(pool_map_fn, chunk_range2))
                    counts = list(map(lambda x: np.stack(x, axis=0), results))
                    # set results, adjust indexing (just in case)
                    self.counts[chunk] = counts

    def process_region_areas(self, pool_size):

        # Initialize keyword arguments
        kwargs_list = []

        for region_id in self.region_ids:
            # get indices of cells from this region
            idx = np.where(region_id == self.cells['TMA_core_id'])[0]
            # get counts if cells are in region
            if len(idx) > 0:
                # fit ellipse from point cloud
                fit_ellipse = FitEllipse()
                idx_fit = fit_ellipse.fit(self.cell_positions[idx][:, [1, 0]], px_to_hull=(100 / self.um_per_px))
                # extract binary mask
                img_tissue_mask = fit_ellipse.img_ellipse
                # down sample for area calculations
                img_tissue_mask_dn = sktran.rescale(img_tissue_mask, self.area_downsample).astype(bool)

                # Create a dictionary for the variables
                kwargs_list.append(
                    (
                        idx,
                        region_id,
                        self.cell_positions[idx][:, [1, 0]] * self.area_downsample,
                        self.arcs_radii,
                        img_tissue_mask_dn,
                        self.arcs_masks
                    )
                )

        # Create a pool of worker processes
        with mp.Pool(processes=pool_size) as pool:
            results = pool.starmap(self.process_cell_areas, kwargs_list)
            print(results)

    def get_counts(self, pool_size=2, save_file=None):
        '''
        get_counts begins the process of identifying the 
        cell counts surrounding each given cell in a 
        dataset
        '''
        self.clear_counts()
        # self.start_pool(pool_size)
        for region_id in tqdm(self.region_ids):
            self.process_region_counts(region_id, pool_size)
        # self.close_pool()

        if save_file is not None:
            column_names = ['%s-%s' % (cell_type, distance) for distance in self.dist_bin_um for cell_type in self.cell_labels.columns.values]
            pd.DataFrame(self.counts.reshape((self.counts.shape[0], -1)), columns=column_names).to_csv(save_file, index=False)

    def get_counts_euc(self, df, dist_bin_px, num_cpus_to_use=2):
        '''
        Another way to get counts, using euclidean distances
        '''
        # Initialize keyword arguments
        images = df['Slide ID'].unique()
        kwargs_list = []

        for image in images:

            df_image = df.loc[df['Slide ID'] == image, :]
            cell_positions = df_image[['Cell X Position', 'Cell Y Position']]
            cell_labels = df_image['Lineage']
            targ_labels = df['Lineage'].unique()
            dist_bin_px = dist_bin_px

            results = self.per_image_cell_counts_euc(image, cell_positions, cell_labels, targ_labels, dist_bin_px)
            print(results.shape)
            kwargs_list.append(
                (
                    image,
                    cell_positions,
                    cell_labels,
                    targ_labels,
                    dist_bin_px
                )
            )

        # Create a pool of worker processes
        # with mp.Pool(processes=num_cpus_to_use) as pool:
        #     results = pool.starmap(self.per_image_cell_counts_euc, kwargs_list)

        return results

    def get_counts_And(self, cpu_pool_size = 8):
        '''
        Andrew's method for getting counts
        '''
        self.counts = self.calculate_density_matrix_for_all_images(cpu_pool_size)

    def get_areas(self, calc_areas, area_threshold, pool_size=2, save_file=None, plots_directory=None):
        '''
        get_areas begins the process of identifying the
        cell areas surrounding each given cell in a dataset
        '''
        self.clear_areas()
        self.cells['area_filter'] = False

        if calc_areas:
            self.start_pool(pool_size)
            self.process_region_areas(pool_size, area_threshold=area_threshold, plots_directory=plots_directory)
            self.close_pool()
        else:
            areas = self.arcs_masks.sum(axis=(0, 1))[np.newaxis, ...]
            self.areas = np.tile(areas, (self.cells.shape[0], 1))

        if save_file is not None:
            pd.DataFrame(self.areas, columns=self.dist_bin_um).to_csv(save_file, index=False)

    def set_train_test(self, n_fit, n_tra, groupby_label = 'TMA_core_id', seed=None, umap_subset_toggle = False):
        '''
        set_test_train() is almost an unecessary method. Ultimately,
        when performing UMAP, we will intend to transform the whole dataset
        and we will never need to actually split the dataset into train and test.

        That said, this method will allow us to identify an even subset of the dataset
        to fit to be fitted to a model quickly, before the rest of the data is 
        transformed based on the UMAP model.

        Args:
            n (int): Minimum number of cells to be included in the training set
            groupby_label (str): Label to group the data by
            seed (int): Random seed for reproducibility

        Returns:
            None
        '''
        # region_ids is a proxy for the collection of images
        self.cells[['umap_train', 'umap_test']] = False

        for region_id, group in self.cells.groupby(groupby_label):
            if group['area_filter'].sum() >= (n_fit * 2):
                idx_train, idx_test, _ = np.split(np.random.default_rng(seed).permutation(group['area_filter'].sum()), [n_fit, n_tra])
                
                self.cells.loc[group.index[group.area_filter][idx_train], 'umap_train'] = True
                if umap_subset_toggle is False:
                    self.cells.loc[group.index[group.area_filter], 'umap_test'] = True
                else:
                    self.cells.loc[group.index[group.area_filter][idx_test], 'umap_test'] = True
        
        print(f'{np.sum(self.cells["umap_train"] == 1)} elements assigned to training data. ~{np.round(100*np.sum(self.cells["umap_train"] == 1)/self.cells.shape[0])}%')
        print(f'{np.sum(self.cells["umap_test"] == 1)} elements assigned to testing data. ~{np.round(100*np.sum(self.cells["umap_test"] == 1)/self.cells.shape[0])}%')

    def calc_densities(self, area_threshold):
        '''
        calculate density base on counts of cells / area of each arc

        Args:
            area_threshold (float): minimum area threshold for a cell to be considered

        Returns:
            None
        '''

        # instantiate our density output matrix
        self.density = np.empty(self.counts.shape)
        # identify those cells that do not have enough other cells around them. Any that
        # do not meet this criteria will be filtered out.
        self.cells['area_filter'] = ((self.areas / self.arcs_masks.sum(axis=(0, 1))[np.newaxis, ...]) > area_threshold).all(axis=1)

        # identify the indices of cells that are pass our filter
        filtIdx = (self.cells['area_filter'] == True)

        # Report the number of cells that pass the filter
        print(f'{np.sum(filtIdx)} cells pass the area filter out of {self.cells.shape[0]}')

        # calculate density (count/area) for filtered cells
        self.density[filtIdx] = self.counts[filtIdx] / self.areas[filtIdx][..., np.newaxis]

    def calc_proportions(self, area_threshold):
        '''
        calculate proportion base on counts of cells / total cells within an arc examine
        '''

        # instantiate our proportion output matrix
        self.proportion = np.empty(self.counts.shape)
        # identify those cells that do not have enough other cells around them. Any that
        # do not meet this criteria will be filtered out.
        self.cells['area_filter'] = ((self.areas / self.arcs_masks.sum(axis=(0, 1))[np.newaxis, ...]) > area_threshold).all(axis=1)

        # identify the indices of cells that are pass our filter
        filtIdx = (self.cells['area_filter'] == True)
        # calculate proportion (count/total_count) for filtered cells
        self.proportion[filtIdx] = self.counts[filtIdx] / self.counts[filtIdx].sum(axis = 2)[..., np.newaxis]

    def mean_measures(self):
        '''
        Setup density values for means
        '''

        dens_umap_test = self.density[self.cells['umap_test'], :, :]

        self.dens_df = pd.DataFrame()
        for clust_label, group in self.df_umap.groupby('clust_label'):

            clust_ind = group.index

            smalldf_D = pd.DataFrame()
            theseDen = dens_umap_test[clust_ind]
            for i, pheno in enumerate(self.phenoLabel):

                theseDen_pheno = theseDen[:,:,i]
                r, c = theseDen_pheno.shape
                theseDen_flat = theseDen_pheno.reshape(-1)

                smalldf_D['dist_bin'] = np.tile(self.dist_bin_um, r)
                smalldf_D['density'] = theseDen_flat
                smalldf_D['phenotype'] = pheno
                smalldf_D['clust_label'] = clust_label

                self.dens_df = pd.concat([self.dens_df, smalldf_D], axis = 0).reset_index(drop=True)

        # Perform Groupby and Mean calculations
        self.dens_df_mean = self.dens_df.groupby(['clust_label', 'phenotype', 'dist_bin'], as_index=False).mean()
        self.dens_df_se   = self.dens_df.groupby(['clust_label', 'phenotype', 'dist_bin'], as_index=False).sem()
        self.dens_df_mean = self.dens_df_mean.rename(columns = {'density': 'density_mean'})
        self.dens_df_se   = self.dens_df_se.rename(columns = {'density': 'density_sem'})
        self.dens_df_mean['density_sem'] = self.dens_df_se['density_sem']

        # Convert to mm^2
        self.dens_df_mean['density_mean'] = self.dens_df_mean['density_mean'] * 1e6
        self.dens_df_mean['density_sem'] = self.dens_df_mean['density_sem'] * 1e6

    def prepare_df_umap_plotting(self, features):
        '''
        Making a simple dataframe for plotting.
        In this case, feature are any and all features that are to be considered
        for plotting downstream of this event. 
        '''

        self.df_umap = pd.DataFrame(data = self.umap_test, columns = ['X', 'Y'])
        self.df_umap['Lineage'] = self.cells['Lineage'].values[self.cells['umap_test']]
        self.df_umap['species_name_short'] = self.cells['species_name_short'].values[self.cells['umap_test']]
        self.df_umap['Cluster'] = self.cells['clust_label'].values[self.cells['umap_test']]

        for feature in features:
            self.df_umap[feature] = self.cells[feature].values[self.cells['umap_test']]

    def makeDummyClinic(self, length):
        '''
        A method for quickly making a clinic dataset if needed 
        to pair with existitng Spatial UMAP methods
        length (int): number of dummy patients (samples) to create
        '''
        sample_number_pat = np.linspace(1, length, length)
        death_5y_pat = np.zeros(length)
        d = {'Sample_number': sample_number_pat, 'Death_5Y': death_5y_pat}
        return pd.DataFrame(data = d)

    def generate_H(self, lineages):
        '''
        Spatial UMAP 2D Density Plots By Lineage and Stratified by 5 Year Survival
        get per specimen density maps
        '''

        # set number of bins and get actual binning points based on whole dataset
        n_bins = 200
        xx = np.linspace(np.min(self.umap_test[:, 0]), np.max(self.umap_test[:, 0]), n_bins + 1)
        yy = np.linspace(np.min(self.umap_test[:, 1]), np.max(self.umap_test[:, 1]), n_bins + 1)
        # initialize holding nd matrix for densities
        n_lineages = len(lineages)
        # last dim is 0:counts, 1:smoothed, density
        H = np.empty([n_bins, n_bins, n_lineages + 1, len(self.patients['Sample_number']), 2])
        for i in range(len(self.patients['Sample_number'])):
            # get cells of this specimen / patient
            idx_pts = self.cells.loc[self.cells['umap_test'], 'Sample_number'] == self.patients['Sample_number'].iloc[i]
            if np.sum(idx_pts) > 0:
                # get counts for lineages
                for j in range(len(lineages)):
                    idx_lineage = self.cells.loc[self.cells['umap_test'], 'Lineage'].isin([lineages[j]])
                    H[:, :, j, i, 0], _, _ = np.histogram2d(self.umap_test[idx_pts & idx_lineage, 0],
                                                            self.umap_test[idx_pts & idx_lineage, 1], bins=[xx, yy])
                # get counts across all lineages
                H[:, :, j + 1, i, 0] = np.nansum(H[:, :, 0:(j + 1), i, 0], axis=2)

                # make smoothed density for lineages
                for j in range(len(lineages)):
                    if np.sum(H[:, :, j, i, 0]) > 0:
                        H[:, :, j, i, 1] = ndi.gaussian_filter(H[:, :, j, i, 0] / np.sum(H[:, :, j, i, 0]), sigma=0.5)
                    else:
                        H[:, :, j, i, 1] = np.nan
                # make smoothed density for all lineages
                if np.sum(H[:, :, j + 1, i, 0]) > 0:
                    H[:, :, j + 1, i, 1] = ndi.gaussian_filter(H[:, :, j + 1, i, 0] / np.sum(H[:, :, j + 1, i, 0]), sigma=0.5)
                else:
                    H[:, :, j + 1, i, 1] = np.nan
            else:
                H[:, :, :, i, :] = np.nan

        return H

class FitEllipse:
    '''FitElipse is a class that defines the concentric circle areas
    used to determine density of cell counts. It has five methods:
    * __init__
    * elipse_function
    * elipse_area
    * draw_ellipse
    * fit

    Parameters:

    Returns:
    class object
    '''
    def __init__(self):
        self.x = None
        self.img_ellipse = None

        self.w = None
        self.h = None
        self.res = None

    @staticmethod
    def ellipse_function(points, x, y, a, b, r):
        t = np.array([np.cos(r), np.sin(r)])
        d = points - np.array([x, y])[np.newaxis, ...]
        return np.square(((t[0] * d[:, 0]) + (t[1] * d[:, 1])) / a) + np.square(((t[1] * d[:, 0]) - (t[0] * d[:, 1])) / b)

    @staticmethod
    def ellipse_area(a, b):
        return np.pi * a * b

    def draw_ellipse(self, x=None):
        assert self.img_ellipse is not None
        _x = x if x is not None else self.x
        xx, yy = skdraw.ellipse(_x[0], _x[1], _x[2], _x[3], self.img_ellipse.shape, _x[4])
        self.img_ellipse[:] = False
        self.img_ellipse[xx, yy] = True

    def fit(self, d, px_to_hull):
        idx_fit = np.ones(d.shape[0], dtype=bool)
        idx_remove = True
        while np.any(idx_remove):
            hull = ConvexHull(d[idx_fit])
            d_h = np.sum(np.square(d[idx_fit][:, np.newaxis, :] - d[idx_fit][hull.vertices][np.newaxis, :, :]), axis=-1)
            idx_remove = np.sum(d_h < np.square(px_to_hull), axis=0) < 5
            idx_fit[np.where(idx_fit)[0][hull.vertices[idx_remove]]] = False
        idx_fit = np.where(idx_fit)[0][np.unique(np.argsort(d_h, axis=0)[:50])]

        self.w, self.h = np.max(d, axis=0).astype(int)
        x_init = np.concatenate([np.array(np.array((self.w, self.h))) / 2, np.log(np.array((self.w, self.h))), [0, ]]).astype(float)
        self.res = optimize.minimize(lambda x: np.mean(np.abs(FitEllipse.ellipse_function(d[idx_fit], x[0], x[1], np.exp(x[2]), np.exp(x[3]), x[4]) - 1)), x_init, method='nelder-mead')
        self.x = self.res.x.copy()
        self.x[2], self.x[3] = np.exp(self.x[[2, 3]])

        self.img_ellipse = np.zeros((self.w, self.h), dtype=bool)
        self.draw_ellipse()

        return idx_fit
