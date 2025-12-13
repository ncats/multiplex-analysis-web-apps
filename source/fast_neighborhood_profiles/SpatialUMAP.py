import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from skimage import draw as skdraw, transform as sktran
import pickle
from time import time
import umap
from multiprocessing import Pool
from functools import partial
from scipy import optimize
import multiprocessing as mp
import fast_neighborhood_profiles.main as fnp_main


import scipy.spatial
def fast_neighbors_counts_for_block2(df_image, image_name, coord_column_names, phenotypes, radii, phenotype_column_name, max_chunk_size_in_mb=200, data_struct="pandas", kdtree_str="kdtree"):
    # A block can be an image, ROI, etc. It's the entity over which it makes sense to calculate the neighbors of centers. Here, we're assuming it's an image, but in the SIT for e.g., we generally want it to refer to a ROI.

    # max_chunk_size_in_mb=200, for a 100K-cell dataset, will yield about 250-row chunks, which will yield about 400 chunks i.e. center KDTrees
    # max_chunk_size_in_mb=500, for a 100K-cell dataset, will yield about 650-row chunks, which will yield about 150 chunks i.e. center KDTrees
    # max_chunk_size_in_mb=1000, for a 100K-cell dataset, will yield about 1300-row chunks, which will yield about 80 chunks i.e. center KDTrees
    # max_chunk_size_in_mb=2000, for a 100K-cell dataset, will yield about 2600-row chunks, which will yield about 40 chunks i.e. center KDTrees
    # max_chunk_size_in_mb=5000, for a 100K-cell dataset, will yield about 6600-row chunks, which will yield about 15 chunks i.e. center KDTrees

    if kdtree_str == "kdtree":
        kdtree_func = scipy.spatial.KDTree
    elif kdtree_str == "ckdtree":
        kdtree_func = scipy.spatial.cKDTree

    # Print the image name
    print(f'Calculating neighbor counts for image {image_name} ({len(df_image)} cells) using the kdtree method {kdtree_str} with a chunk size of {max_chunk_size_in_mb}MB and data structure type {data_struct}...', flush=True)

    # Record the start time
    # start_time = time.time()
    start_time = time()
    elapsed_time_dataframe = 0
    elapsed_time_kdtree = 0

    # Store some properties of the input dataframe
    df_image_index = df_image.index
    num_cells = len(df_image)

    # Get the number of rows per chunk based on the maximum chunk size in MB and assuming every cell will be counted as a neighbor around every center (i.e., use upper limits)
    largest_sublist_size_in_mb = num_cells * 8  / 1024 ** 2  # 8 bytes per int64 --> units = mb / row
    num_rows_per_chunk = max(1, int(max_chunk_size_in_mb / largest_sublist_size_in_mb))

    # Get the corresponding integer indices to index things like the input image dataframe
    start_indices = np.arange(0, num_cells, num_rows_per_chunk, dtype=np.int64)
    stop_indices = start_indices + num_rows_per_chunk

    # Initialize a list to hold the dataframes of neighbor counts for each radius (not each radius range)
    time_before = time()
    if data_struct == "pandas":
        df_counts_holder = [pd.DataFrame(0, index=phenotypes, columns=df_image_index) for _ in radii]
    elif data_struct == "numpy":
        df_counts_holder = [np.zeros((len(phenotypes), len(df_image_index)), dtype=np.int32) for _ in radii]
    elapsed_time_dataframe += time() - time_before

    print(f"data struct: {elapsed_time_dataframe:.2f} seconds, kdtree: {elapsed_time_kdtree:.2f} seconds", flush=True)
    print(f"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA {image_name}", flush=True)
    elapsed_time_dataframe = 0
    elapsed_time_kdtree = 0

    # Convert dataframe columns to numpy arrays for faster access if using numpy data structure.
    if data_struct == "numpy":
        nd_phenotype = df_image[phenotype_column_name].to_numpy()
        nd_coords = df_image[coord_column_names].to_numpy(dtype=np.float64)

    # Pre-calculate the neighbor tree for each phenotype
    neighbor_trees = []
    for neighbor_phenotype in phenotypes:

        # Get the boolean series identifying the current neighbor phenotype
        time_before = time()
        if data_struct == "pandas":
            ser_curr_neighbor_phenotype = df_image[phenotype_column_name] == neighbor_phenotype
        elif data_struct == "numpy":
            ser_curr_neighbor_phenotype = nd_phenotype == neighbor_phenotype
            # nd_curr_neighbor_phenotype = np.where(ser_curr_neighbor_phenotype)[0]
        elapsed_time_dataframe += time() - time_before

        # Construct the KDTree for the current phenotype in the entire current image. This represents the neighbors
        time_before = time()
        if data_struct == "pandas":
            neighbor_trees.append(kdtree_func(df_image.loc[ser_curr_neighbor_phenotype, coord_column_names]))
        elif data_struct == "numpy":
            neighbor_trees.append(kdtree_func(nd_coords[ser_curr_neighbor_phenotype, :]))
        elapsed_time_kdtree += time() - time_before

    print(f"data struct: {elapsed_time_dataframe:.2f} seconds, kdtree: {elapsed_time_kdtree:.2f} seconds", flush=True)
    print(f"BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB {image_name}", flush=True)
    elapsed_time_dataframe = 0
    elapsed_time_kdtree = 0
    
    # For each chunk of centers...
    for start_index, stop_index in zip(start_indices, stop_indices):

        # Construct the KDTree for the current chunk. This represents the centers
        time_before = time()
        if data_struct == "pandas":
            center_tree = kdtree_func(df_image.loc[df_image_index[start_index:stop_index], coord_column_names])
        elif data_struct == "numpy":
            center_tree = kdtree_func(nd_coords[start_index:stop_index, :])
        elapsed_time_kdtree += time() - time_before

        # For each neighbor tree...
        for ineighbor_phenotype, curr_neighbor_tree in enumerate(neighbor_trees):

            # For each radius, which should be monotonically increasing and start with 0...
            for iradius, radius in enumerate(radii):

                # Get the list of lists containing the indices of the neighbors for each center
                time_before = time()
                neighbors_for_radius = center_tree.query_ball_tree(curr_neighbor_tree, radius)
                elapsed_time_kdtree += time() - time_before

                # In the correct dataframe (corresponding to the current radius), set the counts of neighbors (of the current phenotype) for each center
                time_before = time()
                if data_struct == "pandas":
                    df_counts_holder[iradius].iloc[ineighbor_phenotype, start_index:stop_index] = [len(neighbors_for_center) for neighbors_for_center in neighbors_for_radius]
                elif data_struct == "numpy":
                    df_counts_holder[iradius][ineighbor_phenotype, start_index:stop_index] = [len(neighbors_for_center) for neighbors_for_center in neighbors_for_radius]
                    df_counts_holder[iradius][ineighbor_phenotype, start_index:stop_index] = np.fromiter(map(len, neighbors_for_radius), dtype=np.int32)
                elapsed_time_dataframe += time() - time_before

    print(f"data struct: {elapsed_time_dataframe:.2f} seconds, kdtree: {elapsed_time_kdtree:.2f} seconds", flush=True)
    print(f"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC {image_name}", flush=True)
    elapsed_time_dataframe = 0
    elapsed_time_kdtree = 0
    
    # For each annulus, i.e., each radius range...
    df_counts_holder_annulus = []
    if data_struct == "numpy":
        df_counts_annulus_column_holder = []  # "column" not "index" since we transpose
    for iradius in range(len(radii) - 1):

        # Get the counts of neighbors in the current annulus
        time_before = time()
        df_counts_curr_annulus = df_counts_holder[iradius + 1] - df_counts_holder[iradius]
        elapsed_time_dataframe += time() - time_before

        # Rename the index to reflect the current radius range
        radius_range_str = f'({radii[iradius]}, {radii[iradius + 1]}]'
        time_before = time()
        if data_struct == "pandas":
            df_counts_curr_annulus.index = [f'{phenotype} in {radius_range_str}' for phenotype in phenotypes]
        elif data_struct == "numpy":
            df_counts_annulus_column_holder.append([f'{phenotype} in {radius_range_str}' for phenotype in phenotypes])
        elapsed_time_dataframe += time() - time_before

        # Add a transpose of this (so centers are in rows and phenotypes/radii are in columns) to the running list of annulus dataframes
        time_before = time()
        df_counts_holder_annulus.append(df_counts_curr_annulus.T)
        elapsed_time_dataframe += time() - time_before

    print(f"data struct: {elapsed_time_dataframe:.2f} seconds, kdtree: {elapsed_time_kdtree:.2f} seconds", flush=True)
    print(f"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDD {image_name}", flush=True)
    elapsed_time_dataframe = 0
    elapsed_time_kdtree = 0
    
    # Concatenate the annulus dataframes to get the final dataframe of neighbor counts for the current image
    time_before = time()
    if data_struct == "pandas":
        df_curr_counts = pd.concat(df_counts_holder_annulus, axis='columns')
    elif data_struct == "numpy":
        df_curr_counts = np.concatenate(df_counts_holder_annulus, axis=1, dtype=np.int32)
    elapsed_time_dataframe += time() - time_before

    # Convert the dataframe to a more efficient dtype (from int64)
    time_before = time()
    if data_struct == "pandas":
        df_curr_counts = df_curr_counts.astype(np.int32)
    elif data_struct == "numpy":
        df_curr_counts = pd.DataFrame(df_curr_counts, index=df_image_index, columns=[col for sublist in df_counts_annulus_column_holder for col in sublist])
    elapsed_time_dataframe += time() - time_before

    print(f"data struct: {elapsed_time_dataframe:.2f} seconds, kdtree: {elapsed_time_kdtree:.2f} seconds", flush=True)
    
    # Print the time taken to calculate the neighbor counts for the current image
    # print(f'  ...finished calculating neighbor counts for image {image_name} ({len(df_image)} cells) in {time.time() - start_time:.2f} seconds', flush=True)
    print(f'  ...finished calculating neighbor counts for image {image_name} ({len(df_image)} cells) in {time() - start_time:.2f} seconds', flush=True)

    # Return the final dataframe of neighbor counts for the current image
    return df_curr_counts


class SpatialUMAP:

    @staticmethod
    def get_dataframes(results):
        for df in results:
            yield df

    @staticmethod
    def construct_arcs(dist_bin_px):
        # set bool mask of the arcs
        arcs = np.zeros([int(2 * dist_bin_px[-1]) + 1] * 2 + [len(dist_bin_px), ], dtype=bool)
        for i in range(len(dist_bin_px)):
            # circle based on radius
            rr, cc = skdraw.disk(center=(np.array(arcs.shape[:2]) - 1) / 2, radius=dist_bin_px[i] + 1, shape=arcs.shape[:2])
            arcs[rr, cc, i] = True
        # difference logic to produce arcs
        return np.stack([arcs[:, :, 0]] + [arcs[:, :, i] != arcs[:, :, i - 1] for i in range(1, arcs.shape[2])], axis=2)

    @staticmethod
    def process_cell_areas(i, cell_positions, cell_labels, dist_bin_px, img_mask, arcs):
        # true bounds to match arcs
        bounds = np.array([cell_positions[i].astype(int) - dist_bin_px[-1].astype(int), dist_bin_px[-1].astype(int) + 1 + cell_positions[i].astype(int)]).T
        # actual coordinate slices given tissue image
        coords = np.stack([np.maximum(0, bounds[:, 0]), np.array([np.minimum(a, b) for a, b in zip(np.array(img_mask.shape) - 1, bounds[:, 1])])], axis=1)
        # padded extract
        areas = np.pad(img_mask[tuple(map(lambda x: slice(*x), coords))], (bounds - coords) * np.array([-1, 1])[np.newaxis, :], mode='constant', constant_values=0)
        # area in square pixels
        areas = (areas[:, :, np.newaxis] & arcs).sum(axis=(0, 1))
        # return i and areas
        return i, areas

    @staticmethod
    def process_cell_counts(i, cell_positions, cell_labels, dist_bin_px):
        # squared distance
        counts = np.sum(np.square(cell_positions[i][np.newaxis, :] - cell_positions), axis=1)
        # inequalities around arcs
        counts = counts[np.newaxis, :] <= np.square(np.concatenate([[0], dist_bin_px]))[:, np.newaxis]
        # matmul to counts
        counts = np.diff(np.matmul(counts.astype(int), cell_labels.astype(int)), axis=0)
        # return index and counts
        return i, counts

    def __init__(self, dist_bin_um, um_per_px, area_downsample):
        # microns per pixel
        self.um_per_px = um_per_px
        # distance arcs
        self.dist_bin_um = dist_bin_um
        # in pixels
        self.dist_bin_px = self.dist_bin_um / self.um_per_px
        # downsampling factor for area calculations
        self.area_downsample = area_downsample
        self.arcs_radii = (self.dist_bin_px * self.area_downsample).astype(int)
        self.arcs_masks = SpatialUMAP.construct_arcs(self.arcs_radii)

    def clear_counts(self):
        self.counts = np.empty((self.cell_positions.shape[0], len(self.dist_bin_um), self.cell_labels.shape[1]))

    def clear_areas(self):
        self.areas = np.empty((self.cell_positions.shape[0], len(self.dist_bin_um)))

    def start_pool(self, processes, mp_start_method="forkserver"):
        self.pool = mp.get_context(mp_start_method).Pool(processes=processes)

    def close_pool(self):
        self.pool.close()
        del self.pool

    def process_region_counts(self, region_id):
        # get indices of cells from this region
        idx = np.where(region_id == self.cells['TMA_core_id'])[0]
        # get counts if there are cells in region
        if len(idx) > 0:
            # partial for picklable fn for pool for process with data from this region
            args = dict(cell_positions=self.cell_positions[idx], cell_labels=self.cell_labels.values[idx], dist_bin_px=self.dist_bin_px)
            pool_map_fn = partial(SpatialUMAP.process_cell_counts, **args)
            # process
            i, counts = list(map(lambda x: np.stack(x, axis=0), list(zip(*self.pool.map(pool_map_fn, range(len(idx)))))))
            # set results, adjust indexing (just in case)
            self.counts[idx] = counts[i]

    def process_region_areas(self, region_id, area_threshold, plots_directory=None):
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

            # partial for picklable fn for pool for process with data from this region
            args = dict(cell_positions=self.cell_positions[idx][:, [1, 0]] * self.area_downsample, cell_labels=self.cell_labels.values[idx], dist_bin_px=self.arcs_radii, img_mask=img_tissue_mask_dn, arcs=self.arcs_masks)
            pool_map_fn = partial(SpatialUMAP.process_cell_areas, **args)
            # process
            i, areas = list(map(lambda x: np.stack(x, axis=0), list(zip(*self.pool.map(pool_map_fn, range(len(idx)))))))
            # adjust for indexing (just in case)
            areas = areas[i]
            # set filter for cells with adequate area coverage
            filt = ((areas / self.arcs_masks.sum(axis=(0, 1))[np.newaxis, ...]) > area_threshold).all(axis=1)

            # set results
            self.areas[idx] = areas
            self.cells.loc[idx, 'area_filter'] = filt

            if plots_directory is not None:
                plt.ioff()
                f = plt.figure(figsize=(3, 3))
                plt.axes()
                f.axes[0].cla()
                f.axes[0].plot(*self.cell_positions[idx].T, 'k,')
                f.axes[0].plot(*self.cell_positions[idx][idx_fit].T, 'r.', markersize=3, alpha=0.5)
                f.axes[0].plot(*self.cell_positions[idx][filt].T, 'b.', markersize=3, alpha=0.5)
                f.axes[0].imshow(img_tissue_mask, alpha=0.5)
                f.axes[0].axis('off')
                plt.tight_layout(pad=0.1)
                f.savefig('%s/%s.png' % (plots_directory, region_id), format='png')
                plt.close(f)
                del f
                plt.ion()

    def get_counts(self, cpu_pool_size=2, save_file=None):
        self.clear_counts()
        self.start_pool(cpu_pool_size)
        for region_id in tqdm(self.region_ids):
            self.process_region_counts(region_id)
        self.close_pool()

        # self.counts_before_save = self.counts.copy()

        if save_file is not None:
            column_names = ['%s-%s' % (cell_type, distance) for distance in self.dist_bin_um for cell_type in self.cell_labels.columns.values]
            pd.DataFrame(self.counts.reshape((self.counts.shape[0], -1)), columns=column_names).to_csv(save_file, index=False)
            self.counts = pd.read_csv(save_file, sep=',').values.reshape((self.counts.shape[0], self.dist_bin_um.shape[0], self.cell_labels.shape[1]))
            # self.counts_after_load = self.counts.copy()

    def get_counts_And(self, cpu_pool_size = 8, save_file=None, mp_start_method=None):
        '''
        Andrew's method for getting counts
        '''
        self.counts = self.calculate_density_matrix_for_all_images(cpu_pool_size, mp_start_method=mp_start_method)

        # self.counts_before_save = self.counts.copy()

        # May not necessarily be the same label order, though actually probably is (pd.dummies() vs. lazyframe sorting). So commenting out to not imply they're the same.
        # if save_file is not None:
        #     column_names = ['%s-%s' % (cell_type, distance) for distance in self.dist_bin_um for cell_type in self.cell_labels.columns.values]
        #     pd.DataFrame(self.counts.reshape((self.counts.shape[0], -1)), columns=column_names).to_csv(save_file, index=False)
        #     self.counts = pd.read_csv(save_file, sep=',').values.reshape((self.counts.shape[0], self.dist_bin_um.shape[0], self.cell_labels.shape[1]))
        #     # self.counts_after_load = self.counts.copy()

    def calculate_density_matrix_for_all_images(self, cpu_pool_size = 8, mp_start_method=None):
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

        if mp_start_method is None:
            mp_start_method = mp.get_start_method()
        if mp_start_method == 'fork':
            mp_start_method = 'forkserver'
            print(f'Note: We are forcing the multiprocessing module to use the "forkserver" start method instead of the automatically (or manually) chosen "fork" start method.', flush=True)

        df          = self.cells
        phenotypes  = self.species
        radii       = np.concatenate([[0], self.dist_bin_px])

        coord_column_names = ['Xcor', 'Ycor']
        phenotype_column_name = 'Lineage'
        image_column_name     = 'TMA_core_id'
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
                    phenotype_column_name,
                    200,
                    "numpy",
                    "kdtree",
                )
            )

        # Being very explicit with errors because with forkserver/spawn I sometimes get missing output from workers. Actually confirm all workers complete using the following catches.
        print(f"Using start method {mp_start_method} with {cpu_pool_size} CPUs.", flush=True)
        try:
            with mp.get_context(mp_start_method).Pool(processes=cpu_pool_size) as pool:
                # results = pool.starmap(fnp_main.fast_neighbors_counts_for_block2, kwargs_list)
                results = pool.starmap(fast_neighbors_counts_for_block2, kwargs_list)
        except Exception as e:
            # surface hard failures clearly in Streamlit
            raise RuntimeError(f"Parallel run failed: {e}") from e

        assert len(results) == len(kwargs_list), "Got fewer results than tasks."

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

    def get_areas(self, area_threshold, pool_size=2, save_file=None, plots_directory=None):
        self.clear_areas()
        self.cells['area_filter'] = False
        self.start_pool(pool_size)
        for region_id in tqdm(self.region_ids):
            self.process_region_areas(region_id, area_threshold=area_threshold, plots_directory=plots_directory)
        self.close_pool()

        if save_file is not None:
            pd.DataFrame(self.areas, columns=self.dist_bin_um).to_csv(save_file, index=False)

    def set_train_test(self, n, seed=None, train_sample_frac=1.0, test_sample_frac=1.0):
        group_col = 'Sample_number' if 'Sample_number' in self.cells.columns else 'TMA_core_id'
        if group_col == "Sample_number":
            print("NOTE: Using 'Sample_number' to group cells for train/test split, which is the Baras default but is not originally defined in his code. We'd otherwise naturally use 'TMA_core_id'.")
        # region_ids = self.cells['TMA_core_id'].unique()
        self.cells[['umap_train', 'umap_test']] = False
        for region_id, group in self.cells.groupby(group_col):
            if group['area_filter'].sum() >= int((train_sample_frac+test_sample_frac)*n):  # used to be 2n
                idx_train, idx_test, _ = np.split(np.random.default_rng(seed).permutation(group['area_filter'].sum()), [int(train_sample_frac*n), int((train_sample_frac+test_sample_frac)*n)])
                self.cells.loc[group.index[group.area_filter][idx_train], 'umap_train'] = True
                self.cells.loc[group.index[group.area_filter][idx_test], 'umap_test'] = True


class FitEllipse:
    def __init__(self):
        self.x = None
        self.img_ellipse = None

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
