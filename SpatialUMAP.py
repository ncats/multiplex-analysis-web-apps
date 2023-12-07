# Author: Alex Baras, MD, PhD (https://github.com/alexbaras)
# NCATS Maintainer: Dante J Smith, PhD (https://github.com/djsmith17)
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from skimage import draw as skdraw, transform as sktran
from time import time
import multiprocessing as mp
from functools import partial
from scipy import optimize
from scipy import ndimage as ndi
np.seterr(divide='ignore', invalid='ignore')

class SpatialUMAP:
    '''SpatialUMAP is a class that handles the cell count 
    processing and concentric circle area measurements for 
    cell density analysis in SpatialUMAP studies.

    It has the following methods:
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
        # num of expected species. Reset this value in a top level script for your use-case
        self.num_species = 5
        # downsampling factor for area calculations
        self.area_downsample = area_downsample
        self.arcs_radii = (self.dist_bin_px * self.area_downsample).astype(int)
        self.arcs_masks = SpatialUMAP.construct_arcs(self.arcs_radii)

    def clear_counts(self):
        self.counts = np.empty((self.cell_positions.shape[0], len(self.dist_bin_um), self.num_species))

    def clear_areas(self):
        self.areas = np.empty((self.cell_positions.shape[0], len(self.dist_bin_um)))

    def start_pool(self, processes):
        self.pool = mp.get_context(mp.get_start_method()).Pool(processes)

    def close_pool(self):
        self.pool.close()
        self.pool.join()
        del self.pool

    def process_region_counts(self, region_id, pool_size):
        # get indices of cells from this region
        idx = np.where(region_id == self.cells['TMA_core_id'])[0]
        # get counts if there are cells in region
        if len(idx) > 0:
            # partial for picklable fn for pool for process with data from this region
            args = dict(cell_positions=self.cell_positions[idx], 
                        cell_labels=self.cell_labels.values[idx], 
                        dist_bin_px=self.dist_bin_px)
            pool_map_fn = partial(SpatialUMAP.process_cell_counts, **args)
            # process
            i, counts = list(map(lambda x: np.stack(x, axis=0), list(zip(*self.pool.map(pool_map_fn, range(len(idx)))))))
            # set results, adjust indexing (just in case)
            self.counts[idx] = counts[i]

    def process_region_areas(self, region_id, pool_size, area_threshold, plots_directory=None):
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
            args = dict(cell_positions=self.cell_positions[idx][:, [1, 0]] * self.area_downsample, 
                        cell_labels=self.cell_labels.values[idx], 
                        dist_bin_px=self.arcs_radii, img_mask=img_tissue_mask_dn, arcs=self.arcs_masks)
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

    def get_counts(self, pool_size=2, save_file=None):
        self.clear_counts()
        self.start_pool(pool_size)
        for region_id in tqdm(self.region_ids):
            self.process_region_counts(region_id, pool_size)
        self.close_pool()

        if save_file is not None:
            column_names = ['%s-%s' % (cell_type, distance) for distance in self.dist_bin_um for cell_type in self.cell_labels.columns.values]
            pd.DataFrame(self.counts.reshape((self.counts.shape[0], -1)), columns=column_names).to_csv(save_file, index=False)

    def get_areas(self, area_threshold, pool_size=2, save_file=None, plots_directory=None):
        self.clear_areas()
        self.cells['area_filter'] = False
        self.start_pool(pool_size)
        for region_id in tqdm(self.region_ids):
            self.process_region_areas(region_id, pool_size, area_threshold=area_threshold, plots_directory=plots_directory)
        self.close_pool()

        if save_file is not None:
            pd.DataFrame(self.areas, columns=self.dist_bin_um).to_csv(save_file, index=False)

    def set_train_test(self, n, groupByLabel = 'Sample_number', seed=None):
        region_ids = self.cells['TMA_core_id'].unique()
        self.cells[['umap_train', 'umap_test']] = False

        # How many samples do we have in the dataset?
        numSamp = self.cells.groupby(groupByLabel).count().shape[0]
        numSampInc = 0
        for region_id, group in self.cells.groupby(groupByLabel):
            if group['area_filter'].sum() >= (n * 2):
                numSampInc +=1
                idx_train, idx_test, _ = np.split(np.random.default_rng(seed).permutation(group['area_filter'].sum()), [n, n * 2])
                self.cells.loc[group.index[group.area_filter][idx_train], 'umap_train'] = True
                self.cells.loc[group.index[group.area_filter], 'umap_test'] = True
        
        print(f'{numSampInc} Samples used of {numSamp} Samples in dataset')
        print(f'{np.sum(self.cells["umap_train"] == 1)} elements assigned to training data. ~{np.round(100*np.sum(self.cells["umap_train"] == 1)/self.cells.shape[0])}%')
        print(f'{np.sum(self.cells["umap_test"] == 1)} elements assigned to testing data. ~{np.round(100*np.sum(self.cells["umap_test"] == 1)/self.cells.shape[0])}%')
    
    def calc_densities(self, area_threshold):
        '''
        calculate density base on counts of cells / area of each arc examine
        '''

        # instantiate our density output matrix
        self.density = np.empty(self.counts.shape)
        # identify those cells that do not have enough other cells around them. Any that
        # do not meet this criteria will be filtered out. 
        self.cells['area_filter'] = ((self.areas / self.arcs_masks.sum(axis=(0, 1))[np.newaxis, ...]) > area_threshold).all(axis=1)

        # identify the indices of cells that are pass our filter
        filtIdx = (self.cells['area_filter'] == True)
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

        self.dens_df = pd.DataFrame()
        self.prop_df = pd.DataFrame()
        for clust_label, group in self.cells.groupby('clust_label'):
            
            if clust_label != -1:
                ind = group.index

                smalldf_D = pd.DataFrame()
                smalldf_P = pd.DataFrame()
                theseDen = self.density[ind]
                thesePro = self.proportion[ind]
                for i, pheno in enumerate(self.phenoLabel):

                    theseDen_pheno = theseDen[:,:,i]
                    r, c = theseDen_pheno.shape
                    theseDen_flat = theseDen_pheno.reshape(-1)

                    thesePro_pheno = thesePro[:,:,i]
                    r, c = thesePro_pheno.shape
                    thesePro_flat = thesePro_pheno.reshape(-1)
                    
                    smalldf_D['dist_bin'] = np.tile(self.dist_bin_um, r)
                    smalldf_D['density'] = theseDen_flat
                    smalldf_D['phenotype'] = pheno
                    smalldf_D['cluster'] = clust_label

                    smalldf_P['dist_bin'] = np.tile(self.dist_bin_um, r)
                    smalldf_P['density'] = thesePro_flat
                    smalldf_P['phenotype'] = pheno
                    smalldf_P['cluster'] = clust_label

                    self.dens_df = pd.concat([self.dens_df, smalldf_D], 0).reset_index(drop=True)
                    self.prop_df = pd.concat([self.prop_df, smalldf_P], 0).reset_index(drop=True)

        self.dens_df_mean = self.dens_df.groupby(['cluster', 'phenotype', 'dist_bin'], as_index=False).mean()
        self.dens_df_se   = self.dens_df.groupby(['cluster', 'phenotype', 'dist_bin'], as_index=False).sem()
        self.maxdens_df   = 1.05*max(self.dens_df_mean['density'] + self.dens_df_se['density'])
    
    def makeDummyClinic(self, length):
        # A method for quickly making a clinic dataset if needed to pair with existitng Spatial UMAP methods
        # length (int): number of dummy patients (samples) to create
        Sample_numberPat = np.linspace(1, length, length)
        Death_5YPat = np.zeros(length)
        d = {'Sample_number': Sample_numberPat, 'Death_5Y': Death_5YPat}
        return pd.DataFrame(data = d)

    def generate_H(self, lineages):
        ## Spatial UMAP 2D Density Plots By Lineage and Stratified by 5 Year Survival
        # get per specimen density maps

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

    
    def define_phenotype_ordering():
        pass

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
