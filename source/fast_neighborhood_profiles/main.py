import polars as pl
import os
import plotly.express as px
import numpy as np
import pandas as pd
from fast_neighborhood_profiles import SpatialUMAP
import umap
import PlottingTools
import matplotlib.pyplot as plt


# def get_unique_column_values(csv_filename="mawa-unified_datafile-TLS_tissue_SF_-20251112_130129_EST.csv", column_name="Image ID_(standardized)", topdir="."):
#     try:
#         csv_filepath = os.path.join(topdir, "datafiles", csv_filename)
#         lf = pl.scan_csv(csv_filepath)
#         return lf.select(pl.col(column_name).unique()).collect().to_series().to_list()
#     except Exception as e:
#         print(f"An error occurred in function {os.path.basename(__file__)}.{get_unique_column_values.__name__}: {e}")
#         return []


def print_flush(msg):
    print(msg, flush=True)


def get_min_positive_values(pd_df, group_col="TMA_core_id", boolean_column="area_filter"):
    return pd_df.groupby(group_col)[boolean_column].sum().min()


def subset_csv_to_file(csv_filename="mawa-unified_datafile-TLS_tissue_SF_-20251112_130129_EST.csv", handle="two_images", do_filtering=False, filter_column="Image ID_(standardized)", filter_values=["MS_01__cele_1400w", "MS_02__cele_1400w"], topdir=".", subdir="datafiles", file_format="parquet"):
    try:    
        csv_filepath = os.path.join(topdir, subdir, csv_filename)
        filepath = os.path.join(topdir, subdir, handle + "." + file_format)
        lf = pl.scan_csv(csv_filepath)
        if file_format == "parquet":
            write_method = "write_parquet"
        elif file_format == "arrow":
            write_method = "write_ipc"
        elif file_format == "csv":
            write_method = "write_csv"
        else:
            raise ValueError(f"Unsupported file format: {file_format}")
        if do_filtering:
            getattr(lf.filter(pl.col(filter_column).is_in(filter_values)).collect(), write_method)(filepath)
        else:
            getattr(lf.collect(), write_method)(filepath)
        return filepath
    except Exception as e:
        print_flush(f"An error occurred in function {os.path.basename(__file__)}.{subset_csv_to_file.__name__}: {e}")
        return ""


def save_pandas_df_to_file(pd_df, handle="two_images", topdir=".", file_format="parquet", subdir="datafiles"):
    try:    
        filepath = os.path.join(topdir, subdir, handle + "." + file_format)
        pl_df = pl.from_pandas(pd_df)
        if file_format == "parquet":
            pl_df.write_parquet(filepath)
        elif file_format == "arrow":
            pl_df.write_ipc(filepath)
        elif file_format == "csv":
            pl_df.write_csv(filepath)
        else:
            raise ValueError(f"Unsupported file format: {file_format}")
        return filepath
    except Exception as e:
        print_flush(f"An error occurred in function {os.path.basename(__file__)}.{save_pandas_df_to_file.__name__}: {e}")
        return ""


def get_lf(handle, topdir=".", subdir="datafiles", file_format="parquet"):
    try:    
        filepath = os.path.join(topdir, subdir, handle + "." + file_format)
        if file_format == "parquet":
            return pl.scan_parquet(filepath)
        elif file_format == "arrow":
            return pl.scan_ipc(filepath)
        elif file_format == "csv":
            return pl.scan_csv(filepath)
        else:
            raise ValueError(f"Unsupported file format: {file_format}")
    except Exception as e:
        print_flush(f"An error occurred in function {os.path.basename(__file__)}.{get_lf.__name__}: {e}")
        return None


def perform_marker_phenotyping_on_lazyframe(lf, marker_columns, colname_regex_to_replace=r"^Phenotype_\(standardized\)\s+"):

    try:
    
        # Get all column names.
        all_cols = lf.collect_schema().names()

        # Keep only rows that have at least one 1 in marker_columns.
        any_one = pl.any_horizontal([(pl.col(c) == 1) for c in marker_columns])
        lf_filtered = lf.filter(any_one)

        # Get ID columns (all columns that are not marker columns).
        id_cols = [c for c in all_cols if c not in marker_columns]

        # Expand rows: For each row in lf_filtered, create one row per marker column that has a 1.
        marker_phenotyped_lf = (
            lf_filtered.unpivot(
                index=id_cols,
                on=marker_columns,
                variable_name="label_full",
                value_name="value"
            )
            .filter(pl.col("value") == 1)  # Keep only markers that are 1.
            .with_columns(
                pl.col("label_full")
                .str.replace(colname_regex_to_replace, "")
                .alias("label")
            )
            .drop(["label_full", "value"])
        )

        # Validate: Expanded row count must equal total number of 1s across marker columns in original lazyframe.
        total_marker_ones = (
            lf
            .select(pl.sum_horizontal(pl.col(marker_columns)).sum().alias("total_marker_ones"))
            .collect()["total_marker_ones"][0]
        )
        expanded_count = marker_phenotyped_lf.select(pl.len()).collect()["len"][0]

        # Ensure they match.
        assert expanded_count == total_marker_ones, f"Mismatch: num_final_rows={expanded_count}, num_original_ones={total_marker_ones}"
        print_flush(f"Final number of rows: {expanded_count} == total marker 1s: {total_marker_ones}")

        # Return the marker-phenotyped lazyframe.
        return marker_phenotyped_lf
    
    except Exception as e:
        print_flush(f"An error occurred in function {os.path.basename(__file__)}.{perform_marker_phenotyping_on_lazyframe.__name__}: {e}")
        return None


def plot_image_from_frame(
    frame: pl.LazyFrame | pl.DataFrame | pd.DataFrame,
    image_colname: str = "Image ID_(standardized)",
    selected_images: list[str] = ["MS_02__cele_1400w"],
    marker_size: int | float | None = None,  # None -> use Plotly's default
    xcol="Centroid X (µm)_(standardized)",
    ycol="Centroid Y (µm)_(standardized)",
    color_col="label",
    custom_columns=[],
):
    try:

        # Filter and collect to an eager frame; convert to pandas for Plotly Express robustness
        if isinstance(frame, pl.LazyFrame):
            if selected_images:
                df = (
                    frame
                    .filter(pl.col(image_colname).is_in(selected_images))
                    .filter(pl.col(xcol).is_not_null() & pl.col(ycol).is_not_null())
                    .select(custom_columns + [image_colname, xcol, ycol, color_col])
                    .collect()
                )
            else:
                df = (
                    frame
                    .filter(pl.col(xcol).is_not_null() & pl.col(ycol).is_not_null())
                    .select(custom_columns + [image_colname, xcol, ycol, color_col])
                    .collect()
                )
        elif isinstance(frame, pl.DataFrame):
            if selected_images:
                df = (
                    frame
                    .filter(pl.col(image_colname).is_in(selected_images))
                    .filter(pl.col(xcol).is_not_null() & pl.col(ycol).is_not_null())
                    .select(custom_columns + [image_colname, xcol, ycol, color_col])
                )
            else:
                df = (
                    frame
                    .filter(pl.col(xcol).is_not_null() & pl.col(ycol).is_not_null())
                    .select(custom_columns + [image_colname, xcol, ycol, color_col])
                )
        elif isinstance(frame, pd.DataFrame):
            if selected_images:
                mask = frame[image_colname].isin(selected_images)
                mask &= frame[xcol].notna() & frame[ycol].notna()
                df = frame.loc[mask, custom_columns + [image_colname, xcol, ycol, color_col]]
            else:
                mask = frame[xcol].notna() & frame[ycol].notna()
                df = frame.loc[mask, custom_columns + [image_colname, xcol, ycol, color_col]]
        else:
            raise ValueError("Input frame must be a Polars LazyFrame, Polars DataFrame, or Pandas DataFrame.")
        
        hover_data = {
                image_colname: True,
                color_col: True,
                xcol: True,
                ycol: True,
            }
        
        for custom_col in custom_columns:
            hover_data[custom_col] = True

        # Draw the scatter plot.
        fig = px.scatter(
            df,
            x=xcol,
            y=ycol,
            color=color_col,
            title=f"Scatterplot colored by {color_col}",
            hover_data=hover_data,
            render_mode="webgl",
        )

        # Preserve Plotly’s default marker size if marker_size is None.
        if marker_size is not None:
            fig.update_traces(marker=dict(size=marker_size))

        # Ensure 1:1 aspect ratio so spatial distances are faithful.
        fig.update_yaxes(scaleanchor="x", scaleratio=1)
        fig.update_layout(legend_title_text=f"{color_col}")

        # Return the figure.
        return fig

    except Exception as e:
        print_flush(
            f"An error occurred in function "
            f"{os.path.basename(__file__) if '__file__' in globals() else '<interactive>'}."
            f"{plot_image_from_frame.__name__}: {e}"
        )
        return None


def generate_umap(pldf, unique_labels, dist_bin_um_list=[25, 50, 100, 150, 200], area_downsample=0.2, um_per_px=1, cpu_pool_size=None, topdir=".", subdir="results", counts_method="andrew", area_threshold=0.8, custom_areas=True, seed_for_train_test_split=54321, n=2500, keep_images_with_too_little_data=True, train_sample_frac=1.0, test_sample_frac=1.0, de_min_coords=True, mp_start_method=None):
    # Note that cpu_pool_size=None will default to the number of available CPUs.

    # Instantiate the spatial umap object.
    spatial_umap = SpatialUMAP.SpatialUMAP(dist_bin_um=np.array(dist_bin_um_list), um_per_px=um_per_px, area_downsample=area_downsample)

    # Normalize coordinates to start at (0,0) for each image.
    if custom_areas and de_min_coords:
        pldf = pldf.with_columns([
            (pl.col("Xcor") - pl.min("Xcor").over("TMA_core_id")).alias("Xcor"),
            (pl.col("Ycor") - pl.min("Ycor").over("TMA_core_id")).alias("Ycor"),
        ])

    # I don't want to, but convert to pandas for compatibility with SpatialUMAP class.
    spatial_umap.cells = pldf.to_pandas()

    # Set explicitly as numpy array the cell coordinates (x, y).
    spatial_umap.cell_positions = spatial_umap.cells[['Xcor', 'Ycor']].values

    # Set explicitly as one hot data frame the cell labels.
    spatial_umap.cell_labels = pd.get_dummies(spatial_umap.cells['Lineage'])

    # Set the region is to be analyzed (a TMA core is treated similar to a region of a interest).
    spatial_umap.region_ids = spatial_umap.cells.TMA_core_id.unique()

    # Clear metrics.
    spatial_umap.clear_counts()
    spatial_umap.clear_areas()

    # Save the unique labels/lineages/species. This is only used for Andrew's counting method.
    spatial_umap.species = unique_labels

    # Ensure results directory exists.
    os.makedirs(os.path.join(topdir, subdir), exist_ok=True)

    # Determine which counting method to use based on the input parameter. Note I have shown that the methods precisely agree.
    if counts_method == "andrew":
        print_flush("Using Andrew's counts method.")
        spatial_umap.get_counts_And(cpu_pool_size=cpu_pool_size, mp_start_method=mp_start_method)
    else:
        print_flush("Using Baras' counts method.")
        spatial_umap.get_counts(cpu_pool_size=cpu_pool_size)

    # Get the areas of cells and save to pickle file.
    if custom_areas:
        try:
            spatial_umap.get_areas(area_threshold, pool_size=cpu_pool_size, save_file=os.path.join(topdir, subdir, f"areas.csv"), plots_directory=os.path.join(topdir, subdir))  # Sets spatial_umap.cells["area_filter"] and spatial_umap.areas.
        except Exception as e:
            print_flush(f"An error occurred while calculating custom areas, potentially in SpatialUMAP.FitEllipse.fit() in \"hull = ConvexHull(d[idx_fit])\": {e}")
            return spatial_umap, False
    else:
        # Keep in mind areas in the Baras code seem to be in units of pixels squared.
        spatial_umap.cells["area_filter"] = True  # If not using custom areas, set all cells to pass the area filter.
        r0 = np.concatenate(([0], spatial_umap.dist_bin_px))
        # Note not including downsampling at all!
        spatial_umap.areas = (np.pi * (r0[1:] ** 2 - r0[:-1] ** 2)).reshape((1, len(dist_bin_um_list), 1))

    # calculate density base on counts of cells / area of each arc examine
    if custom_areas:
        spatial_umap.density = np.empty(spatial_umap.counts.shape)
        spatial_umap.cells['area_filter'] = ((spatial_umap.areas / spatial_umap.arcs_masks.sum(axis=(0, 1))[np.newaxis, ...]) > area_threshold).all(axis=1)
        spatial_umap.density[spatial_umap.cells['area_filter'].values] = spatial_umap.counts[spatial_umap.cells['area_filter'].values] / spatial_umap.areas[spatial_umap.cells['area_filter'].values][..., np.newaxis]
    else:
        spatial_umap.density = spatial_umap.counts / spatial_umap.areas  # We're not doing this yet (so units are currently in #/px), but we divide by square of um_per_px to get density in units of # per square micron.

    # Output the percentage of the dataset that has been filtered out due to area filtering.
    num_total_cells = spatial_umap.cells.shape[0]
    num_kept_cells = spatial_umap.cells["area_filter"].sum()
    print_flush(f"{100 * (1 - num_kept_cells / num_total_cells):.1f}% of the cells have been filtered out due to area filtering.")

    # Check for dropped images due to insufficient cells passing area filter.
    original_images = set(spatial_umap.region_ids)
    remaining_images = set(spatial_umap.cells[spatial_umap.cells["area_filter"]]["TMA_core_id"].unique())
    dropped_images = original_images - remaining_images
    if dropped_images:
        print_flush(f"WARNING:")
        print_flush(f"  The following images were dropped (i.e., missing scatter plots) due to no cells passing the area filter: {sorted(list(dropped_images))}.")
        print_flush(f"  These images remain: {sorted(list(remaining_images))}.")
        if remaining_images:
            print_flush(f"  Keep in mind that just because some images may not have been dropped, significant numbers of cells in those images may have been dropped.")
        else:
            print_flush(f"  Since no images remain after area filtering, we are aborting UMAP generation.")
            return spatial_umap, False
        remaining_loc = spatial_umap.cells["TMA_core_id"].isin(remaining_images)
        spatial_umap.cells = spatial_umap.cells[remaining_loc]
        spatial_umap.density = spatial_umap.density[remaining_loc.values]

    min_filtered_cells = get_min_positive_values(spatial_umap.cells, group_col="TMA_core_id", boolean_column="area_filter")  # this would be zero if we didn't do the filtering-out line above (spatial_umap.cells = ...)
    print_flush(f"Minimum number of cells passing area filter across all TMA cores: {min_filtered_cells}")

    # Set training and "test" cells for umap training and embedding, respectively. Baras's original code had a hard cutoff of n=2500 so images with fewer than 2*2500 non-filtered-out cells were discarded entirely. n = min(n, min_filtered_cells // 2) allows these images to remain in the analysis with smaller n.
    if keep_images_with_too_little_data:
        n = min(n, min_filtered_cells // 2)

    print_flush(f"Using n_train={int(train_sample_frac*n)} cells per image for UMAP training and n_test={int(test_sample_frac*n)} for testing.")
    spatial_umap.set_train_test(n=n, seed=seed_for_train_test_split, train_sample_frac=train_sample_frac, test_sample_frac=test_sample_frac)

    # Fit umap on training cells.
    spatial_umap.umap_fit = umap.UMAP().fit(spatial_umap.density[spatial_umap.cells['umap_train'].values].reshape((spatial_umap.cells['umap_train'].sum(), -1)))
    print_flush(spatial_umap.umap_fit.embedding_.shape)

    # Apply umap embedding on test cells.
    spatial_umap.umap_test = spatial_umap.umap_fit.transform(spatial_umap.density[spatial_umap.cells['umap_test'].values].reshape((spatial_umap.cells['umap_test'].sum(), -1)))
    print_flush(spatial_umap.umap_test.shape)

    # Save the UMAP coordinates back to the cells dataframe.
    spatial_umap.cells[["umap_1", "umap_2"]] = np.nan
    spatial_umap.cells.loc[spatial_umap.cells['umap_test'].values, ["umap_1", "umap_2"]] = spatial_umap.umap_test

    # # save spatial_umap object as pickle
    # pickle.dump(spatial_umap, open(data_dir + '/pkl/spatial_umap.pkl', 'wb'))

    # Return the spatial UMAP object.
    return spatial_umap, True


def generate_figures(spatial_umap):

    # # load spatial_umap object
    # data_dir = '/home/idies/workspace/Storage/baras/Melanoma'
    # spatial_umap = pickle.load(open(data_dir + '/pkl/spatial_umap.pkl', 'rb'))

    # Spatial UMAP 2D Density Plots By Lineage and with PD-L1 and PD1 overlays

    # set meshgrid / bins for 2d plots based on UMAP x, y distributions
    n_bins = 200
    xx = np.linspace(np.min(spatial_umap.umap_test[:, 0]), np.max(spatial_umap.umap_test[:, 0]), n_bins + 1)
    yy = np.linspace(np.min(spatial_umap.umap_test[:, 1]), np.max(spatial_umap.umap_test[:, 1]), n_bins + 1)
    n_pad = 30

    # set lineages to show and in what order
    lineages_plot = spatial_umap.species
    num_species = len(lineages_plot)

    # get figure and axes
    h = 8
    w = (num_species + 1) / 2 * h
    f, ax = plt.subplots(2, num_species + 1, figsize=(w, h), facecolor='white')

    # color maps
    cmap_viridis = plt.get_cmap('viridis').copy()
    cmap_viridis.set_under('white')
    cmap_magma = plt.get_cmap('magma').copy()
    cmap_magma.set_under('white')
    cmap_bwr = plt.get_cmap('bwr').copy()

    # plot cmaps
    PlottingTools.plt_cmap(ax=ax[1, num_species], cmap=cmap_viridis, extend='max', width=0.01, ylabel='Density')
    # PlottingTools.plt_cmap(ax=ax[2, num_species], cmap=cmap_magma, extend='max', width=0.01, ylabel='PD-L1 MFI')
    # PlottingTools.plt_cmap(ax=ax[3, num_species], cmap=cmap_magma, extend='max', width=0.01, ylabel='PD1 MFI')

    center_ax_col = num_species // 2
    ax_tuples = [(0, i) for i in range(num_species + 1) if i != center_ax_col]

    # clear unneeded axes
    [ax[_].set(visible=False) for _ in ax_tuples]

    # plot 2d denisty in umap of all cells
    PlottingTools.plot_2d_density(spatial_umap.umap_test[:, 0], spatial_umap.umap_test[:, 1], bins=[xx, yy], n_pad=n_pad, ax=ax[0, center_ax_col], cmap=cmap_viridis)
    ax[0, center_ax_col].set(title='All Cells')

    # # get MFI log scaled for PD-L1
    # w = {'PD-L1': spatial_umap.cells['PDL1_520'].values[spatial_umap.cells['umap_test']]}
    # w['PD-L1'] = np.log(0.1 * w['PD-L1'] + 0.1)
    # w['PD-L1'] -= np.min(w['PD-L1'])

    # # get MFI log scaled for PD1
    # w['PD1'] = spatial_umap.cells['PD1_650'].values[spatial_umap.cells['umap_test']]
    # w['PD1'] = np.log(0.1 * w['PD1'] + 0.01)
    # w['PD1'] -= np.min(w['PD1'])
    # w['PD1'] = np.maximum(w['PD1'] - 4, 0)

    for i in range(num_species):
        # cells of lineage(s)
        idx = spatial_umap.cells['Lineage'].values[spatial_umap.cells['umap_test']] == lineages_plot[i]
        ax[1, i].cla()
        ax[1, i].set(title=lineages_plot[i])

        # plot density
        PlottingTools.plot_2d_density(spatial_umap.umap_test[idx, 0], spatial_umap.umap_test[idx, 1], bins=[xx, yy], ax=ax[1, i], cmap=cmap_viridis, vlim=.95)

        # # plot PD-L1 MFI
        # PlottingTools.plot_2d_density(spatial_umap.umap_test[idx, 0], spatial_umap.umap_test[idx, 1], bins=[xx, yy], w=w['PD-L1'][idx], ax=ax[2, i],
        #                 cmap=cmap_magma, vlim=np.array([0, np.quantile(w['PD-L1'], .975)]))
        # # plot PD1 MFI
        # PlottingTools.plot_2d_density(spatial_umap.umap_test[idx, 0], spatial_umap.umap_test[idx, 1], bins=[xx, yy], w=w['PD1'][idx], ax=ax[3, i], cmap=cmap_magma,
        #                 vlim=np.array([0, np.quantile(w['PD1'], .975)]))


    # # Spatial UMAP 2D Density Plots By Lineage and Stratified by 5 Year Survival

    # # get per specimen density maps

    # # set number of bins and get actual binning points based on whole dataset
    # n_bins = 200
    # xx = np.linspace(np.min(spatial_umap.umap_test[:, 0]), np.max(spatial_umap.umap_test[:, 0]), n_bins + 1)
    # yy = np.linspace(np.min(spatial_umap.umap_test[:, 1]), np.max(spatial_umap.umap_test[:, 1]), n_bins + 1)
    # # initialize holding nd matrix for densities
    # n_lineages = len(lineages_plot)
    # # last dim is 0:counts, 1:smoothed, density
    # H = np.empty([n_bins, n_bins, n_lineages + 1, len(spatial_umap.patients['Sample_number']), 2])
    # for i in range(len(spatial_umap.patients['Sample_number'])):
    #     # get cells of this specimen / patient
    #     idx_pts = spatial_umap.cells.loc[spatial_umap.cells['umap_test'], 'Sample_number'] == spatial_umap.patients['Sample_number'].iloc[i]
    #     if np.sum(idx_pts) > 0:
    #         # get counts for lineages
    #         for j in range(len(lineages_plot)):
    #             idx_lineage = spatial_umap.cells.loc[spatial_umap.cells['umap_test'], 'Lineage'].isin([lineages_plot[j]])
    #             H[:, :, j, i, 0], _, _ = np.histogram2d(spatial_umap.umap_test[idx_pts & idx_lineage, 0],
    #                                                     spatial_umap.umap_test[idx_pts & idx_lineage, 1], bins=[xx, yy])
    #         # get counts across all lineages
    #         H[:, :, j + 1, i, 0] = np.nansum(H[:, :, 0:(j + 1), i, 0], axis=2)

    #         # make smoothed density for lineages
    #         for j in range(len(lineages_plot)):
    #             if np.sum(H[:, :, j, i, 0]) > 0:
    #                 H[:, :, j, i, 1] = ndi.gaussian_filter(H[:, :, j, i, 0] / np.sum(H[:, :, j, i, 0]), sigma=0.5)
    #             else:
    #                 H[:, :, j, i, 1] = np.nan
    #         # make smoothed density for all lineages
    #         if np.sum(H[:, :, j + 1, i, 0]) > 0:
    #             H[:, :, j + 1, i, 1] = ndi.gaussian_filter(H[:, :, j + 1, i, 0] / np.sum(H[:, :, j + 1, i, 0]), sigma=0.5)
    #         else:
    #             H[:, :, j + 1, i, 1] = np.nan
    #     else:
    #         H[:, :, :, i, :] = np.nan

    # # specimens with density data across all lineages
    # idx_d = ~np.all(np.all(np.all(np.isnan(H[..., 1]), axis=0), axis=0), axis=0)

    # idx_A = (spatial_umap.patients['Death_5Y'] == 1) & idx_d
    # idx_B = (spatial_umap.patients['Death_5Y'] == 0) & idx_d

    # f, ax = plt.subplots(nrows=3, ncols=7, figsize=(14, 6), facecolor='white')

    # PlottingTools.plt_cmap(ax=ax[0, 6], cmap=cmap_viridis, extend='max', width=0.01, ylabel='Density')
    # PlottingTools.plt_cmap(ax=ax[1, 6], cmap=cmap_viridis, extend='max', width=0.01, ylabel='Density')
    # PlottingTools.plt_cmap(ax=ax[2, 6], cmap=plt.get_cmap('bwr'), extend='both', width=0.01, ylabel='Outcome')

    # d_idx_A = np.nanmean(H[:, :, -1, idx_A, 1], axis=-1).T
    # PlottingTools.plot_2d_density(d_idx_A, bins=[xx, yy], n_pad=30, ax=ax[0, 0], circle_type='bg', cmap=cmap_viridis)
    # ax[0, 0].set(title='All Cells')
    # ax[0, 0].set_ylabel('Survival <= 5 years', rotation='horizontal', ha='right')
    # d_idx_B = np.nanmean(H[:, :, -1, idx_B, 1], axis=-1).T
    # PlottingTools.plot_2d_density(d_idx_B, bins=[xx, yy], n_pad=30, ax=ax[1, 0], circle_type='bg', cmap=cmap_viridis)
    # ax[1, 0].set_ylabel('Survival > 5 years', rotation='horizontal', ha='right')
    # d_diff = d_idx_A - d_idx_B
    # PlottingTools.plot_2d_density(d_diff, bins=[xx, yy], n_pad=30, ax=ax[2, 0], circle_type='arch', cmap=cmap_bwr)
    # ax[2, 0].set_ylabel('Density Differential', rotation='horizontal', ha='right')

    # for i in range(len(lineages_plot)):
    #     d_idx_A = np.nanmean(H[:, :, i, idx_A, 1], axis=-1).T
    #     PlottingTools.plot_2d_density(d_idx_A, bins=[xx, yy], n_pad=30, ax=ax[0, i + 1], circle_type='bg', cmap=cmap_viridis)
    #     ax[0, i + 1].set(title=lineages_plot[i])
    #     d_idx_B = np.nanmean(H[:, :, i, idx_B, 1], axis=-1).T
    #     PlottingTools.plot_2d_density(d_idx_B, bins=[xx, yy], n_pad=30, ax=ax[1, i + 1], circle_type='bg', cmap=cmap_viridis)
    #     d_diff = d_idx_A - d_idx_B
    #     PlottingTools.plot_2d_density(d_diff, bins=[xx, yy], n_pad=30, ax=ax[2, i + 1], circle_type='arch', cmap=cmap_bwr)
