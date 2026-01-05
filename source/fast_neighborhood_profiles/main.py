# Import relevant libraries.
import polars as pl
import os
import plotly.express as px
import numpy as np
import pandas as pd
from fast_neighborhood_profiles import SpatialUMAP
import umap
import PlottingTools
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import framework.utils as framework_utils
import zipfile
import framework.platform_abstraction as pa
import time
import scipy.spatial
import plotly.colors


#### 1. First in load_unified_input_file.py ###############################################################


def get_location_settings():
    return {
        "Available input files": {
            "bucket_name": pa.DATA_OBJECTS_BUCKET_NAME,
            "db_schema": f"{pa.get_user_group(pa.get_current_username())}_group_db.curated_schema",
        },
    }


def get_objects_list(upload_location):
    objects_list = pa.list_objects_in_bucket(
        db_schema=get_location_settings()[upload_location]["db_schema"],
        bucket_name=get_location_settings()[upload_location]["bucket_name"],
    )
    return objects_list


def _subset_csv_to_file(csv_filename="mawa-unified_datafile-TLS_tissue_SF_-20251112_130129_EST.csv", handle="two_images", do_filtering=False, filter_column="Image ID_(standardized)", filter_values=["MS_01__cele_1400w", "MS_02__cele_1400w"], topdir=".", subdir="datafiles", file_format="parquet"):
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
            getattr(lf.filter(pl.col(filter_column).is_in(filter_values)).with_row_index(name="input_index").collect(engine="in-memory"), write_method)(filepath)  # Must specifically be in-memory to get a repeatable index that's consistent with the file's actual rows!
        else:
            getattr(lf.with_row_index(name="input_index").collect(engine="in-memory"), write_method)(filepath)  # Must specifically be in-memory to get a repeatable index that's consistent with the file's actual rows!
        return filepath
    except Exception as e:
        framework_utils.multiprint(f"An error occurred while subsetting a CSV to a file: {e}", (print,))
        raise


def _get_lf(handle, topdir=".", subdir="datafiles", file_format="parquet"):
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
        framework_utils.multiprint(f"An error occurred while getting a lazyframe: {e}", (print,))
        raise


# Load the lazyframe from the specified file using an intermediate file.
def load_unified_input_file_data(file_format, db_schema, bucket_name, object_filename, topdir):

    try:

        # Shortcuts, the first for generalizability.
        full_filenames = [object_filename]
        input_dir = os.path.join(topdir, "input")

        # Download the file from the server.
        pa.download_objects_parallel(
            object_names=full_filenames,
            db_schema=db_schema,
            bucket_name=bucket_name,
            dest_dir=input_dir,
        )

        # Unzip the downloaded file.
        unzipped_paths = []
        for full_filename in full_filenames:
            local_path = os.path.join(input_dir, full_filename)
            lower = full_filename.lower()
            if lower.endswith(".zip"):
                with zipfile.ZipFile(local_path, 'r') as zip_ref:
                    zip_ref.extractall(input_dir)
                base_name = full_filename.removesuffix(".zip")
                os.remove(local_path)
                unzipped_paths.append(os.path.join(input_dir, base_name))
            elif lower.endswith(".gz"):
                # Already gunzipped by download_objects_parallel -> use basename without .gz
                base_name = full_filename.removesuffix(".gz")
                unzipped_paths.append(os.path.join(input_dir, base_name))
            else:
                # Not compressed; use as-is
                unzipped_paths.append(local_path)

        # Generate an intermediate file from which to load the lazyframe.
        for unzipped_path in unzipped_paths:
            filename = os.path.basename(unzipped_path)
            _subset_csv_to_file(csv_filename=filename, handle="unified_input_file", topdir=topdir, subdir="input", file_format=file_format)
            os.remove(unzipped_path)

        # Load the lazyframe.
        lf = _get_lf("unified_input_file", topdir=topdir, subdir="input", file_format=file_format)

        # # Store the local filepath relative to the session directory.
        # local_filepath = filepath.removeprefix(topdir + os.sep)

        # Return the lazyframe.
        return lf
    
    except Exception as e:
        framework_utils.multiprint(f"Unable to load unified input file data (inputs: file_format={file_format}, db_schema={db_schema}, bucket_name={bucket_name}, object_filename={object_filename}): {e}", (print,))
        raise


#### 2. First in phenotype.py ###############################################################


# Get the marker column names from the lazyframe.
def get_marker_columns(lf, prefix="Phenotype_(standardized) "):
    marker_columns = [column for column in lf.collect_schema().names() if column.startswith(prefix)]
    marker_columns_ordered = lf.select(pl.col(marker_columns).sum()).melt(variable_name="column", value_name="sum").sort("sum", descending=True).select("column").collect(engine="streaming").to_series().to_list()
    marker_columns_ordered_no_prefix = [x.removeprefix(prefix) for x in marker_columns_ordered]
    return marker_columns_ordered_no_prefix, marker_columns_ordered


def get_phenotyped_metadata(lf_phenotyped, phenotyping_method):

    # Obtain the resulting labels in decreasing frequency order.
    ordered_labels = (
        lf_phenotyped
        .group_by("label")
        .agg(pl.count().alias("freq"))
        .sort(["freq", "label"], descending=[True, False])
        .select("label")
        .collect(engine="streaming")
        .to_series()
        .to_list()
    )

    colors = px.colors.qualitative.Plotly

    return {
        "num_phenotyped_rows": lf_phenotyped.select(pl.len()).collect(engine="streaming").item(),
        "unique_labels": ordered_labels,
        "unique_image_ids": lf_phenotyped.select(pl.col("Image ID_(standardized)").unique().sort()).collect(engine="streaming").to_series().to_list(),
        "phenotype_color_map": {label: colors[i % len(colors)] for i, label in enumerate(ordered_labels)},
        "phenotyping_method": phenotyping_method,
    }


def perform_marker_phenotyping_on_lazyframe(lf, marker_columns_with_prefix, colname_regex_to_replace=r"^Phenotype_\(standardized\)\s+"):

    try:

        # Keep only rows that have at least one 1 in marker_columns.
        any_one = pl.any_horizontal([(pl.col(c) == 1) for c in marker_columns_with_prefix])
        lf = lf.filter(any_one)

        # Get all column names.
        all_cols = lf.collect_schema().names()

        # Get ID columns (all columns that are not marker columns).
        id_cols = [c for c in all_cols if c not in marker_columns_with_prefix]

        # Expand rows: For each row in lf, create one row per marker column that has a 1.
        marker_phenotyped_lf = (
            lf.unpivot(
                index=id_cols,
                on=marker_columns_with_prefix,
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
            .select(pl.sum_horizontal(pl.col(marker_columns_with_prefix)).sum().alias("total_marker_ones"))
            .collect(engine="streaming")["total_marker_ones"][0]
        )
        expanded_count = marker_phenotyped_lf.select(pl.len()).collect(engine="streaming")["len"][0]

        # Ensure they match.
        assert expanded_count == total_marker_ones, f"Mismatch: num_final_rows={expanded_count}, num_original_ones={total_marker_ones}"
        framework_utils.multiprint(f"Final number of rows: {expanded_count} == total marker 1s: {total_marker_ones}", (print,))

        # Return the marker-phenotyped lazyframe.
        return marker_phenotyped_lf
    
    except Exception as e:
        framework_utils.multiprint(f"An error occurred while performing marker phenotyping on a lazyframe: {e}", (print,))
        raise


def obtain_species_column_from_markers_columns(lf, marker_columns_with_prefix, marker_columns):

    # Keep only rows that have at least one 1 in marker_columns.
    any_one = pl.any_horizontal([(pl.col(c) == 1) for c in marker_columns_with_prefix])
    lf = lf.filter(any_one)

    # We need to map to the short column names below, but sometimes the "to" column name already exists in the dataset (due to how the user might generate the unified datafile). So we first rename only those columns that need renaming, and we assume if they already have the name without the prefix, then it's a duplicate of the prefixed column.
    lf_columns = lf.collect_schema().names()
    column_mapping = {x: y for x, y in zip(marker_columns_with_prefix, marker_columns) if y not in lf_columns}

    # Create a "Species" column that concatenates all marker column names with a 1 in the row.
    lf = lf.rename(column_mapping).with_columns(
        Species = pl.concat_str(
            [
                pl.when(pl.col(c).cast(pl.Int8).fill_null(0) == 1)
                .then(pl.lit(f"{c}+"))
                .otherwise(None)
                for c in marker_columns
            ],
            separator=" ",
            ignore_nulls=True
        )
    )

    # Return the modified lazyframe.
    return lf


def create_species_assignments_table(lf):
    '''
    There is no "collect" in this function and is extremely fast.
    '''

    # Get the species counts.
    lf_species_counts = lf.group_by("Species").agg(pl.count().alias("Count in dataset")).sort("Count in dataset", descending=True)

    # Create a column that we'll edit, and reorder the columns.
    lf_species_counts = lf_species_counts.with_columns(pl.col("Species").alias("Edited name")).select(pl.col(["Species", "Edited name", "Count in dataset"]))

    # Return the assignments lazyframe.
    return lf_species_counts


def perform_species_phenotyping_on_lazyframe(lf, marker_columns_with_prefix, df_species_assignments):

    # For parity with perform_marker_phenotyping_on_lazyframe(), calculate marker_columns, don't send it in.
    marker_columns = [x.removeprefix("Phenotype_(standardized) ") for x in marker_columns_with_prefix]

    # Add the species column to the main lazyframe.
    lf = obtain_species_column_from_markers_columns(lf, marker_columns_with_prefix, marker_columns)

    # Create a polars dataframe out of the needed columns of the assignments table, renaming to the final "label" column.
    assignments_pl = (
        pl.from_pandas(df_species_assignments)
        .select(["Species", "Edited name"])
        .rename({"Edited name": "label"})  # final column name in main LF
        .lazy()
    )

    # Assign the names to the species and drop rows that don't have mappings.
    lf = lf.join(assignments_pl, on="Species", how="inner")

    # Return the phenotyped lazyframe.
    return lf


def plot_image_from_frame(
    frame,
    image_colname: str = "Image ID_(standardized)",
    selected_images: list[str] = ["MS_02__cele_1400w"],
    marker_size=None,
    xcol="Centroid X (µm)_(standardized)",
    ycol="Centroid Y (µm)_(standardized)",
    color_col="label",
    color_map=None,
    custom_columns=[],
    highlight_index_col="sumap_cell_index",
    highlight_indices=(),
    default_marker_size=None,
    default_marker_line_color: str = "rgba(0,0,0,0.3)",
    default_marker_alpha: float = 0.5,
    highlight_marker_size=None,
    highlight_marker_line_color: str = "black",
    highlight_marker_alpha: float = 1.0,
    highlight_show_legend: bool = False,
    plot_faithful_object_sizes: bool = False,
    xmin_col="XMin",
    xmax_col="XMax",
    ymin_col="YMin",
    ymax_col="YMax",
    frame_with_faithful_columns=None,
    common_index="input_index",
    missing_label_value = "Other",
    sort_index_col="",  # If the figure keeps redrawing, it's almost certainly because the streaming engine causes row shuffles so it appears to plotly/etc. that the figure is always changing. E.g., this happens with species phenotyping (not marker phenotyping). So to solve that, we should input a sort_index_col that will sort the dataframe to a consistent order before plotting.
):
    try:

        if plot_faithful_object_sizes:
            if isinstance(frame_with_faithful_columns, (pl.LazyFrame, pl.DataFrame)):
                columns = frame_with_faithful_columns.collect_schema().names()
            elif isinstance(frame_with_faithful_columns, pd.DataFrame):
                columns = frame_with_faithful_columns.columns
            else:
                raise ValueError("Faithful columns frame must be a Polars LazyFrame, Polars DataFrame, or Pandas DataFrame.")
            faithful_columns = [col for col in [xmin_col, xmax_col, ymin_col, ymax_col] if col in columns]
            if len(faithful_columns) < 4:
                plot_faithful_object_sizes = False
                faithful_columns = []
                framework_utils.multiprint(f"Warning: Not all faithful object size columns found in the frame that's thought to contain them. Disabling faithful object size plotting.", (print,))
        else:
            faithful_columns = []

        if plot_faithful_object_sizes:
            if isinstance(frame_with_faithful_columns, (pl.LazyFrame, pl.DataFrame)):
                frame = frame.join(frame_with_faithful_columns.select(pl.col([common_index] + faithful_columns)).unique(subset=[common_index], keep="first"), on=common_index, how="left")
            elif isinstance(frame_with_faithful_columns, pd.DataFrame):
                framework_utils.multiprint("NOTE: WE MAY WANT TO ENFORCE DEDUPLICATION OF THE RHS AS WE DO FOR POLARS ABOVE; NOT DOING THAT FOR PANDAS HERE YET!", (print,))
                frame = frame.merge(frame_with_faithful_columns[[common_index] + faithful_columns], on=common_index, how="left")
            else:
                raise ValueError("Faithful columns frame must be a Polars LazyFrame, Polars DataFrame, or Pandas DataFrame.")

        # Efficiently convert the input frame to a pandas DataFrame with necessary filtering.
        # If there is slowness, we can try keeping as a Polars DataFrame and using Plotly's ability to plot from Polars DataFrames directly, proceeding with polars dataframes in all operations below.
        cols_to_keep = [color_col, image_colname, xcol, ycol] + custom_columns + faithful_columns
        if isinstance(frame, pl.LazyFrame):
            base = frame.filter(pl.col(xcol).is_not_null() & pl.col(ycol).is_not_null())
            if selected_images:
                base = base.filter(pl.col(image_colname).is_in(selected_images))
            lf = base.select(pl.col(cols_to_keep)).with_columns([
                pl.col(xcol).cast(pl.Float32),
                pl.col(ycol).cast(pl.Float32),
            ])
            if sort_index_col and sort_index_col in cols_to_keep:
                df = lf.sort(pl.col(sort_index_col)).collect(engine="streaming").to_pandas()
            else:
                df = lf.collect(engine="streaming").to_pandas()
        elif isinstance(frame, pl.DataFrame):
            base = frame.filter(pl.col(xcol).is_not_null() & pl.col(ycol).is_not_null())
            if selected_images:
                base = base.filter(pl.col(image_colname).is_in(selected_images))
            df = base.select(pl.col(cols_to_keep)).with_columns([
                pl.col(xcol).cast(pl.Float32),
                pl.col(ycol).cast(pl.Float32),
            ]).to_pandas()
        elif isinstance(frame, pd.DataFrame):
            mask = frame[xcol].notna() & frame[ycol].notna()
            if selected_images:
                mask &= frame[image_colname].isin(selected_images)
            df = frame.loc[mask, cols_to_keep].copy()
            df[xcol] = df[xcol].astype(np.float32)
            df[ycol] = df[ycol].astype(np.float32)
        else:
            raise ValueError("Input frame must be a Polars LazyFrame, Polars DataFrame, or Pandas DataFrame.")

        # Prepare hover data columns.
        hover_template = \
        f'<b>{color_col}</b>: %{{customdata[0]}}<br>' + \
        f'<b>{image_colname}</b>: %{{customdata[1]}}<br>' + \
        f'<b>{xcol}</b>: %{{customdata[2]}}<br>' + \
        f'<b>{ycol}</b>: %{{customdata[3]}}' + \
        "".join([f'<br><b>{custom_column}</b>: %{{customdata[{4 + i}]}}' for i, custom_column in enumerate(custom_columns)])
        
        # Determine if we should highlight.
        do_highlight = highlight_index_col in custom_columns and len(highlight_indices) > 0
        if do_highlight:
            mask_high = df[highlight_index_col].isin(highlight_indices)
            df_high = df[mask_high]
            df_other = df[~mask_high]
        else:
            df_high = None
            df_other = df

        # Get unique labels and assign colors using Plotly's default color sequence.
        if not color_map:
            unique_labels = sorted(df[color_col].unique())
            colors = px.colors.qualitative.Plotly
            color_map = {label: colors[i % len(colors)] for i, label in enumerate(unique_labels)}
        else:
            unique_labels = color_map.keys()
        if missing_label_value in unique_labels:
            unique_labels = [missing_label_value] + [label for label in unique_labels if label != missing_label_value]

        # Determine effective marker sizes.
        effective_default_size = default_marker_size if default_marker_size is not None else (marker_size if marker_size is not None else 5)
        effective_highlight_size = highlight_marker_size if highlight_marker_size is not None else (3 * marker_size if marker_size is not None else 15)

        # Create figure and add non-highlighted traces.
        fig = go.Figure()
        for label in unique_labels:
            label_mask = df_other[color_col] == label
            if label_mask.any():
                label_data = df_other[label_mask]
                if not plot_faithful_object_sizes:
                    fig.add_trace(go.Scattergl(
                        x=label_data[xcol],
                        y=label_data[ycol],
                        mode="markers",
                        name=label,
                        marker=dict(
                            size=effective_default_size,
                            color=color_map[label],
                            line=dict(color=default_marker_line_color, width=1),
                        ),
                        opacity=default_marker_alpha,
                        customdata=label_data.values,
                        hovertemplate=hover_template,
                    ))
                else:
                    fig.add_trace(go.Bar(
                        x=((label_data[xmin_col] + label_data[xmax_col]) / 2),
                        y=label_data[ymax_col] - label_data[ymin_col],
                        width=label_data[xmax_col] - label_data[xmin_col],
                        base=label_data[ymin_col],
                        name=label,
                        marker=dict(
                            color=color_map[label],
                            opacity=default_marker_alpha,
                        ),
                        customdata=label_data.values,
                        hovertemplate=hover_template,
                    ))

        # Add highlighted traces if applicable.
        if do_highlight and df_high is not None and not df_high.empty:
            for label in unique_labels:
                label_mask = df_high[color_col] == label
                if label_mask.any():
                    label_data = df_high[label_mask]
                    if not plot_faithful_object_sizes:
                        fig.add_trace(go.Scattergl(
                            x=label_data[xcol],
                            y=label_data[ycol],
                            mode="markers",
                            name=label,
                            showlegend=highlight_show_legend,
                            marker=dict(
                                size=effective_highlight_size,
                                color=color_map[label],
                                line=dict(color=highlight_marker_line_color, width=1),
                            ),
                            opacity=highlight_marker_alpha,
                            customdata=label_data.values,
                            hovertemplate=hover_template,
                        ))
                    else:
                        fig.add_trace(go.Bar(
                            x=((label_data[xmin_col] + label_data[xmax_col]) / 2),
                            y=label_data[ymax_col] - label_data[ymin_col],
                            width=label_data[xmax_col] - label_data[xmin_col],
                            base=label_data[ymin_col],
                            name=label,
                            showlegend=highlight_show_legend,
                            marker=dict(
                                color=color_map[label],
                                opacity=highlight_marker_alpha,
                            ),
                            customdata=label_data.values,
                            hovertemplate=hover_template,
                        ))

        # Update layout.
        fig.update_yaxes(scaleanchor="x", scaleratio=1)
        fig.update_layout(
            title=f"Scatterplot colored by {color_col.lower()}",
            legend_title_text=color_col,
            xaxis_title=xcol,
            yaxis_title=ycol,
        )

        # Return the figure.
        return fig

    except Exception as e:
        framework_utils.multiprint(f"An error occurred while plotting an image from a frame: {e}", (print,))
        raise


#### 3. First in delete_cells.py (none yet) ########################################################


#### 4. First in run_spatial_umap.py ###############################################################


# Wrapper to modify core function for compatibility with job input/output dictionary standards without the need for an intermediate polars dataframe.
def generate_umap_wrapper(**inputs):

    # Modify inputs as needed (standard format, i.e., what's expected of a dictionary as in utils.deserialize_binary_files_to_dictionary() should come in).
    lf = inputs["LAZYFRAMES"]["phenotyped"]["lf"]
    lf = (
        lf
        .rename({"Image ID_(standardized)": "TMA_core_id", "Centroid X (µm)_(standardized)": "Xcor", "Centroid Y (µm)_(standardized)": "Ycor", "label": "Lineage"})
        .select(pl.col(["TMA_core_id", "Xcor", "Ycor", "Lineage", "input_index"]))
        .sort(by="TMA_core_id")
        )
    inputs["lf"] = lf
    del inputs["LAZYFRAMES"]

    # Run the core function.
    outputs = _generate_umap_lf_input(**inputs)  # what comes out of this: dict(spatial_umap=spatial_umap, complete_success=True)

    # Return the properly formatted outputs.
    return outputs


def _get_min_positive_values(pd_df, group_col="TMA_core_id", boolean_column="area_filter"):
    return pd_df.groupby(group_col)[boolean_column].sum().min()


def _generate_umap_lf_input(lf, unique_labels, dist_bin_um_list=[25, 50, 100, 150, 200], area_downsample=0.2, um_per_px=1, cpu_pool_size=None, results_topdir=".", subdir="results", counts_method="andrew", area_threshold=0.8, custom_areas=True, seed_for_train_test_split=54321, n=2500, keep_images_with_too_little_data=True, train_sample_frac=1.0, test_sample_frac=1.0, de_min_coords=True, mp_start_method=None):
    # Note that cpu_pool_size=None will default to the number of available CPUs.

    # Instantiate the spatial umap object.
    spatial_umap = SpatialUMAP.SpatialUMAP(dist_bin_um=np.array(dist_bin_um_list), um_per_px=um_per_px, area_downsample=area_downsample)

    # Normalize coordinates to start at (0,0) for each image.
    if custom_areas and de_min_coords:

        # Print the minima first.
        min_coords_df = lf.group_by("TMA_core_id").agg([
            pl.min("Xcor").alias("min_Xcor"),
            pl.min("Ycor").alias("min_Ycor"),
        ]).collect(engine="streaming")
        framework_utils.multiprint("Minimum coordinates per TMA core:", (print,))
        framework_utils.multiprint(min_coords_df, (print,))

        # Perform normalization.
        lf = lf.with_columns([
            (pl.col("Xcor") - pl.min("Xcor").over("TMA_core_id")).alias("Xcor"),
            (pl.col("Ycor") - pl.min("Ycor").over("TMA_core_id")).alias("Ycor"),
        ])

        #### Note I may need to collect_schema() here potentially to register new columns per my experience needing that with with_row_index(). ####

    # I don't want to, but convert to pandas for compatibility with SpatialUMAP class.
    spatial_umap.cells = lf.collect(engine="streaming").to_pandas()

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

    # Ensure results directory exists. Probably only needed if doing custom areas.
    os.makedirs(os.path.join(results_topdir, subdir), exist_ok=True)

    # Determine which counting method to use based on the input parameter. Note I have shown that the methods precisely agree.
    if counts_method == "andrew":
        framework_utils.multiprint("Using Andrew's counts method.", (print,))
        spatial_umap.get_counts_And(cpu_pool_size=cpu_pool_size, mp_start_method=mp_start_method)
    else:
        framework_utils.multiprint("Using Baras' counts method.", (print,))
        spatial_umap.get_counts(cpu_pool_size=cpu_pool_size)

    # Get the areas of cells and save to pickle file.
    if custom_areas:
        try:
            spatial_umap.get_areas(area_threshold, pool_size=cpu_pool_size, save_file=os.path.join(results_topdir, subdir, f"areas.csv"), plots_directory=os.path.join(results_topdir, subdir))  # Sets spatial_umap.cells["area_filter"] and spatial_umap.areas.
        except Exception as e:
            framework_utils.multiprint(f"An error occurred while calculating custom areas: {e}", (print,))
            raise
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
    framework_utils.multiprint(f"{100 * (1 - num_kept_cells / num_total_cells):.1f}% of the cells have been filtered out due to area filtering.", (print,))

    # Check for dropped images due to insufficient cells passing area filter.
    original_images = set(spatial_umap.region_ids)
    remaining_images = set(spatial_umap.cells[spatial_umap.cells["area_filter"]]["TMA_core_id"].unique())
    dropped_images = original_images - remaining_images
    if dropped_images:
        framework_utils.multiprint(f"WARNING:", (print,))
        framework_utils.multiprint(f"  The following images were dropped (i.e., missing scatter plots) due to no cells passing the area filter: {sorted(list(dropped_images))}.", (print,))
        framework_utils.multiprint(f"  These images remain: {sorted(list(remaining_images))}.", (print,))
        if remaining_images:
            framework_utils.multiprint(f"  Keep in mind that just because some images may not have been dropped, significant numbers of cells in those images may have been dropped.", (print,))
        else:
            framework_utils.multiprint(f"  Since no images remain after area filtering, we are aborting UMAP generation.", (print,))
            return spatial_umap, False
        remaining_loc = spatial_umap.cells["TMA_core_id"].isin(remaining_images)
        spatial_umap.cells = spatial_umap.cells[remaining_loc]
        spatial_umap.density = spatial_umap.density[remaining_loc.values]

    min_filtered_cells = _get_min_positive_values(spatial_umap.cells, group_col="TMA_core_id", boolean_column="area_filter")  # this would be zero if we didn't do the filtering-out line above (spatial_umap.cells = ...)
    framework_utils.multiprint(f"Minimum number of cells passing area filter across all TMA cores: {min_filtered_cells}", (print,))
    # Set training and "test" cells for umap training and embedding, respectively. Baras's original code had a hard cutoff of n=2500 so images with fewer than 2*2500 non-filtered-out cells were discarded entirely. n = min(n, min_filtered_cells // 2) allows these images to remain in the analysis with smaller n.
    if keep_images_with_too_little_data:
        n = min(n, min_filtered_cells // 2)

    framework_utils.multiprint(f"Using n_train={int(train_sample_frac*n)} cells per image for UMAP training and n_test={int(test_sample_frac*n)} for testing.", (print,))
    spatial_umap.set_train_test(n=n, seed=seed_for_train_test_split, train_sample_frac=train_sample_frac, test_sample_frac=test_sample_frac)

    # Fit umap on training cells.
    spatial_umap.umap_fit = umap.UMAP().fit(spatial_umap.density[spatial_umap.cells['umap_train'].values].reshape((spatial_umap.cells['umap_train'].sum(), -1)))
    framework_utils.multiprint(str(spatial_umap.umap_fit.embedding_.shape), (print,))

    # Apply umap embedding on test cells.
    spatial_umap.umap_test = spatial_umap.umap_fit.transform(spatial_umap.density[spatial_umap.cells['umap_test'].values].reshape((spatial_umap.cells['umap_test'].sum(), -1)))
    framework_utils.multiprint(spatial_umap.umap_test.shape, (print,))

    # Save the UMAP coordinates back to the cells dataframe.
    spatial_umap.cells[["umap_1", "umap_2"]] = np.nan
    spatial_umap.cells.loc[spatial_umap.cells['umap_test'].values, ["umap_1", "umap_2"]] = spatial_umap.umap_test

    # # save spatial_umap object as pickle
    # pickle.dump(spatial_umap, open(data_dir + '/pkl/spatial_umap.pkl', 'wb'))

    # Return the result as a dictionary.
    return dict(spatial_umap=spatial_umap, complete_success=True)


# Not currently used.
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


def fast_neighbors_counts_for_block2(df_image, image_name, coord_column_names, phenotypes, radii, phenotype_column_name, max_chunk_size_in_mb=200, data_struct="numpy", kdtree_str="kdtree", num_rows_per_chunk_neighb=100_000, chunk_neighbor_trees=False):
    # Haven't checked neighbor chunking accuracy, though it's probably fine. But check before ever using it in production.

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
    start_time = time.time()
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
    time_before = time.time()
    if data_struct == "pandas":
        df_counts_holder = [pd.DataFrame(0, index=phenotypes, columns=df_image_index) for _ in radii]
    elif data_struct == "numpy":
        df_counts_holder = [np.zeros((len(phenotypes), len(df_image_index)), dtype=np.int32) for _ in radii]
    elapsed_time_dataframe += time.time() - time_before

    print(f"data struct: {elapsed_time_dataframe:.2f} seconds, kdtree: {elapsed_time_kdtree:.2f} seconds", flush=True)
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
        time_before = time.time()
        if data_struct == "pandas":
            ser_curr_neighbor_phenotype = df_image[phenotype_column_name] == neighbor_phenotype
        elif data_struct == "numpy":
            ser_curr_neighbor_phenotype = nd_phenotype == neighbor_phenotype
            if chunk_neighbor_trees:
                nd_curr_neighbor_phenotype = np.where(ser_curr_neighbor_phenotype)[0]
                num_cells_curr_phenotype = len(nd_curr_neighbor_phenotype)
                start_indices_neighb = np.arange(0, num_cells_curr_phenotype, num_rows_per_chunk_neighb, dtype=np.int64)
                stop_indices_neighb = start_indices_neighb + num_rows_per_chunk_neighb
        elapsed_time_dataframe += time.time() - time_before

        # Construct the KDTree for the current phenotype in the entire current image. This represents the neighbors
        time_before = time.time()
        if data_struct == "pandas":
            neighbor_trees.append(kdtree_func(df_image.loc[ser_curr_neighbor_phenotype, coord_column_names]))
        elif data_struct == "numpy":
            if not chunk_neighbor_trees:
                neighbor_trees.append(kdtree_func(nd_coords[ser_curr_neighbor_phenotype, :]))
            else:
                neighbor_trees_curr_phenotype = []
                for start_index_neighb, stop_index_neighb in zip(start_indices_neighb, np.minimum(stop_indices_neighb, num_cells_curr_phenotype)):
                    neighbor_trees_curr_phenotype.append(kdtree_func(nd_coords[nd_curr_neighbor_phenotype[start_index_neighb:stop_index_neighb], :]))
                neighbor_trees.append(neighbor_trees_curr_phenotype)
        elapsed_time_kdtree += time.time() - time_before

    print(f"data struct: {elapsed_time_dataframe:.2f} seconds, kdtree: {elapsed_time_kdtree:.2f} seconds", flush=True)
    elapsed_time_dataframe = 0
    elapsed_time_kdtree = 0
    
    # For each chunk of centers...
    for start_index, stop_index in zip(start_indices, np.minimum(stop_indices, num_cells)):

        # Construct the KDTree for the current chunk. This represents the centers
        time_before = time.time()
        if data_struct == "pandas":
            center_tree = kdtree_func(df_image.loc[df_image_index[start_index:stop_index], coord_column_names])
        elif data_struct == "numpy":
            center_tree = kdtree_func(nd_coords[start_index:stop_index, :])
        elapsed_time_kdtree += time.time() - time_before

        # For each radius, which should be monotonically increasing and start with 0...
        for iradius, radius in enumerate(radii):

            # For each neighbor tree...
            for ineighbor_phenotype, curr_neighbor_tree in enumerate(neighbor_trees):

                # Get the list of lists containing the indices of the neighbors for each center
                time_before = time.time()
                if data_struct == "pandas":
                    neighbors_for_radius = center_tree.query_ball_tree(curr_neighbor_tree, radius)
                elif data_struct == "numpy":
                    if not chunk_neighbor_trees:
                        neighbors_for_radius = center_tree.query_ball_tree(curr_neighbor_tree, radius)
                    else:
                        neighbors_for_radius = []
                        for neighbor_tree_chunk in curr_neighbor_tree:
                            neighbors_for_radius.append(center_tree.query_ball_tree(neighbor_tree_chunk, radius))
                elapsed_time_kdtree += time.time() - time_before

                # In the correct dataframe (corresponding to the current radius), set the counts of neighbors (of the current phenotype) for each center
                time_before = time.time()
                if data_struct == "pandas":
                    df_counts_holder[iradius].iloc[ineighbor_phenotype, start_index:stop_index] = [len(neighbors_for_center) for neighbors_for_center in neighbors_for_radius]
                elif data_struct == "numpy":
                    if not chunk_neighbor_trees:
                        df_counts_holder[iradius][ineighbor_phenotype, start_index:stop_index] = np.fromiter(map(len, neighbors_for_radius), dtype=np.int32)
                    else:
                        counts_per_center = np.zeros(stop_index - start_index, dtype=np.int32)
                        for neighbors_for_radius_chunk in neighbors_for_radius:
                            counts_per_center += np.fromiter(map(len, neighbors_for_radius_chunk), dtype=np.int32)
                        df_counts_holder[iradius][ineighbor_phenotype, start_index:stop_index] = counts_per_center
                elapsed_time_dataframe += time.time() - time_before

    print(f"data struct: {elapsed_time_dataframe:.2f} seconds, kdtree: {elapsed_time_kdtree:.2f} seconds", flush=True)
    elapsed_time_dataframe = 0
    elapsed_time_kdtree = 0
    
    # For each annulus, i.e., each radius range...
    df_counts_holder_annulus = []
    if data_struct == "numpy":
        df_counts_annulus_column_holder = []  # "column" not "index" since we transpose
    for iradius in range(len(radii) - 1):

        # Get the counts of neighbors in the current annulus
        time_before = time.time()
        df_counts_curr_annulus = df_counts_holder[iradius + 1] - df_counts_holder[iradius]
        elapsed_time_dataframe += time.time() - time_before

        # Rename the index to reflect the current radius range
        radius_range_str = f'({radii[iradius]}, {radii[iradius + 1]}]'
        time_before = time.time()
        if data_struct == "pandas":
            df_counts_curr_annulus.index = [f'{phenotype} in {radius_range_str}' for phenotype in phenotypes]
        elif data_struct == "numpy":
            df_counts_annulus_column_holder.append([f'{phenotype} in {radius_range_str}' for phenotype in phenotypes])
        elapsed_time_dataframe += time.time() - time_before

        # Add a transpose of this (so centers are in rows and phenotypes/radii are in columns) to the running list of annulus dataframes
        time_before = time.time()
        df_counts_holder_annulus.append(df_counts_curr_annulus.T)
        elapsed_time_dataframe += time.time() - time_before

    print(f"data struct: {elapsed_time_dataframe:.2f} seconds, kdtree: {elapsed_time_kdtree:.2f} seconds", flush=True)
    elapsed_time_dataframe = 0
    elapsed_time_kdtree = 0
    
    # Concatenate the annulus dataframes to get the final dataframe of neighbor counts for the current image
    time_before = time.time()
    if data_struct == "pandas":
        df_curr_counts = pd.concat(df_counts_holder_annulus, axis='columns')
    elif data_struct == "numpy":
        df_curr_counts = np.concatenate(df_counts_holder_annulus, axis=1, dtype=np.int32)
    elapsed_time_dataframe += time.time() - time_before

    # Convert the dataframe to a more efficient dtype (from int64)
    time_before = time.time()
    if data_struct == "pandas":
        df_curr_counts = df_curr_counts.astype(np.int32)
    elif data_struct == "numpy":
        df_curr_counts = pd.DataFrame(df_curr_counts, index=df_image_index, columns=[col for sublist in df_counts_annulus_column_holder for col in sublist])
    elapsed_time_dataframe += time.time() - time_before

    print(f"data struct: {elapsed_time_dataframe:.2f} seconds, kdtree: {elapsed_time_kdtree:.2f} seconds", flush=True)
    
    # Print the time taken to calculate the neighbor counts for the current image
    print(f'  ...finished calculating neighbor counts for image {image_name} ({len(df_image)} cells) in {time.time() - start_time:.2f} seconds', flush=True)

    # Return the final dataframe of neighbor counts for the current image
    return df_curr_counts


# # Not currently used.
# def test_kdtree_accepts_empty_input():
#     import numpy as np
#     empty = np.empty((0, 2), dtype=float)
#     scipy.spatial.KDTree(empty)  # should succeed; if it raises, you’ll know at test time


def _save_pandas_df_to_file(pd_df, handle="two_images", topdir=".", file_format="parquet", subdir="datafiles", index_column_name="index"):
    try:
        filepath = os.path.join(topdir, subdir, handle + "." + file_format)

        if index_column_name in pd_df.columns:
            framework_utils.multiprint(f"{index_column_name} is already a column in the Pandas DataFrame. This is fine since we want to rebuild objects as much from scratch as possible. Deleting and recreating this column now.", (print,))
            del pd_df[index_column_name]

        # Materialize index into a real column before writing the Pandas dataframe.
        pd_df.insert(0, index_column_name, range(len(pd_df)))

        if file_format == "parquet":
            pd_df.to_parquet(filepath, index=False)
        elif file_format == "arrow":
            pl.from_pandas(pd_df).write_ipc(filepath)
        elif file_format == "csv":
            pd_df.to_csv(filepath, index=False)
        else:
            raise ValueError(f"Unsupported file format: {file_format}")
        return filepath
    except Exception as e:
        framework_utils.multiprint(f"An error occurred while saving a pandas DataFrame to a file: {e}", (print,))
        raise


def save_and_load_pandas_df_to_lf(pd_df, handle, file_format, topdir, index_column_name="index"):
    _save_pandas_df_to_file(pd_df, handle=handle, file_format=file_format, topdir=topdir, subdir="input", index_column_name=index_column_name)
    lf = _get_lf(handle, topdir=topdir, subdir="input", file_format=file_format)
    return lf


def get_true_false_color_map():
    colors = px.colors.qualitative.Plotly
    color_map = {label: colors[i % len(colors)] for i, label in enumerate((True, False))}
    return color_map


#### 5. First in assign_neighborhood_types.py ###############################################################


def plot_neighborhood_profile(data, plot_type, labels_axis_1, labels_axis_2, axis_1_name="Distance bin (µm)", axis_2_name="Phenotype", value_name="Density", color_map=None):
    # Note the axes here refer to the axes of data, not the plot axes.
    # data shape: (axis_0, axis_1, axis_2) where axis_0 is the distribution dimension (e.g., cells), axis_1 is the distance bin, and axis_2 is the phenotype.
    
    # Define color map if not provided.
    if not color_map:
        colors = px.colors.qualitative.Plotly
        unique_labels = sorted(labels_axis_2)
        color_map = {label: colors[i % len(colors)] for i, label in enumerate(unique_labels)}
    
    # Create the figure.
    fig = go.Figure()
    
    # For each phenotype, i.e., series (axis_2)...
    for i2, label_axis_2 in enumerate(labels_axis_2):

        # Get the current color for the series.
        color = color_map[label_axis_2]
        
        # If we want to plot a line plot with shaded area...
        extra_return_info = ""
        if plot_type == "line":

            # Store the quantiles for the 16-84% IQR shading.
            quantiles = np.quantile(data[:, :, i2], [0.16, 0.50, 0.84], axis=0)

            # Get the current color in RGB format.
            rgb = plotly.colors.hex_to_rgb(color)

            # Shaded area between min/max quantiles.
            fig.add_trace(
                go.Scatter(
                    x=list(labels_axis_1) + list(labels_axis_1)[::-1],
                    y=list(quantiles[2, :]) + list(quantiles[0, :])[::-1],
                    fill='toself',
                    fillcolor=f'rgba({rgb[0]}, {rgb[1]}, {rgb[2]}, 0.15)',  # lighter shade
                    line=dict(color='rgba(255,255,255,0)'),
                    hoverinfo="skip",
                    showlegend=False,
                    name=label_axis_2,
                    legendgroup=label_axis_2,
                )
            )

            # Median line.
            fig.add_trace(
                go.Scatter(
                    x=labels_axis_1,
                    y=quantiles[1, :],
                    mode='lines+markers',
                    name=label_axis_2,
                    line=dict(color=color, width=2),
                    legendgroup=label_axis_2,
                )
            )

            # Extra return information.
            extra_return_info = "*Median with 16-84% IQR shaded area."

        # If we want to plot a box plot...
        elif plot_type == 'box':
            for i1, label_axis_1 in enumerate(labels_axis_1):  # For each distance bin (axis_1)...
                fig.add_trace(
                    go.Box(
                        x=[label_axis_1] * data.shape[0],
                        y=data[:, i1, i2],
                        name=label_axis_2,
                        marker_color=color,
                        boxmean=True,
                        showlegend=(i1 == 0),  # Only show legend once per series
                        legendgroup=label_axis_2,
                    )
                )

        # If we want to plot a violin plot...
        elif plot_type == "violin":
            for i1, label_axis_1 in enumerate(labels_axis_1):  # For each distance bin (axis_1)...
                fig.add_trace(
                    go.Violin(
                        x=[label_axis_1] * data.shape[0],
                        y=data[:, i1, i2],
                        name=label_axis_2,
                        legendgroup=label_axis_2,
                        scalegroup=label_axis_2,
                        line_color=color,
                        showlegend=(i1 == 0),  # Only show legend for first distance bin
                        box_visible=False,
                        meanline_visible=True,
                    )
                )

    # Update layout for clarity.
    fig.update_layout(
        xaxis_title=axis_1_name,
        yaxis_title=value_name,
        legend_title=axis_2_name,
    )
    
    # Return the figure.
    return fig, extra_return_info


def add_new_label_column(lf, updates_pd, updates_index_column="sumap_cell_indices", main_index_column="sumap_cell_index", updates_label_column="label", new_label_column="neighborhood_type", keep="last", missing_label_value="Other"):
    # 0) Ensure join key dtype in lf is Int64 (only if needed).
    lf_schema = lf.collect_schema()
    if lf_schema[main_index_column] != pl.Int64:
        lf = lf.with_columns(pl.col(main_index_column).cast(pl.Int64))
 
    # 1) Convert pandas to Polars and explode to make a (index -> label) mapping.
    lf_updates = pl.from_pandas(updates_pd).lazy()
    mapping_lf = (
        lf_updates
        .explode(updates_index_column)                      # each index becomes a row
        .rename({updates_index_column: main_index_column}) # align column name
        .select(
            pl.col(main_index_column).cast(pl.Int64),
            pl.col(updates_label_column).cast(pl.Categorical).alias(new_label_column),
        )
        .unique(subset=[main_index_column], keep=keep)  # resolve duplicates if an index appears in multiple labels
    )

    # 2) Left-join onto main LazyFrame.
    lf2 = (
        lf.join(mapping_lf, on=main_index_column, how="left")
        .with_columns(
            pl.coalesce([pl.col(new_label_column), pl.lit(missing_label_value)])
                .alias(new_label_column)
        )
    )

    # 3) Return updated LazyFrame.
    return lf2


#### 6. First in plot_neighborhood_types.py (none yet) ########################################################
