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


def fast_neighbors_counts_for_block2(df_image, image_name, coord_column_names, phenotypes, radii, phenotype_column_name, max_chunk_size_in_mb=200):
    # A block can be an image, ROI, etc. It's the entity over which it makes sense to calculate the neighbors of centers. Here, we're assuming it's an image, but in the SIT for e.g., we generally want it to refer to a ROI.

    # max_chunk_size_in_mb=200, for a 100K-cell dataset, will yield about 250-row chunks, which will yield about 400 chunks i.e. center KDTrees
    # max_chunk_size_in_mb=500, for a 100K-cell dataset, will yield about 650-row chunks, which will yield about 150 chunks i.e. center KDTrees
    # max_chunk_size_in_mb=1000, for a 100K-cell dataset, will yield about 1300-row chunks, which will yield about 80 chunks i.e. center KDTrees
    # max_chunk_size_in_mb=2000, for a 100K-cell dataset, will yield about 2600-row chunks, which will yield about 40 chunks i.e. center KDTrees
    # max_chunk_size_in_mb=5000, for a 100K-cell dataset, will yield about 6600-row chunks, which will yield about 15 chunks i.e. center KDTrees

    # Print the image name
    print(f'Calculating neighbor counts for image {image_name} ({len(df_image)} cells) using the new kdtree method...')

    # Record the start time
    start_time = time.time()

    # Store some properties of the input dataframe
    df_image_index = df_image.index
    num_cells = len(df_image)

    # Get the number of rows per chunk based on the maximum chunk size in MB and assuming every cell will be counted as a neighbor around every center (i.e., use upper limits)
    largest_sublist_size_in_mb = num_cells * 8  / 1024 ** 2  # 8 bytes per int64 --> units = mb / row
    num_rows_per_chunk = int(max_chunk_size_in_mb / largest_sublist_size_in_mb)

    # Get the corresponding integer indices to index things like the input image dataframe
    start_indices = np.arange(0, num_cells, num_rows_per_chunk)
    stop_indices = start_indices + num_rows_per_chunk

    # Initialize a list to hold the dataframes of neighbor counts for each radius (not each radius range)
    df_counts_holder = [pd.DataFrame(0, index=phenotypes, columns=df_image_index) for _ in radii]

    # Pre-calculate the neighbor tree for each phenotype
    neighbor_trees = []
    for neighbor_phenotype in phenotypes:

        # Get the boolean series identifying the current neighbor phenotype
        ser_curr_neighbor_phenotype = df_image[phenotype_column_name] == neighbor_phenotype

        # Construct the KDTree for the current phenotype in the entire current image. This represents the neighbors
        neighbor_trees.append(scipy.spatial.KDTree(df_image.loc[ser_curr_neighbor_phenotype, coord_column_names]))

    # For each chunk of centers...
    for start_index, stop_index in zip(start_indices, stop_indices):

        # Construct the KDTree for the current chunk. This represents the centers
        center_tree = scipy.spatial.KDTree(df_image.loc[df_image_index[start_index:stop_index], coord_column_names])

        # For each neighbor tree...
        for ineighbor_phenotype, curr_neighbor_tree in enumerate(neighbor_trees):

            # For each radius, which should be monotonically increasing and start with 0...
            for iradius, radius in enumerate(radii):

                # Get the list of lists containing the indices of the neighbors for each center
                neighbors_for_radius = center_tree.query_ball_tree(curr_neighbor_tree, radius)

                # In the correct dataframe (corresponding to the current radius), set the counts of neighbors (of the current phenotype) for each center
                df_counts_holder[iradius].iloc[ineighbor_phenotype, start_index:stop_index] = [len(neighbors_for_center) for neighbors_for_center in neighbors_for_radius]

    # For each annulus, i.e., each radius range...
    df_counts_holder_annulus = []
    for iradius in range(len(radii) - 1):

        # Get the counts of neighbors in the current annulus
        df_counts_curr_annulus = df_counts_holder[iradius + 1] - df_counts_holder[iradius]

        # Rename the index to reflect the current radius range
        radius_range_str = f'({radii[iradius]}, {radii[iradius + 1]}]'
        df_counts_curr_annulus.index = [f'{phenotype} in {radius_range_str}' for phenotype in phenotypes]

        # Add a transpose of this (so centers are in rows and phenotypes/radii are in columns) to the running list of annulus dataframes
        df_counts_holder_annulus.append(df_counts_curr_annulus.T)

    # Concatenate the annulus dataframes to get the final dataframe of neighbor counts for the current image
    df_curr_counts = pd.concat(df_counts_holder_annulus, axis='columns')

    # Convert the dataframe to a more efficient dtype (from int64)
    df_curr_counts = df_curr_counts.astype(np.int32)

    # Print the time taken to calculate the neighbor counts for the current image
    print(f'  ...finished calculating neighbor counts for image {image_name} ({len(df_image)} cells) in {time.time() - start_time:.2f} seconds')

    # Return the final dataframe of neighbor counts for the current image
    return df_curr_counts


# Format the lazyframe, optionally sample it, and convert to a polars dataframe.
def format_lazyframe(lf, sample_size=None, sample_seed=42):

    # Load in cells and patient data.
    lf = (
        lf
        .rename({"Image ID_(standardized)": "TMA_core_id", "Centroid X (µm)_(standardized)": "Xcor", "Centroid Y (µm)_(standardized)": "Ycor", "label": "Lineage"})
        .select(pl.col(["TMA_core_id", "Xcor", "Ycor", "Lineage", "input_index"]))
        )

    # Load in cells and patient data. Sampling will aid in faster testing and development. The sorting after the sampling is crucial to ensure consistent ordering.
    if sample_size is None:
        pldf = (
            lf
            .sort(by="TMA_core_id")
            .collect()
            )
    else:
        pldf = (
            lf
            .collect()
            .sample(n=sample_size, seed=sample_seed)
            .sort(by="TMA_core_id")
            )
        
    return dict(pldf=pldf)


def save_and_load_pandas_df_to_lf(pd_df, handle, file_format, topdir, index_column_name=None):
    save_pandas_df_to_file(pd_df, handle=handle, file_format=file_format, topdir=topdir, subdir="input")
    lf = get_lf(handle, topdir=topdir, subdir="input", file_format=file_format)
    if index_column_name is not None:
        lf = lf.with_row_index(name=index_column_name)
    return lf


def get_true_false_color_map():
    colors = px.colors.qualitative.Plotly
    color_map = {label: colors[i % len(colors)] for i, label in enumerate((True, False))}
    return color_map


# Get the marker column names from the lazyframe.
def get_marker_columns(lf, exclusion_suffix=""):
    if exclusion_suffix:
        marker_columns = [column for column in lf.collect_schema().names() if column.startswith("Phenotype_(standardized) ") and not column.endswith(exclusion_suffix)]
    else:
        marker_columns = [column for column in lf.collect_schema().names() if column.startswith("Phenotype_(standardized) ")]
    return sorted(marker_columns)


def print_flush(msg):
    print(msg, flush=True)


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
            with zipfile.ZipFile(os.path.join(input_dir, full_filename), 'r') as zip_ref:
                zip_ref.extractall(input_dir)
            base_name = full_filename.removesuffix(".zip")
            os.remove(os.path.join(input_dir, full_filename))
            unzipped_paths.append(os.path.join(input_dir, base_name))

        # Generate an intermediate file from which to load the lazyframe.
        for unzipped_path in unzipped_paths:
            filename = os.path.basename(unzipped_path)
            filepath = subset_csv_to_file(csv_filename=filename, handle="unified_input_file", topdir=topdir, subdir="input", file_format=file_format)
            os.remove(unzipped_path)

        # Load the lazyframe.
        lf = get_lf("unified_input_file", topdir=topdir, subdir="input", file_format=file_format)

        # Store the local filepath relative to the session directory.
        local_filepath = filepath.removeprefix(topdir + os.sep)

        # Set the index.
        lf = lf.with_row_index(name="input_index")

        # Return the lazyframe.
        return lf
    
    except Exception as e:
        framework_utils.multiprint(f"Unable to load unified input file data (inputs: file_format={file_format}, db_schema={db_schema}, bucket_name={bucket_name}, object_filename={object_filename}): {e}", (print,))
        raise


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
        raise


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
        raise


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
        raise


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
        raise


def plot_image_from_frame(
    frame: pl.LazyFrame | pl.DataFrame | pd.DataFrame,
    image_colname: str = "Image ID_(standardized)",
    selected_images: list[str] = ["MS_02__cele_1400w"],
    marker_size: int | float | None = None,
    xcol="Centroid X (µm)_(standardized)",
    ycol="Centroid Y (µm)_(standardized)",
    color_col="label",
    color_map: dict | None = None,
    custom_columns=[],
    highlight_indices: list | set = (),
    default_marker_size: int | float | None = None,
    default_marker_line_color: str = "rgba(0,0,0,0.3)",
    default_marker_alpha: float = 0.5,
    highlight_marker_size: int | float | None = None,
    highlight_marker_line_color: str = "black",
    highlight_marker_alpha: float = 1.0,
    highlight_show_legend: bool = False,
):
    try:

        # Efficiently convert the input frame to a pandas DataFrame with necessary filtering.
        # If there is slowness, we can try keeping as a Polars DataFrame and using Plotly's ability to plot from Polars DataFrames directly, proceeding with polars dataframes in all operations below.
        cols_to_keep = [color_col, image_colname, xcol, ycol] + custom_columns
        if isinstance(frame, pl.LazyFrame):
            base = frame.filter(pl.col(xcol).is_not_null() & pl.col(ycol).is_not_null())
            if selected_images:
                base = base.filter(pl.col(image_colname).is_in(selected_images))
            df = base.select(cols_to_keep).collect().to_pandas()
        elif isinstance(frame, pl.DataFrame):
            base = frame.filter(pl.col(xcol).is_not_null() & pl.col(ycol).is_not_null())
            if selected_images:
                base = base.filter(pl.col(image_colname).is_in(selected_images))
            df = base.select(cols_to_keep).to_pandas()
        elif isinstance(frame, pd.DataFrame):
            mask = frame[xcol].notna() & frame[ycol].notna()
            if selected_images:
                mask &= frame[image_colname].isin(selected_images)
            df = frame.loc[mask, cols_to_keep]
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
        index_col = "index"
        do_highlight = index_col in custom_columns and len(highlight_indices) > 0
        if do_highlight:
            mask_high = df[index_col].isin(highlight_indices)
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

        # Determine effective marker sizes.
        effective_default_size = default_marker_size if default_marker_size is not None else (marker_size if marker_size is not None else 5)
        effective_highlight_size = highlight_marker_size if highlight_marker_size is not None else (3 * marker_size if marker_size is not None else 15)

        # Create figure and add non-highlighted traces.
        fig = go.Figure()
        for label in unique_labels:
            label_mask = df_other[color_col] == label
            if label_mask.any():
                label_data = df_other[label_mask]
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

        # Add highlighted traces if applicable.
        if do_highlight and df_high is not None and not df_high.empty:
            for label in unique_labels:
                label_mask = df_high[color_col] == label
                if label_mask.any():
                    label_data = df_high[label_mask]
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
        print_flush(f"An error occurred in function {os.path.basename(__file__) if '__file__' in globals() else '<interactive>'}.{plot_image_from_frame.__name__}: {e}")
        raise


# Function to create a violin plot with distributions from a third axis.
def violin_plot_with_series(data, labels_axis_1, labels_axis_2, axis_1_name="Distance bin (µm)", axis_2_name="Phenotype", value_name="Density", color_map=None):
    # Note the axes here refer to the axes of data, not the plot axes.
    # data shape: (axis_0, axis_1, axis_2) where axis_0 is the distribution dimension (e.g., cells), axis_1 is the distance bin, and axis_2 is the phenotype.
    
    # Define color map if not provided.
    if not color_map:
        colors = px.colors.qualitative.Plotly
        unique_labels = sorted(labels_axis_2)
        color_map = {label: colors[i % len(colors)] for i, label in enumerate(unique_labels)}
    
    # Create the figure.
    fig = go.Figure()
    
    # For each phenotype (axis_2)...
    for i2, label_axis_2 in enumerate(labels_axis_2):
        # Store means for this phenotype to connect with a line
        means = []
        
        # For each distance bin (axis_1)...
        for i1, label_axis_1 in enumerate(labels_axis_1):
            # Extract the distribution across axis_0 for this combination
            distribution_values = data[:, i1, i2]
            
            # Calculate and store mean
            means.append(np.mean(distribution_values))
            
            # Add violin trace
            fig.add_trace(go.Violin(
                x=[label_axis_1] * data.shape[0],
                y=distribution_values,
                name=label_axis_2,
                legendgroup=label_axis_2,
                scalegroup=label_axis_2,
                line_color=color_map[label_axis_2],
                showlegend=(i1 == 0),  # Only show legend for first distance bin
                box_visible=False,
                meanline_visible=True,
            ))
        
        # Add line trace connecting the means
        fig.add_trace(go.Scatter(
            x=labels_axis_1,
            y=means,
            mode='lines+markers',
            line=dict(color=color_map[label_axis_2], width=2),
            marker=dict(size=8, color=color_map[label_axis_2]),
            name=label_axis_2,
            legendgroup=label_axis_2,
            showlegend=False,  # Already shown in violin legend
        ))
    
    # Update layout for clarity.
    fig.update_layout(
        xaxis_title=axis_1_name,
        yaxis_title=value_name,
        legend_title=axis_2_name,
        violinmode='overlay',
    )
    
    # Return the figure.
    return fig


# Function to create a line plot with multiple series.
def line_plot_with_series(data, labels_axis_0, labels_axis_1, axis_0_name="Distance bin (µm)", axis_1_name="Phenotype", value_name="Mean density", color_map=None):
    # Note the axes here refer to the axes of data, not the plot axes.

    # Create a DataFrame from the data.
    data = pl.DataFrame(data, schema=labels_axis_1).with_columns(
        pl.Series(axis_0_name, labels_axis_0)
    )

    # Unpivot the DataFrame for plotting.
    data_unpivoted = data.unpivot(
        index=axis_0_name,
        on=labels_axis_1,
        variable_name=axis_1_name,
        value_name=value_name,
    )

    # Define color map if not provided.
    if not color_map:
        colors = px.colors.qualitative.Plotly
        unique_labels = sorted(labels_axis_1)
        color_map = {label: colors[i % len(colors)] for i, label in enumerate(unique_labels)}

    # Create the line plot.
    fig = px.line(
        data_unpivoted,
        x=axis_0_name,
        y=value_name,
        color=axis_1_name,
        color_discrete_map=color_map,
        markers=True,
    )

    # Update layout for clarity.
    fig.update_layout(
        xaxis_title=axis_0_name,
        yaxis_title=value_name,
        legend_title=axis_1_name,
    )

    # Return the figure.
    return fig


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


# Wrapper to modify core function for compatibility with job input/output dictionary standards without the need for an intermediate polars dataframe.
def generate_umap_wrapper(**inputs):

    # Modify inputs as needed (standard format, i.e., what's expected of a dictionary as in utils.deserialize_binary_files_to_dictionary() should come in).
    lf = inputs["LAZYFRAMES"]["marker_phenotyping"]["lf"]
    lf = (
        lf
        .rename({"Image ID_(standardized)": "TMA_core_id", "Centroid X (µm)_(standardized)": "Xcor", "Centroid Y (µm)_(standardized)": "Ycor", "label": "Lineage"})
        .select(pl.col(["TMA_core_id", "Xcor", "Ycor", "Lineage", "input_index"]))
        .sort(by="TMA_core_id")
        )
    inputs["lf"] = lf
    del inputs["LAZYFRAMES"]

    # Run the core function.
    outputs = generate_umap_lf_input(**inputs)  # what comes out of this: dict(spatial_umap=spatial_umap, complete_success=True)

    # Return the properly formatted outputs.
    return outputs


# def generate_umap(pldf, unique_labels, dist_bin_um_list=[25, 50, 100, 150, 200], area_downsample=0.2, um_per_px=1, cpu_pool_size=None, results_topdir=".", subdir="results", counts_method="andrew", area_threshold=0.8, custom_areas=True, seed_for_train_test_split=54321, n=2500, keep_images_with_too_little_data=True, train_sample_frac=1.0, test_sample_frac=1.0, de_min_coords=True, mp_start_method=None):
#     # Note that cpu_pool_size=None will default to the number of available CPUs.

#     # Instantiate the spatial umap object.
#     spatial_umap = SpatialUMAP.SpatialUMAP(dist_bin_um=np.array(dist_bin_um_list), um_per_px=um_per_px, area_downsample=area_downsample)

#     # Normalize coordinates to start at (0,0) for each image.
#     if custom_areas and de_min_coords:

#         # Print the minima first.
#         min_coords_df = pldf.group_by("TMA_core_id").agg([
#             pl.min("Xcor").alias("min_Xcor"),
#             pl.min("Ycor").alias("min_Ycor"),
#         ])
#         print_flush("Minimum coordinates per TMA core:")
#         print_flush(min_coords_df)

#         # Perform normalization.
#         pldf = pldf.with_columns([
#             (pl.col("Xcor") - pl.min("Xcor").over("TMA_core_id")).alias("Xcor"),
#             (pl.col("Ycor") - pl.min("Ycor").over("TMA_core_id")).alias("Ycor"),
#         ])

#     # I don't want to, but convert to pandas for compatibility with SpatialUMAP class.
#     spatial_umap.cells = pldf.to_pandas()

#     # Set explicitly as numpy array the cell coordinates (x, y).
#     spatial_umap.cell_positions = spatial_umap.cells[['Xcor', 'Ycor']].values

#     # Set explicitly as one hot data frame the cell labels.
#     spatial_umap.cell_labels = pd.get_dummies(spatial_umap.cells['Lineage'])

#     # Set the region is to be analyzed (a TMA core is treated similar to a region of a interest).
#     spatial_umap.region_ids = spatial_umap.cells.TMA_core_id.unique()

#     # Clear metrics.
#     spatial_umap.clear_counts()
#     spatial_umap.clear_areas()

#     # Save the unique labels/lineages/species. This is only used for Andrew's counting method.
#     spatial_umap.species = unique_labels

#     # Ensure results directory exists.
#     os.makedirs(os.path.join(results_topdir, subdir), exist_ok=True)

#     # Determine which counting method to use based on the input parameter. Note I have shown that the methods precisely agree.
#     if counts_method == "andrew":
#         print_flush("Using Andrew's counts method.")
#         spatial_umap.get_counts_And(cpu_pool_size=cpu_pool_size, mp_start_method=mp_start_method)
#     else:
#         print_flush("Using Baras' counts method.")
#         spatial_umap.get_counts(cpu_pool_size=cpu_pool_size)

#     # Get the areas of cells and save to pickle file.
#     if custom_areas:
#         try:
#             spatial_umap.get_areas(area_threshold, pool_size=cpu_pool_size, save_file=os.path.join(results_topdir, subdir, f"areas.csv"), plots_directory=os.path.join(results_topdir, subdir))  # Sets spatial_umap.cells["area_filter"] and spatial_umap.areas.
#         except Exception as e:
#             print_flush(f"An error occurred while calculating custom areas, potentially in SpatialUMAP.FitEllipse.fit() in \"hull = ConvexHull(d[idx_fit])\": {e}")
#             raise
#     else:
#         # Keep in mind areas in the Baras code seem to be in units of pixels squared.
#         spatial_umap.cells["area_filter"] = True  # If not using custom areas, set all cells to pass the area filter.
#         r0 = np.concatenate(([0], spatial_umap.dist_bin_px))
#         # Note not including downsampling at all!
#         spatial_umap.areas = (np.pi * (r0[1:] ** 2 - r0[:-1] ** 2)).reshape((1, len(dist_bin_um_list), 1))

#     # calculate density base on counts of cells / area of each arc examine
#     if custom_areas:
#         spatial_umap.density = np.empty(spatial_umap.counts.shape)
#         spatial_umap.cells['area_filter'] = ((spatial_umap.areas / spatial_umap.arcs_masks.sum(axis=(0, 1))[np.newaxis, ...]) > area_threshold).all(axis=1)
#         spatial_umap.density[spatial_umap.cells['area_filter'].values] = spatial_umap.counts[spatial_umap.cells['area_filter'].values] / spatial_umap.areas[spatial_umap.cells['area_filter'].values][..., np.newaxis]
#     else:
#         spatial_umap.density = spatial_umap.counts / spatial_umap.areas  # We're not doing this yet (so units are currently in #/px), but we divide by square of um_per_px to get density in units of # per square micron.

#     # Output the percentage of the dataset that has been filtered out due to area filtering.
#     num_total_cells = spatial_umap.cells.shape[0]
#     num_kept_cells = spatial_umap.cells["area_filter"].sum()
#     print_flush(f"{100 * (1 - num_kept_cells / num_total_cells):.1f}% of the cells have been filtered out due to area filtering.")

#     # Check for dropped images due to insufficient cells passing area filter.
#     original_images = set(spatial_umap.region_ids)
#     remaining_images = set(spatial_umap.cells[spatial_umap.cells["area_filter"]]["TMA_core_id"].unique())
#     dropped_images = original_images - remaining_images
#     if dropped_images:
#         print_flush(f"WARNING:")
#         print_flush(f"  The following images were dropped (i.e., missing scatter plots) due to no cells passing the area filter: {sorted(list(dropped_images))}.")
#         print_flush(f"  These images remain: {sorted(list(remaining_images))}.")
#         if remaining_images:
#             print_flush(f"  Keep in mind that just because some images may not have been dropped, significant numbers of cells in those images may have been dropped.")
#         else:
#             print_flush(f"  Since no images remain after area filtering, we are aborting UMAP generation.")
#             return spatial_umap, False
#         remaining_loc = spatial_umap.cells["TMA_core_id"].isin(remaining_images)
#         spatial_umap.cells = spatial_umap.cells[remaining_loc]
#         spatial_umap.density = spatial_umap.density[remaining_loc.values]

#     min_filtered_cells = get_min_positive_values(spatial_umap.cells, group_col="TMA_core_id", boolean_column="area_filter")  # this would be zero if we didn't do the filtering-out line above (spatial_umap.cells = ...)
#     print_flush(f"Minimum number of cells passing area filter across all TMA cores: {min_filtered_cells}")

#     # Set training and "test" cells for umap training and embedding, respectively. Baras's original code had a hard cutoff of n=2500 so images with fewer than 2*2500 non-filtered-out cells were discarded entirely. n = min(n, min_filtered_cells // 2) allows these images to remain in the analysis with smaller n.
#     if keep_images_with_too_little_data:
#         n = min(n, min_filtered_cells // 2)

#     print_flush(f"Using n_train={int(train_sample_frac*n)} cells per image for UMAP training and n_test={int(test_sample_frac*n)} for testing.")
#     spatial_umap.set_train_test(n=n, seed=seed_for_train_test_split, train_sample_frac=train_sample_frac, test_sample_frac=test_sample_frac)

#     # Fit umap on training cells.
#     spatial_umap.umap_fit = umap.UMAP().fit(spatial_umap.density[spatial_umap.cells['umap_train'].values].reshape((spatial_umap.cells['umap_train'].sum(), -1)))
#     print_flush(spatial_umap.umap_fit.embedding_.shape)

#     # Apply umap embedding on test cells.
#     spatial_umap.umap_test = spatial_umap.umap_fit.transform(spatial_umap.density[spatial_umap.cells['umap_test'].values].reshape((spatial_umap.cells['umap_test'].sum(), -1)))
#     print_flush(spatial_umap.umap_test.shape)

#     # Save the UMAP coordinates back to the cells dataframe.
#     spatial_umap.cells[["umap_1", "umap_2"]] = np.nan
#     spatial_umap.cells.loc[spatial_umap.cells['umap_test'].values, ["umap_1", "umap_2"]] = spatial_umap.umap_test

#     # # save spatial_umap object as pickle
#     # pickle.dump(spatial_umap, open(data_dir + '/pkl/spatial_umap.pkl', 'wb'))

#     # Return the result as a dictionary.
#     return dict(spatial_umap=spatial_umap, complete_success=True)


def generate_umap_lf_input(lf, unique_labels, dist_bin_um_list=[25, 50, 100, 150, 200], area_downsample=0.2, um_per_px=1, cpu_pool_size=None, results_topdir=".", subdir="results", counts_method="andrew", area_threshold=0.8, custom_areas=True, seed_for_train_test_split=54321, n=2500, keep_images_with_too_little_data=True, train_sample_frac=1.0, test_sample_frac=1.0, de_min_coords=True, mp_start_method=None):
    # Note that cpu_pool_size=None will default to the number of available CPUs.

    # Instantiate the spatial umap object.
    spatial_umap = SpatialUMAP.SpatialUMAP(dist_bin_um=np.array(dist_bin_um_list), um_per_px=um_per_px, area_downsample=area_downsample)

    # Normalize coordinates to start at (0,0) for each image.
    if custom_areas and de_min_coords:

        # Print the minima first.
        min_coords_df = lf.group_by("TMA_core_id").agg([
            pl.min("Xcor").alias("min_Xcor"),
            pl.min("Ycor").alias("min_Ycor"),
        ]).collect()
        print_flush("Minimum coordinates per TMA core:")
        print_flush(min_coords_df)

        # Perform normalization.
        lf = lf.with_columns([
            (pl.col("Xcor") - pl.min("Xcor").over("TMA_core_id")).alias("Xcor"),
            (pl.col("Ycor") - pl.min("Ycor").over("TMA_core_id")).alias("Ycor"),
        ])

        #### Note I may need to collect_schema() here potentially to register new columns per my experience needing that with with_row_index(). ####

    # I don't want to, but convert to pandas for compatibility with SpatialUMAP class.
    spatial_umap.cells = lf.collect().to_pandas()

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
        print_flush("Using Andrew's counts method.")
        spatial_umap.get_counts_And(cpu_pool_size=cpu_pool_size, mp_start_method=mp_start_method)
    else:
        print_flush("Using Baras' counts method.")
        spatial_umap.get_counts(cpu_pool_size=cpu_pool_size)

    # Get the areas of cells and save to pickle file.
    if custom_areas:
        try:
            spatial_umap.get_areas(area_threshold, pool_size=cpu_pool_size, save_file=os.path.join(results_topdir, subdir, f"areas.csv"), plots_directory=os.path.join(results_topdir, subdir))  # Sets spatial_umap.cells["area_filter"] and spatial_umap.areas.
        except Exception as e:
            print_flush(f"An error occurred while calculating custom areas, potentially in SpatialUMAP.FitEllipse.fit() in \"hull = ConvexHull(d[idx_fit])\": {e}")
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

    # Return the result as a dictionary.
    return dict(spatial_umap=spatial_umap, complete_success=True)


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
