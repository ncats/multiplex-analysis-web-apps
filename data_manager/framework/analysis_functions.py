import time
import math
import os
import basic_phenotyper_lib as bpl  # Useful functions for phenotyping collections of cells


def run_analysis_job(function_name, inputs, job_dir):
    try:
        outputs_dir = os.path.join(job_dir, "outputs")  # This demonstrates that for a potentially asynchronous job that generates files, you should place the results in /tmp/multiplex_analysis_web_apps/job_data/<JOB_ID>/outputs specifically so the results are stored together with the worker output results in memory.
        if function_name == "find_primes_up_to":
            function_to_run = find_primes_up_to
        elif function_name == "init_spatial_umap":
            function_to_run = init_spatial_umap
        outputs = function_to_run(**inputs, results_topdir=outputs_dir)
        return outputs
    except Exception as e:
        print(f"Error occurred while running analysis job {function_name}: {e}")
        return None


def find_primes_up_to(limit, results_subdir, results_topdir):
    """
    Find all prime numbers up to a given limit using trial division.
    Returns the list of primes and timing information.

    On my laptop this takes about 6-8 seconds: primes, duration = find_primes_up_to(4000000).
    """
    start_time = time.time()

    if limit < 2:
        return [], 0

    primes = []

    for num in range(2, limit + 1):
        is_prime = True

        # Check if num is prime by testing divisibility
        for i in range(2, int(math.sqrt(num)) + 1):
            if num % i == 0:
                is_prime = False
                break

        if is_prime:
            primes.append(num)

    end_time = time.time()
    duration = end_time - start_time

    results_dir = os.path.join(results_topdir, results_subdir)

    # Create results directory if it doesn't exist.
    os.makedirs(results_dir, exist_ok=True)

    # Save results to text files
    with open(os.path.join(results_dir, "primes.txt"), "w") as f:
        f.write(f"Found {len(primes)} primes up to {limit} in {duration:.2f} seconds\n")
        f.write(f"First 10 primes: {primes[:10]}\n")
        f.write(f"Last 10 primes: {primes[-10:]}\n")

    return {"primes": primes, "duration": duration}


def init_spatial_umap(calc_unique_areas_toggle, area_filter_per, df, marker_multi_sel, phenoOrder, datafile_min_img_size, cpu_pool_size, results_topdir):
    '''
    Initalizing the spatial_umap object
    '''

    # Andrew moved this to its own button in Neighborhood_Profiles.py; search for "Initialize Neighborhood Profiles" there.
    # Reset the settings required for Neighborhood Analysis
    # st.session_state = ndl.reset_neigh_profile_settings(st.session_state)

    if not calc_unique_areas_toggle:
        area_filter = 0
    else:
        area_filter = area_filter_per

    # st.session_state.bc.startTimer()
    # with st.spinner('Calculating Cell Counts and Areas', show_time=True):
    spatial_umap = bpl.setup_Spatial_UMAP(df = df,
                                                            marker_names = marker_multi_sel,
                                                            pheno_order = phenoOrder,
                                                            smallest_image_size = datafile_min_img_size)

    spatial_umap = bpl.perform_density_calc(spatial_umap,
                                                            #  st.session_state.bc,
                                                                calc_unique_areas_toggle,
                                                                cpu_pool_size,
                                                                area_threshold = area_filter)

    # Record time elapsed
    # st.session_state.bc.set_value_df('time_to_run_counts', st.session_state.bc.elapsedTime())

    return {"spatial_umap": spatial_umap, "density_completed": True}

    # Save checkpoint for Neighborhood Profile structure
    # save_neipro_struct()
