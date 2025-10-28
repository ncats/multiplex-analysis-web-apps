import time
import math
import os
import sys


def run_analysis_job(function_name, inputs, job_dir):
    try:
        outputs_dir = os.path.join(job_dir, "outputs")  # This demonstrates that for a potentially asynchronous job that generates files, you should place the results in /tmp/multiplex_analysis_web_apps/job_data/<JOB_ID>/outputs specifically so the results are stored together with the worker output results in memory.
        if function_name == "find_primes_up_to":
            function_to_run = find_primes_up_to
            outputs = function_to_run(**inputs, results_topdir=outputs_dir)
        elif function_name == "init_spatial_umap_analysis":
            function_to_run = init_spatial_umap_analysis
            # Add results_topdir to inputs for this function
            inputs_with_topdir = inputs.copy()
            inputs_with_topdir["results_topdir"] = outputs_dir
            outputs = function_to_run(**inputs_with_topdir)
        else:
            raise ValueError(f"Unknown function name: {function_name}")
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


def init_spatial_umap_analysis(df, marker_names, pheno_order, smallest_image_size, 
                              calc_unique_areas_toggle, cpu_pool_size, area_threshold, 
                              results_subdir, results_topdir):
    """
    Initialize spatial UMAP analysis with cell counts and density calculations.
    This function wraps the spatial UMAP initialization for asynchronous execution.
    
    Args:
        df: Input dataframe with cell data
        marker_names: List of marker names to use
        pheno_order: Phenotype order
        smallest_image_size: Smallest image size in the dataset
        calc_unique_areas_toggle: Whether to calculate unique areas
        cpu_pool_size: Number of CPUs to use for parallel processing
        area_threshold: Area filter threshold
        results_subdir: Subdirectory for results
        results_topdir: Top directory for results
        
    Returns:
        Dictionary containing the spatial_umap object and timing information
    """
    # Import modules inside the function to avoid Streamlit context issues
    import sys
    import os
    
    # Add the source directory to Python path
    source_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if source_dir not in sys.path:
        sys.path.insert(0, source_dir)
    
    # Import the required modules
    import basic_phenotyper_lib as bpl
    from benchmark_collector import benchmark_collector
    
    print("DEBUG: init_spatial_umap_analysis function called!")
    print(f"DEBUG: Function arguments received")
    
    start_time = time.time()
    
    # Create a benchmark collector for timing
    print("DEBUG: About to create benchmark_collector")
    bc = benchmark_collector()
    print("DEBUG: benchmark_collector created successfully")
    
    try:
        # Debug: Check all parameters before calling functions
        print(f"DEBUG: df type = {type(df)}, shape = {df.shape if hasattr(df, 'shape') else 'No shape'}")
        print(f"DEBUG: marker_names = {marker_names}, type = {type(marker_names)}")
        print(f"DEBUG: pheno_order = {pheno_order}, type = {type(pheno_order)}")
        print(f"DEBUG: smallest_image_size = {smallest_image_size}, type = {type(smallest_image_size)}")
        print(f"DEBUG: calc_unique_areas_toggle = {calc_unique_areas_toggle}, type = {type(calc_unique_areas_toggle)}")
        print(f"DEBUG: cpu_pool_size = {cpu_pool_size}, type = {type(cpu_pool_size)}")
        print(f"DEBUG: area_threshold = {area_threshold}, type = {type(area_threshold)}")
        print(f"DEBUG: bc = {bc}, type = {type(bc)}")
        
        # Setup Spatial UMAP object
        print("DEBUG: About to call bpl.setup_Spatial_UMAP")
        bc.startTimer()
        spatial_umap = bpl.setup_Spatial_UMAP(df=df,
                                            marker_names=marker_names,
                                            pheno_order=pheno_order,
                                            smallest_image_size=smallest_image_size)
        print("DEBUG: setup_Spatial_UMAP completed successfully")
        
        # Perform density calculations
        print("DEBUG: About to call bpl.perform_density_calc")
        spatial_umap = bpl.perform_density_calc(spatial_umap,
                                              bc,
                                              calc_unique_areas_toggle,
                                              cpu_pool_size,
                                              area_threshold=area_threshold)
        print("DEBUG: perform_density_calc completed successfully")
        
        elapsed_time = bc.elapsedTime()
        
        # Create results directory if it doesn't exist
        results_dir = os.path.join(results_topdir, results_subdir)
        os.makedirs(results_dir, exist_ok=True)
        
        # Save timing information to a file
        with open(os.path.join(results_dir, "spatial_umap_timing.txt"), "w") as f:
            f.write(f"Spatial UMAP initialization completed in {elapsed_time:.2f} seconds\n")
            f.write(f"Number of cells processed: {len(spatial_umap.cells)}\n")
            f.write(f"Calc unique areas: {calc_unique_areas_toggle}\n")
            f.write(f"CPU pool size: {cpu_pool_size}\n")
        
        end_time = time.time()
        total_duration = end_time - start_time
        
        return {
            "spatial_umap": spatial_umap,
            "elapsed_time": elapsed_time,
            "total_duration": total_duration
        }
        
    except Exception as e:
        print(f"Error in init_spatial_umap_analysis: {e}")
        raise e
