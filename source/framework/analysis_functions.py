import time
import math
import os
import pickle
import dill


def run_analysis_job(function_name, inputs, job_dir):
    try:
        outputs_dir = os.path.join(job_dir, "outputs")  # This demonstrates that for a potentially asynchronous job that generates files, you should place the results in /tmp/multiplex_analysis_web_apps/job_data/<JOB_ID>/outputs specifically so the results are stored together with the worker output results in memory.
        if function_name == "find_primes_up_to":
            function_to_run = find_primes_up_to
        elif function_name == "run_phenograph_clustering":
            function_to_run = run_phenograph_clustering
        elif function_name == "run_neighb_clustering":
            function_to_run = run_neighb_clustering
        else:
            raise ValueError(f"Unknown function name: {function_name}")
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


def run_phenograph_clustering(adata_object_id, n_neighbors, clustering_algo, min_cluster_size, 
                             primary_metric, resolution_parameter, nn_method, random_seed, 
                             n_principal_components, n_jobs, n_iterations, fast, results_subdir, results_topdir):
    """
    Run phenograph clustering asynchronously.
    
    Parameters match those from RunPhenographClust function in Pheno_Cluster_a.py
    """
    import sys
    import os
    
    # Add the source directory to Python path to ensure pages2 can be imported
    source_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if source_dir not in sys.path:
        sys.path.insert(0, source_dir)
    
    from pages2.Pheno_Cluster_a import RunPhenographClust
    import framework.platform_abstraction as pa
    
    start_time = time.time()
    
    results_dir = os.path.join(results_topdir, results_subdir)
    os.makedirs(results_dir, exist_ok=True)
    
    # Download adata from object storage
    try:
        bucket_name = os.getenv('DATA_OBJECTS_BUCKET_NAME', 'objects')
        adata_buffer = pa.download_object_data(bucket_name, adata_object_id)
        
        # Load adata from buffer
        adata = pickle.loads(adata_buffer)
        
        # Clean up the object from storage after loading
        try:
            pa.delete_object_data(bucket_name, adata_object_id)
        except:
            pass  # Don't fail if cleanup fails
            
    except Exception as e:
        print(f"Failed to download adata from object storage: {e}")
        return None
    
    # Run the clustering
    adata_result = RunPhenographClust(
        adata=adata,
        n_neighbors=n_neighbors,
        clustering_algo=clustering_algo,
        min_cluster_size=min_cluster_size,
        primary_metric=primary_metric,
        resolution_parameter=resolution_parameter,
        nn_method=nn_method,
        random_seed=random_seed,
        n_principal_components=n_principal_components,
        n_jobs=n_jobs,
        n_iterations=n_iterations,
        fast=fast
    )
    
    end_time = time.time()
    duration = end_time - start_time
    
    # Save the result
    result_path = os.path.join(results_dir, "clustering_result.pkl")
    with open(result_path, 'wb') as f:
        pickle.dump(adata_result, f)
    
    # Save summary
    with open(os.path.join(results_dir, "clustering_summary.txt"), "w") as f:
        f.write(f"Phenograph clustering completed in {duration:.2f} seconds\n")
        f.write(f"Number of clusters: {len(adata_result.obs['Cluster'].unique())}\n")
        f.write(f"Parameters used:\n")
        f.write(f"  n_neighbors: {n_neighbors}\n")
        f.write(f"  clustering_algo: {clustering_algo}\n")
        f.write(f"  resolution: {resolution_parameter}\n")
    
    return {"adata_result": adata_result, "duration": duration, "result_path": result_path}


def run_neighb_clustering(adata_object_id, n_neighbors, metric, resolution, random_state, 
                         n_principal_components, n_jobs, n_iterations, fast, transformer, 
                         results_subdir, results_topdir):
    """
    Run scanpy neighbor clustering asynchronously.
    
    Parameters match those from RunNeighbClust function in Pheno_Cluster_a.py
    """
    import sys
    import os
    
    # Add the source directory to Python path to ensure pages2 can be imported
    source_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if source_dir not in sys.path:
        sys.path.insert(0, source_dir)
    
    from pages2.Pheno_Cluster_a import RunNeighbClust
    import framework.platform_abstraction as pa
    
    start_time = time.time()
    
    results_dir = os.path.join(results_topdir, results_subdir)
    os.makedirs(results_dir, exist_ok=True)
    
    # Download adata from object storage
    try:
        bucket_name = os.getenv('DATA_OBJECTS_BUCKET_NAME', 'objects')
        adata_buffer = pa.download_object_data(bucket_name, adata_object_id)
        
        # Load adata from buffer
        adata = pickle.loads(adata_buffer)
        
        # Clean up the object from storage after loading
        try:
            pa.delete_object_data(bucket_name, adata_object_id)
        except:
            pass  # Don't fail if cleanup fails
            
    except Exception as e:
        print(f"Failed to download adata from object storage: {e}")
        return None
    
    # Run the clustering
    adata_result = RunNeighbClust(
        adata=adata,
        n_neighbors=n_neighbors,
        metric=metric,
        resolution=resolution,
        random_state=random_state,
        n_principal_components=n_principal_components,
        n_jobs=n_jobs,
        n_iterations=n_iterations,
        fast=fast,
        transformer=transformer
    )
    
    end_time = time.time()
    duration = end_time - start_time
    
    # Save the result
    result_path = os.path.join(results_dir, "clustering_result.pkl")
    with open(result_path, 'wb') as f:
        pickle.dump(adata_result, f)
    
    # Save summary
    with open(os.path.join(results_dir, "clustering_summary.txt"), "w") as f:
        f.write(f"Scanpy clustering completed in {duration:.2f} seconds\n")
        f.write(f"Number of clusters: {len(adata_result.obs['Cluster'].unique())}\n")
        f.write(f"Parameters used:\n")
        f.write(f"  n_neighbors: {n_neighbors}\n")
        f.write(f"  metric: {metric}\n")
        f.write(f"  resolution: {resolution}\n")
    
    return {"adata_result": adata_result, "duration": duration, "result_path": result_path}
