# Async Clustering Implementation Summary

## Overview
Added asynchronous execution capability for phenograph and scanpy clustering methods in the Pheno Cluster page. This allows long-running clustering operations to execute in background worker containers while keeping the UI responsive.

## Files Modified

### 1. `/source/framework/analysis_functions.py`
**Purpose**: Added worker functions for async clustering execution

**Changes**:
- Added `run_phenograph_clustering()` function - executes phenograph clustering in worker container
- Added `run_neighb_clustering()` function - executes scanpy neighbor clustering in worker container  
- Updated `run_analysis_job()` to recognize the new clustering function names
- Both functions handle object storage download/upload and cleanup of temporary data

**Key Features**:
- Downloads adata object from MinIO object storage
- Calls the actual clustering functions from Pheno_Cluster_a.py
- Saves results and timing information
- Cleans up temporary objects from storage

### 2. `/source/framework/platform_abstraction.py`
**Purpose**: Added object storage utility functions for sharing data between containers

**Changes**:
- Added `upload_object_data()` - uploads raw object data to MinIO storage
- Added `download_object_data()` - downloads raw object data from MinIO storage
- Added `delete_object_data()` - removes objects from MinIO storage

**Why New Functions Were Needed**:

The existing platform abstraction already had object storage functions, but they were designed for specific use cases:
- **Existing functions** (`upload_zip_object_data`, `download_zip_object_data`): Handle ZIP archives for session state persistence and file-based data exchange
- **`generate_results.py` compatibility**: The primes generation function only used simple data types (integers, strings) that could be passed directly as job parameters and returned as basic Python types

**The Clustering Data Challenge**:

*Why `generate_results.py` worked without new functionality:*
- **Input data**: Simple integers (`limit`, `results_subdir`) - easily passed as job parameters
- **Processing**: Pure computational work with no complex data structures  
- **Output**: Basic Python types (list of integers, float) - directly returnable

*Why clustering functions needed new object storage:*
- **Input data**: Complex AnnData objects containing scipy sparse matrices, pandas DataFrames, and nested metadata - too large and complex for direct parameter passing
- **Container isolation**: Frontend and worker containers are separate processes that cannot share memory
- **Data serialization**: AnnData objects require pickle serialization and binary storage

The new functions enable passing complex, serialized Python objects between containers via shared object storage, which was not needed for the simple arithmetic operations in `generate_results.py`.

**Why Needed**: 
- Worker containers run in isolation and can't access frontend container's filesystem
- Need shared storage (MinIO) to pass adata objects between frontend and worker containers

### 3. `/source/pages2/Pheno_Cluster_a.py` (MOVED BACK from main directory)
**Purpose**: Clustering page with UI and async support, returned to pages2/ subdirectory

**Changes**:
- Moved back from `/source/Pheno_Cluster_a.py` to `/source/pages2/Pheno_Cluster_a.py` 
- Retained all async functionality including `submit_clustering_job()` function
- Contains both sync and async execution modes for clustering algorithms

### 4. `/source/app.py`
**Purpose**: Updated import to reference clustering page in pages2/ subdirectory

**Changes**:
- Changed `import Pheno_Cluster_a` back to `from pages2 import Pheno_Cluster_a`
- This restores the original file organization while maintaining async functionality

### 5. `/source/framework/analysis_functions.py` 
**Purpose**: Updated worker functions to handle pages2/ import path

**Changes**:
- Added dynamic Python path modification in both clustering functions
- Uses `sys.path.insert()` to add source directory to Python path at runtime
- Imports from `pages2.Pheno_Cluster_a` instead of direct import
- This enables worker containers to find clustering functions in subdirectory

### 6. `/docker-compose.yml`
**Purpose**: Fixed async job execution and enabled live code development

**Changes**:
- Fixed `WORKER_IMAGE` environment variable to include proper tag (`frontend:${IMAGE_TAG}`)
- Added volume mount for source code (`./source:/app/source`) to enable live development
- This resolved the original "500 Server Error" when trying to run async jobs

## How It Works

### Sync Mode (Original Behavior)
1. User selects clustering method and parameters
2. Clicks "Run Clustering" 
3. Clustering executes directly in frontend container
4. UI blocks until completion
5. Results displayed immediately

### Async Mode (New Feature)
1. User selects clustering method and parameters
2. Toggles "Run [method] clustering asynchronously" 
3. Clicks "Run [method] clustering"
4. Frontend saves adata object to MinIO storage
5. Job submitted to docker orchestrator
6. Worker container launched to execute clustering
7. Worker downloads adata, runs clustering, saves results
8. Frontend polls for completion and loads results
9. UI remains responsive during execution

## Technical Architecture

```
Frontend Container                    Worker Container
     |                                      |
     | 1. Save adata to MinIO              |
     |-------------------------------->    |
     | 2. Submit job                       |
     |-------------------------------->    |
     |                                     | 3. Download adata
     |                                     | 4. Run clustering
     |                                     | 5. Save results
     |                                     |
     | 6. Poll for completion         <----|
     | 7. Load results                     |
```

### Import Path Resolution for Subdirectories

**The Challenge**: Worker containers execute in a different Python path context than the frontend container. When clustering functions are in `pages2/` subdirectory, worker containers cannot automatically resolve the import path.

**The Solution**: Dynamic Python path modification in worker functions:

```python
# Add the source directory to Python path to ensure pages2 can be imported
source_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if source_dir not in sys.path:
    sys.path.insert(0, source_dir)

from pages2.Pheno_Cluster_a import RunPhenographClust
```

**How it works**:
1. `__file__` points to `analysis_functions.py` in `framework/` subdirectory
2. `os.path.dirname(os.path.dirname())` goes up two levels to reach `source/` directory  
3. Adding `source/` to `sys.path` allows Python to resolve `pages2.Pheno_Cluster_a` imports
4. This happens at runtime in each worker container, ensuring imports work regardless of working directory

**Why this approach**:
- Maintains original file organization with pages in subdirectories
- Enables async functionality without restructuring the entire codebase
- Self-contained solution that doesn't require external configuration changes
- Works consistently across different container execution contexts

## Benefits
- **Responsiveness**: UI doesn't freeze during long clustering operations
- **Scalability**: Multiple clustering jobs can run simultaneously  
- **Resource Isolation**: Heavy computations run in dedicated worker containers
- **Live Development**: Code changes reflected immediately without container rebuilds

## Testing
- Async functionality works for both phenograph and scanpy clustering methods
- PARC and UTAG clustering remain synchronous (not yet implemented for async)
- Volume mounting enables immediate code updates during development
- Fixed original docker orchestrator image resolution issues

## Future Enhancements
- Could extend async support to PARC and UTAG clustering methods
- Could add progress indicators for running jobs
- Could implement job cancellation functionality
