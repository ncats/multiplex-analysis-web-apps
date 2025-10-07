import importlib
import sys
from typing import Dict, Any

def import_and_show_versions(package_names: list) -> Dict[str, Any]:
    """Import packages and display their versions"""
    imported_packages = {}
    version_info = []
    
    for pkg_name in package_names:
        try:
            # Import the package
            module = importlib.import_module(pkg_name)
            imported_packages[pkg_name] = module
            
            # Try to get version
            version = getattr(module, '__version__', 'Unknown')
            if version == 'Unknown':
                # Some packages store version differently
                try:
                    version = getattr(module, 'version', 'Unknown')
                except:
                    version = 'Unknown'
            
            # Convert version to string to handle cases where it's not a string
            version_str = str(version) if version != 'Unknown' else 'Unknown'
            
            version_info.append((pkg_name, version_str, '✅'))
            
        except ImportError as e:
            version_info.append((pkg_name, str(e), '❌'))
        except Exception as e:
            version_info.append((pkg_name, f"Error: {str(e)}", '⚠️'))
    
    # Display results in a nice table format
    print(f"{'Package':<20} {'Version':<15} {'Status'}")
    print("=" * 50)
    for name, version, status in version_info:
        print(f"{name:<20} {version:<15} {status}")
    
    # Add imported packages to global namespace
    globals().update(imported_packages)
    
    return imported_packages

# Your package list
# Insert output from e.g. "grep "=" source/environment.yml | awk -v FS="- " '{print $2}' | awk -v FS== '{gsub("-", "_", $1); printf("\"%s\", ", $1)}' | less"
packages = [
    "streamlit", "psycopg2", "yaml", "dill", "git", "requests", "minio", "polars", "snowflake.snowpark", "matplotlib", "natsort", "numpy", "palantir", "pandas", "plotly", "scipy", "seaborn", "skimage", "sklearn", "split_file_reader", "streamlit_extras", "tqdm", "umap", "pympler", "objsize", "phenograph", "parmap", "setuptools_scm", "pynndescent", "plotnine", "shapely", "hnswlib", "spatialdata", "dask", "OpenSSL", "numba", "pip", "st_pages", "streamlit_javascript", "parc", "sklearn_ann", "anndata", "annoy", "boto3", "squidpy"
]

# Import all packages and show versions
imported = import_and_show_versions(packages)

# Now all packages are available in the global namespace
# e.g., you can use numpy as np if you want
np = imported.get('numpy')
pd = imported.get('pandas')
if 'matplotlib' in imported:
    plt = imported['matplotlib'].pyplot