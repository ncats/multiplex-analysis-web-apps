import subprocess
import importlib

def is_mamba_installed():
    try:
        # Run the 'mamba --version' command
        result = subprocess.run(['mamba', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        # Check if the command was successful
        if result.returncode == 0:
            print("&&&& Mamba is installed.")
            return True
        else:
            print("&&&& Mamba is not installed.")
            return False
    except FileNotFoundError:
        # The command was not found
        print("&&&& Mamba is not installed.")
        return False

def install_with_mamba(packages):
    print(f"&&&& Attempting to install {', '.join(packages)} with mamba.")
    try:
        # Run the 'mamba install <packages>' command
        result = subprocess.run(['mamba', 'install', '-y'] + packages, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        # Check if the command was successful
        if result.returncode == 0:
            print(f"&&&& {', '.join(packages)} have been installed successfully with mamba.")
            print(result.stdout)
        else:
            print(f"&&&& Failed to install {', '.join(packages)} with mamba.")
            print(result.stderr)
    except Exception as e:
        print(f"&&&& An error occurred while trying to install {', '.join(packages)} with mamba: {e}")

def install_with_conda(packages):
    print(f"&&&& Attempting to install {', '.join(packages)} with conda.")
    try:
        # Run the 'conda install <packages>' command
        result = subprocess.run(['conda', 'install', '-y'] + packages, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        # Check if the command was successful
        if result.returncode == 0:
            print(f"&&&& {', '.join(packages)} have been installed successfully with conda.")
            print(result.stdout)
        else:
            print(f"&&&& Failed to install {', '.join(packages)} with conda.")
            print(result.stderr)
    except Exception as e:
        print(f"&&&& An error occurred while trying to install {', '.join(packages)} with conda: {e}")

def install_with_pip(packages):
    print(f"&&&& Attempting to install {', '.join(packages)} with pip.")
    try:
        # Run the 'pip install <packages>' command
        result = subprocess.run(['pip', 'install'] + packages, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        # Check if the command was successful
        if result.returncode == 0:
            print(f"&&&& {', '.join(packages)} have been installed successfully with pip.")
            print(result.stdout)
        else:
            print(f"&&&& Failed to install {', '.join(packages)} with pip.")
            print(result.stderr)
    except Exception as e:
        print(f"&&&& An error occurred while trying to install {', '.join(packages)} with pip: {e}")

def live_package_installation():
    packages_to_install = ['hnswlib', 'parc', 'sklearn_ann', 'annoy', 'pyNNDescent']  # last two probably only needed for published dashboards
    installers_to_use = ['mamba', 'conda', 'pip']

    for package in packages_to_install:
        try:
            importlib.import_module(package.lower())
            print(f"&&&& {package} is already installed.")
        except ImportError:
            print(f"&&&& {package} is not installed.")
            for installer in installers_to_use:
                print(f"&&&& Trying to install {package} using {installer}.")
                if installer == 'mamba':
                    if is_mamba_installed():
                        install_with_mamba([package])
                    else:
                        print(f"&&&& mamba is not installed. Trying the next installer.")
                        continue
                elif installer == 'conda':
                    install_with_conda([package])
                elif installer == 'pip':
                    install_with_pip([package])
                
                try:
                    importlib.import_module(package.lower())
                    print(f"&&&& {package} has been successfully installed using {installer}.")
                    break
                except ImportError:
                    print(f"&&&& {package} was not successfully installed with {installer}.")
            else:
                print(f"&&&& {package} could not be installed after trying all installers.")
