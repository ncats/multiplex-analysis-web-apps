'''
A set of methods for installing required packages
that are not previously installed
'''

import subprocess
import importlib

def is_mamba_installed():
    '''
    Function to check if mamba is installed

    Returns:
        Boolean: Is mamba installed or not
    '''
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
    '''
    Install list of packages using mamba
    '''
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
    '''
    Install list of packages using conda
    '''

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
    '''
    Install list of packages using pip
    '''
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
    '''
    Function to check if packages are installed
    '''
    
    # last two probably only needed for published dashboards
    packages_to_install = ['hnswlib', 'parc', 'sklearn_ann', 'annoy', 'pyNNDescent']
    installers_to_use = ['mamba', 'pip']

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
                        print("&&&& mamba is not installed. Trying the next installer.")
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
