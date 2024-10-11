import streamlit as st
import os
import tempfile
import shutil
import subprocess
import platform_io
import numpy as np


def create_unique_user_directory(base_path: str, debug: bool=False) -> str:
    """
    Create a unique temporary directory within the base path.

    Parameters:
    base_path (str): The base directory where the unique directory will be created.

    Returns:
    str: The path to the created unique directory.
    """
    try:
        unique_user_directory = tempfile.mkdtemp(dir=base_path)
        if debug:
            st.success(f'Directory {unique_user_directory} created.')
        return unique_user_directory
    except Exception as e:
        if debug:
            st.error(f"Error creating directory: {e}")
        raise


def delete_directory_recursively(directory_path: str, debug: bool=False) -> None:
    """
    Recursively delete a directory and all its contents.

    Parameters:
    directory_path (str): The path to the directory to be deleted.
    """
    try:
        if os.path.exists(directory_path):
            shutil.rmtree(directory_path)
            if debug:
                st.success(f"Directory {directory_path} has been deleted.")
        else:
            if debug:
                st.warning(f"Directory {directory_path} does not exist.")
    except Exception as e:
        if debug:
            st.error(f"Error deleting directory: {e}")
        raise


def set_up_output_directory(output_user_dir: str, debug: bool=False) -> None:
    os.makedirs(os.path.join(output_user_dir, 'checkpoints', 'neighborhood_profiles'), exist_ok=True)
    if not os.path.exists(os.path.join(output_user_dir, 'MAWA_Suite_Benchmarking.csv')):
        with open(os.path.join(output_user_dir, 'MAWA_Suite_Benchmarking.csv'), 'w') as f:
            f.write('id,on_NIDAP,file,nSlides,nCells,CellsxSlide,time_load_data,time_to_run_counts,time_to_run_UMAP,time_to_run_cluster\n')
        if debug:
            st.success(f"File {os.path.join(output_user_dir, 'MAWA_Suite_Benchmarking.csv')} created.")
    else:
        if debug:
            st.warning(f"File {os.path.join(output_user_dir, 'MAWA_Suite_Benchmarking.csv')} already exists")


def platform_is_nidap():
    '''
    Check if the Streamlit application is operating on NIDAP
    '''
    return np.any(['nidap.nih.gov' in x for x in subprocess.run('conda config --show channels', shell=True, capture_output=True).stdout.decode().split('\n')[1:-1]])


def platform_is_streamlit_community_cloud():
    '''
    Check if the Streamlit application is operating on the Streamlit Community Cloud
    '''
    return os.getenv('STREAMLIT_SERVER_HEADLESS') == 'true'


def check_for_platform(session_state):
    '''
    Set the platform parameters based on the platform the Streamlit app is running on
    '''
    # Initialize the platform object
    if 'platform' not in session_state:
        session_state['platform'] = platform_io.Platform(platform=('nidap' if platform_is_nidap() else 'local'))
    return session_state


def get_user_dirs_for_platform(input_top_dir: str, output_top_dir: str, saved_states_top_dir: str) -> tuple:
    if platform_is_streamlit_community_cloud():
        input_user_dir = create_unique_user_directory(input_top_dir)
        basename = os.path.basename(input_user_dir)
        output_user_dir = os.path.join(output_top_dir, basename)
        saved_states_user_dir = os.path.join(saved_states_top_dir, basename)
        return (input_user_dir, output_user_dir, saved_states_user_dir)
    else:
        return (input_top_dir, output_top_dir, saved_states_top_dir)
        

def set_up_user_directories(input_top_dir: str, output_top_dir: str, saved_states_top_dir: str, debug: bool=False) -> None:
    try:
        st.session_state['input_user_dir'], st.session_state['output_user_dir'], st.session_state['saved_states_user_dir'] = get_user_dirs_for_platform(input_top_dir, output_top_dir, saved_states_top_dir)

        # Create directories if they do not exist
        if not os.path.exists(st.session_state['output_user_dir']):
            os.makedirs(st.session_state['output_user_dir'])
            if debug:
                st.success(f"Output directory {st.session_state['output_user_dir']} created.")
        else:
            if debug:
                st.warning(f"Output directory {st.session_state['output_user_dir']} already exists.")

        set_up_output_directory(st.session_state['output_user_dir'])
        
        if not os.path.exists(st.session_state['saved_states_user_dir']):
            os.makedirs(st.session_state['saved_states_user_dir'])
            if debug:
                st.success(f"Saved states directory {st.session_state['saved_states_user_dir']} created.")
        else:
            if debug:
                st.warning(f"Saved states directory {st.session_state['saved_states_user_dir']} already exists.")
    except Exception as e:
        if debug:
            st.error(f"Error during initialization: {e}")
        raise


def delete_user_directories(debug: bool=False) -> None:
    # Remove directories from session state and delete them from the filesystem
    if 'input_user_dir' in st.session_state:
        delete_directory_recursively(st.session_state.pop('input_user_dir'))
    else:
        if debug:
            st.warning('Input user directory does not exist in session state.')
    
    if 'output_user_dir' in st.session_state:
        delete_directory_recursively(st.session_state.pop('output_user_dir'))
    else:
        if debug:
            st.warning('Output user directory does not exist in session state.')
    
    if 'saved_states_user_dir' in st.session_state:
        delete_directory_recursively(st.session_state.pop('saved_states_user_dir'))
    else:
        if debug:
            st.warning('Saved states user directory does not exist in session state.')


def get_num_cpus_for_sit():
    if platform_is_nidap():
        return 7
    elif platform_is_streamlit_community_cloud():
        return 2
    else:  # local
        return 4


def main() -> None:

    input_top_dir = os.path.join('.', 'input')
    output_top_dir = os.path.join('.', 'output')
    saved_states_top_dir = os.path.join('.', 'saved_streamlit_session_states')

    if 'input_user_dir' not in st.session_state:
        set_up_user_directories(input_top_dir, output_top_dir, saved_states_top_dir)

    if 'num_cpus_for_sit' not in st.session_state:
        st.session_state['num_cpus_for_sit'] = get_num_cpus_for_sit()

    st.write(f"Input user directory: {st.session_state['input_user_dir']}")
    st.write(f"Output user directory: {st.session_state['output_user_dir']}")
    st.write(f"Saved states user directory: {st.session_state['saved_states_user_dir']}")

    # Output the contents of these directories
    st.write(f"Input user directory contents: {os.listdir(st.session_state['input_user_dir'])}")
    st.write(f"Output user directory contents: {os.listdir(st.session_state['output_user_dir'])}")
    st.write(f"Saved states user directory contents: {os.listdir(st.session_state['saved_states_user_dir'])}")

    if st.button('Delete user directories'):
        delete_user_directories()
        

if __name__ == '__main__':
    main()
