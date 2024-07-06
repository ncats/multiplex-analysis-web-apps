'''
Set of functions for managing the MAWA platform
'''

# Import relevant libraries
import os
import time
import shutil
import pandas as pd
import streamlit as st
import streamlit_dataframe_editor as sde
import utils
from pages2 import memory_analyzer

# Constant
local_input_dir = os.path.join('.', 'input')
local_output_dir = os.path.join('.', 'output')

# Write a dataframe from a file listing with columns for selection, filename, # of files inside (for directories), and modification time, sorted descending by modification time
# Note this is primarily for local listings, not remote listings
def make_complex_dataframe_from_file_listing(dirpath, item_names, df_session_state_key_basename=None, editable=True):
    import time
    num_contents = [len(os.listdir(os.path.join(dirpath, x))) if os.path.isdir(os.path.join(dirpath, x)) else None for x in item_names]
    modification_times = [os.path.getmtime(os.path.join(dirpath, x)) for x in item_names]
    selecteds = [False for _ in item_names]
    df = pd.DataFrame({'Selected': selecteds, 'File or directory name': item_names, '# of files within': num_contents, 'Modification time': [time.ctime(x) for x in modification_times], 'mod_time_sec': modification_times}).sort_values('mod_time_sec', ascending=False).reset_index(drop=True)
    if editable:
        ss_de_key_name = 'loader__de_' + df_session_state_key_basename
        ss_df_key_name = 'loader__df_' + df_session_state_key_basename
        # st.session_state[ss_df_key_name] = st.data_editor(df.iloc[:, :-1], key=(ss_df_key_name + '_input__do_not_persist'))
        if ss_de_key_name in st.session_state:
            if set(df['File or directory name']) != set(st.session_state[ss_de_key_name].reconstruct_edited_dataframe()['File or directory name']):
                del st.session_state[ss_de_key_name]
        if ss_de_key_name not in st.session_state:
            st.session_state[ss_de_key_name] = sde.DataframeEditor(df_name=ss_df_key_name, default_df_contents=df.iloc[:, :-1])
        st.session_state[ss_de_key_name].dataframe_editor(reset_data_editor_button_text='Reset file selections')
    else:
        st.dataframe(df.iloc[:, 1:-1])
        if df_session_state_key_basename is not None:
            st.warning('Session state key {} is not being assigned since editable=False was selected in call to make_complex_dataframe_from_file_listing()'.format(ss_df_key_name))

# Write an editable dataframe (simple, having only a selection column and filenames with .zip removed) for the available files, also saving the filenames with the possible .zip extensions to a separate Series
def make_simple_dataframe_from_file_listing(available_files, df_session_state_key_basename=None, streamlit_key_for_available_filenames_srs=None, editable=True):

    # Save (to Streamlit, analogous to how it's done for the data editor, below), the full filenames of the available files
    if streamlit_key_for_available_filenames_srs is not None:
        st.session_state[streamlit_key_for_available_filenames_srs] = pd.Series(available_files, name='doesnt_matter')

    # Make a copy of the available files with the .zip extensions removed, if present
    available_files = [filename.split('.zip')[0] for filename in available_files]

    # Create a simple dataframe of the available files, with a selection column and stripped of any .zip extensions
    df = pd.DataFrame({'Selected': [False for _ in available_files], 'File or directory name': available_files})

    # Display an editable dataframe version of this
    if editable:
        ss_de_key_name = 'loader__de_' + df_session_state_key_basename
        ss_df_key_name = 'loader__df_' + df_session_state_key_basename
        # st.session_state[ss_df_key_name] = st.data_editor(df, key=(ss_df_key_name + '_input__do_not_persist'))
        if ss_de_key_name in st.session_state:
            if set(df['File or directory name']) != set(st.session_state[ss_de_key_name].reconstruct_edited_dataframe()['File or directory name']):
                del st.session_state[ss_de_key_name]
        if ss_de_key_name not in st.session_state:
            st.session_state[ss_de_key_name] = sde.DataframeEditor(df_name=ss_df_key_name, default_df_contents=df)
        st.session_state[ss_de_key_name].dataframe_editor(reset_data_editor_button_text='Reset file selections')
    else:
        st.dataframe(df)
        if df_session_state_key_basename is not None:
            st.warning('Session state key {} is not being assigned since editable=False was selected in call to make_simple_dataframe_from_file_listing()'.format(ss_df_key_name))

# Delete selected files/dirs from a directory
def delete_selected_files_and_dirs(directory, selected_files):
    for curr_file in selected_files:
        curr_path = os.path.join(directory, curr_file)
        if os.path.isfile(curr_path):
            os.remove(curr_path)
        elif os.path.isdir(curr_path):
            shutil.rmtree(curr_path)

# Write the current settings to disk, including a timestamp, hostname, and git commit
def write_current_tool_parameters_to_disk(output_dir):
    import socket
    import subprocess
    import yaml
    import streamlit_utils
    import sys
    print('Writing the current tool parameters to disk...')
    settings_yaml_filename = 'settings_as_of_{}.yml'.format(utils.get_timestamp())
    pathname = os.path.join(output_dir, settings_yaml_filename)
    if not os.path.exists(pathname):
        with open(pathname, mode='wt') as file:
            file.write('# Timestamp: {}\n'.format(utils.get_timestamp(pretty=True)))
            file.write('# Hostname: {}\n'.format(socket.gethostname()))
            if os.path.exists('.git'):
                file.write('# Git commit: {}\n'.format(subprocess.run('git rev-parse HEAD', shell=True, capture_output=True).stdout.decode().split('\n')[:-1][0]))
            else:
                file.write('# Git commit: Currently unknown because there is no .git directory present\n')
            file.write('# Python version (may conflict with environment.yml, showing strange system setup): {}\n'.format(sys.version.split('\n')[0]))
            file.write('\n')
            yaml.dump(streamlit_utils.get_current_settings(), file, sort_keys=False)
            if 'sit__used_settings' in st.session_state:
                file.write('\n')
                file.write('Actual settings used in the Spatial Interaction Tool:\n')
                file.write('\n')
                yaml.dump(st.session_state['sit__used_settings'], file, sort_keys=False)
    else:
        st.warning('File {} already exists; not overwriting it'.format(pathname))

# Write the current conda/pip environment to disk
def write_current_environment_to_disk(output_dir):
    import subprocess
    import os
    print('Writing the current conda/pip environment to disk...')
    environment_yaml_filename = 'environment_as_of_{}.yml'.format(utils.get_timestamp())
    pathname = os.path.join(output_dir, environment_yaml_filename)
    if not os.path.exists(pathname):
        subprocess.run('conda env export > {}'.format(pathname), shell=True, capture_output=True)
    else:
        st.warning('File {} already exists; not overwriting it'.format(pathname))

# Create an empty output archive directory
def create_empty_output_archive(new_output_archive_name, output_dir):
    archive_dirname = 'output_archive-{}-{}'.format(utils.get_timestamp(), new_output_archive_name)
    archive_path = os.path.join(output_dir, archive_dirname)
    if not os.path.exists(archive_path):
        os.mkdir(archive_path)
        st.info('Directory {} has been created'.format(archive_path))
    else:
        st.warning('Directory {} already exists; not creating it'.format(archive_path))
    return archive_path

# Copy all contents from the output directory to a new output archive
# Currently only applies to local
def copy_output_dir_contents_to_output_archive(new_output_archive_name, output_dir):
    print(f'Copying all contents from {output_dir} to a new output archive with descriptive name "{new_output_archive_name}"...')
    archive_dirname = 'output_archive-{}-{}'.format(utils.get_timestamp(), new_output_archive_name)
    archive_path = os.path.join(output_dir, archive_dirname)
    if not os.path.exists(archive_path):
        shutil.copytree(output_dir, archive_path, ignore=shutil.ignore_patterns('output_archive-*'))
        st.info('Data has been saved to {}'.format(archive_path))
    else:
        st.warning('Archive directory {} already exists; not saving'.format(archive_path))

# Because there's no clean and efficient way in Python (or on NIDAP, which doesn't have the zip bash function) to create a zip file with ignores without doing a full crawl using the zipfile library--which in hindsight is likely better than the following and is what is indeed implemented below in create_zipfile_from_files_in_dir() and should be used instead of this function!--create a function to create a zipfile while ignoring some files/dirs based on a prefix
def create_zipfile_with_ignores(zipfile_dirpath, basename_suffix_for_zipfile, prefix_to_ignore, local_output_dir, local_tmp_dir):

    # Import relevant libraries
    import subprocess

    # Ultimate zipfile path without the .zip extension
    zipfile_basename = os.path.join(zipfile_dirpath, 'output_archive-{}-{}'.format(utils.get_timestamp(), basename_suffix_for_zipfile))

    # Determine if items exist in the local output directory that we want to ignore; otherwise, don't jump through any hoops!
    items_to_ignore_exist = len([x for x in os.listdir(local_output_dir) if x.startswith(prefix_to_ignore)]) > 0

    # If there are items to ignore...
    if items_to_ignore_exist:

        # Ensure a temporary directory exists and is empty
        ensure_empty_directory(local_tmp_dir)
        
        # Move all directories and files starting with "output_archive-" (or the generic prefix to ignore) in the local output directory to the new temporary directory
        subprocess.run('mv {}{}{}* {}'.format(local_output_dir, os.path.sep, prefix_to_ignore, local_tmp_dir), shell=True)

    # Zip the entire local output directory, which now (or already) contains everything but the "output_archive-" files/dirs, into the desired zip file
    shutil.make_archive(zipfile_basename, 'zip', local_output_dir + os.path.sep)

    # If there are items to ignore...
    if items_to_ignore_exist:

        # Move everything in the temporary directory back into the local output directory
        subprocess.run('mv {}/* {}'.format(local_tmp_dir, local_output_dir), shell=True)

        # Delete the temporary directory
        shutil.rmtree(local_tmp_dir)

    # Return the full path of the created zipfile
    return zipfile_basename + '.zip'

# Ensure a directory exists but is empty
def ensure_empty_directory(dirpath):
    import shutil
    import os
    if os.path.exists(dirpath):
        shutil.rmtree(dirpath)  # if it exists, then delete it
        print('Tree {} deleted'.format(dirpath))
    os.mkdir(dirpath)  # create an empty destination directory
    print('New directory {} created'.format(dirpath))

# Logic in case currently selected value is no long present in directory listing as might happen when it has just been deleted
def account_for_stale_streamlit_values(session_state_key, possible_values):
    if session_state_key not in st.session_state:
        if len(possible_values) > 0:
            st.session_state[session_state_key] = possible_values[0]
        else:
            st.session_state[session_state_key] = None
    if len(possible_values) > 0:
        if st.session_state[session_state_key] not in possible_values:
            st.session_state[session_state_key] = possible_values[0]
    else:
        st.session_state[session_state_key] = None

# Create a class that takes the platform type (e.g., local, nidap) as input and creates corresponding methods, with consisting naming, for performing the same function on different platforms
class Platform:

    # Object instantiation
    def __init__(self, platform='local'):
        self.platform = platform
        self.available_inputs = None
        self.available_archives = None

    # Get a list of the files available to import locally to serve as input files for a workflow
    def get_available_inputs_listing(self):
        # Potentially slow

        # When running locally, this is irrelevant as we can easily place all input files into the same ./input directory, as opposed to having to read them in from a remote
        if self.platform == 'local':
            available_inputs = []

        # On NIDAP, load the metadata for the "input" unstructured dataset
        elif self.platform == 'nidap':
            import nidap_io
            dataset = nidap_io.get_foundry_dataset(alias='input')
            dataset_file_objects = nidap_io.get_file_objects_from_dataset(dataset)  # slow
            available_inputs = nidap_io.list_files_in_dataset(dataset_file_objects)
            self.dataset_file_objects_for_available_inputs = dataset_file_objects  # save the values *that we'll need later* that result from the long calculation (1-2 sec) as properties of the object so they're stored rather than discarded

        # Save the values *that we'll need later* that result from the long calculation (1-2 sec) as properties of the object so they're stored rather than discarded
        self.available_inputs = sorted(available_inputs)
    
    # Write a dataframe of the available inputs on the remote
    def display_available_inputs_df(self):

        st.subheader('Input Data on NIDAP :open_file_folder:')
        
        # Again, irrelevant for local
        if self.platform == 'local':
            pass

        # If on NIDAP...
        elif self.platform == 'nidap':

            # If we've never determined the inputs available on the remote (e.g., when the script first starts), do so now
            if self.available_inputs is None:
                self.get_available_inputs_listing()

            # Get a shortcut to the available input list
            available_inputs = self.available_inputs

            # Create a simple editable dataframe of the available input filenames
            make_simple_dataframe_from_file_listing(available_files=available_inputs, df_session_state_key_basename='available_inputs', streamlit_key_for_available_filenames_srs='srs_available_input_filenames', editable=True)
            
    # Add a button to re-read the available input files on the remote
    def add_refresh_available_inputs_button(self):

        # Irrelevant for local
        if self.platform == 'local':
            pass

        # If on NIDAP, create a button to simply update the available inputs
        elif self.platform == 'nidap':
            if st.button(':arrows_clockwise: Refresh available input data'):
                self.get_available_inputs_listing()
                st.rerun()  # rerun since this potentially changes outputs... rule of thumb for rerunning the page should probably be that if this method changes outputs, will those possibly changed outputs definitely get redrawn? If not, do a rerun! Consider where this method falls in the top-down rerun of the calling script, are the outputs before or after the method is called?
    
    # Load any selected available inputs on the remote to the local machine
    def load_selected_inputs(self):

        st.subheader('Load into MAWA :arrow_forward:')
        load_button = st.button('Load NIDAP input data into MAWA :arrow_right:')
        # Irrelevant for local
        if self.platform == 'local':
            pass

        # If on NIDAP...
        elif self.platform == 'nidap':

            # If a load button is clicked...
            if load_button:

                # Import relevant libraries
                import nidap_io
                import sys

                # # Get "shortcuts" to the object properties
                # dataset_file_objects = self.dataset_file_objects_for_available_inputs

                # Get the selected filenames
                df_available_inputs = st.session_state['loader__de_available_inputs'].reconstruct_edited_dataframe()
                srs_available_input_filenames = st.session_state['srs_available_input_filenames']
                selected_input_filenames = srs_available_input_filenames[df_available_inputs['Selected']].tolist()

                # Download the selected files
                all_downloaded_files = nidap_io.download_files_from_dataset(nidap_io.get_foundry_dataset(alias='input'), dataset_filter_func=lambda f: f.path in selected_input_filenames, limit=15)

                # For each downloaded file, move it to the local input directory
                for selected_input_filename, local_download_path in all_downloaded_files.items():

                # # For each selected available input file...
                # for selected_input_filename in selected_input_filenames:

                    # # Download the file and get its local download path
                    # dataset_file_object = nidap_io.get_dataset_file_object(dataset_file_objects, selected_filename=selected_input_filename)
                    # local_download_path = nidap_io.download_file_from_dataset(dataset_file_object)  # slow

                    # Populate the local input directory with the current selection
                    # Rules:
                    #   * If it's not a .zip file, just copy it over
                    #   * If it's a .zip file, there can be either 1 or 2 periods in the full filename including extension
                    #   * If there's just a single period (e.g., asdf.zip), it must be a zipped directory with name following either DIRNAME.zip or DIRNAME--bleh.zip
                    #   * If there are two periods (e.g., asdf.csv.zip), it must be a zipped datafile with name following asdf.csv
                    if selected_input_filename.endswith('.zip'):
                        splitted = selected_input_filename.split('.')  # should be of length 2 or 3 (for, e.g., asdf.csv.zip)
                        num_periods = len(splitted) - 1  # should be 1 or 2
                        if (num_periods < 1) or (num_periods > 2):
                            st.error('Available .zip input filename {} has a bad number of periods ({}... it should have 1-2 periods); please fix this.'.format(selected_input_filename, num_periods))
                            sys.exit()
                        if num_periods == 1:  # it's a zipped directory, by specification
                            if '--' not in selected_input_filename:
                                dirpath = os.path.join(local_input_dir, selected_input_filename.rstrip('.zip'))
                            else:
                                dirpath = os.path.join(local_input_dir, selected_input_filename.split('--')[0])
                            ensure_empty_directory(dirpath)
                            shutil.unpack_archive(local_download_path, dirpath)
                        elif num_periods == 2:  # it's a zipped datafile
                            shutil.unpack_archive(local_download_path, local_input_dir)
                    else:
                        shutil.copy(local_download_path, local_input_dir)
    
    # Save a MAWA-unified datafile to NIDAP
    def save_selected_input(self):

        # If working on NIDAP...
        if self.platform == 'nidap':

            # Write a header
            st.subheader(':tractor: Save MAWA-unified datafile to NIDAP')

            # Create a list of the CSV files having a "mawa-unified_datafile-" prefix and ".csv" suffix in the local input directory
            mawa_unified_datafiles = [x for x in os.listdir(local_input_dir) if x.startswith('mawa-unified_datafile-') and x.endswith('.csv')]

            # Create a dictionary of the stripped filenames and their corresponding full filenames
            mawa_unified_datafiles_dict = {x.split('mawa-unified_datafile-')[1].split('.csv')[0]: x for x in mawa_unified_datafiles}
            keys = list(mawa_unified_datafiles_dict.keys())

            # If the session state key doesn't exist, create it and set it to the first key in the list (if it exists)
            if ('loader__mawa_unified_datafile_to_save' not in st.session_state) or (st.session_state['loader__mawa_unified_datafile_to_save'] not in keys):
                st.session_state['loader__mawa_unified_datafile_to_save'] = keys[0] if keys else None
            st.selectbox('Select MAWA-unified datafile to save:', keys, key='loader__mawa_unified_datafile_to_save')

            # Create a button to zip the selected file and save it to NIDAP
            if st.button('Save selected (above) MAWA-unified datafile to NIDAP :arrow_left:', help='This will zip the selected file and save it to NIDAP. We generally don\'t want to save a file **generated** in the app to the **`input`** dataset on NIDAP on principle, but this is a reasonable exception so that the file can be used again or in other use cases.', disabled=st.session_state['loader__mawa_unified_datafile_to_save'] is None):

                # Create a spinner to indicate that the zipping and saving is in progress
                with st.spinner('Zipping and saving...'):

                    # Zip the selected file
                    selected_mawa_unified_datafile = mawa_unified_datafiles_dict[st.session_state['loader__mawa_unified_datafile_to_save']]
                    shutil.make_archive(os.path.join(local_input_dir, selected_mawa_unified_datafile), 'zip', local_input_dir, selected_mawa_unified_datafile)

                    # Transfer the zipped file to NIDAP
                    import nidap_io
                    dataset = nidap_io.get_foundry_dataset(alias='input')
                    upload_single_file_to_dataset((dataset, local_input_dir, selected_mawa_unified_datafile + '.zip'))
                    # nidap_io.upload_file_to_dataset(dataset, selected_filepath=os.path.join(local_input_dir, selected_mawa_unified_datafile + '.zip'))

                    # Delete the zipped file from the local input directory
                    os.remove(os.path.join(local_input_dir, selected_mawa_unified_datafile + '.zip'))

    # Get a listing of the files/dirs in the local input directory, which is platform-independent because it's local
    def get_local_inputs_listing(self):
        return sorted([x for x in os.listdir(local_input_dir) if not x.endswith('.zip')])  # ignore zip files, which can appear locally only on a local platform, since for a remote platform such as NIDAP, per above, all zip files get unzipped
    
    # Write a dataframe of the local input files, which we don't want to be editable because we don't want to mess with the local inputs (for now), even though they're basically a local copy
    def display_local_inputs_df(self):
        st.subheader('Input Data in MAWA :open_file_folder:')
        local_inputs = self.get_local_inputs_listing()
        if self.platform == 'local':  # not editable locally because deletion is disabled anyway so there'd be nothing to do with selected files
            make_complex_dataframe_from_file_listing(dirpath=local_input_dir, item_names=local_inputs, editable=False)
        elif self.platform == 'nidap':  # editable on NIDAP because deletion is enabled since it's safe to delete loaded input files since they're backed up to NIDAP
            make_complex_dataframe_from_file_listing(dirpath=local_input_dir, item_names=local_inputs, df_session_state_key_basename='local_inputs', editable=True)

    # Possibly allow for deletion of loaded input files
    def add_delete_local_inputs_button(self):

        # search for other calls to make_complex_dataframe_from_file_listing() and search for the key to ensure I'm consistent with how I'm grabbing the selected values from the edited dataframe!!

        # Don't delete local input files because there are ostensibly no backups
        if self.platform == 'local':
            pass

        # Delete local input files on NIDAP which is safe because they are backed up to NIDAP and deletion here will only be for the *loaded* input files
        elif self.platform == 'nidap':
            if st.button(':x: Delete selected (above) loaded input files'):
                df_local_inputs = st.session_state['loader__de_local_inputs'].reconstruct_edited_dataframe()
                local_input_files_to_delete = df_local_inputs[df_local_inputs['Selected']]['File or directory name']
                delete_selected_files_and_dirs(local_input_dir, local_input_files_to_delete)
                st.rerun()
    
    # List the results archives on the remote
    def get_archives_listing(self):

        # List the output_archive-* folders
        if self.platform == 'local':
            available_archives = [x for x in os.listdir(local_output_dir) if (x.startswith('output_archive-') and (not x.endswith('.zip')))]  # locally, archives shouldn't be zipped (rather in directories), though for setup/transfer-to-nidap purposes, there may exist corresponding zip files
            available_archives_trimmed = available_archives

        # List the contents of the output unstructured dataset (there should only be output_archive-*.zip files)
        elif self.platform == 'nidap':
            import nidap_io
            dataset = nidap_io.get_foundry_dataset(alias='output')
            dataset_file_objects = nidap_io.get_file_objects_from_dataset(dataset)  # slow
            available_archives = [x for x in nidap_io.list_files_in_dataset(dataset_file_objects) if (x.startswith('output_archive-') and ('.zip' in x))]
            available_archives_trimmed = []
            for archive_basename in set([x.split('.zip')[0] for x in available_archives]):
                curr_parts_files = [x for x in available_archives if x.startswith(archive_basename + '.zip.')]
                num_parts_files = len(curr_parts_files)
                if num_parts_files > 0:  # it's in parts
                    num_expected_parts = int(curr_parts_files[0].split('_')[-1])
                    if num_parts_files == num_expected_parts:
                        print('{} is a complete set of zip parts; adding it to the list'.format(archive_basename))
                        available_archives_trimmed.append(archive_basename + '.zip.')
                    else:
                        print('WARNING: {} is not a complete set of zip parts: {} files expected, found {}. Not adding it to the list'.format(archive_basename, num_expected_parts, num_parts_files))
                else:  # it's a single normal zip file
                    print('{} is a regular, single zip file; adding it to the list'.format(archive_basename))
                    available_archives_trimmed.append(archive_basename + '.zip')
            self.dataset_file_objects_for_available_archives = dataset_file_objects  # save the results of the long calculations that we'll need later as object properties

        # Save the results of the long calculations that we'll need later as object properties
        self.available_archives = sorted(available_archives_trimmed)
    
    # Write a dataframe of the available archives
    def display_archives_df(self):

        # Always get the listing locally because that's fast
        if self.platform == 'local':
            st.subheader(':open_file_folder: Available results archives (i.e., saved results)')
            self.get_archives_listing()
            make_complex_dataframe_from_file_listing(dirpath=local_output_dir, item_names=self.available_archives, df_session_state_key_basename='available_archives', editable=True)

        # Only get the listing on NIDAP when it's not already loaded (or when the refresh button is hit, below) because that's "slow"
        elif self.platform == 'nidap':
            st.subheader(':open_file_folder: Available results archives (i.e., saved results) on NIDAP')
            if self.available_archives is None:
                self.get_archives_listing()
            make_simple_dataframe_from_file_listing(available_files=self.available_archives, editable=False)
    
    # Add a button for deleting available archives
    def add_delete_archives_button(self):

        # If working locally...
        if self.platform == 'local':

            # If the button is hit...
            if st.button(':x: Delete selected (above) results archives'):

                # Get the editable dataframe of the available archives from Streamlit
                df_available_archives = st.session_state['loader__de_available_archives'].reconstruct_edited_dataframe()

                # Get the names of just the selected directories
                dirs_to_delete = df_available_archives[df_available_archives['Selected']]['File or directory name']

                # Delete them from the output results directory
                delete_selected_files_and_dirs(local_output_dir, dirs_to_delete)

                # Rerun since deleting archives changes outputs
                st.rerun()

        # Don't do this on NIDAP, if it were even possible
        elif self.platform == 'nidap':
            pass

    # Add a button to refresh the results archives listing
    def add_refresh_archives_button(self):

        # No need to do this if working locally because the listing is always refreshed when the dataframe is drawn
        if self.platform == 'local':
            pass

        # Do this on NIDAP manually because it's "slow"
        elif self.platform == 'nidap':
            if st.button(':arrows_clockwise: Refresh available results archives'):
                self.get_archives_listing()
                st.rerun()  # this may change outputs so refresh
    
    # Load the selected results archives so calculations can be resumed or results can be visualized
    # def load_selected_archive(self, nworkers_for_data_transfer=8):
    def load_selected_archive(self):

        st.subheader(':tractor: Load results')

        # If working locally...
        if self.platform == 'local':

            # Account for possibly stale values of the session state key
            account_for_stale_streamlit_values('archive_to_load', self.available_archives)

            # Make a dropdown of currently available results archives
            st.selectbox('Select available results archive to load:', self.available_archives, key='archive_to_load')

            # If the user wants to load the selected archive...
            if st.button('Load selected (above) results archive :arrow_right:', help='WARNING: This will copy the contents of the selected archive to the results directory and will overwrite currently loaded results; please ensure they are backed up (you can just use the functions on this page)!'):

                # First delete everything in currently in the output results directory (i.e., all currently loaded data) that's not an output archive
                delete_selected_files_and_dirs(local_output_dir, self.get_local_results_listing())

                # Copy everything from the selected output archive to the output directory
                shutil.copytree(os.path.join(local_output_dir, st.session_state['archive_to_load']), local_output_dir, dirs_exist_ok=True)

                # Mimic (sort of) callback behavior, especially because we want to see updated available session states
                st.rerun()

        # If working on NIDAP...
        elif self.platform == 'nidap':

            # Account for possibly stale values of the session state key
            current_options = [x.split('.zip')[0] for x in self.available_archives]
            account_for_stale_streamlit_values('archive_to_load', current_options)

            # Make a dropdown of currently available results archives
            st.selectbox('Select available results archive to load:', current_options, key='archive_to_load')

            # If the user wants to load the selected archive...
            if st.button('Load selected (above) results archive :arrow_right:', help='WARNING: This will copy the contents of the selected archive to the results directory and will overwrite currently loaded results; please ensure they are backed up (you can just use the functions on this page)!'):

                # Import relevant libraries
                import nidap_io
                # import utils
                # import multiprocessing

                # Delete all files currently present in the output results directory
                delete_selected_files_and_dirs(local_output_dir, self.get_local_results_listing())
                
                # Obtain the full filename corresponding to the selected archive to load and run a check
                list_of_len_1 = [x for x in self.available_archives if x.startswith(st.session_state['archive_to_load'])]
                if len(list_of_len_1) != 1:
                    import sys
                    print('ERROR: More than one available archive found ({}) for selected archive to load ({})'.format(list_of_len_1, st.session_state['archive_to_load']))
                    sys.exit()
                selected_archive_with_proper_extension = list_of_len_1[0]

                # If it's a single normal zip file, download and unzip it
                if not selected_archive_with_proper_extension.endswith('.'):

                    # dataset_file_object = nidap_io.get_dataset_file_object(self.dataset_file_objects_for_available_archives, selected_filename=selected_archive_with_proper_extension)

                    start_time = time.time()
                    # local_download_path = nidap_io.download_file_from_dataset(dataset_file_object)
                    all_downloaded_files = nidap_io.download_files_from_dataset(nidap_io.get_foundry_dataset(alias='output'), dataset_filter_func=lambda f: f.path == selected_archive_with_proper_extension, limit=15)
                    local_download_path = all_downloaded_files[selected_archive_with_proper_extension]
                    filesize = os.path.getsize(local_download_path) / 1024 ** 2
                    duration = time.time() - start_time
                    print('  Download of {} ({:5.3f} MB) from Compass to Workspaces took {:3.1f} seconds --> {:3.1f} MB/s'.format(selected_archive_with_proper_extension, filesize, duration, filesize / duration))

                    extract_zipfile_to_directory(zipfile_name=local_download_path, extraction_path=local_output_dir)

                # If it corresponds to a chunked set of zip files...
                else:

                    # Obtain the corresponding chunked set of zip files
                    matching_archives_files = sorted([x for x in nidap_io.list_files_in_dataset(self.dataset_file_objects_for_available_archives) if x.startswith(selected_archive_with_proper_extension)])  # there must be at least one

                    # Download the files from the dataset in parallel
                    all_downloaded_files = nidap_io.download_files_from_dataset(nidap_io.get_foundry_dataset(alias='output'), dataset_filter_func=lambda f: f.path in matching_archives_files, limit=15)
                    local_download_paths = [all_downloaded_files[zip_file_chunk] for zip_file_chunk in matching_archives_files]

                    # Extract all downloaded parts
                    extract_zipfile_to_directory(filepaths=local_download_paths, extraction_path=local_output_dir)
    
                # Mimic (sort of) callback behavior, especially because we want to see updated available session states
                st.rerun()

    # Save the currently loaded results to a results archive
    # def save_results_to_archive(self, nworkers_for_data_transfer=1):
    def save_results_to_archive(self):

        st.subheader(':tractor: Save results')

        # Allow the user to add a custom basename for the new archive
        st.text_input('Suffix for the basename of the new results archive to create:', key='basename_suffix_for_new_results_archive', help='The name will go after "output_archive" and a timestamp, e.g., "output_archive-20230920_003801-project_xxxx_panel_07".')

        # If the user is ready to save the current results to a new results archive...
        if st.button(':arrow_left: Save current results to a new archive', help='This will copy all current results to a new archive, including job settings and environment information.'):

            # Copy a YAML file of the current tool settings to the current/loaded results
            write_current_tool_parameters_to_disk(local_output_dir)

            # Save the current environment to the current/loaded results
            write_current_environment_to_disk(local_output_dir)

            # Delete any files in the local output directory that start with "streamlit_session_state-" because we're about to create a current one and we don't want to back up more than one as they're generally large
            delete_selected_files_and_dirs(local_output_dir, [x for x in os.listdir(local_output_dir) if x.startswith('streamlit_session_state-')])

            # Save the current session state to the current/loaded results
            memory_analyzer.save_session_state(local_output_dir)

            # If working locally...
            if self.platform == 'local':

                # Copy everything in the local output directory except for files/dirs like ^output_archive- to a new archive directory
                copy_output_dir_contents_to_output_archive(st.session_state['basename_suffix_for_new_results_archive'], local_output_dir)

            # If working on NIDAP...
            elif self.platform == 'nidap':

                # Back up everything in the local output directory to the output dataset on NIDAP
                back_up_results_to_nidap(local_output_dir, st.session_state['basename_suffix_for_new_results_archive'])

            # Rerun since this potentially changes outputs
            st.rerun()
            
    # List all currently loaded results that aren't output archives, which is platform-independent
    def get_local_results_listing(self):
        return sorted([x for x in os.listdir(local_output_dir) if not x.startswith('output_archive-')])  # only locally will there exist files/dirs that start with output_archive- but it doesn't hurt to keep this here
    
    # Write a dataframe of the results in the local output directory, also obviously platform-independent
    def display_local_results_df(self):
        st.subheader(':open_file_folder: Results loaded in the tool')
        make_complex_dataframe_from_file_listing(local_output_dir, self.get_local_results_listing(), df_session_state_key_basename='local_results', editable=True)

    # Delete selected items from the output results directory
    def add_delete_local_results_button(self):

        # If the user wants to delete selected items...
        if st.button(':x: Delete selected (above) results files or directories'):

            # Store the output results dataframe
            df_local_results = st.session_state['loader__de_local_results'].reconstruct_edited_dataframe()

            # Get just the file or directory names they want to delete
            selected_items_to_delete = df_local_results[df_local_results['Selected']]['File or directory name']

            # Delete them
            delete_selected_files_and_dirs(local_output_dir, selected_items_to_delete)
    
            # Rerun since this potentially changes outputs
            st.rerun()
        
    # Write a YAML file of the current tool parameters to the loaded results directory
    def write_settings_to_local_results(self):
        st.subheader(':tractor: Write current tool parameters to loaded results')
        if st.button(':pencil2: Write current tool settings to the results directory', help='Note you can subsequently load these parameters from the "Tool parameter selection" tab at left'):
            write_current_tool_parameters_to_disk(local_output_dir)
            st.rerun()  # rerun since this potentially changes outputs

    # Write a YAML file of the current environment to the loaded results directory
    def write_environment_to_local_results(self):
        st.subheader(':tractor: Write current conda/pip environment to loaded results')
        if st.button(':pencil2: Write current environment to the results directory'):
            write_current_environment_to_disk(local_output_dir)
            st.rerun()  # rerun since this potentially changes outputs
    
    # Create a local empty results archive
    def add_create_empty_archive_button(self):

        st.subheader(':tractor: Create empty results archive')

        # Allow the user to add a custom basename to the new local directory name
        st.text_input('Suffix for the basename of the local archive directory to create:', key='basename_suffix_for_new_local_archive_dir', help='The name will go after "output_archive-" and a timestamp, e.g., "output_archive-20230920_003801-project_xxxx_panel_07".')

        # Create a new local directory with that basename suffix
        if st.button(':pencil2: Create empty results archive directory'):
            _ = create_empty_output_archive(st.session_state['basename_suffix_for_new_local_archive_dir'], local_output_dir)
            st.rerun()  # rerun since this potentially changes outputs
    
    # Delete local empty results output archive directories
    def add_delete_empty_archives_button(self):

        st.subheader(':tractor: Delete empty results archives')

        # If the user wants to delete all empty local results archives...
        if st.button(':x: Delete empty local results archive directories'):
        
            # Store the local results dataframe
            df_local_results = st.session_state['loader__de_local_results'].reconstruct_edited_dataframe()

            # Obtain the directories to delete as the empty ones that start with "output_archive-"
            dirs_to_delete = df_local_results[(df_local_results['# of files within'] == 0) & (df_local_results['File or directory name'].apply(lambda x: x.startswith('output_archive-')))]['File or directory name']

            # Delete the empty archive directories
            delete_selected_files_and_dirs(local_output_dir, dirs_to_delete)

            # Rerun since this potentially changes outputs
            st.rerun()

    # Create a snapshot of the session state in the "output" directory
    def create_session_state_snapshot(self):
        st.subheader(':tractor: Create snapshot of session state')
        st.button(':camera: Create snapshot', on_click=memory_analyzer.save_session_state, args=(local_output_dir))

# Determine whether a full string contains any of a tuple of substrings
def multi_contains(full_str, substrs):
    return (True if sum([substr in full_str for substr in substrs]) > 0 else False)

# List files/links within a top directory
def get_recursive_file_listing_of_directory(topdir, dirpath_prefixes_to_exclude=(), dirpath_suffixes_to_exclude=(), dirpath_substrs_to_exclude=(), filename_prefixes_to_exclude=(), filename_suffixes_to_exclude=(), filename_substrs_to_exclude=()):
    # Sample usage: platform_io.get_recursive_file_listing_of_directory(os.path.join('.', 'config'))

    # Import relevant library
    import os

    # Initialize a list holding the file listing
    file_listing = []

    # For every directory within topdir...
    for dirpath, _, filenames in os.walk(topdir):  # equivalent of "find . -type f,l"

        # For every file or link in the current directory...
        for filename in filenames:

            # If the current directory or filename does not include any strings that we want to exclude...
            if (not dirpath.startswith(dirpath_prefixes_to_exclude)) and (not dirpath.endswith(dirpath_suffixes_to_exclude)) and (not multi_contains(dirpath, dirpath_substrs_to_exclude)) and (not filename.startswith(filename_prefixes_to_exclude)) and (not filename.endswith(filename_suffixes_to_exclude)) and (not multi_contains(filename, filename_substrs_to_exclude)):

                # Add it to the running file listing
                file_listing.append(os.path.join(dirpath, filename))

    # Return the full file listing
    return file_listing

# From a file string that's either a path or just the filename, get the directory name and pure filename
def get_dirname_and_basename_from_file(file_str):
    import os
    if os.path.sep in file_str:
        file_dirname = os.path.dirname(file_str)
        file_basename = os.path.basename(file_str)
    else:
        file_dirname = '.'
        file_basename = file_str
    return file_dirname, file_basename

# Append the number of files in a group of files to all members of the group
def append_group_size_to_all_files_in_group(file_path):
    # zipfile_path can be e.g. (1) os.path.join('..', 'my_zipfile.zip') or (2) 'my_zipfile.zip'

    # Import relevant library
    import os

    # Get the directory and basename of the file
    file_dirname, file_basename = get_dirname_and_basename_from_file(file_path)

    # Get the paths to all the matching zipfile parts
    file_parts_paths = sorted([os.path.join(file_dirname, x) for x in os.listdir(file_dirname) if x.startswith(file_basename)])

    # Determine the number of matches
    num_parts = len(file_parts_paths)

    # For every matching path...
    for curr_path in file_parts_paths:

        # Append the number of files in the group to the end of the file
        os.rename(curr_path, '{}_of_{:03d}'.format(curr_path, num_parts))

    # Print what we did
    print('{} files renamed'.format(num_parts))

    # Return the number of files in the group
    return(num_parts)

# Create a zip file from an iterable of filenames
def zipfile_creation_from_filenames(file_or_buffer, filenames):

    # Import relevant library
    import zipfile
    
    # Open the zipfile for writing using the zipfile library
    with zipfile.ZipFile(file=file_or_buffer, mode='x', compression=zipfile.ZIP_DEFLATED) as myzip:

        # For every file to write to the archive...
        for filename in filenames:

            # Print what we're doing for the current file
            print('Writing {} to zip'.format(filename))

            # Add it to the archive
            myzip.write(filename=filename)

# Extract a zip file to a directory
def zipfile_extraction_to_a_directory(file_or_buffer, extraction_path):

    # Import relevant library
    import zipfile

    # Open the zipfile for extraction using the zipfile library
    with zipfile.ZipFile(file=file_or_buffer, mode='r') as myzip:

        # Extract all contents of the zip file to the desired directory
        myzip.extractall(path=extraction_path)

# Create a zip file from the files present in a directory of interest, with possible exclusions
def create_zipfile_from_files_in_dir(zipfile_name, topdir, chunksize_in_mb=None, dirpath_prefixes_to_exclude=(), dirpath_suffixes_to_exclude=(), dirpath_substrs_to_exclude=(), filename_prefixes_to_exclude=(), filename_suffixes_to_exclude=(), filename_substrs_to_exclude=()):
    # Sample usage:
    #   platform_io.create_zipfile_from_files_in_dir('/home/weismanal/projects/spatial-interaction-tool/app-dev/repo/config.zip', os.path.join('output', 'output_archive-probably_good_recent_lci_results_from_original_dataset-20230921_020450'), filename_prefixes_to_exclude=('original_gmb_phenotype_ids',))
    #   platform_io.create_zipfile_from_files_in_dir('../dude2.zip', 'output/output_archive-probably_good_recent_lci_results_from_original_dataset-20230921_020450', chunksize_in_mb=250)

    # Import relevant libraries
    import os
    import split_file_reader.split_file_writer

    # Get the absolute path of the zip file
    zipfile_name = os.path.abspath(zipfile_name)

    # If the zip file(s) to create does not already exist...
    if chunksize_in_mb is None:
        zipfiles_exist = os.path.exists(zipfile_name)
    else:
        zipfile_dirname, zipfile_basename = get_dirname_and_basename_from_file(zipfile_name)
        zipfiles_exist = len([x for x in os.listdir(zipfile_dirname) if x.startswith(zipfile_basename)]) != 0
    if not zipfiles_exist:

        # Enter the directory to zip and update the "top directory" accordingly
        orig_cwd = os.getcwd()
        os.chdir(topdir)
        topdir = '.'

        # Recursively list the files present, with possible exclusions
        filenames = get_recursive_file_listing_of_directory(topdir, dirpath_prefixes_to_exclude=dirpath_prefixes_to_exclude, dirpath_suffixes_to_exclude=dirpath_suffixes_to_exclude, dirpath_substrs_to_exclude=dirpath_substrs_to_exclude, filename_prefixes_to_exclude=filename_prefixes_to_exclude, filename_suffixes_to_exclude=filename_suffixes_to_exclude, filename_substrs_to_exclude=filename_substrs_to_exclude)

        # If we don't want to do chunking, create a single zipfile without the splitting library
        if chunksize_in_mb is None:
            zipfile_creation_from_filenames(file_or_buffer=zipfile_name, filenames=filenames)
            num_parts = 1
            suffix = ''

        # If we do want to do chunking, create a series of zipfiles using the splitting library
        else:
            with split_file_reader.split_file_writer.SplitFileWriter(zipfile_name + '.', (chunksize_in_mb * 1024**2)) as sfw:
                zipfile_creation_from_filenames(file_or_buffer=sfw, filenames=filenames)
            
            # Append the number of files to each created file
            num_parts = append_group_size_to_all_files_in_group(zipfile_name)
            suffix = '.*'

        # Print what we just did
        print('{} zip file(s) {}{} created from {} files'.format(num_parts, zipfile_name, suffix, len(filenames)))

        # Change back to the original directory
        os.chdir(orig_cwd)

    # If the zip file to create already exists, say so
    else:
        print('Zip file(s) {} NOT created because the file(s) already exists'.format(zipfile_name))

# Extract a zip file to a directory
def extract_zipfile_to_directory(zipfile_name='', extraction_path='', filepaths=None):
    # Sample usage:
    #   platform_io.extract_zipfile_to_directory(zipfile_name='/home/weismanal/projects/spatial-interaction-tool/app-dev/repo/config.zip', extraction_path=os.path.join('/home/weismanal/windows_home/Downloads/test'))
    #   platform_io.extract_zipfile_to_directory(zipfile_name='../dude2.zip', extraction_path='./tmp2')
    #   Use either zipfile_name or filepaths!

    # Import relevant library
    import split_file_reader.split_file_reader
    import os

    # Get the list of zip files matching zipfile_name
    if filepaths is None:
        zipfile_dirname, zipfile_basename = get_dirname_and_basename_from_file(zipfile_name)
        filepaths = sorted([os.path.join(zipfile_dirname, x) for x in os.listdir(zipfile_dirname) if x.startswith(zipfile_basename)])
    else:
        zipfile_name = os.path.basename(filepaths[0]).split('.zip')[0] + '.zip'

    # Detect whether we're dealing with the standard, single-file zip or the chunked, multiple-file zip
    standard_zip = False
    for filepath in filepaths:
        if filepath.endswith('.zip'):
            standard_zip = True
            break

    # Run some checks
    if standard_zip:
        if len(filepaths) != 1:
            import sys
            print('ERROR: At least one of the detected files ends with ".zip" alone, but more than one file was detected! Detected files: {}'.format(filepaths))
            sys.exit()
    else:
        if sum([not x.endswith('.zip') for x in filepaths]) != len(filepaths):
            import sys
            print('ERROR: Not all detected files don\'t end purely with ".zip"! Detected files: {}'.format(filepaths))
            sys.exit()

    # If we're doing the standard unzipping of a single file, do so
    if standard_zip:
        zipfile_extraction_to_a_directory(zipfile_name, extraction_path)
        num_parts = 1
        suffix=''

    # Otherwise, we're unzipping a series of zip file parts
    else:
        with split_file_reader.split_file_reader.SplitFileReader(filepaths) as sfr:
            zipfile_extraction_to_a_directory(sfr, extraction_path)
        num_parts = len(filepaths)
        suffix='.*'

    # Print what we just did
    print('{} zip file(s) {}{} extracted to {}'.format(num_parts, zipfile_name, suffix, extraction_path))

def upload_single_file_to_dataset(args_as_single_tuple):
    import os
    import time
    import nidap_io
    dataset, filedir, filename = args_as_single_tuple
    print('Uploading zip file chunk {}...'.format(filename))
    filesize = os.path.getsize(os.path.join(filedir, filename)) / 1024 ** 2
    start_time = time.time()
    nidap_io.upload_file_to_dataset(dataset, selected_filepath=os.path.join(filedir, filename))
    duration = time.time() - start_time
    print('  Upload of {} ({:5.3f} MB) from Workspaces to Compass took {:3.1f} seconds --> {:3.1f} MB/s'.format(filename, filesize, duration, filesize / duration))

def back_up_results_to_nidap(local_output_dir, basename_suffix_for_new_results_archive, chunksize_in_mb=200):

    # Import relevant library
    import nidap_io

    # Create a temporary transfer directory to hold the generated zip files
    local_transfer_dir = os.path.join('.', 'transfer')
    os.makedirs(local_transfer_dir, exist_ok=True)

    # Get the file like ../transfer/myfile.zip for which parts will be created like ../transfer/myfile.zip.000_of_007 etc.
    zipfile_part_prefix = os.path.join(local_transfer_dir, 'output_archive-{}-{}.zip'.format(utils.get_timestamp(), basename_suffix_for_new_results_archive))

    # Zip all loaded results to a temporary directory
    print(f'Creating {chunksize_in_mb} MB zip file parts for backup to NIDAP...')
    create_zipfile_from_files_in_dir(
        zipfile_part_prefix,
        local_output_dir,
        chunksize_in_mb=chunksize_in_mb,  # note that 750 MB is the largest chunk filesize that can be reliably transferred without timeout issues on NIDAP's end. Making this much smaller (to 250 MB) though since I'm about to implement parallelization of file transfers between NIDAP and Workspaces and probably the more files, the better. Note I recently saw the timeout issue with 750 MB so I'm splitting the difference and calling it 500 MB for now. Note that on 3/[9-10]/24, I no longer such a network-related limit, and instead a more filesystem-y one at 2000 MB, so we could potentially increase this. However, decreasing it to 200 MB to create more files and utilize more parallelism.
        dirpath_prefixes_to_exclude=('output_archive-',)
        )

    # Initialize the output dataset to which to transfer the zipfile parts
    dataset = nidap_io.get_foundry_dataset(alias='output')

    # Recursively delete any files or directories in local_transfer_dir that do not start with os.path.basename(zipfile_part_prefix)
    for filename in os.listdir(local_transfer_dir):
        if not filename.startswith(os.path.basename(zipfile_part_prefix)):
            path = os.path.join(local_transfer_dir, filename)
            if os.path.isfile(path):
                os.remove(path)
            elif os.path.isdir(path):
                shutil.rmtree(path)

    # This will upload the files in parallel
    print(f'Uploading {len(os.listdir(local_transfer_dir))} zip file parts in parallel to NIDAP...')
    nidap_io.upload_dir_to_dataset(dataset, path_to_dir_to_upload=local_transfer_dir)

    # Delete the temporary transfer directory and its contents
    shutil.rmtree(local_transfer_dir)

    # Print what we just did
    print('All results backed up to NIDAP!')
