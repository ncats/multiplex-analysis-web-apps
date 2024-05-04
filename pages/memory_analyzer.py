# Import relevant libraries
import streamlit as st
import app_top_of_page as top
import streamlit_dataframe_editor as sde
import pandas as pd
from pympler import asizeof
import numpy as np
import pickle
import dill
import os
import utils
import time

# For each custom class, add a key-value pair where the class is the key and the value is a list of picklable attributes of that class. Only do this if the size of that attribute can be larger than 1 MB, which you can assess by using this app. See possible classes (at least as of 5/1/24) in the get_object_class function below, which is not used right now
picklable_attributes_per_class = {
    'dataset_formats.Standardized': ['data']
    }

# Constant
bytes_to_mb = 1024 ** 2


def get_object_class(value):

    # No longer actually used as of 20240502_1706. However, it still lists all the custom classes that are used/present in MAWA. Note that isinstance() doesn't work reliably, so we really need to convert the dtype to a string, which we do elsewhere in this script.

    import neighbors_counts_for_neighborhood_profiles as custom_module__neighbors_counts_for_neighborhood_profiles
    import benchmark_collector as custom_module__benchmark_collector
    import dataset_formats as custom_module__dataset_formats
    import platform_io as custom_module__platform_io
    import time_cell_interaction_lib as custom_module__time_cell_interaction_lib
    import SpatialUMAP as custom_module__SpatialUMAP
    import foundry_IO_lib as custom_module__foundry_IO_lib
    import streamlit_dataframe_editor as custom_module__streamlit_dataframe_editor
    import neighborhood_profiles as custom_module__neighborhood_profiles

    if isinstance(value, custom_module__streamlit_dataframe_editor.DataframeEditor):
        class_str = 'streamlit_dataframe_editor.DataframeEditor'
    elif isinstance(value, custom_module__neighbors_counts_for_neighborhood_profiles.dummySessionState):
        class_str = 'neighbors_counts_for_neighborhood_profiles.dummySessionState'
    elif isinstance(value, custom_module__benchmark_collector.benchmark_collector):
        class_str = 'benchmark_collector.benchmark_collector'
    elif isinstance(value, custom_module__dataset_formats.Standardized):
        class_str = 'dataset_formats.Standardized'
    elif isinstance(value, custom_module__platform_io.Platform):
        class_str = 'platform_io.Platform'
    elif isinstance(value, custom_module__time_cell_interaction_lib.TIMECellInteraction):
        class_str = 'time_cell_interaction_lib.TIMECellInteraction'
    elif isinstance(value, custom_module__SpatialUMAP.SpatialUMAP):
        class_str = 'SpatialUMAP.SpatialUMAP'
    elif isinstance(value, custom_module__SpatialUMAP.FitEllipse):
        class_str = 'SpatialUMAP.FitEllipse'
    elif isinstance(value, custom_module__foundry_IO_lib.foundry_IO_lib):
        class_str = 'foundry_IO_lib.foundry_IO_lib'
    elif isinstance(value, custom_module__streamlit_dataframe_editor.DataframeEditor):
        class_str = 'streamlit_dataframe_editor.DataframeEditor'
    elif isinstance(value, custom_module__neighborhood_profiles.NeighborhoodProfiles):
        class_str = 'neighborhood_profiles.NeighborhoodProfiles'
    elif isinstance(value, custom_module__neighborhood_profiles.UMAPDensityProcessing):
        class_str = 'neighborhood_profiles.UMAPDensityProcessing'
    else:
        class_str = '<"NATIVE">'

    return class_str


def deserialize_file_to_dict(filepath, serialization_lib, output_func=print, calculate_size_in_mem=False):
    # This is as fast as can be as long as calculate_size_in_mem == False
    with open(filepath, 'rb') as f:
        dict_to_load = serialization_lib.load(f)
    if calculate_size_in_mem:
        predicted_size_in_mb = asizeof.asizeof(dict_to_load) / bytes_to_mb
    else:
        predicted_size_in_mb = 0
    actual_size_in_mb = os.path.getsize(filepath) / bytes_to_mb
    output_func(f'Loading dict_to_load from {filepath} (which is {actual_size_in_mb} MB) of predicted size {predicted_size_in_mb:.2f} MB using serialization library `{serialization_lib}`')
    return dict_to_load


def serialize_dict_to_file(dict_to_save, filepath, serialization_lib, output_func=print, calculate_size_in_mem=False):
    # This is as fast as can be as long as calculate_size_in_mem == False and dict_to_save contains objects of types that are quickly dumped depending on the serialization library, which is the point of this module
    if calculate_size_in_mem:
        predicted_size_in_mb = asizeof.asizeof(dict_to_save) / bytes_to_mb
    else:
        predicted_size_in_mb = 0
    output_func(f'Saving dict_to_save ({len(dict_to_save)} objects) to {filepath} using serialization library `{serialization_lib}`. This should take around {predicted_size_in_mb:.2f} MB...', end='', flush=True)
    with open(filepath, 'wb') as f:
        serialization_lib.dump(dict_to_save, f)
    actual_size_in_mb = os.path.getsize(filepath) / bytes_to_mb
    output_func(f' {actual_size_in_mb:.2f} MB saved, which is {actual_size_in_mb - predicted_size_in_mb:.2f} MB larger than the predicted size.')


def load_session_state_from_disk(saved_streamlit_session_states_dir, saved_streamlit_session_state_prefix='streamlit_session_state-', saved_streamlit_session_state_key='session_selection', selected_session=None):

    # This essentially replaces streamlit_session_state_management.load_session_state()

    # This is fast as can be

    # Print what we're doing
    print('Loading session state...')

    # Get the selected session basename to load
    if selected_session is None:
        selected_session = st.session_state[saved_streamlit_session_state_key]

    # If no session file was explicitly input and if no session files exist, do nothing
    if selected_session is None:
        st.warning(f'{utils.get_timestamp(pretty=True)}: No session state files exist so none were loaded')
        return

    # Choose one of the selected sessions... if not a manually input one (nominally the most recent), then the one selected in the session state file selection dropdown
    filepath_without_extension = os.path.join(saved_streamlit_session_states_dir, saved_streamlit_session_state_prefix + selected_session)

    # Delete every key in the current session state except for the selected session
    for key in st.session_state.keys():
        if key != saved_streamlit_session_state_key:
            del st.session_state[key]

    # Load the state (as a dictionary) from the binary files
    pickle_dict = deserialize_file_to_dict(filepath_without_extension + '.pickle', pickle)
    dill_dict = deserialize_file_to_dict(filepath_without_extension + '.dill', dill)

    # Check that there is no key overlap between the pickle and dill dictionaries since we're about to combine them
    assert not set(pickle_dict.keys()) & set(dill_dict.keys()), f'There is a key overlap of keys between the pickle and dill dictionaries: {set(pickle_dict.keys()) & set(dill_dict.keys())}'

    # Combine the two dictionaries into one
    session_dict = {**pickle_dict, **dill_dict}

    # Load each key-value pair individually into session_state
    for key, value in session_dict.items():
        print(f'Loading {key} of type {type(value)}')
        st.session_state[key] = value

    # Output a success message
    print(f'{utils.get_timestamp(pretty=True)}: State loaded from {filepath_without_extension}.pickle/.dill, which is {os.path.getsize(filepath_without_extension) / 1024 ** 2:.2f} MB')


def write_session_state_to_disk(ser_serialization_lib, saved_streamlit_session_states_dir, saved_streamlit_session_state_prefix='streamlit_session_state-'):

    # This essentially replaces streamlit_session_state_management.save_session_state()... well probably use entire button block in main() function below

    # Print what we're doing
    print('Saving session state...')
    
    # Create the output directory for saving session state if it doesn't exist
    os.makedirs(saved_streamlit_session_states_dir, exist_ok=True)

    # Generate the basename of the files to write with the save date and time
    filepath_without_extension = os.path.join(saved_streamlit_session_states_dir, saved_streamlit_session_state_prefix + utils.get_timestamp())

    # Save large, picklable objects to disk using pickle, which is faster than dill for large objects. At this point, all large objects have necessarily been made picklable using the functionality in this script
    pickle_dict = {key: st.session_state[key] for key in ser_serialization_lib[ser_serialization_lib == 'pickle'].index}
    serialize_dict_to_file(pickle_dict, filepath_without_extension + '.pickle', pickle)

    # Save small objects to disk using dill, which can serialize custom objects but is slow for large objects. At this point, all custom objects have been made small enough to be saved efficiently using dill using the functionality in this script
    dill_dict = {key: st.session_state[key] for key in ser_serialization_lib[ser_serialization_lib == 'dill'].index}
    serialize_dict_to_file(dill_dict, filepath_without_extension + '.dill', dill)

    # Output a success message
    print(f'{utils.get_timestamp(pretty=True)}: State saved to {filepath_without_extension}.pickle/.dill, which is {os.path.getsize(filepath_without_extension) / 1024 ** 2:.2f} MB')


def recombine_picklable_attributes_with_custom_object(ser_memory_usage_in_mb, update_memory_usage=True):

    # This is fast except when asizeof is called, which should be infrequent.

    # Initialize a list of main objects to which we will set the picklable attributes
    if update_memory_usage:
        main_objects = []

    # For every item in ser_memory_usage_in_mb...
    for key in ser_memory_usage_in_mb.index:

        # If the key is a separately-saved data attribute, then recombine it with the corresponding object
        if key.startswith('memory_analyzer__'):

            # Get the key of the object with which we want to recombine the picklable attribute
            key_of_object = key.removeprefix('memory_analyzer__key_').split('_class_')[0]

            # Get the name of the picklable attribute
            attribute_name = key.split('_attribute_')[-1]

            # Make sure the key of the object exists in the session state
            assert key_of_object in st.session_state, f'Key {key_of_object} does not exist in the session state, so we cannot recombine the {attribute_name} attribute with its corresponding object.'

            # Recombine the data attribute with the object
            setattr(st.session_state[key_of_object], attribute_name, st.session_state[key])

            # Delete the separately-saved data attribute
            del st.session_state[key]

            # Only if we want to update the memory usage...
            if update_memory_usage:
            
                # Delete the memory usage of the separately-saved data attribute
                ser_memory_usage_in_mb.drop(key, inplace=True)

                # Store the current main object
                main_objects.append(key_of_object)

    # Only if we want to update the memory usage...
    if update_memory_usage:
    
        # For every main object...
        for main_object in list(set(main_objects)):
                
            # Store the memory usage of the current object
            ser_memory_usage_in_mb[main_object] = asizeof.asizeof(st.session_state[main_object]) / bytes_to_mb

        # Return the updated memory usage series
        return ser_memory_usage_in_mb


def split_off_picklable_attributes_from_custom_object(ser_memory_usage_in_mb, output_func=print):

    # See comments elsewhere in this file but point is that for keys that are: (1) in st.session_state, (2) large, and (3) custom, we have problems.
    
    # This is only done when the object is "large" (> 1 MB).
    
    # This is fast except when asizeof is called, which should be infrequent.

    # For every item in ser_memory_usage_in_mb...
    for key in ser_memory_usage_in_mb.index:

        # Store the type of the current object
        type_str = str(type(st.session_state[key]))

        # If the current object in the session state is large...
        if (ser_memory_usage_in_mb[key] > 1):
            
            # For every custom class containing potentially large picklable attributes...
            for class_str in picklable_attributes_per_class.keys():

                # If the current object is of the current custom class...
                if type_str == f"<class '{class_str}'>":

                    # For every picklable attribute of the current object...
                    for picklable_attribute in picklable_attributes_per_class[class_str]:

                        # Write what we're doing
                        output_func(f'Key {key} in the session state has a value that is large in size (> 1 MB) and is of custom format {class_str}. Storing the "{picklable_attribute}" attribute separately (this may or may not contribute to the large size)...')

                        # Store the picklable attribute separately in the session state and store its memory usage
                        attribute_key = 'memory_analyzer__key_' + key + f'_class_{class_str}_attribute_{picklable_attribute}'
                        st.session_state[attribute_key] = getattr(st.session_state[key], picklable_attribute)

                        # Delete the picklable attribute from the current object
                        delattr(st.session_state[key], picklable_attribute)

                        # Store the memory usage of the standalone picklable attribute
                        ser_memory_usage_in_mb[attribute_key] = asizeof.asizeof(st.session_state[attribute_key]) / bytes_to_mb  # note this is slow

                    # Store the new memory usage of the current object
                    ser_memory_usage_in_mb[key] = asizeof.asizeof(st.session_state[key]) / bytes_to_mb  # note this is slow

    # Return the updated memory usage series
    return ser_memory_usage_in_mb


def assess_whether_same_object(df):

    # String to indicate that the previous is not the same as the current
    different_string = 'Nothing is the same as above'

    # Concatenate together all the values (which are lists) of the picklable_attributes_per_class dictionary, storing the result in attribute_list
    attribute_list = list(set([item for attribute_list_per_class in picklable_attributes_per_class.values() for item in attribute_list_per_class]))

    # Initialize the comparison series
    ser_is_same_as_above = pd.Series(None, name='is_same_as_above', index=df.index)
    ser_comparison_string = pd.Series(None, name='comparison_string', index=df.index)

    # For all but the first row of the dataframe...
    for i in range(1, len(df)):

        # If the object is large enough to be saved using pickle, i.e., is large...
        if df.iloc[i]['size_mb'] > 1:

            # Get the current and previous data objects
            current_data = st.session_state.get(df.iloc[i]['key'], 'current_data')
            previous_data = st.session_state.get(df.iloc[i-1]['key'], 'previous_data')

            # Initialize comparison data for the current row
            ser_is_same_as_above.iloc[i] = False
            ser_comparison_string.iloc[i] = different_string

            # Get the items and their names to compare
            current_items = [current_data] + [getattr(current_data, attribute, f'current_{attribute}_attribute') for attribute in attribute_list]
            previous_items = [previous_data] + [getattr(previous_data, attribute, f'previous_{attribute}_attribute') for attribute in attribute_list]
            current_item_names = ['current top-level object'] + [f'current {attribute} attribute' for attribute in attribute_list]
            previous_item_names = ['previous top-level object'] + [f'previous {attribute} attribute' for attribute in attribute_list]

            # For each item, compare it to the corresponding item from the previous object
            for current_item, current_item_name in zip(current_items, current_item_names):
                for previous_item, previous_item_name in zip(previous_items, previous_item_names):
                    if current_item is previous_item:
                        ser_is_same_as_above.iloc[i] = True
                        same_string = f'{current_item_name} is the same as {previous_item_name}'
                        if ser_comparison_string.iloc[i] == different_string:
                            ser_comparison_string.iloc[i] = same_string
                        else:  # if there's more than one match, which shouldn't happen, then don't overwrite the first match, just add to a running string of matches
                            ser_comparison_string.iloc[i] = f'{ser_comparison_string.iloc[i]} & {same_string}'

    # Return both series
    return ser_is_same_as_above, ser_comparison_string


def get_session_state_object_info(ser_memory_usage_in_mb, return_val=None, write_dataframe=False):

    # This is generally slow because asizeof is slow for large dataframes. However, if the elements in ser_memory_usage_in_mb are not None, then that function is not called and therefore this function is much faster.

    # If we don't want to return anything from this function, then we must mean we want to write the dataframe to the screen (and note that therefore we need to calculate everything in the dataframe)
    if not return_val:
        write_dataframe = True
    
    # Initialize some of the dataframe columns
    key_holder = []
    if write_dataframe:
        type_holder2 = []
    size_holder = []

    # For every item in ser_memory_usage_in_mb, save the key, type, and size of the current object
    for key in ser_memory_usage_in_mb.index:
        key_holder.append(key)
        if write_dataframe:
            type_holder2.append(str(type(st.session_state[key])))
        if not np.isnan(ser_memory_usage_in_mb[key]):
            size_holder.append(ser_memory_usage_in_mb[key])
        else:
            size_holder.append(asizeof.asizeof(st.session_state[key]) / bytes_to_mb)  # note this is slow

    # Create a dataframe from these data, sorting by decreasing size
    if write_dataframe:
        df_ss_object_info = pd.DataFrame({'key': key_holder, 'type': type_holder2, 'size_mb': size_holder}).sort_values(by='size_mb', ascending=False)
    else:
        df_ss_object_info = pd.DataFrame({'key': key_holder, 'size_mb': size_holder}).sort_values(by='size_mb', ascending=False)

    # Assess whether the current object is the same as the previous object
    if write_dataframe:
        ser_is_same_as_above, ser_comparison_string = assess_whether_same_object(df_ss_object_info)

    # Determine whether the key should be saved with pickle (large objects of "native" dtypes) or dill (small objects of any dtype)
    ser_pickle_or_dill = df_ss_object_info['size_mb'].apply(lambda x: 'pickle' if x > 1 else 'dill')  # this takes trivial time so we'll always do it even though we don't need to if return_val == 'memory'
    ser_pickle_or_dill.name = 'pickle_or_dill'

    # Add the resulting three columns to the informational dataframe
    if write_dataframe:
        df_ss_object_info = pd.concat([df_ss_object_info, ser_is_same_as_above, ser_comparison_string, ser_pickle_or_dill], axis='columns')
    else:
        df_ss_object_info = pd.concat([df_ss_object_info, ser_pickle_or_dill], axis='columns')

    # Check that no rows were added during the concatenation
    if write_dataframe:
        assert len(df_ss_object_info) == len(ser_is_same_as_above) == len(ser_comparison_string) == len(ser_pickle_or_dill), f'Lengths of dataframes (df_to_write: {len(df_ss_object_info)}, ser_is_same_as_above: {len(ser_is_same_as_above)}, ser_comparison_string: {len(ser_comparison_string)}, ser_pickle_or_dill: {len(ser_pickle_or_dill)}) do not match.'
    else:
        assert len(df_ss_object_info) == len(ser_pickle_or_dill), f'Lengths of dataframes (df_to_write: {len(df_ss_object_info)}, ser_pickle_or_dill: {len(ser_pickle_or_dill)}) do not match.'

    # Write the dataframe to the screen if desired
    if write_dataframe:
        st.dataframe(df_ss_object_info)

    # Return the desired value, if any
    if return_val == 'memory':
        return df_ss_object_info.set_index('key')['size_mb']
    elif return_val == 'serialization library':
        return df_ss_object_info.set_index('key')['pickle_or_dill']


def initialize_memory_usage_series(saved_streamlit_session_state_key='session_selection'):

    # This is fast

    # Initialize a holder for the selected keys in the session state
    key_holder = []

    # For every item in the session state...
    for key in st.session_state:

        # If the current key is like '__do_not_persist' or 'FormSubmitter:' or if we're looking at the saved_streamlit_session_state_key, then skip. Otherwise, store the key in the key_holder list
        if (not key.endswith('__do_not_persist')) and (not key.startswith('FormSubmitter:')) and (key != saved_streamlit_session_state_key):
            key_holder.append(key)

    # Create a series with the keys as the index and the default value of None, which gets converted to np.nan when the series is created
    ser_memory_usage_in_mb = pd.Series(None, index=key_holder)

    # Return the series
    return ser_memory_usage_in_mb


def main():
    """
    Main function for the memory analyzer page.
    """

    # The idea is that pickle cannot generally pickle objects from custom classes so we want to use dill for that. However, dill is much slower for large data structures such as large pandas dataframes, which are supported by pickle. Therefore, we want to use pickle for large "native" objects and use dill for everything else. Therefore we need to know both the type and size of each object and save the small ones using dill and save the large ones using pickle. However, sometimes there are large objects that are of custom classes, which we can't save using pickle. Therefore, we need to save the data attribute of such objects separately and then delete the data attribute from the object before saving the object using pickle. Upon loading the object, we need to check for the existence of the separately-saved data attribute and re-combine it with the object.

    # Parameters
    saved_streamlit_session_state_key = 'session_selection'
    saved_streamlit_session_states_dir = 'saved_streamlit_session_states'
    saved_streamlit_session_state_prefix = 'streamlit_session_state-'

    if st.button('Reset screen'):
        pass
    
    if st.button('Analyze memory usage of relevant objects in the session state'):
    
        # Initialize the memory usage series so this is the only function that iterates through the keys in the session state; the rest iterate over the index in ser_memory_usage_in_mb
        ser_memory_usage_in_mb = initialize_memory_usage_series(saved_streamlit_session_state_key=saved_streamlit_session_state_key)

        # Calculate the memory used by every object and write the session state object information to screen
        st.write(f'Initial session state object information ({len(ser_memory_usage_in_mb)} relevant objects):')
        ser_memory_usage_in_mb = get_session_state_object_info(ser_memory_usage_in_mb, return_val='memory', write_dataframe=True)

        # Split off the data attribute from the dataset_formats objects
        ser_memory_usage_in_mb = split_off_picklable_attributes_from_custom_object(ser_memory_usage_in_mb, output_func=st.write)

        # Write the session state object information to screen
        st.write(f'Session state object information after splitting off picklable attributes from large custom objects ({len(ser_memory_usage_in_mb)} relevant objects):')
        get_session_state_object_info(ser_memory_usage_in_mb, return_val=None, write_dataframe=True)

        # Recombine the data attribute with the dataset_formats objects
        ser_memory_usage_in_mb = recombine_picklable_attributes_with_custom_object(ser_memory_usage_in_mb)

        # Write the session state object information to screen
        st.write(f'Session state object information after recombining picklable attributes with custom objects ({len(ser_memory_usage_in_mb)} relevant objects):')
        get_session_state_object_info(ser_memory_usage_in_mb, return_val=None, write_dataframe=True)

    if st.button('Save objects in the session state to disk using pickle and dill'):

        # This whole block is as fast as it can be

        start_time = time.time()


        # Initialize the memory usage series so this is the only function that iterates through the keys in the session state; the rest iterate over the index in ser_memory_usage_in_mb
        ser_memory_usage_in_mb = initialize_memory_usage_series(saved_streamlit_session_state_key=saved_streamlit_session_state_key)  # fast

        # Calculate the memory used by every object
        ser_memory_usage_in_mb = get_session_state_object_info(ser_memory_usage_in_mb, return_val='memory', write_dataframe=False)  # as fast as it can be

        # Split off the data attribute from the dataset_formats objects
        ser_memory_usage_in_mb = split_off_picklable_attributes_from_custom_object(ser_memory_usage_in_mb)  # as fast as it can be

        # Get the series specifying the serialization library to use for each relevant object in the session state
        ser_serialization_lib = get_session_state_object_info(ser_memory_usage_in_mb, return_val='serialization library', write_dataframe=False)  # fast

        # Save the session state to disk using pickle and dill
        write_session_state_to_disk(ser_serialization_lib, saved_streamlit_session_states_dir, saved_streamlit_session_state_prefix=saved_streamlit_session_state_prefix)  # as fast as it can be

        # Recombine the data attribute with the dataset_formats objects
        recombine_picklable_attributes_with_custom_object(ser_memory_usage_in_mb, update_memory_usage=False)  # fast


        st.write(f'Time to save objects in the session state to disk using pickle and dill: {time.time() - start_time:.2f} seconds')

    if st.button('Load objects to the session state from pickle+dill files on disk'):
        pass


# Run the main function
if __name__ == '__main__':

    # Set page settings
    page_name = 'Memory Analyzer'
    st.set_page_config(layout='wide', page_title=page_name)
    st.title(page_name)

    # Run streamlit-dataframe-editor library initialization tasks at the top of the page
    st.session_state = sde.initialize_session_state(st.session_state)

    # Run Top of Page (TOP) functions
    st.session_state = top.top_of_page_reqs(st.session_state)

    # Call the main function
    main()

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)
