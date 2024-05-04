# Import relevant libraries
import streamlit as st
import app_top_of_page as top
import streamlit_dataframe_editor as sde
import pandas as pd
from pympler import asizeof

# For each custom class, add a key-value pair where the class is the key and the value is a list of picklable attributes of that class. Only do this is the size of that attribute can be larger than 1 MB, which you can assess by using this app
picklable_attributes_per_class = {'dataset_formats.Standardized': ['data']}  # see possible classes (at least as of 5/1/24 in the get_object_class function, which is not used right now)

# Constant
bytes_to_mb = 1024 ** 2


def get_object_class(value):

    # No longer actually used as of 20240502_1706

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


def initialize_memory_usage_series(saved_streamlit_session_state_key='session_selection'):

    # Initialize a holder for the selected keys in the session state
    key_holder = []

    # For every item in the session state...
    for key in st.session_state:

        # If the current key is like '__do_not_persist' or 'FormSubmitter:' or if we're looking at the saved_streamlit_session_state_key, then skip. Otherwise, store the key in the key_holder list
        if (not key.endswith('__do_not_persist')) and (not key.startswith('FormSubmitter:')) and (key != saved_streamlit_session_state_key):
            key_holder.append(key)

    # Create a series with the keys as the index and the default value of None
    memory_usage_in_mb = pd.Series(None, index=key_holder)

    # Return the series
    return memory_usage_in_mb


def split_off_picklable_attributes_from_custom_object(memory_usage_in_mb):

    # See comments elsewhere in this file but point is that for keys: (1) in st.session_state, (2) are large, and (3) are custom, we have problems
    
    # This is only done when the object is "large" (> 1 MB).
    
    # This is fast except when asizeof is called, which should be infrequent.

    # For every item in memory_usage_in_mb...
    for key in memory_usage_in_mb.index:

        # Store the type of the current object
        type_str = str(type(st.session_state[key]))

        # If the current object in the session state is large...
        if (memory_usage_in_mb[key] > 1):
            
            # For every custom class containing potentially large picklable attributes...
            for class_str in picklable_attributes_per_class.keys():

                # If the current object is of the current custom class...
                if type_str == f"<class '{class_str}'>":

                    # For every picklable attribute of the current object...
                    for picklable_attribute in picklable_attributes_per_class[class_str]:

                        # Write what we're doing
                        st.write(f'Key {key} in the session state has a value that is large in size (> 1 MB) and is of custom format {class_str}. Storing the "{picklable_attribute}" attribute separately (this may or may not contribute to the large size)...')

                        # Store the picklable attribute separately in the session state and store its memory usage
                        attribute_key = 'memory_analyzer__key_' + key + f'_class_{class_str}_attribute_{picklable_attribute}'
                        st.session_state[attribute_key] = getattr(st.session_state[key], picklable_attribute)

                        # Delete the picklable attribute from the current object
                        delattr(st.session_state[key], picklable_attribute)

                        # Store the memory usage of the standalone picklable attribute
                        memory_usage_in_mb[attribute_key] = asizeof.asizeof(st.session_state[attribute_key]) / bytes_to_mb  # note this is slow

                    # Store the new memory usage of the current object
                    memory_usage_in_mb[key] = asizeof.asizeof(st.session_state[key]) / bytes_to_mb  # note this is slow

    # Return the updated memory usage series
    return memory_usage_in_mb


def recombine_picklable_attributes_with_custom_object(memory_usage_in_mb):

    # This is fast except when asizeof is called, which should be infrequent.

    # Initialize a list of main objects to which we will set the picklable attributes
    main_objects = []

    # For every item in memory_usage_in_mb...
    for key in memory_usage_in_mb.index:

        # If the key is a separately-saved data attribute, then recombine it with the corresponding object
        if key.startswith('memory__analyzer'):

            # Get the key of the object with which we want to recombine the picklable attribute
            key_of_object = key.removeprefix('memory__analyzer__key_').split('_class_')[0]

            # Get the name of the picklable attribute
            attribute_name = key.split('_attribute_')[-1]

            # Make sure the key of the object exists in the session state
            assert key_of_object in st.session_state, f'Key {key_of_object} does not exist in the session state, so we cannot recombine the {attribute_name} attribute with its corresponding object.'

            # Recombine the data attribute with the object
            setattr(st.session_state[key_of_object], attribute_name, st.session_state[key])

            # Delete the separately-saved data attribute
            del st.session_state[key]

            # Delete the memory usage of the separately-saved data attribute
            memory_usage_in_mb.drop(key, inplace=True)

            # Store the current main object
            main_objects.append(key_of_object)

    # For every main object...
    for main_object in list(set(main_objects)):
            
        # Store the memory usage of the current object
        memory_usage_in_mb[main_object] = asizeof.asizeof(st.session_state[main_object]) / bytes_to_mb

    # Return the updated memory usage series
    return memory_usage_in_mb


def get_session_state_object_info(memory_usage_in_mb, return_val=None, write_dataframe=False):

    # This is generally slow because asizeof is slow for large dataframes. However, if memory_usage_in_mb is not None, then that function is not called and therefore this function is much faster.

    # If we don't want to return anything from this function, then we must mean we want to write the dataframe to the screen (and note that therefore we need to calculate everything in the dataframe)
    if not return_val:
        write_dataframe = True
    
    # Initialize some of the dataframe columns
    key_holder = []
    if write_dataframe:
        type_holder2 = []
    size_holder = []

    # For every item in memory_usage_in_mb, save the key, type, and size of the current object
    for key in memory_usage_in_mb.index:
        key_holder.append(key)
        if write_dataframe:
            type_holder2.append(str(type(st.session_state[key])))
        if memory_usage_in_mb[key]:
            size_holder.append(memory_usage_in_mb[key])
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


def assess_whether_same_object(df):
    # Note that if everything is commented out except for `return df`, the same number of session state keys exists, but instead all are "keys that would actually be saved" and there are none that "would not actually be saved." This is true for calls below like `df_do_not_save = assess_whether_same_object(df_do_not_save)` and `df_do_save = assess_whether_same_object(df_do_save)`. Makse no sense to me. Same goes if I merely instead say `df_do_not_save = df_do_not_save` and `df_do_save = df_do_save`.
    # Maybe this won't happen anymore since I'm just creating a separate series so I'm not modifying df in any way

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


def write_session_state_to_disk():

    # Split off the data attribute from the dataset_formats objects in the session state
    memory_usage_in_mb = get_session_state_object_info(return_val='memory', write_dataframe=False)
    memory_usage_in_mb = split_off_picklable_attributes_from_custom_object(memory_usage_in_mb)

    # Determine whether to write each object in the session state to disk using pickle or dill
    serialization_lib = get_session_state_object_info(write_dataframe=False, return_val='serialization library')

    # Recombine the data attribute with the dataset_formats objects
    recombine_picklable_attributes_with_custom_object()

    # Write the session state object information to screen
    st.write('Session state object information after recombining data attributes with dataset_formats objects:')
    get_session_state_object_info(write_dataframe=False)


def main():
    """
    Main function for the memory analyzer page.
    """

    # The idea is that pickle cannot generally pickle objects from custom classes so we want to use dill for that. However, dill is much slower for large data structures such as large pandas dataframes, which are supported by pickle. Therefore, we want to use pickle for large "native" objects and use dill for everything else. Therefore we need to know both the type and size of each object and save the small ones using dill and save the large ones using pickle. However, sometimes there are large objects that are of custom classes, which we can't save using pickle. Therefore, we need to save the data attribute of such objects separately and then delete the data attribute from the object before saving the object using pickle. Upon loading the object, we need to check for the existence of the separately-saved data attribute and re-combine it with the object.

    memory_usage_in_mb = initialize_memory_usage_series()

    # Write the session state object information to screen
    st.write('Initial session state object information:')
    memory_usage_in_mb = get_session_state_object_info(return_val='memory', write_dataframe=True)

    # Split off the data attribute from the dataset_formats objects
    memory_usage_in_mb = split_off_picklable_attributes_from_custom_object(memory_usage_in_mb)

    # Write the session state object information to screen
    st.write('Session state object information after splitting off data attributes from dataset_formats objects:')
    get_session_state_object_info(write_dataframe=True)

    # Recombine the data attribute with the dataset_formats objects
    recombine_picklable_attributes_with_custom_object()

    # Write the session state object information to screen
    st.write('Session state object information after recombining data attributes with dataset_formats objects:')
    get_session_state_object_info(write_dataframe=True)


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
