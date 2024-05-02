# Import relevant libraries
import streamlit as st
import app_top_of_page as top
import streamlit_dataframe_editor as sde
import pandas as pd
from pympler import asizeof


bytes_to_mb = 1024 ** 2


def get_object_class(value):

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


def split_off_dataset_formats_data_attribute(memory_usage_in_mb, saved_streamlit_session_state_key='session_selection'):

    # For every item in the session state...
    for key, value in st.session_state.items():

        # If we don't want to save the key (__do_not_persist and FormSubmitter keys) or if we're looking at the saved_streamlit_session_state_key, then skip
        if (not key.endswith('__do_not_persist')) and (not key.startswith('FormSubmitter:')) and (key != saved_streamlit_session_state_key):

            # Get the type and size of the current object
            # type1 = type(st.session_state[key])
            type2 = get_object_class(value)
            # predicted_size_in_mb = asizeof.asizeof(value) / bytes_to_mb  # note this is slow
            predicted_size_in_mb = memory_usage_in_mb[key]

            # If in st.session_state, are large, and are custom, we have problems, see overall comment above
            if (predicted_size_in_mb > 1) and (type2 == 'dataset_formats.Standardized'):

                # Write what we're doing
                st.write(f'Key {key} in the session state has a value that is large in size (> 1 MB) and is of custom format {type2}. Storing the "data" attribute separately...')

                # Store the data attribute (which is a pandas.DataFrame) separately in the session state
                st.session_state[key + '_dataset_formats_data_attribute'] = value.data

                # Delete the data attribute from the current object
                del st.session_state[key].data


def recombine_data_attribute_with_dataset_formats_object():

    # For every item in the session state...
    for key, value in st.session_state.items():

        # If the key is a separately-saved data attribute, then recombine it with the corresponding object
        if key.endswith('_dataset_formats_data_attribute'):

            # Get the key of the object that we want to recombine the data attribute with
            key_of_object = key.replace('_dataset_formats_data_attribute', '')

            # Make sure the key of the object exists in the session state
            assert key_of_object in st.session_state, f'Key {key_of_object} does not exist in the session state, so we cannot recombine the data attribute {key} with its corresponding object.'

            # Recombine the data attribute with the object
            st.session_state[key_of_object].data = value

            # Delete the separately-saved data attribute
            del st.session_state[key]


def write_dataframe_info_and_get_memory_usage(saved_streamlit_session_state_key='session_selection', return_val=None):
    key_holder = []
    type_holder1 = []
    type_holder2 = []
    size_holder = []
    for key, value in st.session_state.items():
        if (not key.endswith('__do_not_persist')) and (not key.startswith('FormSubmitter:')) and (key != saved_streamlit_session_state_key):
            type1 = type(st.session_state[key])
            type2 = get_object_class(value)
            predicted_size_in_mb = asizeof.asizeof(value) / bytes_to_mb  # note this is slow
            key_holder.append(key)
            type_holder1.append(type1)
            type_holder2.append(type2)
            size_holder.append(predicted_size_in_mb)
    df_to_write = pd.DataFrame({'key': key_holder, 'type1': type_holder1, 'type2': type_holder2, 'size_mb': size_holder}).sort_values(by='size_mb', ascending=False)
    ser_is_same_as_above = assess_whether_same_object(df_to_write)
    ser_pickle_or_dill = df_to_write['size_mb'].apply(lambda x: 'pickle' if x > 1 else 'dill')
    ser_pickle_or_dill.name = 'pickle_or_dill'
    df_to_write = pd.concat([df_to_write, ser_is_same_as_above, ser_pickle_or_dill], axis='columns')
    assert len(df_to_write) == len(ser_is_same_as_above) == len(ser_pickle_or_dill), f'Lengths of dataframes (df_to_write: {len(df_to_write)}, ser_is_same_as_above: {len(ser_is_same_as_above)}, ser_pickle_or_dill: {len(ser_pickle_or_dill)}) do not match.'
    st.dataframe(df_to_write)

    if return_val == 'memory':
        return df_to_write.set_index('key')['size_mb']
    elif return_val == 'dataframe':
        return df_to_write


def assess_whether_same_object(df):
    # Note that if everything is commented out except for `return df`, the same number of session state keys exists, but instead all are "keys that would actually be saved" and there are none that "would not actually be saved." This is true for calls below like `df_do_not_save = assess_whether_same_object(df_do_not_save)` and `df_do_save = assess_whether_same_object(df_do_save)`. Makse no sense to me. Same goes if I merely instead say `df_do_not_save = df_do_not_save` and `df_do_save = df_do_save`.
    # Maybe this won't happen anymore since I'm just creating a separate series so I'm not modifying df in any way
    ser_is_same_as_above = pd.Series(None, name='is_same_as_above', index=df.index)
    for i in range(1, len(df)):
        if df.iloc[i]['size_mb'] > 1:
            current_key = df.iloc[i]['key']
            previous_key = df.iloc[i-1]['key']
            current_data = st.session_state.get(current_key, 'current_data')
            previous_data = st.session_state.get(previous_key, 'previous_data')
            current_data_attributes = getattr(current_data, 'data', 'current_data_attributes')
            previous_data_attributes = getattr(previous_data, 'data', 'previous_data_attributes')
            if (current_data is previous_data) or (current_data is previous_data_attributes) or (previous_data is current_data_attributes) or (current_data_attributes is previous_data_attributes):
                ser_is_same_as_above.iloc[i] = True
            else:
                ser_is_same_as_above.iloc[i] = False
    return ser_is_same_as_above


def main():
    """
    Main function for the page.
    """

    # The idea is that pickle cannot generally pickle objects from custom classes so we want to use dill for that. However, dill is much slower for large data structures such as large pandas dataframes, which are supported by pickle. Therefore, we want to use pickle for large "native" objects and use dill for everything else. Therefore we need to know both the type and size of each object and save the small ones using dill and save the large ones using pickle. However, sometimes there are large objects that are of custom classes, which we can't save using pickle. Therefore, we need to save the data attribute of such objects separately and then delete the data attribute from the object before saving the object using pickle. Upon loading the object, we need to check for the existence of the separately-saved data attribute and re-combine it with the object.

    # Define the key of the session state object that we want to save
    saved_streamlit_session_state_key = 'session_selection'

    # Write the session state object information to screen
    st.write('Initial session state object information:')
    memory_usage_in_mb = write_dataframe_info_and_get_memory_usage(saved_streamlit_session_state_key=saved_streamlit_session_state_key, return_val='memory')

    # Split off the data attribute from the dataset_formats objects
    split_off_dataset_formats_data_attribute(memory_usage_in_mb, saved_streamlit_session_state_key=saved_streamlit_session_state_key)

    # Write the session state object information to screen
    st.write('Session state object information after splitting off data attributes from dataset_formats objects:')
    write_dataframe_info_and_get_memory_usage(saved_streamlit_session_state_key=saved_streamlit_session_state_key)

    # Recombine the data attribute with the dataset_formats objects
    recombine_data_attribute_with_dataset_formats_object()

    # Write the session state object information to screen
    st.write('Session state object information after recombining data attributes with dataset_formats objects:')
    write_dataframe_info_and_get_memory_usage(saved_streamlit_session_state_key=saved_streamlit_session_state_key)


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
