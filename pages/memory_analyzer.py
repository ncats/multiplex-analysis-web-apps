# Import relevant libraries
import streamlit as st
import app_top_of_page as top
import streamlit_dataframe_editor as sde
import pandas as pd
import streamlit_dataframe_editor
from pympler import asizeof


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


def assess_whether_same_object(df):
    # Note that if everything is commented out except for `return df`, the same number of session state keys exists, but instead all are "keys that would actually be saved" and there are none that "would not actually be saved." This is true for calls below like `df_do_not_save = assess_whether_same_object(df_do_not_save)` and `df_do_save = assess_whether_same_object(df_do_save)`. Make no sense to me. Same goes if I merely instead say `df_do_not_save = df_do_not_save` and `df_do_save = df_do_save`.
    df = df.reset_index(drop=True)
    df['is_same_as_above'] = None
    for i in range(1, len(df)):
        if df.loc[i, 'size_mb'] < 1:
            df.loc[i, 'is_same_as_above'] = None
        else:
            current_key = df.loc[i, 'key']
            previous_key = df.loc[i-1, 'key']
            current_data = st.session_state.get(current_key, 'current_data')
            previous_data = st.session_state.get(previous_key, 'previous_data')
            current_data_attributes = getattr(current_data, 'data', 'current_data_attributes')
            previous_data_attributes = getattr(previous_data, 'data', 'previous_data_attributes')
            if (current_data is previous_data) or (current_data is previous_data_attributes) or (previous_data is current_data_attributes) or (current_data_attributes is previous_data_attributes):
                df.loc[i, 'is_same_as_above'] = True
            else:
                df.loc[i, 'is_same_as_above'] = False
    return df


def analyze_memory_usage(saved_streamlit_session_state_key='session_selection'):

    # This function is largely copied from streamlit_session_state_management.save_session_state()
    
    key_holder = []
    type_holder1 = []
    type_holder2 = []
    size_holder = []
    bytes_to_mb = 1024 ** 2
    for key, value in st.session_state.items():
        if (not key.endswith('__do_not_persist')) and (not key.startswith('FormSubmitter:')) and (key != saved_streamlit_session_state_key):
            type1 = type(st.session_state[key])
            type2 = get_object_class(value)
            predicted_size_in_mb = asizeof.asizeof(value) / bytes_to_mb  # note this is slow
            key_holder.append(key)
            type_holder1.append(type1)
            type_holder2.append(type2)
            size_holder.append(predicted_size_in_mb)

            # If in st.session_state, are large, and are custom, we have problems
            if (predicted_size_in_mb > 1) and (type2 == 'dataset_formats.Standardized'):
                st.write(f'Key {key} in the session state has a value that is large in size (> 1 MB) and is of custom format {type2}')
                # WRITE CODE HERE THAT WILL STORE THE DATA ATTRIBUTE SEPARATELY IN THE SESSION STATE AND DELETE IT FROM THE CURRENT OBJECT
                # THEN, rewrite dataframe showing size of each object to confirm the following will be correct
                # Save all eligible objects in the session state that are large (which should now be all native) using pickle, and save the rest of the eligible objects in the session state using dill
                # When loading in the resulting pickle and dill files, check for existence of the separately-saved data attribute and re-combine it with the corresponding object

    st.dataframe(pd.DataFrame({'key': key_holder, 'type1': type_holder1, 'type2': type_holder2, 'size_mb': size_holder}))

    st.write(f'Size of session state: {len(st.session_state)}')

    # Create a dictionary of most items in the session state
    session_dict = {}
    keys_to_exclude = []
    for key, value in st.session_state.items():
        if (not key.endswith('__do_not_persist')) and (not key.startswith('FormSubmitter:')) and (key != saved_streamlit_session_state_key):
            if isinstance(value, (sde.DataframeEditor, streamlit_dataframe_editor.DataframeEditor)):  # if this still doesn't seem to catch all DataframeEditor objects, try converting the type to a string and then checking if it contains 'DataframeEditor' or something like that
                print(f'Analyzing components of dataframe editor {key}')
                dataframe_editor_components = {
                    'df_name': value.df_name,
                    'default_df_contents': value.default_df_contents,
                    'edited_dataframe': value.reconstruct_edited_dataframe()
                }
                dataframe_editor_components_name = 'dataframe_editor_components__' + key
                session_dict[dataframe_editor_components_name] = dataframe_editor_components
                keys_to_exclude.append(value.df_name)
                keys_to_exclude.append(value.df_name + '_changes_dict')
                keys_to_exclude.append(value.df_name + '_key')
            else:
                print(f'Would save {key} of type {type(value)}')
                session_dict[key] = value

    # For each key to exclude, delete it from session_dict, i.e., the thing being saved to disk. This should allow a new key to be assigned to the st.data_editor() object within any dataframe editor objects that are later initialized in the load_session_state() function. Otherwise, the key will always be the same so upon loading a session state, the st.data_editor() object will not be redrawn. The solution to that is to force the st.data_editor() object to be redrawn by forcing its key to change.
    bytes_to_mb = 1024 ** 2
    tot_size_in_memory_do_not_save = 0
    keys_holder_do_not_save = []
    size_holder_do_not_save = []
    for key in keys_to_exclude:
        if key in session_dict:
            value = session_dict[key]
            predicted_size_in_mb = asizeof.asizeof(value) / bytes_to_mb
            print(f'Wouldn\'t actually save {key} of type {type(value)}, which would take around {predicted_size_in_mb:.2f} MB...')
            tot_size_in_memory_do_not_save += predicted_size_in_mb
            keys_holder_do_not_save.append(key)
            size_holder_do_not_save.append(predicted_size_in_mb)
            del session_dict[key]

    # Print the total size in memory
    print(f'Total size in memory (i.e., predicted) of keys that would not actually be saved: {tot_size_in_memory_do_not_save:.2f} MB')

    # Save the dictionary to the pickle file. Note this no longer randomly crashes with "PicklingError: Can't pickle <class 'streamlit_dataframe_editor.DataframeEditor'>: it's not the same object as streamlit_dataframe_editor.DataframeEditor" because we're no longer saving the DataframeEditor object itself, but rather the initialization data and the current contents. Note that using dill did also solve the problem, which if we were to use dill, we could try saving the entire session at once (instead of individual objects) and also thereby include difficult items such as st.form objects. NOTE: If we start getting the error again, try either using dill or probably better yet, excluding other custom object types from being saved in the first place, e.g., class 'platform_io.Platform'. Such exclusion would be done in keys_to_exclude as above.
    tot_size_in_memory_mb_do_save = 0
    keys_holder_do_save = []
    size_holder_do_save = []
    for key, value in session_dict.items():
        predicted_size_in_mb = asizeof.asizeof(value) / bytes_to_mb
        print(f'Would save {key} of type {type(value)} to a pickle file, which would take around {predicted_size_in_mb:.2f} MB...')
        tot_size_in_memory_mb_do_save += predicted_size_in_mb
        keys_holder_do_save.append(key)
        size_holder_do_save.append(predicted_size_in_mb)

    # Print the total size in memory
    print(f'Total size in memory (i.e., predicted) of keys that would indeed be saved: {tot_size_in_memory_mb_do_save:.2f} MB')

    # Write a dataframe of the keys and sizes that will not be saved, sorted in descending order of size
    df_do_not_save = pd.DataFrame({'key': keys_holder_do_not_save, 'size_mb': size_holder_do_not_save})
    df_do_not_save = df_do_not_save.sort_values(by='size_mb', ascending=False)
    # df_do_not_save = df_do_not_save.copy()  # this does different stuff than df_do_not_save = df_do_not_save a bit strangely
    df_do_not_save = assess_whether_same_object(df_do_not_save.copy())
    st.write(f'Keys that would not actually be saved (total {tot_size_in_memory_do_not_save:.2f} MB):')
    st.dataframe(df_do_not_save)

    # Write a dataframe of the keys and sizes that will be saved, sorted in descending order of size
    df_do_save = pd.DataFrame({'key': keys_holder_do_save, 'size_mb': size_holder_do_save})
    df_do_save = df_do_save.sort_values(by='size_mb', ascending=False)
    # df_do_save = df_do_save.copy()  # this does different stuff than df_do_save = df_do_save a bit strangely
    df_do_save = assess_whether_same_object(df_do_save.copy())
    st.write(f'Keys that would actually be saved (total {tot_size_in_memory_mb_do_save:.2f} MB):')
    st.dataframe(df_do_save)


def main():
    """
    Main function for the page.
    """

    analyze_memory_usage()


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
