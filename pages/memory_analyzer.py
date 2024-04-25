# Import relevant libraries
import streamlit as st
import app_top_of_page as top
import streamlit_dataframe_editor as sde
import pandas as pd
import streamlit_dataframe_editor
from pympler import asizeof


def analyze_memory_usage(saved_streamlit_session_state_key='session_selection'):

    # This function is largely copied from streamlit_session_state_management.save_session_state()

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
    st.write(f'Keys that would not actually be saved (total {tot_size_in_memory_do_not_save:.2f} MB):')
    st.dataframe(df_do_not_save)

    # Write a dataframe of the keys and sizes that will be saved, sorted in descending order of size
    df_do_save = pd.DataFrame({'key': keys_holder_do_save, 'size_mb': size_holder_do_save})
    df_do_save = df_do_save.sort_values(by='size_mb', ascending=False)
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
