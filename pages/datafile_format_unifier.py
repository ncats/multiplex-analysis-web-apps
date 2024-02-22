import os
import streamlit as st
import pandas as pd
import app_top_of_page as top
import streamlit_dataframe_editor as sde

def list_files(directory, extensions):
    files = []
    for file in sorted(os.listdir(directory)):
        if file.endswith(extensions):
            files.append(file)
    return files

def main():

    # Set page settings
    st.set_page_config(layout='wide', page_title='File Detector')
    st.title('File Detector')

    # Run streamlit-dataframe-editor library initialization tasks at the top of the page
    st.session_state = sde.initialize_session_state(st.session_state)

    # Run Top of Page (TOP) functions
    st.session_state = top.top_of_page_reqs(st.session_state)

    st.write('Detecting .csv and .tsv files in the \'input\' directory')

    directory = 'input'
    extensions = ('.csv', '.tsv')

    files = list_files(directory, extensions)

    if len(files) == 0:
        st.write('No files found.')
    else:
        st.write('Found files:')
        df_input_files = pd.DataFrame(files, columns=['Filename'])
        df_input_files.insert(0, 'Selected', False)

        if 'unifier__de_datafile_selection' not in st.session_state:
            st.session_state['unifier__de_datafile_selection'] = sde.DataframeEditor(df_name='unifier__df_datafile_selection', default_df_contents=df_input_files)

        st.session_state['unifier__de_datafile_selection'].dataframe_editor(reset_data_editor_button_text='Reset file selections')

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)

if __name__ == '__main__':
    main()
