'''
This is the python script which produces the PHENOTYPING PAGE
'''
import time
import streamlit as st
from streamlit_javascript import st_javascript
from streamlit_extras.add_vertical_space import add_vertical_space 

# Import relevant libraries
import nidap_dashboard_lib as ndl   # Useful functions for dashboards connected to NIDAP
import basic_phenotyper_lib as bpl  # Useful functions for phenotyping collections of cells

def main():
    '''
    Main function for running the page
    '''

    # Use the whole page width
    st.set_page_config(page_title="Data Import and Export",
                       layout="wide")

    for key, val in st.session_state.items():
        if (not key.endswith('__do_not_persist')) and (not key.startswith('FormSubmitter:')):
            st.session_state[key] = val

    data_header = st.columns([2, 3])
    with data_header[0]:
        st.header('Data Import and Export\nNCATS-NCI-DMAP')
    with data_header[1]:
        add_vertical_space(6)
        st.markdown('### Project Name: OMAL/reec-UMAP')

    ### SIDE BAR ORGANIZATION ###
    with st.sidebar:
        url = st_javascript("await fetch('').then(r => window.parent.location.href)")
        st.write(f'''[Open app in new Tab]({url})\n (MS Edge/ Google Chrome)''')

    ### SELECT NIDAP DATASET ###
    # Use select box to choose a project
    select_proj_u = st.session_state.usDatasetPaths[0]
    st.markdown('## NIDAP Input Dataset:')

    # Identify the available NIDAP 'files'
    files_dict = st.session_state.files_dict[select_proj_u]

    # Use select box to choose from list of NIDAP 'files'
    st.session_state.datafileU = st.selectbox(
            'Select your DataFrame',
            (list(files_dict.keys())))
    ############################

    errorStatus = None
    ### LOAD NIDAP BUTTON ###
    if st.button('Load NIDAP Data', help='Load the selected data file', key = 'LoadUDataButton__do_not_persist'):
        tLoadSt = time.time() # Setup some timing
        # Load the data into memory and perform df checks for req features
        df_NIDAP    = ndl.load_dataset(st.session_state.fiol, select_proj_u, files_dict, st.session_state.datafileU, loadCompass=True)
        up_file_rdy = ndl.check_upload_df(df_NIDAP, st.session_state.reqFeatures, st.session_state.marker_pre)

        if not up_file_rdy: # ERROR!!!
            errorStatus = True
            err_msg_inputs = st.session_state.errmsg_wrongCol
        else:               # SUCCESS!!!
            errorStatus = False
            err_msg_inputs = st.session_state.errmsg_def2row
            st.session_state = ndl.loadDataButton(st.session_state, df_NIDAP, select_proj_u, st.session_state.datafileU[:-4])
        load_elapsed = time.time() - tLoadSt

        # Save answers to benchmarking dataframe
        st.session_state.bc.set_value_df('file', st.session_state.datafileU)
        st.session_state.bc.set_value_df('nRows', st.session_state.df.shape[0])
        st.session_state.bc.set_value_df('data_import_loc', 'Compass_Unstructured')
        st.session_state.bc.set_value_df('time_load_data', load_elapsed)
    ############################

    if errorStatus is None:
        add_vertical_space(6)
    elif errorStatus is True:
        st.error(err_msg_inputs)
    elif errorStatus is False:
        st.success('Datafile succesfully imported')

    # Exported Files Path Selection
    selectProjOutCSV_U = st.session_state.OutputCSVPaths_U[0]
    
    st.markdown('## NIDAP Output Dataset:')
    ### Create containers for Button and checkbox
    export_cols = st.columns([1,1,2,1])
    
    with export_cols[0]:
        st.markdown('### File Export List')
    with export_cols[1]:
        exportButton = st.button('Export phenotyping files')
    with export_cols[2]:
        if exportButton:
            numFiles = st.session_state.files_to_export.shape[0]
            my_bar = st.progress(0)
            for index, row in st.session_state.files_to_export.iterrows():
                my_bar.progress((0 + index/numFiles), f'Uploading {row["File Name"]}')
                ndl.export_results_dataset(st.session_state.fiol, 
                                           st.session_state.spec_summ, # Change this
                                           selectProjOutCSV_U, 
                                           row['File Name'], 
                                           saveCompass=True, 
                                           type='U')

            my_bar.progress(100, 'Files have been exported!')

    st.dataframe(data=st.session_state.files_to_export, hide_index=False)

if __name__ == '__main__':
    main()
