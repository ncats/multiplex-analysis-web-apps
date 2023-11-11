'''
This is the python script which produces the PHENOTYPING PAGE
'''
import time
import streamlit as st
from streamlit_javascript import st_javascript

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

    st.header('Data Import and Export\nNCATS-NCI-DMAP')

    ### SIDE BAR ORGANIZATION ###
    with st.sidebar:
        url = st_javascript("await fetch('').then(r => window.parent.location.href)")
        st.write(f'''[Open app in new Tab]({url})\n (MS Edge/ Google Chrome)''')

    NIDAP_Load_U, Local_Upload = st.tabs(['NIDAP Load [Unstuctured]',
                                          'Upload Local Files'])

    ### LOCAL UPLOAD TAB
    with Local_Upload:
        up_file_rdy = False
        uploaded_file = st.file_uploader('Select a .csv file from your computer',
                                         type = 'csv', accept_multiple_files=False,
                                         key = 'file_uploader__do_not_persist')
        if uploaded_file is not None:
            tLoadSt = time.time() # Setup some timing
            # Load the data into memory and perform df checks for req features
            df_upload   = ndl.load_dataset(st.session_state.fiol, uploaded_file,
                                           files_dict=None, file_path=None, loadCompass=False)
            up_file_rdy = ndl.check_upload_df(df_upload,
                                              st.session_state.reqFeatures,
                                              st.session_state.marker_pre)

            if not up_file_rdy: # ERROR!!!
                err_msg_inputs = st.session_state.errmsg_wrongCol
            else:               # SUCCESS!!!
                err_msg_inputs = st.session_state.errmsg_def2row
            load_elapsed = time.time() - tLoadSt
            st.markdown(err_msg_inputs)

        if st.button('Apply Upload Data to Phenotyper', disabled = (not up_file_rdy), 
                     key = 'File_upload_button__do_not_persist'):
            st.session_state = ndl.loadDataButton(st.session_state, df_upload, 'Local Upload', uploaded_file.name[:-4])

            # Save answers to benchmarking dataframe
            st.session_state.bc.set_value_df('file', uploaded_file.name)
            st.session_state.bc.set_value_df('nRows', st.session_state.df.shape[0])
            st.session_state.bc.set_value_df('data_import_loc', 'Local Upload')
            st.session_state.bc.set_value_df('time_load_data', load_elapsed)

    ### NIDAP UNSTRUCTURED LOAD TAB ###
    with NIDAP_Load_U:
        ### SELECT NIDAP DATASET ###
        # Use select box to choose a project
        select_proj_u = st.session_state.usDatasetPaths[0]
        st.markdown('## NIDAP Dataset Import:')
        st.markdown(f'#### {select_proj_u}')

        # Identify the available NIDAP 'files'
        files_dict = st.session_state.files_dict[select_proj_u]

        # Use select box to choose from list of NIDAP 'files'
        st.session_state.datafileU = st.selectbox(
                'Select your DataFrame',
                (list(files_dict.keys())))
        ############################

        ### LOAD NIDAP BUTTON ###
        if st.button('Load NIDAP Data', help='Load the selected data file', key = 'LoadUDataButton__do_not_persist'):
            tLoadSt = time.time() # Setup some timing
            # Load the data into memory and perform df checks for req features
            df_NIDAP    = ndl.load_dataset(st.session_state.fiol, select_proj_u, files_dict, st.session_state.datafileU, loadCompass=True)
            up_file_rdy = ndl.check_upload_df(df_NIDAP, st.session_state.reqFeatures, st.session_state.marker_pre)

            if not up_file_rdy: # ERROR!!!
                err_msg_inputs = st.session_state.errmsg_wrongCol
            else:               # SUCCESS!!!
                err_msg_inputs = st.session_state.errmsg_def2row
                st.session_state = ndl.loadDataButton(st.session_state, df_NIDAP, select_proj_u, st.session_state.datafileU[:-4])
            load_elapsed = time.time() - tLoadSt

            # Save answers to benchmarking dataframe
            st.session_state.bc.set_value_df('file', st.session_state.datafileU)
            st.session_state.bc.set_value_df('nRows', st.session_state.df.shape[0])
            st.session_state.bc.set_value_df('data_import_loc', 'Compass_Unstructured')
            st.session_state.bc.set_value_df('time_load_data', load_elapsed)

            st.markdown(err_msg_inputs)
        ############################

    # Exported Files Path Selection
    selectProjOutCSV_U = st.session_state.OutputCSVPaths_U[0]
    
    st.markdown('## NIDAP Dataset Export:')
    st.markdown(f'#### {selectProjOutCSV_U}')
    
    # Text-input box for phenotype_summary file name
    st.text_input('Phenotype Summary File Name',
                    key = 'pheno_assign_filename_U')
    
    # Text-input box for updated dataset file name
    st.text_input('Updated Dataset File Name',
                    key = 'df_update_filename_U')

    ### Create containers for Button and checkbox
    filtOutU1, filtOutU2, filtOutU3 = st.columns([3,2,1])

    with filtOutU2:
        exportFiltU = st.checkbox('Export Filtered Dataset', key = 'exportFiltU')
    with filtOutU1:
        # Every form must have a submit button.
        submitU = st.button('Export phenotyping files')
        if submitU:
            # Export Assigned Phenotypes output
            ndl.export_results_dataset(st.session_state.fiol, st.session_state.spec_summ, selectProjOutCSV_U, st.session_state.pheno_assign_filename_U, saveCompass=True, type='U')
            # Export Updated dataset output
            ndl.export_results_dataset(st.session_state.fiol, st.session_state.df, selectProjOutCSV_U, st.session_state.df_update_filename_U, saveCompass=True, type='U')
            # Check to export filtered-phenotyped-datasets
            if exportFiltU:
                fileFileName = st.session_state.df_update_filename_U + '_filtered'
                ndl.export_results_dataset(st.session_state.fiol, st.session_state.df_filt, selectProjOutCSV_U, fileFileName, saveCompass=True, type='U')
            # Add a message that this file writing is complete
            st.write('Files have been exported!')

if __name__ == '__main__':
    main()
