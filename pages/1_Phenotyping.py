'''
This is the python script which produces the PHENOTYPING PAGE
'''
import time
import streamlit as st
from streamlit_javascript import st_javascript

# Import relevant libraries
import nidap_dashboard_lib as ndl   # Useful functions for dashboards connected to NIDAP
import basic_phenotyper_lib as bpl  # Useful functions for phenotyping collections of cells

def data_editor_change_callback():
    '''
    data_editor_change_callback is a callback function for the streamlit data_editor widget
    which updates the saved value of the user-created changes after every instance of the 
    data_editor on_change method. This ensures the dashboard can remake the edited data_editor
    when the user navigates to a different page.
    '''
    st.session_state['saved_dataeditor_values'] = st.session_state['dataeditor__do_not_persist']

def update_input_data_editor():
    '''
    update_input_data_editor is a function that remakes the input dataframe to the streamlit
    data_editor and still allows for the changes to persist when the user navigates to a 
    different page. Take note that there are global variables used here and may need to be 
    adjusted for your specific use case.
    '''
    for key, value in st.session_state['saved_dataeditor_values']['edited_rows'].items():
        for key2, value2 in value.items():
            st.session_state.spec_summ_dataeditor.loc[key, key2] = value2

def main():
    '''
    Main function for running the page
    '''

    # Use the whole page width
    st.set_page_config(page_title="Phenotyping",
                       layout="wide")

    for key, val in st.session_state.items():
        if (not key.endswith('__do_not_persist')) and (not key.startswith('FormSubmitter:')):
            st.session_state[key] = val

    st.header('Phenotyper\nNCATS-NCI-DMAP')

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

    ### NIDAP STRUCTURED LOAD TAB ###
    # with NIDAP_Load_S:
    #     selectProj_s = st.selectbox(
    #         'Select your NIDAP Project folder',
    #         (st.session_state.projectPaths))
    #     # Use select box to choose from list of NIDAP 'files'
    #     st.session_state.datafileS = st.text_input('Type the name of your Structured Dataset',
    #                                             value=st.session_state.datafileS)
    #     ### LOAD NIDAP BUTTON ###
    #     if st.button('Load NIDAP Data', help='Load the selected data file', key = 'LoadSDataButton'):
    #         tLoadSt = time.time() # Setup some timing
    #         # Load the data into memory and perform df checks for req features
    #         df_NIDAP    = ndl.load_dataset(st.session_state.fiol, selectProj_s, files_dict=None, file_path=st.session_state.datafileS, loadCompass=True)
    #         up_file_rdy = ndl.check_upload_df(df_NIDAP, st.session_state.reqFeatures, st.session_state.marker_pre)

    #         if (not up_file_rdy): # ERROR!!!
    #             err_msg_inputs = st.session_state.errmsg_wrongCol
    #         else:                 # SUCCESS!!!
    #             err_msg_inputs = st.session_state.errmsg_def2row
    #             st.session_state = ndl.loadDataButton(st.session_state, df_NIDAP, selectProj_s, st.session_state.datafileS)
    #         loadElapsed = time.time() - tLoadSt

    #         # Save answers to benchmarking dataframe
    #         st.session_state.bc.set_value_df('file', st.session_state.datafileS)
    #         st.session_state.bc.set_value_df('nRows', st.session_state.df.shape[0])
    #         st.session_state.bc.set_value_df('data_import_loc', 'Compass_Structured')
    #         st.session_state.bc.set_value_df('time_load_data', loadElapsed)

    #         st.markdown(err_msg_inputs)
    #     ############################

    ### NIDAP UNSTRUCTURED LOAD TAB ###
    with NIDAP_Load_U:
        ### SELECT NIDAP DATASET ###
        # Use select box to choose a project
        select_proj_u = st.selectbox('Select your NIDAP Unstructured Dataset',
                                    (st.session_state.usDatasetPaths))

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

    ### Data Phenotyping Container ###
    with st.expander('Phenotying Options'):
        with st.form('Analysis Levers'):

            phenotyping_labels = ('Species', 'Marker', 'Custom')
            st.radio("Choose a Phenotyping Method",
                     phenotyping_labels,
                     horizontal= True,
                     key = 'phenoMeth')

            # Every form must have a submit button.
            submitted = st.form_submit_button('Apply Phenotyping Method')
            if submitted:
                st.session_state = ndl.updatePhenotyping(st.session_state)
                st.session_state.pointstSliderVal_Sel = st.session_state.calcSliderVal

    ### Data Filters Container ###
    with st.expander('Data Filters'):
        with st.form('Filter Levers'):
            filt_col_1, filt_col_2, filt_col_3 = st.columns(3)
            with filt_col_1:
                # Select Box Features
                for feat in st.session_state.SEL_feat:
                    st.selectbox(feat,
                                (st.session_state.df_raw[feat].unique()),
                                key = 'sel' + feat)

            with filt_col_2:
                # Check Box Features
                for feat in st.session_state.CHK_feat:
                    st.checkbox(feat,
                                key = 'sel' + feat)

            submitted = st.form_submit_button('Apply Filters')
            if submitted:
                # Filtered dataset
                st.session_state.df_filt = ndl.perform_filtering(st.session_state)

                # Update and reset Figure Objects
                st.session_state = ndl.setFigureObjs(st.session_state)
                st.session_state.pointstSliderVal_Sel = st.session_state.calcSliderVal

    ### SIDE BAR ORGANIZATION ###
    with st.sidebar:
        url = st_javascript("await fetch('').then(r => window.parent.location.href)")
        st.write(f'''[Open app in new Tab]({url})\n (MS Edge/ Google Chrome)''')

    ## In-App Instructions
    if st.session_state.data_loaded is False:
        st.warning('Data not loaded (above)', icon="⚠️")
    elif st.session_state.selected_phenoMeth == st.session_state.noPhenoOpt:
        st.warning('No phenotyping method applied (above)', icon="⚠️")
    else:
        st.markdown(f'## Current Phenotyping Method: {st.session_state.selected_phenoMeth}')

    ### Define Visualization Columns
    vizCol1, vizCol2 = st.columns(2)

    # Second column on the page
    with vizCol2:

        ### PHENOTYPE SUMMARY TABLE ###
        st.markdown('## Phenotype Summary')
        if st.session_state.selected_phenoMeth == 'Custom':

            if 'saved_dataeditor_values' in st.session_state:
                update_input_data_editor()

            # Allow the user to edit the value counts dataframe
            st.write('Replace the "unassigned" values in the "phenotype_name" column below with a descriptive name, such as "M2 macrophage". Don\'t change the values in any other column!')
            st.session_state.spec_summ = st.data_editor(st.session_state.spec_summ_dataeditor,
                                                        key='dataeditor__do_not_persist',
                                                        use_container_width=True,
                                                        on_change=data_editor_change_callback)

            # Update the Dataset with the Species Summary changes
            st.session_state.df = bpl.update_df_phenotype(st.session_state.df, st.session_state.spec_summ)
            st.session_state.df_filt = ndl.perform_filtering(st.session_state)

            ## Assign Special spec_sum based on current spec_sum
            st.session_state.spec_summ_load = st.session_state.spec_summ.copy()

            # Update the Assigned Phenotypes dataframe with Species Summary changes
            st.session_state.assign_pheno = bpl.init_assign_pheno(st.session_state.df)

            # Set Figure Objects
            st.session_state = ndl.setFigureObjs(st.session_state, st.session_state.pointstSliderVal_Sel)

        else:
            st.dataframe(st.session_state.spec_summ, use_container_width=True)

        ### ASSIGNED PHENOTYPES TABLE ###
        st.write('## Assigned Phenotypes')
        st.write('The following phenotypes will update as the table above is modified. Double-click a cell to see all its contents at once.')
        st.dataframe(st.session_state.assign_pheno, use_container_width=True)

        ### EXPORT PHENOTYPES CONTAINER ###
        st.markdown('## Export Phenotypes')
        st.markdown('''Click the button below to export:  
            1. A .csv file of the Phenotype Summary table  
            2. An updated .csv file of the original dataset including a column with the assigned phenotypes.''')
        
        # Prepare for Exporting
        st.session_state.df_update = st.session_state.df.copy().drop(['mark_bits', 'species_name_long', 'species_name_short'], axis=1)

        # NIDAP_Export_S, NIDAP_Export_U = st.tabs(['NIDAP Export [Stuctured]',
        #                                         'NIDAP Export [Unstuctured]'])
        # with NIDAP_Export_S:
        #     with st.form('Structured Export'):
        #         # Exported Files Path Selection
        #         selectProjOutCSV_S = st.selectbox(
        #                                 'Select your Project',
        #                                 (st.session_state.OutputCSVPaths_S),
        #                                 key = 'ProjOutSel_S')
        #         # Text Boxes for naming your files
        #         st.session_state.pheno_assign_filename_S = st.text_input('Phenotype Summary File Name', st.session_state.pheno_assign_filename_S, key = 'phenoAssignText_S')
        #         st.session_state.df_update_filename_S = st.text_input('Updated Dataset File Name', st.session_state.datafile + '_updated', key = 'updateDataText_S')
        #         ### Create containers for Button and checkbox
        #         filtOutS1, filtOutS2, filtOutS3 = st.columns([3,2,1])

        #         with filtOutS2:
        #             exportFiltS = st.checkbox('Export Filtered Dataset', key = 'exportFiltS')
        #         with filtOutS1:
        #             # Every form must have a submit button.
        #             submitS = st.form_submit_button('Export phenotyping files')
        #             if submitS:
        #                 # Export Assigned Phenotypes output
        #                 ndl.export_results_dataset(st.session_state.fiol, st.session_state.spec_summ, selectProjOutCSV_S, st.session_state.pheno_assign_filename_S, saveCompass=True, type='S')
        #                 # Export Updated dataset output
        #                 ndl.export_results_dataset(st.session_state.fiol, st.session_state.df, selectProjOutCSV_S, st.session_state.df_update_filename_S, saveCompass=True, type='S')
        #                 # Check to export filtered-phenotyped-datasets
        #                 if exportFiltS:
        #                     fileFileName = st.session_state.df_update_filename_S + '_filtered'
        #                     ndl.export_results_dataset(st.session_state.fiol, st.session_state.df_filt, selectProjOutCSV_S, fileFileName, saveCompass=True, type='S')
        #                 # Add a message that this file writing is complete
        #                 st.write('Files have been exported!')

        # with NIDAP_Export_U:
        with st.form('Unstructured Export'):
            # Exported Files Path Selection
            selectProjOutCSV_U = st.selectbox('Select your Unstructured Dataset',
                                              (st.session_state.OutputCSVPaths_U),
                                              key = 'ProjOutSel_U')
            
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
                submitU = st.form_submit_button('Export phenotyping files')
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

    # First column on the page
    with vizCol1:
        # Print a column header
        st.header('Phenotype Plot')

        plotSlide1, plotSlide2 = st.columns(2)
        with plotSlide1:
            with st.form('Plotting Num'):
                st.slider('How many points to plot (%)', 0, 100, key = 'pointstSliderVal_Sel')
                update_pixels_button = st.form_submit_button('Update Scatterplot')
                if update_pixels_button:
                    st.session_state = ndl.setFigureObjs(st.session_state,
                                                         st.session_state.pointstSliderVal_Sel)

        with plotSlide2:
            if st.session_state.calcSliderVal < 100:
                st.warning('Not plotting full scatterplot', icon="⚠️")
            else:
                st.markdown('### Plotting full scatterplot')

            st.write(f'Drawing {st.session_state.drawnPoints} points')

        st.session_state.bc.startTimer()
        st.pyplot(st.session_state.phenoFig)
        st.session_state.bc.printElapsedTime(f'Seaborn plotting of {st.session_state.drawnPoints} points')

        ### Save Figure Button ###
        st.markdown('## Save Figure')

        selectProjOutPNG = st.selectbox(
                            'Select your Unstructured Dataset',
                            (st.session_state.OutputPNGPaths))
        st.text('File name saves as the current DateTime.\nProvide a suffix below to append the end of the .png file name (Optional)')

        # png saving settings columns
        set_suff_col, save_png_col = st.columns([3, 1])

        # set suffix col
        with set_suff_col:
            st.session_state.imgFileSuffixText = st.text_input('.png file suffix (Optional)',
                                                               st.session_state.imgFileSuffixText,
                                                               label_visibility = 'collapsed')
        # png Save Button
        with save_png_col:
            if st.button('Save Figure'):

                # Define the output file name
                art = ''
                if st.session_state.imgFileSuffixText != '':
                    art = '_'
                timestr = time.strftime("%Y%m%d-%H%M%S")
                png_file_name = timestr + art + st.session_state.imgFileSuffixText

                # Save the png to NIDAP
                ndl.save_png_dataset(st.session_state.fiol, selectProjOutPNG,
                                     png_file_name, st.session_state.phenoFig)

if __name__ == '__main__':
    main()
