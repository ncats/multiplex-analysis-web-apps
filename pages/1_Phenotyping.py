'''
This is the python script which produces the PHENOTYPING PAGE
'''
import os
import time
import streamlit as st
from st_pages import show_pages_from_config, add_indentation
from streamlit_extras.add_vertical_space import add_vertical_space 
from streamlit_extras.app_logo import add_logo
import dataset_formats

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

def setFiltering_features(file_format):
    if file_format == 'REEC':
        st.session_state.SEL_feat = ['tNt']
        st.session_state.CHK_feat = ['GOODNUC']
    elif file_format == 'QuPath':
        st.session_state.SEL_feat = ['Slide_ID']
        st.session_state.CHK_feat = []
    elif file_format == 'OMAL':
        st.session_state.SEL_feat = ['Slide_ID']
        st.session_state.CHK_feat = []
    elif file_format == 'GMBSecondGeneration':
        st.session_state.SEL_feat = ['Slide_ID']
        st.session_state.CHK_feat = []
    elif file_format == 'Native':
        st.session_state.SEL_feat = ['Slide_ID']
        st.session_state.CHK_feat = []
    else:
        st.session_state.SEL_feat = []
        st.session_state.CHK_feat = []

def main():
    '''
    Main function for running the page
    '''

    # Set a wide layout
    st.set_page_config(page_title="Phenotyping",
                       layout="wide")

    # Remove key values from session_state that should not persist
    for key, val in st.session_state.items():
        if (not key.endswith('__do_not_persist')) and (not key.startswith('FormSubmitter:')):
            st.session_state[key] = val

    # Apply pages order and indentation
    add_indentation()
    show_pages_from_config()

    # Sidebar organization
    with st.sidebar:
        st.write('**:book: [Documentation](https://ncats.github.io/multiplex-analysis-web-apps)**')

    # Add logo to page
    add_logo('app_images/mawa_logo-width315.png', height=150)

    if 'init' not in st.session_state:
        settings_yaml_file = 'config_files/OMAL_REEC.yml'
        # Initialize session_state values for streamlit processing
        st.session_state = ndl.init_session_state(st.session_state, settings_yaml_file)

    st.header('Phenotyper\nNCATS-NCI-DMAP')

    input_directory = os.path.join('.', 'input')
    output_directory = os.path.join('.', 'output')
    options_for_input_datafiles = [x for x in os.listdir(input_directory) if x.endswith(('.csv', '.tsv'))]
    phenoFileOptions = [x for x in os.listdir(output_directory) if (x.startswith('phenotype_summary')) and (x.endswith(('.csv', '.tsv')))]

    dataLoadedCols = st.columns([2,2,2])
    with dataLoadedCols[0]:
        st.selectbox(label = 'Choose a datafile', options = options_for_input_datafiles, key = 'datafileU')
        st.number_input('x-y coordinate units (microns):', min_value=0.0, key='phenotyping_micron_coordinate_units', help='E.g., if the coordinates in the input datafile were pixels, this number would be a conversion to microns in units of microns/pixel.', format='%.4f', step=0.0001)

        if (st.button('Load Data')) and (st.session_state.datafileU is not None):
            start = time.time()
            input_datafile = os.path.join('input', st.session_state.datafileU)
            _, _, _, _, file_format, _ = dataset_formats.extract_datafile_metadata(input_datafile)

            if file_format == 'HALO':
                file_format = 'OMAL'

            dataset_class = getattr(dataset_formats, file_format)  # done this way so that the format (e.g., “REEC”) can be select programmatically
            setFiltering_features(file_format)
            dataset_obj = dataset_class(input_datafile, 
                                        coord_units_in_microns = st.session_state.phenotyping_micron_coordinate_units, 
                                        extra_cols_to_keep=['tNt', 'GOODNUC', 'HYPOXIC', 'NORMOXIC', 'NucArea', 'RelOrientation'])
            dataset_obj.process_dataset()
            stop_dataload = time.time()
            elapsed_dataload = round(stop_dataload - start, 2)
            print(f'{input_datafile} tool {elapsed_dataload}s to load into memory')
            
            st.session_state = ndl.loadDataButton(st.session_state, dataset_obj.data, 'Input', st.session_state.datafileU[:-4])
            stop_phenotyping = time.time()
            elapsed_phenotyping = round(stop_phenotyping-stop_dataload, 2)
            print(f'{input_datafile} took {elapsed_phenotyping}s to perform phenotyping')

    with dataLoadedCols[1]:
        st.selectbox(label = 'Choose a previous phenotyping file', options = phenoFileOptions, key = 'phenoFileSelect')
        if (st.button('Load Phenotyping File')) and (st.session_state.phenoFileSelect is not None):
            phenotype_file = os.path.join('output', st.session_state.phenoFileSelect)
            st.session_state.spec_summ_load = bpl.load_previous_species_summary(phenotype_file)
            st.session_state.phenoMeth = 'Custom'
            st.session_state = ndl.updatePhenotyping(st.session_state)
            st.session_state.pointstSliderVal_Sel = st.session_state.calcSliderVal

    with dataLoadedCols[2]:
    ### Data Phenotyping Container ###
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

    #
    if st.session_state.selected_phenoMeth != 'Not Selected':
        st.session_state.phenotyping_completed = True

    ### Data Filters Container ###
    with st.expander('Data Filters'):
        with st.form('Filter Levers'):
            filt_col = st.columns([1, 2])
            with filt_col[0]:
                # Select Box Features
                for feat in st.session_state.SEL_feat:
                    st.selectbox(feat,
                                (st.session_state.df_raw[feat].unique()),
                                key = 'sel' + feat)

            with filt_col[1]:
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

        # Prepare for Exporting
        st.session_state.df_update = st.session_state.df.copy().drop(['mark_bits', 'species_name_long', 'species_name_short'], axis=1)

        phen_summ_cols = st.columns([2, 1])
        with phen_summ_cols[0]:
            st.text_input('Phenotype Summary File Name', key = 'pheno_assign_filename_U')
        with phen_summ_cols[1]:
            add_vertical_space(2)
            if st.button('Append Export List', key = 'appendexportbutton_phenotypesummary__do_not_persist'):
                ndl.save_csv(st.session_state.spec_summ, st.session_state.pheno_assign_filename_U)
                st.toast(f'Added {st.session_state.pheno_assign_filename_U} to export list ')

        updated_df_cols = st.columns([2, 1])
        with updated_df_cols[0]:
            st.text_input('Updated Dataset File Name', key = 'df_update_filename_U')
        with updated_df_cols[1]:
            add_vertical_space(2)
            if st.button('Append Export List', key = 'appendexportbutton_updateddf__do_not_persist'):
                ndl.save_csv(st.session_state.df_update, st.session_state.df_update_filename_U)
                st.toast(f'Added {st.session_state.df_update_filename_U} to export list ')

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

        pheno_scat_upload = st.columns([2, 1])
        with pheno_scat_upload[0]:
            st.text_input('.png file suffix (Optional)', key = 'imgFileSuffixText')
        with pheno_scat_upload[1]:
            add_vertical_space(2)
            if st.button('Append Export List', key = 'appendexportbutton_phenotypescatter__do_not_persist'):
                ndl.save_png(st.session_state.phenoFig, 'Phenotype Scatterplot', st.session_state.imgFileSuffixText)
                st.toast(f'Added {st.session_state.imgFileSuffixText} to export list ')

if __name__ == '__main__':
    main()
