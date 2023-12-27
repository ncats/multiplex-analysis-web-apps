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
import app_top_of_page as top
import streamlit_dataframe_editor as sde

def data_editor_change_callback():
    '''
    data_editor_change_callback is a callback function for the streamlit data_editor widget
    which updates the saved value of the user-created changes after every instance of the 
    data_editor on_change method. This ensures the dashboard can remake the edited data_editor
    when the user navigates to a different page.
    '''

    st.session_state.df = bpl.assign_phenotype_custom(st.session_state.df, st.session_state['pheno__de_phenotype_assignments'].reconstruct_edited_dataframe())

    # Create Phenotypes Summary Table based on 'phenotype' column in df
    st.session_state.pheno_summ = bpl.init_pheno_summ(st.session_state.df)

    # Perform filtering
    st.session_state.df_filt = ndl.perform_filtering(st.session_state)

    # Set Figure Objects based on updated df
    st.session_state = ndl.setFigureObjs(st.session_state, st.session_state.pointstSliderVal_Sel)

def slide_id_prog_left_callback():
    '''
    callback function when the left Cell_ID progression button is clicked
    '''
    if st.session_state['idxSlide ID'] > 0:
        st.session_state['idxSlide ID'] -=1
        st.session_state['selSlide ID'] = st.session_state['uniSlide ID'][st.session_state['idxSlide ID']]
        st.session_state['selSlide ID_short'] = st.session_state['uniSlide ID_short'][st.session_state['idxSlide ID']]
        filter_and_plot()

def slide_id_prog_right_callback():
    '''
    callback function when the right Cell_ID progression button is clicked
    '''
    if st.session_state['idxSlide ID'] < st.session_state['numSlide ID']-1:
        st.session_state['idxSlide ID'] +=1
        st.session_state['selSlide ID'] = st.session_state['uniSlide ID'][st.session_state['idxSlide ID']]
        st.session_state['selSlide ID_short'] = st.session_state['uniSlide ID_short'][st.session_state['idxSlide ID']]
        filter_and_plot()

def slide_id_callback():
    '''
    callback function when the Cell_ID select box changes
    '''
    st.session_state['idxSlide ID'] = st.session_state['uniSlide ID_short'].index(st.session_state['selSlide ID_short'])
    st.session_state['selSlide ID'] = st.session_state['uniSlide ID'][st.session_state['idxSlide ID']]
    filter_and_plot()

def filter_and_plot():
    '''
    function to update the filtering and the figure plotting
    '''
    st.session_state.prog_left_disabeled  = False
    st.session_state.prog_right_disabeled = False

    if st.session_state['idxSlide ID'] == 0:
        st.session_state.prog_left_disabeled = True

    if st.session_state['idxSlide ID'] == st.session_state['numSlide ID']-1:
        st.session_state.prog_right_disabeled = True

    # Filtered dataset
    st.session_state.df_filt = ndl.perform_filtering(st.session_state)

    # Update and reset Figure Objects
    st.session_state = ndl.setFigureObjs(st.session_state)

def marker_multiselect_callback():
    st.session_state.marker_names = st.session_state.marker_multi_sel
    st.session_state = ndl.set_phenotyping_elements(st.session_state, st.session_state.df_raw)

def main():
    '''
    Main function for running the page
    '''

    # Set a wide layout
    st.set_page_config(page_title="Phenotyping",
                       layout="wide")

    # Run streamlit-dataframe-editor library initialization tasks at the top of the page
    st.session_state = sde.initialize_session_state(st.session_state)

    # Apply pages order and indentation
    add_indentation()
    show_pages_from_config()

    if 'init' not in st.session_state:
        settings_yaml_file = 'config_files/OMAL_REEC.yml'
        # Initialize session_state values for streamlit processing
        st.session_state = ndl.init_session_state(st.session_state, settings_yaml_file)

    # Sidebar organization
    with st.sidebar:
        st.write('**:book: [Documentation](https://ncats.github.io/multiplex-analysis-web-apps)**')

        # st.button('Record Benchmarking', on_click=st.session_state.bc.save_run_to_csv)

    # Add logo to page
    add_logo('app_images/mawa_logo-width315.png', height=150)

    # Run Top of Page (TOP) functions
    st.session_state = top.check_for_platform(st.session_state)

    st.header('Phenotyper\nNCATS-NCI-DMAP')

    input_directory = os.path.join('.', 'input')
    output_directory = os.path.join('.', 'output')
    options_for_input_datafiles = [x for x in os.listdir(input_directory) if x.endswith(('.csv', '.tsv'))]
    phenoFileOptions = [x for x in os.listdir(output_directory) if (x.startswith('phenotype_summary')) and (x.endswith(('.csv', '.tsv')))]

    dataLoadedCols = st.columns([2,2,2])
    with dataLoadedCols[0]:
        st.selectbox(label = 'Choose a datafile', options = options_for_input_datafiles, key = 'datafileU')
        st.number_input('x-y coordinate units (microns):', min_value=0.0, key='phenotyping_micron_coordinate_units', help='E.g., if the coordinates in the input datafile were pixels, this number would be a conversion to microns in units of microns/pixel.', format='%.4f', step=0.0001)

        databuttonCols = st.columns([1, 2])
        with databuttonCols[0]:
            if (st.button('Load Data')) and (st.session_state.datafileU is not None):
                st.session_state.bc.startTimer()
                input_datafile = os.path.join('input', st.session_state.datafileU)
                _, _, _, _, file_format, _ = dataset_formats.extract_datafile_metadata(input_datafile)

                if file_format == 'HALO':
                    file_format = 'OMAL'

                st.session_state.file_format = file_format
                dataset_class = getattr(dataset_formats, file_format)  # done this way so that the format (e.g., “REEC”) can be select programmatically
                dataset_obj = dataset_class(input_datafile, 
                                            coord_units_in_microns = st.session_state.phenotyping_micron_coordinate_units, 
                                            extra_cols_to_keep=['tNt', 'GOODNUC', 'HYPOXIC', 'NORMOXIC', 'NucArea', 'RelOrientation'])
                dataset_obj.process_dataset(do_calculate_minimum_coordinate_spacing_per_roi=False)
                st.session_state.bc.printElapsedTime(msg = f'Loading {input_datafile} into memory')
                
                st.session_state.bc.startTimer()
                st.session_state = ndl.loadDataButton(st.session_state, dataset_obj.data, 'Input', st.session_state.datafileU[:-4])
                st.session_state.bc.printElapsedTime(msg = f'Performing Phenotyping on {input_datafile}')

        with databuttonCols[1]:
            if (st.button('Load Multi-axial Gating Data')) & ('mg__df' in st.session_state):
                st.session_state.bc.startTimer()
                st.session_state = ndl.loadDataButton(st.session_state, st.session_state['mg__df'], 'Mutli-axial Gating', st.session_state.mg__input_datafile_filename[:-4])
                st.session_state.bc.printElapsedTime(msg = f'Performing Phenotyping')
        st.session_state.bc.set_value_df('time_load_data', st.session_state.bc.elapsedTime())

    with dataLoadedCols[1]:
        st.selectbox(label = 'Choose a previous phenotyping file', options = phenoFileOptions, key = 'phenoFileSelect', help='Loaded .csv files populate here when the file name begins with "phenotype_summary"')
        if (st.button('Load Phenotyping File')) and (st.session_state.phenoFileSelect is not None):
            phenotype_file = os.path.join('output', st.session_state.phenoFileSelect)
            st.session_state.spec_summ_load = bpl.load_previous_species_summary(phenotype_file)
            st.session_state.phenoMeth = 'Custom'
            st.session_state.selected_phenoMeth = 'Custom'
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
                st.session_state.selected_phenoMeth = st.session_state.phenoMeth
                st.session_state = ndl.updatePhenotyping(st.session_state)
                st.session_state.pointstSliderVal_Sel = st.session_state.calcSliderVal

    #
    if st.session_state.selected_phenoMeth != 'Not Selected':
        st.session_state.phenotyping_completed = True

    midCol = st.columns(2)
    with midCol[0]:
    ### Data Filters Container ###
        with st.expander('Data Filters'):
            with st.form('Filter Levers'):
                filt_col = st.columns([1, 2])
                with filt_col[0]:
                    # Select Box Features
                    for feat in st.session_state.SEL_feat_widg:
                        st.selectbox(feat,
                                    (st.session_state.df_raw[feat].unique()),
                                    key = 'sel' + feat)

                with filt_col[1]:
                    # Check Box Features
                    for feat in st.session_state.CHK_feat_widg:
                        st.checkbox(feat,
                                    key = 'sel' + feat)

                submitted = st.form_submit_button('Apply Filters')
                if submitted:
                    filter_and_plot()
                    st.session_state.pointstSliderVal_Sel = st.session_state.calcSliderVal
    with midCol[1]:
         with st.expander('Choose Markers to include'):
            st.multiselect('Markers', options = st.session_state.loaded_marker_names,
                                      key = 'marker_multi_sel',
                                      on_change=marker_multiselect_callback)

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

        ### PHENOTYPE ASSIGNMENTS TABLE ###
        st.markdown('## Phenotype Assignments')
        if st.session_state.selected_phenoMeth == 'Custom':

            # Allow the user to edit the value counts dataframe
            st.write('Replace the "unassigned" values in the "phenotype_name" column below with a descriptive name, such as "M2 macrophage". Don\'t change the values in any other column!')
            if 'pheno__de_phenotype_assignments' not in st.session_state:
                st.session_state['pheno__de_phenotype_assignments'] = sde.DataframeEditor(df_name='pheno__df_phenotype_assignments', default_df_contents=bpl.init_pheno_assign(st.session_state.df))
            st.session_state['pheno__de_phenotype_assignments'].dataframe_editor(on_change=data_editor_change_callback, reset_data_editor_button_text='Reset phenotype assignments')  # note there is no return variable
        else:
            st.dataframe(st.session_state.spec_summ, use_container_width=True)

        ### PHENOTYPE SUMMARY TABLE ###
        st.write('## Phenotype Summary')
        st.write('The following phenotypes will update as the table above is modified. Double-click a cell to see all its contents at once.')
        st.dataframe(st.session_state.pheno_summ, use_container_width=True)

        # Prepare for Exporting
        st.session_state.df_update = st.session_state.df.copy().drop(['mark_bits', 'species_name_long', 'species_name_short'], axis=1)

        phen_summ_cols = st.columns([2, 1])
        with phen_summ_cols[0]:
            st.text_input('Phenotype Summary File Name', key = 'pheno_assign_filename_U')
        with phen_summ_cols[1]:
            add_vertical_space(2)
            if st.button('Append Export List', key = 'appendexportbutton_phenotypesummary__do_not_persist'):
                if st.session_state.selected_phenoMeth == 'Custom':
                    ndl.save_csv(st.session_state['pheno__de_phenotype_assignments'].reconstruct_edited_dataframe(), st.session_state.pheno_assign_filename_U)  # use dataframe editor
                else:
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
            st.checkbox('Omit drawing cells with all negative markers',
                        key = 'selhas_pos_mark',
                        on_change=filter_and_plot)

        imageProgCol = st.columns([3, 1, 1, 2])
        with imageProgCol[0]:
            st.selectbox('Slide ID',
                         (st.session_state['uniSlide ID_short']),
                         key = 'selSlide ID_short',
                         on_change=slide_id_callback)
        with imageProgCol[1]:
            add_vertical_space(2)
            st.button('←', on_click=slide_id_prog_left_callback, disabled=st.session_state.prog_left_disabeled)
        with imageProgCol[2]:
            add_vertical_space(2)
            st.button('→', on_click=slide_id_prog_right_callback, disabled=st.session_state.prog_right_disabeled)
        with imageProgCol[3]:
            add_vertical_space(2)
            st.write(f'Image {st.session_state["idxSlide ID"]+1} of {st.session_state["numSlide ID"]}')

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

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)

if __name__ == '__main__':
    main()
