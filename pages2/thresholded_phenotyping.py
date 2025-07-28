'''
This is the python script which produces the PHENOTYPING PAGE
'''
import os
import pandas as pd
import streamlit as st
from streamlit_extras.add_vertical_space import add_vertical_space
import nidap_dashboard_lib as ndl   # Useful functions for dashboards connected to NIDAP
import basic_phenotyper_lib as bpl  # Useful functions for phenotyping collections of cells
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

def slide_id_prog_left_callback():
    '''
    callback function when the left Cell_ID progression button is clicked
    '''
    if st.session_state['idxSlide ID'] > 0:
        st.session_state['idxSlide ID'] -=1
        st.session_state['selSlide ID'] = st.session_state['uniSlide ID'][st.session_state['idxSlide ID']]
        st.session_state['selSlide ID_short'] = st.session_state['uniSlide ID_short'][st.session_state['idxSlide ID']]

def slide_id_prog_right_callback():
    '''
    callback function when the right Cell_ID progression button is clicked
    '''
    if st.session_state['idxSlide ID'] < st.session_state['numSlide ID']-1:
        st.session_state['idxSlide ID'] +=1
        st.session_state['selSlide ID'] = st.session_state['uniSlide ID'][st.session_state['idxSlide ID']]
        st.session_state['selSlide ID_short'] = st.session_state['uniSlide ID_short'][st.session_state['idxSlide ID']]

def slide_id_callback():
    '''
    callback function when the Cell_ID select box changes
    '''
    st.session_state['idxSlide ID'] = st.session_state['uniSlide ID_short'].index(st.session_state['selSlide ID_short'])
    st.session_state['selSlide ID'] = st.session_state['uniSlide ID'][st.session_state['idxSlide ID']]

def filter_and_plot(plot_by_slider = False):
    '''
    function to update the filtering and the figure plotting
    '''
    st.session_state.prog_left_disabled = False
    st.session_state.prog_right_disabled = False

    if st.session_state['idxSlide ID'] == 0:
        st.session_state.prog_left_disabled = True

    if st.session_state['idxSlide ID'] == st.session_state['numSlide ID']-1:
        st.session_state.prog_right_disabled = True

    if plot_by_slider:
        slider_val = st.session_state.point_slider_val
    else:
        slider_val = None

    # Filtered dataset
    df_filt = ndl.perform_filtering(st.session_state)

    st.session_state.pheno_summ_filt = bpl.init_pheno_summ(df_filt)

    # Redefine the Phenotyping Order
    st.session_state.phenoOrder = list(st.session_state.pheno_summ.loc[st.session_state.pheno_summ['phenotype_count'].index, 'phenotype'])

    # Update and reset Figure Objects
    st.session_state = ndl.set_figure_objs(session_state = st.session_state,
                                           df_plot = df_filt,
                                           slider_val = slider_val)

def marker_multiselect_callback():
    '''
    Callback function for selecting which markers to include in phenotyping
    '''
    st.session_state.marker_names = st.session_state.marker_multi_sel
    st.session_state = ndl.set_phenotyping_elements(st.session_state, st.session_state.df_raw)

def make_sample_phenotype_assignments(csv_filename):
    '''
    Allow sample phenotype assignments to be made for quick testing
    '''
    # NEXT:
    #   * Make button below like in the Gater
    #   * make_sample_phenotype_assignments(csv_filename='sample_phenotype_assignments.csv')

    # Get the current version of the phenotype assignments dataframe
    df_to_which_to_update = st.session_state['pheno__de_phenotype_assignments'].reconstruct_edited_dataframe()

    # Load in a dataframe containing the translations to make, from species_name_short to phenotype
    df_to_assign = pd.read_csv(os.path.join('.', 'sample_phenotyping', csv_filename))

    # Create a dictionary of these translations
    assignments = dict(zip(df_to_assign['species_name_short'], df_to_assign['phenotype']))

    # Use the translations dictionary to perform actual translation
    df_to_which_to_update['phenotype'] = df_to_which_to_update['species_name_short'].apply(lambda species_name_short: assignments[species_name_short])

    # Update the official phenotype assignments dataframe with the dataframe to which to update it
    st.session_state['pheno__de_phenotype_assignments'].update_editor_contents(df_to_which_to_update, reset_key=False, additional_callback=data_editor_change_callback)

def main():
    '''
    Main function for running the page
    '''

    # If 'input_dataset' isn't in the session state, print an error message and return
    if 'input_dataset' not in st.session_state:
        st.error('An input dataset has not yet been opened. Please do so using the "Open File" page in the sidebar.')
        return

    output_directory = os.path.join('.', 'output')
    phenoFileOptions = [x for x in os.listdir(output_directory) if (x.startswith('phenotype_summary')) and (x.endswith(('.csv', '.tsv')))]

    data_load_cols = st.columns([2,2,2])
    with data_load_cols[0]:
        if st.button('Load Data'):
            st.session_state.bc.startTimer()
            dataset_obj = st.session_state['input_dataset']
            st.session_state.bc.printElapsedTime(msg = f'Loading {st.session_state["input_metadata"]["datafile_path"]} into memory, Datafile Unifier')
            st.session_state.bc.startTimer()
            if st.session_state['input_metadata']['datafile_path'] is not None:
                datafile_name = os.path.splitext(os.path.basename(st.session_state['input_metadata']['datafile_path']))[0]
            else:
                datafile_name = 'from_memory'
            st.session_state = ndl.loadDataButton(st.session_state, dataset_obj.data, 'Input', datafile_name)
            st.session_state.bc.printElapsedTime(msg = f'Performing Phenotyping on {st.session_state["input_metadata"]["datafile_path"]}')
        st.session_state.bc.set_value_df('time_load_data', st.session_state.bc.elapsedTime())

    with data_load_cols[1]:
        st.selectbox(label = 'Choose a previous phenotyping file', options = phenoFileOptions, key = 'phenoFileSelect', help='Loaded .csv files populate here when the file name begins with "phenotype_summary"')
        if (st.button('Load Phenotyping File')) and (st.session_state.phenoFileSelect is not None):
            phenotype_file = os.path.join('output', st.session_state.phenoFileSelect)
            st.session_state.spec_summ_load = bpl.load_previous_species_summary(phenotype_file)
            st.session_state.phenoMeth = 'Custom'
            st.session_state.selected_phenoMeth = 'Custom'
            st.session_state = ndl.updatePhenotyping(st.session_state)
            st.session_state.pointstSliderVal_Sel = st.session_state.calcSliderVal

    with data_load_cols[2]:
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
                if 'pheno__de_phenotype_assignments' in st.session_state:
                    del st.session_state['pheno__de_phenotype_assignments']

    if st.session_state.selected_phenoMeth != 'Not Selected':
        st.session_state.phenotyping_completed = True

    mid_col = st.columns(2)
    with mid_col[1]:
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
    viz_col = st.columns(2)

    # Second column on the page
    with viz_col[1]:

        ### PHENOTYPE ASSIGNMENTS TABLE ###
        st.markdown('## Phenotype Assignments')
        if st.session_state.selected_phenoMeth == 'Custom':

            # Allow the user to edit the value counts dataframe
            st.write('Replace the "unassigned" values in the "phenotype_name" column below with a descriptive name, such as "M2 macrophage". Don\'t change the values in any other column!')
            if 'pheno__de_phenotype_assignments' not in st.session_state:
                st.session_state['pheno__de_phenotype_assignments'] = sde.DataframeEditor(df_name='pheno__df_phenotype_assignments', default_df_contents=bpl.init_pheno_assign(st.session_state.df))
            st.session_state['pheno__de_phenotype_assignments'].dataframe_editor(on_change=data_editor_change_callback, reset_data_editor_button_text='Reset phenotype assignments')  # note there is no return variable

            # # Allow a sample gating table to be loaded
            # st.button('Load sample gating table', on_click=make_sample_phenotype_assignments, kwargs={'csv_filename': 'sample_phenotype_assignments.csv'})

        else:
            st.dataframe(st.session_state.spec_summ, use_container_width=True)

        ### PHENOTYPE SUMMARY TABLE ###
        st.write('## Phenotype Summary')
        pheno_summ_tabs = st.tabs(['Summary Table', 'Summary Figure'])
        with pheno_summ_tabs[0]:
            st.write('The following phenotypes will update as the table above is modified. Double-click a cell to see all its contents at once.')
            st.dataframe(st.session_state.pheno_summ, use_container_width=True)

            # Prepare for Exporting
            df_update = st.session_state.df.drop(['mark_bits', 'species_name_long', 'species_name_short'], axis=1)

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
                    ndl.save_csv(df_update, st.session_state.df_update_filename_U)
                    st.toast(f'Added {st.session_state.df_update_filename_U} to export list ')

        with pheno_summ_tabs[1]:
            st.radio('Bar plot view', options=['Full data', 'Slide data'], key='pheno_bar_plot_view', horizontal=True)
            if st.session_state.pheno_bar_plot_view == 'Full data':
                st.markdown('### Full Dataset')
                pheno_df = st.session_state.pheno_summ
            else:
                st.markdown(f'### {st.session_state["selSlide ID_short"]} Dataset')
                pheno_df = st.session_state.pheno_summ_filt
            st.session_state.pheno_summ_bar_fig = bpl.draw_pheno_summ_bar_fig(pheno_df, st.session_state.selhas_pos_mark)
            st.plotly_chart(st.session_state.pheno_summ_bar_fig)

    # First column on the page
    with viz_col[0]:
        # Print a column header
        st.header('Phenotype Plot')

        if 'pheno_flip_yaxis' not in st.session_state:
            st.session_state.pheno_flip_yaxis = False

        filter_and_plot(plot_by_slider=True)

        plot_slide = st.columns(2)
        with plot_slide[0]:
            st.slider('How many points to plot (%)', 0, 100, key = 'point_slider_val')

        with plot_slide[1]:
            if st.session_state.point_slider_val < 100:
                st.warning('Not plotting full scatterplot', icon="⚠️")
            else:
                st.markdown('### Plotting full scatterplot')

            st.write(f'Drawing {st.session_state.drawnPoints} points')
            st.toggle('Omit cells with all negative markers',
                        key = 'selhas_pos_mark')
            st.toggle('Flip Y-axis', key='pheno_flip_yaxis')

        image_prog_col = st.columns([3, 1, 1, 2])
        with image_prog_col[0]:
            st.selectbox('Slide ID',
                         (st.session_state['uniSlide ID_short']),
                         key = 'selSlide ID_short',
                         on_change=slide_id_callback)
        with image_prog_col[1]:
            add_vertical_space(2)
            st.button('←', on_click=slide_id_prog_left_callback, disabled=st.session_state.prog_left_disabled)
        with image_prog_col[2]:
            add_vertical_space(2)
            st.button('→', on_click=slide_id_prog_right_callback, disabled=st.session_state.prog_right_disabled)
        with image_prog_col[3]:
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

if __name__ == '__main__':
    main()
