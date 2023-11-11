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

        # Prepare for Exporting
        st.session_state.df_update = st.session_state.df.copy().drop(['mark_bits', 'species_name_long', 'species_name_short'], axis=1)

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

if __name__ == '__main__':
    main()
