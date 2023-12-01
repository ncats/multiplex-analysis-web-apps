# Import relevant libraries
import os
import streamlit as st
import utils as utils
from streamlit_extras.app_logo import add_logo
from st_pages import show_pages_from_config, add_indentation

import app_top_of_page as top

def main():

    # Create functions to update correlated widgets
    def update_roi_index(roi_names):
        st.session_state['roi_index_to_visualize'] = roi_names.index(st.session_state['roi_name_to_visualize'])
    def update_roi_name(roi_names):
        st.session_state['roi_name_to_visualize'] = roi_names[st.session_state['roi_index_to_visualize']]

    # Set a wide layout
    st.set_page_config(layout="wide")

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

    # Run Top of Page (TOP) functions
    st.session_state = top.check_for_platform(st.session_state)

    # Display page heading
    st.title('Individual ROI heatmaps')

    if os.path.exists(os.path.join('.', 'output', 'images', 'single_roi_outlines_on_whole_slides')) and os.path.exists(os.path.join('.', 'output', 'images', 'roi_plots')) and os.path.exists(os.path.join('.', 'output', 'images', 'dens_pvals_per_roi')):

        # Create an expander to hide some optional widgets
        with st.expander('Optional: Image path extraction', expanded=False):

            # Get the possible analysis radii by analyzing the subdirectories present in the results directory
            # analysis_dir_listing = os.listdir(os.path.join(os.getcwd(), '..', 'results', 'webpage'))
            # analysis_dir_listing = os.listdir(os.path.join('.', 'output', 'images'))
            # st.selectbox('Select analysis radius in microns:', [int(x.lstrip('slices_1x')) for x in analysis_dir_listing], key='analysis_radius_in_microns')
            # st.write('Used to be here: "Select analysis radius in microns:"')

            # Create a button to optionally re-run the image path collection, which is done by default if it hasn't been done already
            collect_image_paths = st.button('Re-run image path extraction')

        # Run the image path collection, assigning the result to the session state
        if ('df_paths_per_roi' not in st.session_state) or collect_image_paths:
            st.session_state['df_paths_per_roi'] = utils.get_paths_for_rois()
            image_paths_extracted = True
        else:
            image_paths_extracted = False

        # Print a message when the image paths have been extracted from disk
        if image_paths_extracted:
            st.info('Image paths extracted from disk')

        # Get shortcut variables to the paths dataframe and a list version of its index
        df_paths_per_roi = st.session_state['df_paths_per_roi']
        roi_names = list(df_paths_per_roi.index)

        # Allow user to select by name the ROI data to display
        st.selectbox('Select ROI name to visualize...', roi_names, key='roi_name_to_visualize', on_change=update_roi_index, args=(roi_names,))

        # Allow user to select by index the ROI data to display
        st.number_input('...OR, select ROI index to visualize:', min_value=0, max_value=(len(roi_names) - 1), key='roi_index_to_visualize', on_change=update_roi_name, args=(roi_names,))

        # Split up the window into two columns
        display_col1, display_col2 = st.columns(2)

        # Display the three images for the currently selected ROI
        with display_col1:
            st.image(df_paths_per_roi.loc[st.session_state['roi_name_to_visualize'], 'roi'])
            image_path_entry = df_paths_per_roi.loc[st.session_state['roi_name_to_visualize'], 'heatmap']
            # if isinstance(image_path_entry, str):
            if image_path_entry != '':
                st.image(image_path_entry)
            else:
                st.write('No heatmap data available')
        with display_col2:
            st.image(df_paths_per_roi.loc[st.session_state['roi_name_to_visualize'], 'outline'])

    else:
        st.warning('At least one of the three sets of per-ROI plots does not exist; please run all per-ROI components of the workflow on the "Run workflow" page', icon='⚠️')

if __name__ == '__main__':
    main()