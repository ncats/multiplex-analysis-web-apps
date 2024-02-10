# Import relevant libraries
import streamlit as st
import utils as utils
import os

import app_top_of_page as top
import streamlit_dataframe_editor as sde

def main():

    # Create functions to update correlated widgets
    def update_slide_index(slide_names):
        st.session_state['slide_index_to_visualize'] = slide_names.index(st.session_state['slide_name_to_visualize'])
    def update_slide_name(slide_names):
        st.session_state['slide_name_to_visualize'] = slide_names[st.session_state['slide_index_to_visualize']]
    def update_center_species_index(center_species_names):
        st.session_state['center_species_index_to_visualize'] = center_species_names.index(st.session_state['center_species_name_to_visualize'])
    def update_center_species_name(center_species_names):
        st.session_state['center_species_name_to_visualize'] = center_species_names[st.session_state['center_species_index_to_visualize']]
    def update_neighbor_species_index(neighbor_species_names):
        st.session_state['neighbor_species_index_to_visualize'] = neighbor_species_names.index(st.session_state['neighbor_species_name_to_visualize'])
    def update_neighbor_species_name(neighbor_species_names):
        st.session_state['neighbor_species_name_to_visualize'] = neighbor_species_names[st.session_state['neighbor_species_index_to_visualize']]

    # Set a wide layout
    st.set_page_config(layout="wide")

    # Run streamlit-dataframe-editor library initialization tasks at the top of the page
    st.session_state = sde.initialize_session_state(st.session_state)

    # Run Top of Page (TOP) functions
    st.session_state = top.top_of_page_reqs(st.session_state)

    # Display page heading
    st.title('ROI P values overlaid on slides')

    if os.path.exists(os.path.join('.', 'output', 'images', 'density_pvals_over_slide_spatial_plot')):

        # Create an expander to hide some optional widgets
        with st.expander('Optional: Overlay information extraction', expanded=False):

            # # Get the possible analysis radii by analyzing the subdirectories present in the results directory
            # analysis_dir_listing = os.listdir(os.path.join(os.getcwd(), '..', 'results', 'webpage'))
            # st.selectbox('Select analysis radius in microns:', [int(x.lstrip('slices_1x')) for x in analysis_dir_listing], key='analysis_radius_in_microns')

            # Create a button to optionally re-run the overlay information extraction, which is done by default if it hasn't been done already
            collect_overlay_info = st.button('Re-run overlay information extraction')

        # Run the overlay information extraction, assigning the result to the session state
        if ('overlay_info' not in st.session_state) or collect_overlay_info:
            st.session_state['overlay_info'] = utils.get_overlay_info()
            overlay_info_extracted = True
            st.session_state['slide_name_to_visualize'] = st.session_state['overlay_info']['slide_names'][0]
            update_slide_index(st.session_state['overlay_info']['slide_names'])
            st.session_state['center_species_name_to_visualize'] = st.session_state['overlay_info']['center_species'][0]
            update_center_species_index(st.session_state['overlay_info']['center_species'])
            st.session_state['neighbor_species_name_to_visualize'] = st.session_state['overlay_info']['neighbor_species'][0]
            update_neighbor_species_index(st.session_state['overlay_info']['neighbor_species'])
        else:
            overlay_info_extracted = False

        # Print a message when the overlay information has been extracted from disk
        if overlay_info_extracted:
            st.info('Overlay information extracted from disk')

        # Get shortcut variables to the information contained in the overlay dictionary
        overlay_dir = st.session_state['overlay_info']['overlay_dir']
        slide_names = st.session_state['overlay_info']['slide_names']
        center_species = st.session_state['overlay_info']['center_species']
        neighbor_species = st.session_state['overlay_info']['neighbor_species']

        # Create three columns for the filename component selection
        filename_component_col1, filename_component_col2, filename_component_col3 = st.columns(3)

        # Allow user to select by name or by index the slide data to display
        with filename_component_col1:
            st.selectbox('Select slide name to visualize...', slide_names, key='slide_name_to_visualize', on_change=update_slide_index, args=(slide_names,))
            st.number_input('...OR, select slide index to visualize:', min_value=0, max_value=(len(slide_names) - 1), key='slide_index_to_visualize', on_change=update_slide_name, args=(slide_names,))

        # Allow user to select by name or by index the center species data to display
        with filename_component_col2:
            st.selectbox('Select center species name to visualize...', center_species, key='center_species_name_to_visualize', on_change=update_center_species_index, args=(center_species,))
            st.number_input('...OR, select center species index to visualize:', min_value=0, max_value=(len(center_species) - 1), key='center_species_index_to_visualize', on_change=update_center_species_name, args=(center_species,))

        # Allow user to select by name or by index the neighbor species data to display
        with filename_component_col3:
            st.selectbox('Select neighbor species name to visualize...', neighbor_species, key='neighbor_species_name_to_visualize', on_change=update_neighbor_species_index, args=(neighbor_species,))
            st.number_input('...OR, select neighbor species index to visualize:', min_value=0, max_value=(len(neighbor_species) - 1), key='neighbor_species_index_to_visualize', on_change=update_neighbor_species_name, args=(neighbor_species,))

        # Split up the window into two columns
        display_col1, display_col2 = st.columns(2)

        # Display the left and right P values for the current analysis selections
        with display_col1:
            st.subheader('Left P value')
            st.image(os.path.join(overlay_dir, '{}-with_log_dens_pvals_per_roi__center_{}__neighbor_{}__left_pvals.png'.format(st.session_state['slide_name_to_visualize'], st.session_state['center_species_name_to_visualize'], st.session_state['neighbor_species_name_to_visualize'])))
        with display_col2:
            st.subheader('Right P value')
            st.image(os.path.join(overlay_dir, '{}-with_log_dens_pvals_per_roi__center_{}__neighbor_{}__right_pvals.png'.format(st.session_state['slide_name_to_visualize'], st.session_state['center_species_name_to_visualize'], st.session_state['neighbor_species_name_to_visualize'])))

    else:
        st.warning('The component "Plot density P values for each ROI over slide spatial plot" of the workflow does not appear to have been run; please select it on the "Run workflow" page', icon='⚠️')

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)

if __name__ == '__main__':
    main()
