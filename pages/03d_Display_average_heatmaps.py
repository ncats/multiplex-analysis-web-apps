# Import relevant libraries
import streamlit as st
from st_pages import show_pages_from_config, add_indentation
import utils as utils
from streamlit_extras.app_logo import add_logo
import os

def main():

    # Create functions to update correlated widgets
    def update_slide_index(slide_names):
        st.session_state['slide_index_to_visualize'] = slide_names.index(st.session_state['slide_name_to_visualize'])
    def update_slide_name(slide_names):
        st.session_state['slide_name_to_visualize'] = slide_names[st.session_state['slide_index_to_visualize']]

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

    # Display page heading
    st.title('Average heatmaps per slide')

    if os.path.exists(os.path.join('.', 'output', 'images', 'whole_slide_patches')) and os.path.exists(os.path.join('.', 'output', 'images', 'dens_pvals_per_slide')):

        # Create an expander to hide some optional widgets
        with st.expander('Optional: Image path extraction', expanded=False):

            # # Get the possible analysis radii by analyzing the subdirectories present in the results directory
            # analysis_dir_listing = os.listdir(os.path.join(os.getcwd(), '..', 'results', 'webpage'))
            # st.selectbox('Select analysis radius in microns:', [int(x.lstrip('slices_1x')) for x in analysis_dir_listing], key='analysis_radius_in_microns')

            # Create a button to optionally re-run the image path collection, which is done by default if it hasn't been done already
            collect_image_paths = st.button('Re-run image path extraction')

        # Run the image path collection, assigning the result to the session state
        if ('df_paths_per_slide' not in st.session_state) or collect_image_paths:
            st.session_state['df_paths_per_slide'] = utils.get_paths_for_slides()
            image_paths_extracted = True
        else:
            image_paths_extracted = False

        # Print a message when the image paths have been extracted from disk
        if image_paths_extracted:
            st.info('Image paths extracted from disk')

        # Get shortcut variables to the paths dataframe and a list version of its index
        df_paths_per_slide = st.session_state['df_paths_per_slide']
        slide_names = list(df_paths_per_slide.index)

        # Allow user to select by name the slide data to display
        st.selectbox('Select slide name to visualize...', slide_names, key='slide_name_to_visualize', on_change=update_slide_index, args=(slide_names,))

        # Allow user to select by index the slide data to display
        st.number_input('...OR, select slide index to visualize:', min_value=0, max_value=(len(slide_names) - 1), key='slide_index_to_visualize', on_change=update_slide_name, args=(slide_names,))

        # Split up the window into two columns
        display_col1, display_col2 = st.columns(2)

        # Display the three images for the currently selected slide
        with display_col1:
            st.image(df_paths_per_slide.loc[st.session_state['slide_name_to_visualize'], 'heatmap'])
            st.radio('Display slide patching at right?', ['not patched', 'patched'], key='display_slide_patching')
        with display_col2:
            slide_suffix = ('' if st.session_state['display_slide_patching'] == 'not patched' else '_patched')
            st.image(df_paths_per_slide.loc[st.session_state['slide_name_to_visualize'], 'slide{}'.format(slide_suffix)])

    else:
        st.warning('At least one of the two sets of per-slide plots does not exist; please run all per-slide components of the workflow on the "Run workflow" page', icon='⚠️')

if __name__ == '__main__':
    main()
