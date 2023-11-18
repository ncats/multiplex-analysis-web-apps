# Import relevant libraries
import streamlit as st
from st_pages import show_pages_from_config, add_indentation
import annotations
import os
from streamlit_extras.app_logo import add_logo
import streamlit_utils

def main():

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
    st.title('Annotation plots')

    if os.path.exists(os.path.join('.', 'output', 'images', 'raw_weights_check')) and os.path.exists(os.path.join('.', 'output', 'images', 'all_annotation_data.png')) and os.path.exists(os.path.join('.', 'output', 'images', 'weight_heatmaps_on_annot')) and os.path.exists(os.path.join('.', 'output', 'images', 'pixel_plot')) and os.path.exists(os.path.join('.', 'output', 'images', 'analysis_overlaid_on_annotation')) and os.path.exists(os.path.join('.', 'output', 'images', 'dens_pvals_per_annotation')):

        # Constant: annotation-dependent plot types (there are two more that are annotation-independent; see the "Optional plots" section)
        annotation_dependent_plot_types = ['Average density P value heatmaps', 'Per-pixel "images" of annotation data', 'Annotation weights heatmaps overlaid on annotation data', 'Analysis data overlaid on annotation data']

        # Set some default values in the session state
        streamlit_utils.assign_default_values_in_session_state('plots_to_display_for_each_annot', annotation_dependent_plot_types)
        streamlit_utils.assign_default_values_in_session_state('plot_all_annotation_data', False)
        streamlit_utils.assign_default_values_in_session_state('plot_raw_weights_scatter_plots', False)
        streamlit_utils.assign_default_values_in_session_state('selected_image_id', '')

        # Determine the paths to all annotation-related images
        # df_paths = annotations.get_annotation_plots_paths(top_plot_dir='../results/webpage/slices_1x{}/real'.format(st.session_state['thickness']))
        df_paths = annotations.get_annotation_plots_paths(top_plot_dir=os.path.join('.', 'output', 'images'))

        # Get arrays of the unique parameters on which the annotation plots depend
        image_ids = [x for x in df_paths['image_id'].unique() if x is not None]
        region_types = [x for x in df_paths['region_type'].unique() if x is not None]
        weight_column_prefixes = [x[:-1] for x in df_paths['weight_column_prefix'].unique() if x is not None]
        # do_log_transforms = [x for x in df_paths['do_log_transform'].unique() if x is not None]

        # Set the default values of the main plot parameter widgets
        if st.session_state['selected_image_id'] not in image_ids:
            st.session_state['selected_image_id'] = image_ids[0]
        streamlit_utils.assign_default_values_in_session_state('selected_image_id', image_ids[0])
        streamlit_utils.update_button_activities(image_ids, selected_option='selected_image_id', button_name='image_button')
        streamlit_utils.assign_default_values_in_session_state('selected_weight_column_prefix', 'footprint_integer_area')
        streamlit_utils.update_button_activities(weight_column_prefixes, selected_option='selected_weight_column_prefix', button_name='weight_type_button')
        streamlit_utils.assign_default_values_in_session_state('selected_do_log_transform', True)

        # Top row of user selection options (widgets)
        cols_top = st.columns(3)
        with cols_top[0]:
            st.selectbox('Image:', image_ids, key='selected_image_id', on_change=streamlit_utils.update_button_activities, args=(image_ids,), kwargs={'selected_option': 'selected_image_id', 'button_name': 'image_button'})
        with cols_top[1]:
            st.selectbox('ROI weight type:', weight_column_prefixes, key='selected_weight_column_prefix')
        with cols_top[2]:
            st.multiselect('Plots to display for each annotation region type:', annotation_dependent_plot_types, key='plots_to_display_for_each_annot')

        # Bottom row of user selection options
        cols_bottom = st.columns(6)
        with cols_bottom[0]:
            st.button('Previous image', use_container_width=True, on_click=streamlit_utils.previous_option, disabled=st.session_state['disable_previous_image_button'], args=(image_ids,), kwargs={'selected_option': 'selected_image_id', 'button_name': 'image_button'})
        with cols_bottom[1]:
            st.button('Next image', use_container_width=True, on_click=streamlit_utils.next_option, disabled=st.session_state['disable_next_image_button'], args=(image_ids,), kwargs={'selected_option': 'selected_image_id', 'button_name': 'image_button'})
        with cols_bottom[2]:
            st.button('Previous weight type', use_container_width=True, on_click=streamlit_utils.previous_option, disabled=st.session_state['disable_previous_weight_type_button'], args=(weight_column_prefixes,), kwargs={'selected_option': 'selected_weight_column_prefix', 'button_name': 'weight_type_button'})
        with cols_bottom[3]:
            st.button('Next weight type', use_container_width=True, on_click=streamlit_utils.next_option, disabled=st.session_state['disable_next_weight_type_button'], args=(weight_column_prefixes,), kwargs={'selected_option': 'selected_weight_column_prefix', 'button_name': 'weight_type_button'})
        with cols_bottom[4]:
            st.toggle('Log transform the weights', key='selected_do_log_transform')

        # Initialize the same number of columns as annotation region types
        annot_cols = st.columns(len(region_types))

        # For each region type / annotation column...
        for annot_col, region_type in zip(annot_cols, region_types):

            # In the current column...
            with annot_col:

                # Display the region type
                st.subheader(region_type)

                # Plot the average density P value heatmaps corresponding to the current region type and user selections
                if 'Average density P value heatmaps' in st.session_state['plots_to_display_for_each_annot']:
                    df_result = df_paths[(df_paths['plot_type'] == 'average_p_value_heatmap') &
                                        (df_paths['region_type'] == region_type) &
                                        (df_paths['image_id'] == st.session_state['selected_image_id']) &
                                        (df_paths['weight_column_prefix'] == (st.session_state['selected_weight_column_prefix'] + '_')) &
                                        (df_paths['do_log_transform'] == st.session_state['selected_do_log_transform'])]
                    if len(df_result) == 1:
                        st.image(df_result.iloc[0]['path'])
                    else:
                        st.warning('Image of average density P value heatmap is not present', icon="⚠️")

                # Plot the annotation weights heatmaps overlaid on the annotation data
                if 'Annotation weights heatmaps overlaid on annotation data' in st.session_state['plots_to_display_for_each_annot']:
                    df_result = df_paths[(df_paths['plot_type'] == 'weights_heatmap_overlay') &
                                        (df_paths['region_type'] == region_type) &
                                        (df_paths['image_id'] == st.session_state['selected_image_id']) &
                                        (df_paths['weight_column_prefix'] == (st.session_state['selected_weight_column_prefix'] + '_')) &
                                        (df_paths['do_log_transform'] == st.session_state['selected_do_log_transform'])]
                    if len(df_result) == 1:
                        st.image(df_result.iloc[0]['path'])
                    else:
                        st.warning('Image of annotation weights heatmap overlaid on annotation data is not present', icon="⚠️")

                # Plot per-pixel "images" of the annotation data
                if 'Per-pixel "images" of annotation data' in st.session_state['plots_to_display_for_each_annot']:
                    df_result = df_paths[(df_paths['plot_type'] == 'pixel_plot') &
                                        (df_paths['region_type'] == region_type) &
                                        (df_paths['image_id'] == st.session_state['selected_image_id'])]
                    if len(df_result) == 1:
                        st.image(df_result.iloc[0]['path'])
                    else:
                        st.warning('Per-pixel "image" of annotation data is not present', icon="⚠️")

                # Plot the analysis data overlaid on the annotation data
                if 'Analysis data overlaid on annotation data' in st.session_state['plots_to_display_for_each_annot']:
                    df_result = df_paths[(df_paths['plot_type'] == 'analysis_overlay') &
                                        (df_paths['region_type'] == region_type) &
                                        (df_paths['image_id'] == st.session_state['selected_image_id'])]
                    if len(df_result) == 1:
                        st.image(df_result.iloc[0]['path'])
                    else:
                        st.warning('Image of analysis data overlaid on annotation data is not present', icon="⚠️")

        # Add a new section for the two sets of optional plots
        st.divider()
        st.header('Optional plots')

        # Allow the user to select whether the optional plots should be displayed since they may slow down page redraws (albeit probably trivially)
        cols_checkboxes = st.columns(6)
        with cols_checkboxes[1]:
            st.checkbox('Plot all annotation data', key='plot_all_annotation_data')
        with cols_checkboxes[4]:
            st.checkbox('Plot scatter plots of the raw weights', key='plot_raw_weights_scatter_plots')

        # Plot the one image of all the annotation data for all images and annotation types
        cols_other_images = st.columns(2)
        with cols_other_images[0]:
            if st.session_state['plot_all_annotation_data']:
                df_result = df_paths[(df_paths['plot_type'] == 'all_annotation_data')]
                if len(df_result) == 1:
                    st.image(df_result.iloc[0]['path'])
                else:
                    st.warning('Image of all annotation data is not present', icon="⚠️")

        # Plot one of two images showing the scatter plots of the raw weights for all images and region types
        with cols_other_images[1]:
            if st.session_state['plot_raw_weights_scatter_plots']:
                df_result = df_paths[(df_paths['plot_type'] == 'raw_weights_scatter_plots') &
                                    (df_paths['do_log_transform'] == st.session_state['selected_do_log_transform'])]
                if len(df_result) == 1:
                    st.image(df_result.iloc[0]['path'])
                else:
                    st.warning('Image of the scatter plots of the raw weights is not present', icon="⚠️")

    else:
        st.warning('The component "Average density P values over ROIs for each annotation region type" of the workflow does not appear to have been run; please select it on the "Run workflow" page', icon='⚠️')

if __name__ == '__main__':
    main()
