'''
Script for creating custom features for your dataset
'''

import streamlit as st

def create_new_column():
    '''
    Create a new column in the input dataset.
    '''
    if st.session_state['fs_new_column_name']:
        if len(st.session_state['fs_columns_to_concatenate']) == 2:

            # Concatenate the selected columns (assume two to start)
            new_column = st.session_state['fs_columns_to_concatenate'][0] + \
                         st.session_state['input_dataset'].data[st.session_state['fs_columns_to_concatenate'][0]].astype(str) + \
                         '__' + \
                         st.session_state['fs_columns_to_concatenate'][1] + \
                         st.session_state['input_dataset'].data[st.session_state['fs_columns_to_concatenate'][1]].astype(str)

            st.session_state['input_dataset'].data[st.session_state['fs_new_column_name']] = new_column
            st.success(f"Added new column: {st.session_state['fs_new_column_name']}")
            update_list_columns()

            st.session_state['fc_columns_select'].append(st.session_state['fs_new_column_name'])

            # Reset the input field
            st.session_state['fs_new_column_name'] = ''
        else:
            st.error("Please select exactly two columns to concatenate.")
    else:
        st.error("Please enter a valid column name.")

def update_list_columns():
    '''
    Function to update the list of columns in the session state.
    '''

    st.session_state['fc_input_df_columns'] = st.session_state['input_dataset'].data.columns.tolist()

def update_image_list():
    '''
    Function to update the list of images in the session state.
    '''

    # Check if the column 'Image ID_(standardized)' exists
    if 'Image ID_(standardized)' in st.session_state['input_dataset'].data.columns:
        st.session_state.fc_image_select_disabled = False
        st.session_state['fc_input_df_images'] = st.session_state['input_dataset'].data['Image ID_(standardized)'].unique().tolist()

        # Append 'All Images' to the front of this list
        st.session_state['fc_input_df_images'].insert(0, 'All Images')
    else:
        st.session_state.fc_image_select_disabled = True

def main():
    '''
    Main function for creating custom features
    '''

    if 'input_dataset' not in st.session_state:
        st.error('An input dataset has not yet been opened. Please do so using the "Open File" page in the sidebar.')
        return

    # Show Columns
    st.subheader('Available Columns')

    if 'fc_input_df_columns' not in st.session_state:
        update_list_columns()
        st.session_state['fc_columns_select'] = st.session_state['fc_input_df_columns']

    update_image_list()

    with st.expander("Select columns to display", expanded=True):
        selected_columns = st.multiselect(
            'Select columns:',
            options=st.session_state['fc_input_df_columns'],
            key='fc_columns_select'
        )

    # Add a new column in a third-width column
    cols = st.columns([1, 1, 1])
    with cols[0]:
        st.text_input("New column name:", key='fs_new_column_name')
        st.button("Add Column", on_click=create_new_column)
    with cols[1]:
        st.selectbox(
            "Type of operation:",
            options=['Concatenate'],
            key='fs_operation_type'
        )
    with cols[2]:
        st.multiselect(
            "Select columns to concatenate:",
            options=st.session_state['fc_columns_select'],
            key='fs_columns_to_concatenate'
        )

    if selected_columns:
        st.subheader(f"Number of cells: {st.session_state['input_dataset'].data.shape[0]}")
        st.dataframe(st.session_state['input_dataset'].data[selected_columns])
    else:
        st.info("No columns selected.")

    # Add a dataframe which displays the number of unique values in selected columns
    cols = st.columns([1, 1, 1])
    with cols[1]:
        if selected_columns:
            unique_counts = st.session_state['input_dataset'].data[selected_columns].nunique()
            unique_counts = unique_counts.rename('Unique Values')
            st.dataframe(unique_counts)
        else:
            st.info("No columns selected.")

    # Draw a bar chat of the count of each unique value for a given column
    with cols[2]:

        sub_cols = st.columns([1, 1])
        with sub_cols[0]:
            unique_bar_col = st.selectbox(
                "Select column for unique value counts:",
                options=st.session_state['fc_input_df_columns'],
                key='fs_unique_bar_col'
            )

        with sub_cols[1]:
            if not st.session_state.fc_image_select_disabled:
                image_to_view = st.selectbox(
                    "Select image to view:",
                    options=st.session_state['fc_input_df_images'],
                    key='fs_image_to_view'
                )
            else:
                image_to_view = 'All Images'
                st.write("No Column named 'Image ID_(standardized)'. Please complete the Datafile Unification")

        if unique_bar_col:

            # Find max value from unique counts across all images
            max_unique_count = st.session_state['input_dataset'].data[unique_bar_col].value_counts().max()

            if image_to_view == 'All Images':
                bar_chart_data = st.session_state['input_dataset'].data[unique_bar_col].value_counts()
            else:
                bar_chart_data = st.session_state['input_dataset'].data[unique_bar_col][st.session_state['input_dataset'].data['Image ID_(standardized)'] == image_to_view].value_counts()

            st.subheader(f"Feature: {unique_bar_col}")
            st.write(f'''**Image**: {image_to_view}  
                         **Number of cells**: {bar_chart_data.sum()}''')
            st.bar_chart(bar_chart_data)
        else:
            st.info("No columns selected.")

if __name__ == '__main__':
    main()
