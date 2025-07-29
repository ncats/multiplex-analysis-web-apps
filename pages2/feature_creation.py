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
            st.error("Please select at two columns to concatenate.")
    else:
        st.error("Please enter a valid column name.")

def update_list_columns():
    '''
    Function to update the list of columns in the session state.
    '''

    st.session_state['fc_input_df_columns'] = st.session_state['input_dataset'].data.columns.tolist()

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
        st.dataframe(st.session_state['input_dataset'].data[selected_columns])
    else:
        st.info("No columns selected.")
if __name__ == '__main__':
    main()
