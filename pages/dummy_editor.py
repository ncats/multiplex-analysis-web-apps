import streamlit as st
import streamlit_dataframe_editor as sde

def data_editor_change_callback():
    st.session_state['saved_dataeditor_values'] = st.session_state['dataeditor__do_not_persist']

def update_input_data_editor():
    for key, value in st.session_state['saved_dataeditor_values']['edited_rows'].items():
        for key2, value2 in value.items():
            st.session_state.df_dataeditor_input.loc[key, key2] = value2

def main():

    # Run streamlit-dataframe-editor library initialization tasks at the top of the page
    st.session_state = sde.initialize_session_state(st.session_state)

    if 'saved_dataeditor_values' in st.session_state:
        update_input_data_editor()

    df_edited = st.data_editor(st.session_state.df_dataeditor_input,
                               key='dataeditor__do_not_persist',
                               on_change=data_editor_change_callback)

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)

if __name__ == '__main__':
    main()
