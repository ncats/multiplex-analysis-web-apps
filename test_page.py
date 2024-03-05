# Import relevant libraries
import streamlit as st
import pandas as pd
import app_top_of_page as top
import streamlit_dataframe_editor as sde

def main():
    """
    Main function.
    """

    # Set page settings
    st.set_page_config(layout='wide', page_title='Test Page')
    st.title('Test Page')

    # Run streamlit-dataframe-editor library initialization tasks at the top of the page
    st.session_state = sde.initialize_session_state(st.session_state)

    # Run Top of Page (TOP) functions
    st.session_state = top.top_of_page_reqs(st.session_state)

    if 'dataframe_editor' not in st.session_state:
        st.session_state['dataframe_editor'] = sde.DataframeEditor(df_name='my_dataframe', default_df_contents=pd.DataFrame({'a': [1, 2, 3], 'b': [4, 5, 6]}))
    st.session_state['dataframe_editor'].dataframe_editor()
    st.write(st.session_state['dataframe_editor'].reconstruct_edited_dataframe())

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)

# Run the main function
if __name__ == '__main__':
    main()
