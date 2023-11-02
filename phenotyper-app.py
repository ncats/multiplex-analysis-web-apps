'''
Python script for performing simple phenotyping
'''
import streamlit as st
from streamlit_javascript import st_javascript

# Import relevant libraries
import nidap_dashboard_lib as ndl   # Useful functions for dashboards connected to NIDAP

def main():
    '''
    Main function for running the app
    '''
    # Use the whole page width
    st.set_page_config(page_title="Welcome",
                       layout="wide")

    for key, val in st.session_state.items():
        if (not key.endswith('__do_not_persist')) and (not key.startswith('FormSubmitter:')):
            st.session_state[key] = val

    # Initalize environment variables and session-state variables
    if 'init' not in st.session_state:
        settings_yaml_file = 'config_files/TaIS.yml'
        # Initialize session_state values for streamlit processing
        st.session_state = ndl.init_session_state(st.session_state, settings_yaml_file)

    ### SIDE BAR ORGANIZATION ###
    with st.sidebar:
        url = st_javascript("await fetch('').then(r => window.parent.location.href)")
        st.write(f'''[Open app in new Tab]({url})\n (MS Edge/ Google Chrome)''')

    intro_markdown = ndl.read_markdown_file('markdown/Phenotyping_TutorialWelcome.md')
    st.markdown(intro_markdown, unsafe_allow_html=True)

if __name__ == '__main__':
    main()
