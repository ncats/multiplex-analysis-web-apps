import numpy as np
import platform_io
import subprocess
from streamlit_javascript import st_javascript
from streamlit_extras.app_logo import add_logo

def top_of_page_reqs(session_state):
    '''
    Top of the page requirements. These are the commands
    that are run at the top of the page that help maintain 
    '''

    # # Remove key values from session_state that should not persist
    # for key, val in session_state.items():
    #     if (not key.endswith('__do_not_persist')) and (not key.startswith('FormSubmitter:')):
    #         session_state[key] = val
    # Add Logo
    # add_logo('app_images/mawa_logo-width315.png', height=150)
    
    # Get the URL of the current page, per https://discuss.streamlit.io/t/what-is-the-current-page-in-use-multipage-app/41898, which we've used before but keeping the reference here anyway. Remember something strange happens here, with the script running twice or the like and not picking up the session state fully or vice versa, so sometimes we see strange behavior as a result though it's usually not a bit deal. In this case I've taken care of it anyway below by resetting the index on default_df_contents2
    session_state.curr_url = st_javascript("await fetch('').then(r => window.parent.location.href)")

    # Check the platform
    session_state = check_for_platform(session_state)

    return session_state

def platform_is_nidap():
    return np.any(['nidap.nih.gov' in x for x in subprocess.run('conda config --show channels', shell=True, capture_output=True).stdout.decode().split('\n')[1:-1]])

def check_for_platform(session_state):
    # Initialize the platform object
    if 'platform' not in session_state:
        session_state['platform'] = platform_io.Platform(platform=('nidap' if platform_is_nidap() else 'local'))
    return session_state