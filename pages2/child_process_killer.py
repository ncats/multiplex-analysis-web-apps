# Import relevant libraries
import streamlit as st
import app_top_of_page as top
import streamlit_dataframe_editor as sde
import os
import subprocess
import psutil
import pandas as pd
import time

def get_system_info():
    # Run the top command and get its output
    output = subprocess.check_output(
        ["top", "-b", "-n", "1"], 
        universal_newlines=True
    )

    # Split the output into lines
    lines = output.split("\n")

    # Filter out the lines that contain the process information
    process_lines = [line for line in lines if line.startswith(" ")]

    # Parse the process information into a list of dictionaries
    processes = []
    for line in process_lines:
        parts = line.split()
        processes.append({
            'pid': parts[0],
            'name': parts[11],
            'cpu_percent': parts[8],
            'memory_percent': parts[9],
        })

    # Convert the list of processes to a DataFrame
    df = pd.DataFrame(processes)

    # Calculate the total CPU and memory usage
    cpu_total = df['cpu_percent'].sum()
    memory_total = df['memory_percent'].sum()

    # Add the total usage to the DataFrame
    new_row = pd.DataFrame([{'pid': 'Total', 'name': '', 'cpu_percent': cpu_total, 'memory_percent': memory_total}])
    df = pd.concat([df, new_row], ignore_index=True)

    return df

    
def kill_child_processes(dry_run=False):
    parent_pid = os.getpid()  # Get the process ID of the current process

    if dry_run:
        # Get a list of child processes
        parent = psutil.Process(parent_pid)
        children = parent.children(recursive=True)
        for child in children:
            print(f'Would kill process ID {child.pid}, name {child.name()}')
    else:
        for child in psutil.Process(parent_pid).children(recursive=True):
            child.kill()  # Send a SIGKILL signal


def main():
    """
    Main function for the page.
    """

    if st.button('Show system info'):
        df = get_system_info()
        st.dataframe(df)
    
    if st.button('Show what child processes *would* be killed if the following button is clicked'):
        kill_child_processes(dry_run=True)
    
    if st.button('Kill child processes'):
        kill_child_processes()


# Run the main function
if __name__ == '__main__':

    # Set page settings
    page_name = 'Child Process Killer'
    st.set_page_config(layout='wide', page_title=page_name)
    st.title(page_name)

    # Run streamlit-dataframe-editor library initialization tasks at the top of the page
    st.session_state = sde.initialize_session_state(st.session_state)

    # Run Top of Page (TOP) functions
    st.session_state = top.top_of_page_reqs(st.session_state)

    # Call the main function
    main()

    # Run streamlit-dataframe-editor library finalization tasks at the bottom of the page
    st.session_state = sde.finalize_session_state(st.session_state)
