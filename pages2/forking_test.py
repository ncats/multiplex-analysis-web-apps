import streamlit as st
import utils
import time


def sleep_task(tuple_of_args):
    sleep_time_sec = tuple_of_args[0]
    print(f'Single task running, waiting for {sleep_time_sec} seconds')
    time.sleep(sleep_time_sec)


def main():

    if 'nworkers' not in st.session_state:
        st.session_state['nworkers'] = 1
    st.number_input('Number of workers', min_value=1, max_value=10, key='nworkers')

    if 'num_tasks' not in st.session_state:
        st.session_state['num_tasks'] = 10
    st.number_input('Number of tasks', min_value=1, max_value=100, key='num_tasks')

    start_time = time.time()
    if st.button('Run forking test'):

        st.write('Forking test button pressed')

        utils.execute_data_parallelism_potentially(sleep_task, [(1,)] * st.session_state['num_tasks'], nworkers=st.session_state['nworkers'], task_description='Sleeping tasks')

        st.write('Done')

    st.write(f'Total time: {time.time() - start_time:.2f} seconds')


if __name__ == '__main__':
    main()
