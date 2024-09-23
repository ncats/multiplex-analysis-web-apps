import streamlit as st
import time
import multiprocessing
import os

def workder(n):
    print(f'Worker {n} started in process {os.getpid()}')
    time.sleep(0.5)
    print(f'Worker {n} finished in process {os.getpid()}')


def main():
    if 'nworkers' not in st.session_state:
        st.session_state['nworkers'] = 1
    st.number_input('Number of workers', min_value=1, max_value=10, key='nworkers')

    if 'num_tasks' not in st.session_state:
        st.session_state['num_tasks'] = 10
    st.number_input('Number of tasks', min_value=1, max_value=100, key='num_tasks')

    if st.button('Run forking test'):
        start_time = time.time()
        st.write('Forking test button pressed')
        
        multiprocessing.set_start_method('forkserver', force=True)
        with multiprocessing.Pool(st.session_state['nworkers']) as pool:
            pool.map(workder, range(st.session_state['num_tasks']))

        st.write('Done')

        st.write(f'Total time: {time.time() - start_time:.2f} seconds')


if __name__ == '__main__':
    multiprocessing.set_start_method('forkserver', force=True)
    main()
