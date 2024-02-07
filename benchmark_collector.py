'''
Class for benchmarking the analyses used within MAWA.
Specifically to mark the time it takes for functions to run
and save their values in a spreadsheet (if wanted)
'''

import time
import numpy as np
import pandas as pd
from datetime import datetime

class benchmark_collector:
    '''
    benchmark_collector (bc) class. Used for timing the execution
    of functions withing MAWA.
    '''
    def __init__(self, fiol = None):
        '''
        Initialize the benchmarking class
        '''

        self.fiol = fiol
        if self.fiol is None:
            on_nidap = False
        else:
            on_nidap = self.fiol.onNIDAP

        self.benchmark_csv = "C:/Users/smithdaj/OneDrive - National Institutes of Health/Documents - NCATS-NCI-DMAP/MAWA/MAWA_Suite_Benchmarking.csv"
        self.benchmark_project_path = '/NIH/Data Management & Analysis Program (DMAP)/benchmarking/'
        self.benchmark_dataset     = 'Neighborhood-Profiles-Benchmarks'

        d = {'id': [datetime.now()],
             'on_NIDAP': [on_nidap],
             'file': [None],
             'nSlides': [None],
             'nCells': [None],
             'CellsxSlide': [None],
             'time_load_data': [None],
             'time_to_run_counts': [None],
             'time_to_run_UMAP': [None],
             'time_to_run_cluster': [None]}
        self.benchmarkDF = pd.DataFrame(data = d)

        self.stTimer = None
        self.spTimer = None

    def startTimer(self):
        '''
        Set the Start time to the current date-time
        '''
        self.stTimer = time.time()

    def stopTimer(self):
        '''
        Set the Stop time to the current date-time
        '''
        self.spTimer = time.time()

    def elapsedTime(self):
        '''
        Calculate the elapsed time from the spTimer and the stTimer
        '''
        if self.stTimer is not None:
            self.stopTimer()
            elapsed_time =  np.round(self.spTimer - self.stTimer, 3)
        else:
            elapsed_time = None
        return elapsed_time

    def print_elapsed_time(self, msg):
        '''
        Print the current value of elapsed time 
        '''
        print(f'{msg} took {self.elapsedTime()} s')

    def check_df(self):
        '''
        Check the current head of benchmark dataframe
        '''
        print(self.benchmarkDF.head())

    def set_value_df(self, field, value):
        '''
        Add a field/value combo to the dataframe
        '''
        self.benchmarkDF[field] = value

    def save_run_to_csv(self):
        '''
        Saves the benchmarking datafile to a csv
        '''
        self.benchmarkDF.to_csv(self.benchmark_csv, mode='a', index=False, header=False)
        print('Saved run to Benchmarking csv')

    def send_to_nidap(self):
        '''
        Export the benchmarking dataframe to NIDAP
        '''
        self.fiol.export_results_dataset(self.benchmarkDF,
                                         path = self.benchmark_project_path,
                                         filename = self.benchmark_dataset,
                                         saveCompass = True,
                                         type = 'S',
                                         create = False)
