'''
Class for benchmarking the analyses used within MAWA.
Specifically to mark the time it takes for functions to run
and save their values in a spreadsheet (if wanted)
'''

import os
import time
import numpy as np
import pandas as pd
from datetime import datetime

class benchmark_collector:
    '''
    benchmark_collector (bc) class. Used for timing the execution
    of functions within MAWA.
    '''
    def __init__(self, fiol = None):
        '''
        Initialize the benchmarking class
        '''

        column_names = ['id',
                'on_NIDAP', 
                'file',
                'nSlides', 
                'nCells',
                'CellsxSlide',
                'time_load_data',
                'time_to_run_counts',
                'time_to_run_UMAP',
                'time_to_run_cluster']
        self.benchmarkDF = pd.DataFrame(columns = column_names)

        self.fiol = fiol
        if self.fiol is None:
            on_nidap = False
        else:
            on_nidap = self.fiol.onNIDAP

        if on_nidap:
            localdir = './input'
        else:
            localdir = "C:/Users/smithdaj/OneDrive - National Institutes of Health/Documents - NCATS-NCI-DMAP/MAWA/"

        self.benchmark_csv = os.path.join(localdir, 'MAWA_Suite_Benchmarking.csv')

        if os.path.exists(self.benchmark_csv) is False:
            print('Could not find benchmarking file, creating new one')
            self.create_new_csv()

        self.benchmark_project_path = '/NIH/Data Management & Analysis Program (DMAP)/benchmarking/'
        self.benchmark_dataset      = 'Neighborhood-Profiles-Benchmarks'

        self.benchmarkDF.loc[0, 'id']       = datetime.now()
        self.benchmarkDF.loc[0, 'on_NIDAP'] = on_nidap
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

    def printElapsedTime(self, msg):
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

    def create_new_csv(self):
        '''
        Create a new csv file for benchmarking
        '''
        self.benchmarkDF.to_csv(self.benchmark_csv, mode='w', index=False)
        print('Created new Benchmarking csv')

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
