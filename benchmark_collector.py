import time
import numpy as np
import pandas as pd
from datetime import datetime

class benchmark_collector:
    def __init__(self, fiol = None):

        self.fiol = fiol
        if self.fiol is None:
            onNIDAP = False
        else:
            onNIDAP = self.fiol.onNIDAP
            
        self.benchmark_csv = "C:/Users/smithdaj/OneDrive - National Institutes of Health/General/neighborhood-profiles/Neighborhood-Profiles-Benchmarks.csv"
        self.benchmark_project_path = '/NIH/Data Management & Analysis Program (DMAP)/benchmarking/'
        self.benchmark_dataset     = 'Neighborhood-Profiles-Benchmarks'

        d = {'id': [datetime.now()],
             'on_NIDAP': [onNIDAP],
             'NUMBA_NUM_THREADS': ['NaN'],
             'file': ['NaN'],
             'nRows': ['NaN'],
             'data_import_loc': ['NaN'],
             'time_load_data': ['NaN'],
             'time_finish_draw': ['NaN'],
             'time_to_run_counts': ['NaN'],
             'time_to_run_UMAP': ['NaN'],
             'time_to_run_cluster': ['NaN'],
             'counts_multiprocess': [False],
             'cpu_pool_size': [1],
             'new_pooling': [True]}
        self.benchmarkDF = pd.DataFrame(data = d)

        self.stTimer = 'None'
        self.spTimer = 'None'
    
    def startTimer(self):
        self.stTimer = time.time()

    def stopTimer(self):
        self.spTimer = time.time()

    def printElapsedTime(self, msg):
        if self.stTimer is not None:
            self.stopTimer()
            elapsedTime = np.round(self.spTimer - self.stTimer, 3)
            print(f'{msg} took {elapsedTime} s')



    def check_df(self):
        print(self.benchmarkDF.head())

    def set_value_df(self, field, value):
        self.benchmarkDF[field] = value

    def save_run_toCsv(self):
        self.benchmarkDF.to_csv(self.benchmark_csv, mode='a', index=False, header=False)
        print('Saved run to Benchmarking csv')

    def send_to_NIDAP(self):
        self.fiol.export_results_dataset(self.benchmarkDF, 
                                         path = self.benchmark_project_path, 
                                         filename = self.benchmark_dataset, 
                                         saveCompass = True,
                                         type = 'S',
                                         create = False)