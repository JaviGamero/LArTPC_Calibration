import os 
import sys 

# path to python main folder in this project
libraries = os.path.abspath(os.path.join(os.path.abspath(os.path.dirname(__file__)), 
                                         os.pardir)) 
sys.path.append(libraries) 

import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
from time import time

################################################################################
################################################################################
print('\nThis script checks how many of the total series accomplish the ',
      'differents criteria of this project.\n')

################################################################################
################################################################################

PATH = os.path.abspath(os.path.join(os.getcwd())) # main path

t_path = os.path.join(PATH, 'data_preproc/LightSignal_t.csv')
t = pd.read_csv(t_path, sep=';', header=None) # t[0] is nonsensen, remove it
t.set_index(0, inplace=True)

t0 = 150 # (ns), moment to start considering the slow component, EXPERIMENTAL
# t_idx = np.where(t>t0)[1] 
t = np.array(t.iloc[0, :]).reshape(-1) #1D array from t0 and on

e_GT_path = os.path.join(PATH, 'data_preproc/LightSignal_decomp_e.csv')
e_signals = pd.read_csv(e_GT_path, sep=';', header=0, index_col=0)
# e_signals.set_index(0, inplace=True)

count_total = e_signals.shape[0]

################################################################################
################################################################################

count_criteria_1 = e_signals.shape[0] # electron energy > 0.001
count_criteria_2 = 0 # peak of the electron signal in a moment > t0
count_criteria_3 = 0 # height of the peak of the electron > 5
count_both = 0 # both criterias 2 and 3
for idx in e_signals.index: 
    max_idx = np.argmax(e_signals.loc[idx, :])
    max_peak = np.max(e_signals.loc[idx, :])
    
    if t[max_idx]>t0:
        count_criteria_2+=1
        
    if max_peak > 5: 
        count_criteria_3+=1 
        
    if (t[max_idx]>t0) and (max_peak > 5): 
        count_both+=1 
 
print('Total number of series:', count_total)   
print('Number of series that accomplish criteria 1:', count_criteria_1)
print('Number of series that accomplish criteria 2:', count_criteria_2)
print('Number of series that accomplish criteria 3:', count_criteria_3)
print('Number of series that accomplish all criterias:', count_both)

print()

print('Ratio of series that accomplish criteria 1:', count_criteria_1/count_total*100, '%')
print('Ratio of series that accomplish criteria 2:', 100 - count_criteria_2/count_total*100, '%')
print('Ratio of series that accomplish criteria 3:', 100 - count_criteria_3/count_total*100, '%')
print('Ratio of series that accomplish all criterias:', 100 - count_both/count_total*100, '%')

print()
# print(t[0], t[1])