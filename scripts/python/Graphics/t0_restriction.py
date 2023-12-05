import os 

import numpy as np 
import pandas as pd

PATH = os.getcwd() # main path
print(PATH)

e_signals = pd.read_csv(os.path.join(PATH, 'data_preproc/LightSignal_decomp_e.csv'), 
                        header=0)
mu_signals = pd.read_csv(os.path.join(PATH, 'data_preproc/LightSignal_decomp_mu.csv'))

print(e_signals.head())