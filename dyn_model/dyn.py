"""
Dyn model
"""

import pickle
import numpy as np

from models import DataModel
from funcs import processDynamic

# parameters
cellID = 'A123'     # cell identifier
numpoles = 1        # number of resistor-capacitor pairs in final model
temps = [-25, -15, -5, 5, 15, 25, 35, 45]   # temperatures
mags = [10, 10, 30, 45, 45, 50, 50, 50]     # A123

# read model OCV file, previously computed by runProcessOCV
modelocv = pickle.load(open('../ocv_model/modelocv.pickle', 'rb'))

# initialize array to store battery cell data
data = np.zeros(len(mags), dtype=object)

# load battery cell data for each temperature as objects then store in data array
# note that data files are in the dyn_data folder
for idx, temp in enumerate(temps):
    mag = mags[idx]
    if temp < 0:
        tempfmt = f'{abs(temp):02}'
        files = [f'../dyn_data/{cellID}_DYN_{mag}_N{tempfmt}_s1.csv',
                 f'../dyn_data/{cellID}_DYN_{mag}_N{tempfmt}_s2.csv',
                 f'../dyn_data/{cellID}_DYN_{mag}_N{tempfmt}_s3.csv']
        data[idx] = DataModel(temp, files)
        print('Loaded files', *files, sep='\n')
    else:
        tempfmt = f'{abs(temp):02}'
        files = [f'../dyn_data/{cellID}_DYN_{mag}_P{tempfmt}_s1.csv',
                 f'../dyn_data/{cellID}_DYN_{mag}_P{tempfmt}_s2.csv',
                 f'../dyn_data/{cellID}_DYN_{mag}_P{tempfmt}_s3.csv']
        data[idx] = DataModel(temp, files)
        print('Loaded files', *files, sep='\n')

processDynamic(data, modelocv, numpoles, 2)

