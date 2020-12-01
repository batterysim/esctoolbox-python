"""
Dyn model
"""

import numpy as np
import pickle

from funcs import processDynamic
from models import DataModel, ModelOcv

from pathlib import Path

# parameters
cellID = 'A123'     # cell identifier
numpoles = 1        # number of resistor-capacitor pairs in final model
temps = [-25, -15, -5, 5, 15, 25, 35, 45]   # temperatures
mags = [10, 10, 30, 45, 45, 50, 50, 50]     # A123

# read model OCV file, previously computed by runProcessOCV
modelocv = ModelOcv.load(Path(f'./ocv_model/modelocv.pickle'))

# initialize array to store battery cell data
data = np.zeros(len(mags), dtype=object)

# load battery cell data for each temperature as objects then store in data array
# note that data files are in the dyn_data folder
print('Load files')
for idx, temp in enumerate(temps):
    mag = mags[idx]
    if temp < 0:
        tempfmt = f'{abs(temp):02}'
        files = [Path(f'./dyn_data/{cellID}_DYN_{mag}_N{tempfmt}_s1.csv'),
                 Path(f'./dyn_data/{cellID}_DYN_{mag}_N{tempfmt}_s2.csv'),
                 Path(f'./dyn_data/{cellID}_DYN_{mag}_N{tempfmt}_s3.csv')]
        data[idx] = DataModel(temp, files)
        print(*files, sep='\n')
    else:
        tempfmt = f'{abs(temp):02}'
        files = [Path(f'./dyn_data/{cellID}_DYN_{mag}_P{tempfmt}_s1.csv'),
                 Path(f'./dyn_data/{cellID}_DYN_{mag}_P{tempfmt}_s2.csv'),
                 Path(f'./dyn_data/{cellID}_DYN_{mag}_P{tempfmt}_s3.csv')]
        data[idx] = DataModel(temp, files)
        print(*files, sep='\n')

modeldyn = processDynamic(data, modelocv, numpoles, 1)

# store ocv and dyn results in model object then save to disk
pickle.dump(modeldyn, open('modeldyn.pickle', 'wb'))