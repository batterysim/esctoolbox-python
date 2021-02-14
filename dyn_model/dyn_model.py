"""
Dyn model
"""

import numpy as np
import json

from funcs import processDynamic
from models import DataModel, ModelOcv

from pathlib import Path

# parameters
cellID = 'A123'     # cell identifier
numpoles = 3        # number of resistor-capacitor pairs in final model
temps = [-25, -15, -5, 5, 15, 25, 35, 45]   # temperatures
mags = [10, 10, 30, 45, 45, 50, 50, 50]     # A123
doHyst = 1          # 1 "find M, M0 and G params" or 0 "make hys params 0" 

# read model OCV file, previously computed by runProcessOCV
modelocv = ModelOcv.load(Path(f'./modelocv.json'))

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

modeldyn = processDynamic(data, modelocv, numpoles, doHyst)

# convert ocv and dyn results model object to dict, then save in JSON to disk 
modeldyn = {k:v.tolist() if isinstance(v, np.ndarray) else v for k,v in modeldyn.__dict__.items()}
if True:
    if doHyst:
        with open('modeldyn.json', 'w') as json_file:
            json.dump(modeldyn, json_file, indent=4)
    else:
        with open('modeldyn-no-hys.json', 'w') as json_file:
            json.dump(modeldyn,json_file, indent=4)