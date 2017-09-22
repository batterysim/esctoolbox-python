"""
Plot voltages recorded for the three scripts related to the A123 cell. Data is
read from csv files located in the dyn_data/ folder. The csv files were created
from the mat files located in the Matlab ESCtoolbox at DYN/A123_DYN. 
"""

import matplotlib.pyplot as plt
import pandas as pd


# list of csv files based on magnitude and temperature
magtemps = ['10_N15', '10_N25', '30_N05', '45_P05', '45_P15', '50_P25', '50_P35', '50_P45']

# choose which group of files to plot, index should be a number from 0-7
mt = magtemps[0]

# plot data from the csv files
plt.ion()
plt.close('all')

for s in ['s1', 's2', 's3']:
    name = f'A123_DYN_{mt}_{s}'
    nfile = f'../dyn_data/{name}.csv'
    df = pd.read_csv(nfile, sep=', ', engine='python')
    volt = df['voltage'].values
    time = df['time'].values
    t = (time - time[0])/3600

    plt.figure()
    plt.plot(t, volt, lw=2)
    plt.xlabel('Time (min)')
    plt.ylabel('Voltage (V)')
    plt.title(name)


