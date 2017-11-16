# ESCtoolbox (Python version)

This is a Python version of Gregory Plett's enhanced self-correcting (ESC)
battery cell model. The original code is written in Matlab which is available
in the ESCtoolbox at
[mocha-java.uccs.edu/BMS1/](http://mocha-java.uccs.edu/BMS1/).

## OCV model

The open-circuit voltage (OCV) files are located in the `ocv_model/` folder
where the OCV model is the `ocv.py` file. The results and plots generated from
this model should be similar to the Matlab `runProcessOCV.m` plots. The
`funcs.py` file contains the `OCVfromSOCtemp` function while the `models.py`
file contains model objects used by the OCV model. The `data.py` file plots the
experimental data from the csv files located in the `ocv_data/`.

Battery test data for the A123 cell is available in the `ocv_data/` folder as
csv files. The data files were exported from the Excel spreadsheets in the
`OCV_Files/A123_OCV` directory of the Matlab ESCtoolbox. Plots of the battery
test data are generated with the `data.py` script. Change the tc variable to
view plots from the other temperature tests. For example, change the `tc`
string to `N05` to create plots for the CSV files named A123_OCV_N05_S1,
A123_OCV_N05_S2, A123_OCV_N05_S3, and A123_OCV_N05_S4. The figures produced
from the script should be similar to the graphs shown in the `A123_OCV_N05_S1`,
`A123_OCV_N05_S2`, `A123_OCV_N05_S3`, and `A123_OCV_N05_S4` Excel spreadsheets.

See the comments in each file for more information.

## DYN model

The dynamic model files are located in the `dyn_model/` folder where the
`dyn_model.py` is the main file to run. The data from the dynamic battery tests
are located in the `dyn_data/` folder.

See the comments in each file for more information.

## Installation

Requires Python 3.6, Matplotlib, NumPy, and Pandas. The preferred method to
install Python 3 and associated packages is with the Anaconda or Miniconda
distribtion available at
[continuum.io/downloads](https://www.continuum.io/downloads).

## Usage

Clone or download the files to your local machine. Start iPython from within
the `ocv_model/` directory then type `run ocv.py` to run the OCV model.

