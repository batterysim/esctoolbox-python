"""
Models for the DYN calculations
"""

import pandas as pd


class DataModel:
    """
    Data from battery script tests.
    """

    def __init__(self, temp, csvfiles):
        """
        Initialize from script data.
        """
        self.temp = temp
        self.s1 = Script(csvfiles[0])
        self.s2 = Script(csvfiles[1])
        self.s3 = Script(csvfiles[2])


class Script:
    """
    Object to represent script data.
    """

    def __init__(self, csvfile):
        """
        Initialize with data from csv file.
        """
        df = pd.read_csv(csvfile, sep=', ', engine='python')
        time = df['time'].values
        step = df['step'].values
        current = df['current'].values
        voltage = df['voltage'].values
        chgAh = df['chgAh'].values
        disAh = df['disAh'].values

        self.time = time
        self.step = step
        self.current = current
        self.voltage = voltage
        self.chgAh = chgAh
        self.disAh = disAh


