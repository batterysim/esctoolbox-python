"""
Models for the OCV calculations
"""

import pandas as pd


class BatteryScript:
    """
    Script or experiment performed on the battery cell.
    """

    def __init__(self, csvdata):
        """
        Initialize the script measurements.
        """
        columns = ['Discharge_Capacity(Ah)', 'Charge_Capacity(Ah)', 'Step_Index', 'Voltage(V)'] 
        df = pd.read_csv(csvdata, usecols=columns)
        self.disAh = df['Discharge_Capacity(Ah)'].values
        self.chgAh = df['Charge_Capacity(Ah)'].values
        self.step = df['Step_Index'].values
        self.voltage = df['Voltage(V)'].values


class BatteryData:
    """
    Object to store battery measurements from script or experiment for a
    certain temperature.
    """

    def __init__(self, csvfiles):
        """
        Initialize with list of CSV data files.
        """
        self.s1 = BatteryScript(csvfiles[0])
        self.s2 = BatteryScript(csvfiles[1])
        self.s3 = BatteryScript(csvfiles[2])
        self.s4 = BatteryScript(csvfiles[3])


class FileData:
    """
    Calculated data from file.
    """

    def __init__(self, disV, disZ, chgV, chgZ, rawocv, temp):
        self.disV = disV
        self.disZ = disZ
        self.chgV = chgV
        self.chgZ = chgZ
        self.rawocv = rawocv
        self.temp = temp


class ModelOcv:
    """
    Model representing OCV results.
    """

    def __init__(self, OCV0, OCVrel, SOC, OCV, SOC0, SOCrel, OCVeta, OCVQ):
        self.OCV0 = OCV0
        self.OCVrel = OCVrel
        self.SOC = SOC
        self.OCV = OCV
        self.SOC0 = SOC0
        self.SOCrel = SOCrel
        self.OCVeta = OCVeta
        self.OCVQ = OCVQ


