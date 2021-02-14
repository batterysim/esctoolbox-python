"""
Class objects for the dyn_model calculations
"""

# Modules
# ------------------------------------------------------------------------------

import json
import pandas as pd
import numpy as np

# Class Objects
# ------------------------------------------------------------------------------

class DataModel:
    """
    Data from battery script tests. Requires the Script class which reads the
    csv file and assigns the data to class attributes.
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
        df = pd.read_csv(csvfile)
        time = df['time'].values
        step = df[' step'].values
        current = df[' current'].values
        voltage = df[' voltage'].values
        chgAh = df[' chgAh'].values
        disAh = df[' disAh'].values

        self.time = time
        self.step = step
        self.current = current
        self.voltage = voltage
        self.chgAh = chgAh
        self.disAh = disAh


class ModelOcv:
    """
    Model representing OCV results.
    """
    # pylint: disable=too-many-instance-attributes

    def __init__(self, OCV0, OCVrel, SOC, OCV, SOC0, SOCrel, OCVeta, OCVQ):
        self.OCV0 = np.array(OCV0)
        self.OCVrel = np.array(OCVrel)
        self.SOC = np.array(SOC)
        self.OCV = np.array(OCV)
        self.SOC0 = np.array(SOC0)
        self.SOCrel = np.array(SOCrel)
        self.OCVeta = np.array(OCVeta)
        self.OCVQ = np.array(OCVQ)

    @classmethod
    def load(cls, pfile):
        """
        Load attributes from json file where pfile is string representing
        path to the json file.
        """
        ocv = json.load(open(pfile, 'r'))
        return cls(ocv['OCV0'], ocv['OCVrel'], ocv['SOC'], ocv['OCV'], ocv['SOC0'], ocv['SOCrel'], ocv['OCVeta'], ocv['OCVQ'])


class ModelDyn:
    """
    Model representing results from the dynamic calculations.
    """
    # pylint: disable=too-many-instance-attributes

    def __init__(self):
        self.temps = None
        self.etaParam = None
        self.QParam = None
        self.GParam = None
        self.M0Param = None
        self.MParam = None
        self.R0Param = None
        self.RCParam = None
        self.RParam = None
        self.SOC = None
        self.OCV0 = None
        self.OCVrel = None


