"""
Class objects for the dyn_model calculations
"""

# Modules
# ------------------------------------------------------------------------------

import pickle
import pandas as pd

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
        self.OCV0 = OCV0
        self.OCVrel = OCVrel
        self.SOC = SOC
        self.OCV = OCV
        self.SOC0 = SOC0
        self.SOCrel = SOCrel
        self.OCVeta = OCVeta
        self.OCVQ = OCVQ

    @classmethod
    def load(cls, pfile):
        """
        Load attributes from pickle file where pfile is string representing
        path to the pickle file.
        """
        ocv = pickle.load(open(pfile, 'rb'))
        return cls(ocv.OCV0, ocv.OCVrel, ocv.SOC, ocv.OCV, ocv.SOC0, ocv.SOCrel, ocv.OCVeta, ocv.OCVQ)


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


