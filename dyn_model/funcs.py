"""
Functions for DYN model
"""

import ipdb
import numpy as np
from scipy.optimize import fminbound


def OCVfromSOCtemp(soc, temp, model):
    """ OCV function """
    soccol = soc             # force soc to be column vector
    SOC = model.SOC          # force to be column vector
    OCV0 = model.OCV0        # force to be column vector
    OCVrel = model.OCVrel    # force to be column vector

    tempcol = temp * np.ones(len(soccol))

    diffSOC = SOC[1] - SOC[0]           # spacing between SOC points - assume uniform
    ocv = np.zeros(len(soccol))         # initialize output to zero
    I1, = np.where(soccol <= SOC[0])    # indices of socs below model-stored data
    I2, = np.where(soccol >= SOC[-1])   # and of socs above model-stored data
    I3, = np.where((soccol > SOC[0]) & (soccol < SOC[-1]))   # the rest of them
    I6 = np.isnan(soccol)               # if input is "not a number" for any locations

    # for voltages less than lowest stored soc datapoint, extrapolate off
    # low end of table
    if len(I1) != 0:
        dv = (OCV0[1] + tempcol*OCVrel[1]) - (OCV0[0] + tempcol*OCVrel[0])
        ocv[I1] = (soccol[I1] - SOC[0])*dv[I1]/diffSOC + OCV0[0] + tempcol[I1]*OCVrel[0]

    # for voltages greater than highest stored soc datapoint, extrapolate off
    # high end of table
    if len(I2) != 0:
        dv = (OCV0[-1] + tempcol*OCVrel[-1]) - (OCV0[-2] + tempcol*OCVrel[-2])
        ocv[I2] = (soccol[I2] - SOC[-1])*dv[I2]/diffSOC + OCV0[-1] + tempcol[I2]*OCVrel[-1]

    # for normal soc range, manually interpolate (10x faster than "interp1")
    I4 = (soccol[I3] - SOC[0])/diffSOC  # using linear interpolation
    I5 = np.floor(I4)
    I5 = I5.astype(int)
    I45 = I4 - I5
    omI45 = 1 - I45
    ocv[I3] = OCV0[I5]*omI45 + OCV0[I5+1]*I45
    ocv[I3] = ocv[I3] + tempcol[I3]*(OCVrel[I5]*omI45 + OCVrel[I5+1]*I45)
    ocv[I6] = 0     # replace NaN SOCs with zero voltage
    return ocv


def SISOsubid(y, u, n):
    """
    siso function
    """
    ny = len(y)
    i = 2*n
    twoi = 4*n

    # Determine the number of columns in the Hankel matrices
    j = ny - twoi + 1

    # Make Hankel matrices Y and U
    Y = np.zeros((twoi, j))
    U = np.zeros((twoi, j))

    for k in range(2*i):
        Y[k] = y[k:k+j]
        U[k] = u[k:k+j]

    # Compute the R factor
    # TODO continue line 353 in processDynamic.m

    ipdb.set_trace()

    return A


def optfn(x, data, model, theTemp):
    """
    optfn function
    """

    idx, = np.where(np.array(model.temps) == theTemp)
    model.GParam[idx] = abs(x)

    cost, _ = minfn(data, model, theTemp)
    return cost


def minfn(data, model, theTemp):
    """
    minfn function
    """

    alltemps = [d.temp for d in data]
    ind, = np.where(np.array(alltemps) == theTemp)[0]

    G = abs(model.GParam[ind])

    Q = abs(model.QParam[ind])
    eta = abs(model.etaParam[ind])
    RC = abs(model.RCParam[ind])
    numpoles = 1

    ik = data[ind].s1.current.copy()
    vk = data[ind].s1.voltage.copy()
    tk = np.arange(len(vk))
    etaik = ik.copy()
    etaik[ik < 0] = etaik[ik < 0] * eta

    hh = 0*ik
    sik = 0*ik
    fac = np.exp(-abs(G * etaik/(3600*Q)))

    for k in range(1, len(ik)):
        hh[k] = (fac[k-1]*hh[k-1]) - ((1-fac[k-1])*np.sign(ik[k-1]))
        sik[k] = np.sign(ik[k])
        if abs(ik[k]) < Q/100:
            sik[k] = sik[k-1]

    # First modeling step: Compute error with model = OCV only
    vest1 = data[ind].OCV
    verr = vk - vest1

    ipdb.set_trace()
    # Second modeling step: Compute time constants in "A" matrix
    A = SISOsubid(-np.diff(verr), np.diff(etaik), numpoles)
    
    return cost, model


def processDynamic(data, modelocv, numpoles, doHyst):
    """
    DYN function
    """

    class Model:
        """ Model containing results from this function """
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


    # Step 1: Compute capacity and coulombic efficiency for every test
    # ------------------------------------------------------------------

    alltemps = [d.temp for d in data]
    alletas = np.zeros(len(alltemps))
    allQs = np.zeros(len(alltemps))

    ind25, = np.where(np.array(alltemps) == 25)[0]
    not25, = np.where(np.array(alltemps) != 25)

    k = ind25

    totDisAh = data[k].s1.disAh[-1] + data[k].s2.disAh[-1] + data[k].s3.disAh[-1]
    totChgAh = data[k].s1.chgAh[-1] + data[k].s2.chgAh[-1] + data[k].s3.chgAh[-1]
    eta25 = totDisAh/totChgAh
    data[k].eta = eta25
    alletas[k] = eta25
    data[k].s1.chgAh = data[k].s1.chgAh * eta25
    data[k].s2.chgAh = data[k].s2.chgAh * eta25
    data[k].s3.chgAh = data[k].s3.chgAh * eta25

    Q25 = data[k].s1.disAh[-1] + data[k].s2.disAh[-1] - data[k].s1.chgAh[-1] - data[k].s2.chgAh[-1]
    data[k].Q = Q25
    allQs[k] = Q25

    eta25 = np.mean(alletas[ind25])

    for k in not25:
        data[k].s2.chgAh = data[k].s2.chgAh*eta25
        data[k].s3.chgAh = data[k].s3.chgAh*eta25
        eta = (data[k].s1.disAh[-1] + data[k].s2.disAh[-1] + data[k].s3.disAh[-1] - data[k].s2.chgAh[-1] - data[k].s3.chgAh[-1])/data[k].s1.chgAh[-1]

        data[k].s1.chgAh = eta*data[k].s1.chgAh
        data[k].eta = eta
        alletas[k] = eta

        Q = data[k].s1.disAh[-1] + data[k].s2.disAh[-1] - data[k].s1.chgAh[-1] - data[k].s2.chgAh[-1]
        data[k].Q = Q
        allQs[k] = Q

    model = Model()
    model.temps = alltemps
    model.etaParam = alletas
    model.QParam = allQs

    # Step 2: Compute OCV for "discharge portion" of test
    # ------------------------------------------------------------------

    for k, _ in enumerate(data):
        etaParam = model.etaParam[k]
        etaik = data[k].s1.current.copy()
        etaik[etaik < 0] = etaParam*etaik[etaik < 0]
        data[k].Z = 1 - np.cumsum(etaik) * 1/(data[k].Q * 3600)
        data[k].OCV = OCVfromSOCtemp(data[k].Z, alltemps[k], modelocv)

    # Step 3: Now, optimize!)
    # ------------------------------------------------------------------

    model.GParam = np.zeros(len(model.temps))   # gamma hysteresis parameter
    model.M0Param = np.zeros(len(model.temps))  # M0 hysteresis parameter
    model.MParam = np.zeros(len(model.temps))   # M hysteresis parameter
    model.R0Param = np.zeros(len(model.temps))  # R0 ohmic resistance parameter
    model.RCParam = np.zeros(len(model.temps))  # time constant
    model.RParam = np.zeros(len(model.temps))   # Rk

    tempidx = 0
    theTemp = model.temps[tempidx]
    print('Processing temperature', theTemp, 'C')

    g = abs(fminbound(optfn, 1, 250, args=(data, model, theTemp)))


