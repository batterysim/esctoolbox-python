"""
Functions used by the dyn_model
"""

# Modules
# ------------------------------------------------------------------------------

import ipdb
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fminbound, nnls
from models import ModelDyn

# Functions
# ------------------------------------------------------------------------------

def OCVfromSOCtemp(soc, temp, model):
    """ OCV function """
    SOC = model.SOC          # force to be column vector
    OCV0 = model.OCV0        # force to be column vector
    OCVrel = model.OCVrel    # force to be column vector

    # if soc is scalar then make it a vector
    soccol = np.asarray(soc)
    if soccol.ndim == 0:
        soccol = soccol[None]

    tempcol = temp * np.ones(np.size(soccol))

    diffSOC = SOC[1] - SOC[0]           # spacing between SOC points - assume uniform
    ocv = np.zeros(np.size(soccol))     # initialize output to zero
    I1, = np.where(soccol <= SOC[0])    # indices of socs below model-stored data
    I2, = np.where(soccol >= SOC[-1])   # and of socs above model-stored data
    I3, = np.where((soccol > SOC[0]) & (soccol < SOC[-1]))   # the rest of them
    I6 = np.isnan(soccol)               # if input is "not a number" for any locations

    # for voltages less than lowest stored soc datapoint, extrapolate off
    # low end of table
    if I1.any():
        dv = (OCV0[1] + tempcol*OCVrel[1]) - (OCV0[0] + tempcol*OCVrel[0])
        ocv[I1] = (soccol[I1] - SOC[0])*dv[I1]/diffSOC + OCV0[0] + tempcol[I1]*OCVrel[0]

    # for voltages greater than highest stored soc datapoint, extrapolate off
    # high end of table
    if I2.any():
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
    Identify state-space "A" matrix from input-output data.
    y: vector of measured outputs
    u: vector of measured inputs
    n: number of poles in solution
    A: discrete-time state-space state-transition matrix.

    Theory from "Subspace Identification for Linear Systems Theory - Implementation
    - Applications" Peter Van Overschee / Bart De Moor (VODM) Kluwer Academic
      Publishers, 1996. Combined algorithm: Figure 4.8 page 131 (robust). Robust
      implementation: Figure 6.1 page 169.

    Code adapted from "subid.m" in "Subspace Identification for Linear Systems"
    toolbox on MATLAB CENTRAL file exchange, originally by Peter Van Overschee,
    Dec. 1995
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
    UY = np.concatenate((U, Y))     # combine U and Y into one array
    _, r = np.linalg.qr(UY.T)       # QR decomposition
    R = r.T                         # transpose of upper triangle

    # STEP 1: Calculate oblique and orthogonal projections
    # ------------------------------------------------------------------

    Rf = R[-i:]                                 # future outputs
    Rp = np.concatenate((R[:i], R[2*i:3*i]))    # past inputs and outputs
    Ru = R[i:twoi, :twoi]                       # future inputs

    RfRu = np.linalg.lstsq(Ru.T, Rf[:, :twoi].T, rcond=None)[0].T
    RfRuRu = RfRu.dot(Ru)
    tm1 = Rf[:, :twoi] - RfRuRu
    tm2 = Rf[:, twoi:4*i]
    Rfp = np.concatenate((tm1, tm2), axis=1)    # perpendicular future outputs

    RpRu = np.linalg.lstsq(Ru.T, Rp[:, :twoi].T, rcond=None)[0].T
    RpRuRu = RpRu.dot(Ru)
    tm3 = Rp[:, :twoi] - RpRuRu
    tm4 = Rp[:, twoi:4*i]
    Rpp = np.concatenate((tm3, tm4), axis=1)    # perpendicular past inputs and outputs

    # The oblique projection is computed as (6.1) in VODM, page 166.
    # obl/Ufp = Yf/Ufp * pinv(Wp/Ufp) * (Wp/Ufp)
    # The extra projection on Ufp (Uf perpendicular) tends to give
    # better numerical conditioning (see algo on VODM page 131)

    # Funny rank check (SVD takes too long)
    # This check is needed to avoid rank deficiency warnings

    nmRpp = np.linalg.norm(Rpp[:, 3*i-3:-i], ord='fro')
    if nmRpp < 1e-10:
        # oblique projection as (Rfp*pinv(Rpp')') * Rp
        Ob = Rfp.dot(np.linalg.pinv(Rpp.T).T).dot(Rp)
    else:
        # oblique projection as (Rfp/Rpp) * Rp
        Ob = (np.linalg.lstsq(Rpp.T, Rfp.T, rcond=None)[0].T).dot(Rp)

    # STEP 2: Compute weighted oblique projection and its SVD
    #         Extra projection of Ob on Uf perpendicular
    # ------------------------------------------------------------------

    ObRu = np.linalg.lstsq(Ru.T, Ob[:, :twoi].T, rcond=None)[0].T
    ObRuRu = ObRu.dot(Ru)
    tm5 = Ob[:, :twoi] - ObRuRu
    tm6 = Ob[:, twoi:4*i]
    WOW = np.concatenate((tm5, tm6), axis=1)

    U, S, _ = np.linalg.svd(WOW, full_matrices=False)
    ss = np.diag(S)

    # STEP 3: Partitioning U into U1 and U2 (the latter is not used)
    # ------------------------------------------------------------------

    U1 = U[:, :n]       # determine U1

    # STEP 4: Determine gam = Gamma(i) and gamm = Gamma(i-1)
    # ------------------------------------------------------------------

    gam = U1 * np.diag(np.sqrt(ss[:n]))
    gamm = gam[0, i-2]
    gam_inv = np.linalg.pinv(gam)[0]            # pseudo inverse of gam
    gamm2 = np.array([[gamm], [gamm]])
    gamm_inv = np.linalg.pinv(gamm2)[0][0]*2    # pseudo inverse of gamm

    # STEP 5: Determine A matrix (also C, which is not used)
    # ------------------------------------------------------------------

    tm7 = np.concatenate((gam_inv.dot(R[-i:, :-i]), np.zeros(n)))
    tm8 = R[i:twoi, 0:3*i+1]
    Rhs = np.vstack((tm7, tm8))
    Lhs = np.vstack((gamm_inv*R[-i+1, :-i+1], R[-i, :-i+1]))
    sol = np.linalg.lstsq(Rhs.T, Lhs.T, rcond=None)[0].T    # solve least squares for [A; C]
    A = sol[n-1, n-1]                           # extract A

    return A


def minfn(data, model, theTemp):
    """
    Using an assumed value for gamma (already stored in the model), find optimum
    values for remaining cell parameters, and compute the RMS error between true
    and predicted cell voltage
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

    # Second modeling step: Compute time constants in "A" matrix
    y = -np.diff(verr)
    u = np.diff(etaik)
    A = SISOsubid(y, u, numpoles)
    # eigA = A

    RCfact = A
    RC = -1/np.log(A)

    # Simulate the R-C filters to find R-C currents
    # assume no control-system toolbox
    vrcRaw = np.zeros(len(etaik))
    for vrcK in range(len(etaik)-1):
        vrcRaw[vrcK+1] = RCfact*vrcRaw[vrcK] + (1-RCfact)*etaik[vrcK]

    # Third modeling step: Hysteresis parameters
    H = np.vstack((hh, sik, -etaik, -vrcRaw)).T
    W = nnls(H, verr)
    M = W[0][0]
    M0 = W[0][1]
    R0 = W[0][2]
    Rfact = W[0][3]

    idx, = np.where(np.array(model.temps) == data[ind].temp)[0]
    model.R0Param[idx] = R0
    model.M0Param[idx] = M0
    model.MParam[idx] = M
    model.RCParam[idx] = RC
    model.RParam[idx] = Rfact

    vest2 = vest1 + M*hh + M0*sik - R0*etaik - vrcRaw*Rfact
    verr = vk - vest2

    # plot voltages
    plt.figure(1)
    plt.plot(tk[::10]/60, vk[::10], label='voltage')
    plt.plot(tk[::10]/60, vest1[::10], label='vest1 (OCV)')
    plt.plot(tk[::10]/60, vest2[::10], label='vest2 (DYN)')
    plt.xlabel('Time (min)')
    plt.ylabel('Voltage (V)')
    plt.title(f'Voltage and estimates at T = {data[ind].temp} C')
    plt.legend(loc='best', numpoints=1)
    #plt.show()

    # plot modeling errors
    plt.figure(2)
    plt.plot(tk[::10]/60, verr[::10], label='verr')
    plt.xlabel('Time (min)')
    plt.ylabel('Error (V)')
    plt.title(f'Modeling error at T = {data[ind].temp} C')
    #plt.show()

    # Compute RMS error only on data roughly in 5% to 95% SOC
    v1 = OCVfromSOCtemp(0.95, data[ind].temp, model)[0]
    v2 = OCVfromSOCtemp(0.05, data[ind].temp, model)[0]
    N1 = np.where(vk < v1)[0][0]
    N2 = np.where(vk < v2)[0][0]

    rmserr = np.sqrt(np.mean(verr[N1:N2]**2))
    cost = np.sum(rmserr)
    print(f'RMS error = {cost*1000:.2f} mV')

    return cost, model


def optfn(x, data, model, theTemp):
    """
    This minfn works for the enhanced self-correcting cell model
    """

    idx, = np.where(np.array(model.temps) == theTemp)
    model.GParam[idx] = abs(x)

    cost, _ = minfn(data, model, theTemp)
    return cost


def processDynamic(data, modelocv, numpoles, doHyst):
    """
    Technical note: PROCESSDYNAMIC assumes that specific Arbin test scripts have
    been executed to generate the input files.  "makeMATfiles.m" converts the raw
    Excel data files into "MAT" format where the MAT files have fields for time,
    step, current, voltage, chgAh, and disAh for each script run.

    The results from three scripts are required at every temperature.
    The steps in each script file are assumed to be:
    Script 1 (thermal chamber set to test temperature):
        Step 1: Rest @ 100% SOC to acclimatize to test temperature
        Step 2: Discharge @ 1C to reach ca. 90% SOC
        Step 3: Repeatedly execute dynamic profiles (and possibly intermediate
        rests) until SOC is around 10%
    Script 2 (thermal chamber set to 25 degC):
        Step 1: Rest ca. 10% SOC to acclimatize to 25 degC
        Step 2: Discharge to min voltage (ca. C/3)
        Step 3: Rest
        Step 4: Constant voltage at vmin until current small (ca. C/30)
        Steps 5-7: Dither around vmin
        Step 8: Rest
    Script 3 (thermal chamber set to 25 degC):
        Step 2: Charge @ 1C to max voltage
        Step 3: Rest
        Step 4: Constant voltage at vmax until current small (ca. C/30)
        Steps 5-7: Dither around vmax
        Step 8: Rest

    All other steps (if present) are ignored by PROCESSDYNAMIC. The time step
    between data samples must be uniform -- we assume a 1s sample period in this
    code.

    The inputs:
    - data: An array, with one entry per temperature to be processed.
        One of the array entries must be at 25 degC. The fields of "data" are:
        temp (the test temperature), script1, script 2, and script 3, where the
        latter comprise data collected from each script.  The sub-fields of
        these script structures that are used by PROCESSDYNAMIC are the
        vectors: current, voltage, chgAh, and disAh
    - model: The output from processOCV, comprising the OCV model
    - numpoles: The number of R-C pairs in the model
    - doHyst: 0 if no hysteresis model desired; 1 if hysteresis desired

    The output:
    - model: A modified model, which now contains the dynamic fields filled in.
    """

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

    modeldyn = ModelDyn()
    modeldyn.temps = alltemps
    modeldyn.etaParam = alletas
    modeldyn.QParam = allQs

    # Step 2: Compute OCV for "discharge portion" of test
    # ------------------------------------------------------------------

    for k, _ in enumerate(data):
        etaParam = modeldyn.etaParam[k]
        etaik = data[k].s1.current.copy()
        etaik[etaik < 0] = etaParam*etaik[etaik < 0]
        data[k].Z = 1 - np.cumsum(etaik) * 1/(data[k].Q * 3600)
        data[k].OCV = OCVfromSOCtemp(data[k].Z, alltemps[k], modelocv)

    # Step 3: Now, optimize!
    # ------------------------------------------------------------------

    modeldyn.GParam = np.zeros(len(modeldyn.temps))   # gamma hysteresis parameter
    modeldyn.M0Param = np.zeros(len(modeldyn.temps))  # M0 hysteresis parameter
    modeldyn.MParam = np.zeros(len(modeldyn.temps))   # M hysteresis parameter
    modeldyn.R0Param = np.zeros(len(modeldyn.temps))  # R0 ohmic resistance parameter
    modeldyn.RCParam = np.zeros(len(modeldyn.temps))  # time constant
    modeldyn.RParam = np.zeros(len(modeldyn.temps))   # Rk

    modeldyn.SOC = modelocv.SOC        # copy SOC values from OCV model
    modeldyn.OCV0 = modelocv.OCV0      # copy OCV0 values from OCV model
    modeldyn.OCVrel = modelocv.OCVrel  # copy OCVrel values from OCV model

    for theTemp in modeldyn.temps:
        print('Processing temperature', theTemp, 'C')

        if doHyst:
            g = abs(fminbound(optfn, 1, 250, args=(data, modeldyn, theTemp)))
            print('g =', g)

        else:
            g = 0 
    return modeldyn   
