"""
Python version of the runProcessOCV Matlab file for A123_OCV battery cell.
"""

import matplotlib.pyplot as plt
import numpy as np

from models import BatteryScript, BatteryData, FileData, ModelOcv
from funcs import OCVfromSOCtemp

# Parameters and Data
# ------------------------------------------------------------------------------

# temperatures for cell experiments
temps = np.array([-25, -15, -5, 5, 15, 25, 35, 45])

minV = 2.00     # minimum cell voltage, used for plotting results
maxV = 3.75     # maximum cell voltage, used for plotting results

SOC = np.arange(0, 1+0.005, 0.005)  # range for state of charge

# initialize variables to store calculations
eta = np.zeros(len(temps))  # coulombic efficiency
Q = np.zeros(len(temps))    # apparent total capacity

# initialize array to store battery cell data
data = np.zeros(len(temps), dtype=object)     

# load battery cell data for each temperature as objects then store in data array
for idx, temp in enumerate(temps):
    if temp < 0:
        tempfmt = f'{abs(temp):02}'
        files = [f'../a123_data/A123_OCV_N{tempfmt}_S1.csv', f'../a123_data/A123_OCV_N{tempfmt}_S2.csv',
                 f'../a123_data/A123_OCV_N{tempfmt}_S3.csv', f'../a123_data/A123_OCV_N{tempfmt}_S4.csv']
        data[idx] = BatteryData(files)
    else:
        tempfmt = f'{abs(temp):02}'
        files = [f'../a123_data/A123_OCV_P{tempfmt}_S1.csv', f'../a123_data/A123_OCV_P{tempfmt}_S2.csv',
                 f'../a123_data/A123_OCV_P{tempfmt}_S3.csv', f'../a123_data/A123_OCV_P{tempfmt}_S4.csv']
        data[idx] = BatteryData(files)

# initial array to store calculated data
filedata = np.zeros(len(temps), dtype=object)

# Process 25 degC data to find raw OCV relationship and eta25
# ------------------------------------------------------------------------------

k, = np.where(temps == 25)[0]   # index where temperature is 25 degC
p25 = data[k]

# compute total discharge in ampere hours, Ah
totDisAh = p25.s1.disAh[-1] + p25.s2.disAh[-1] + p25.s3.disAh[-1] + p25.s4.disAh[-1]

# compute total charge in ampere hours, Ah
totChgAh = p25.s1.chgAh[-1] + p25.s2.chgAh[-1] + p25.s3.chgAh[-1] + p25.s4.chgAh[-1]

# the 25 degC coulombic efficiency
eta25 = totDisAh/totChgAh
eta[k] = eta25

# adjust charge Ah in all scripts per eta25
p25.s1.chgAh = p25.s1.chgAh * eta25
p25.s2.chgAh = p25.s2.chgAh * eta25
p25.s3.chgAh = p25.s3.chgAh * eta25
p25.s4.chgAh = p25.s4.chgAh * eta25

# compute cell capacity at 25 degC, should be essentially same at
# all temps, but we're computing them individually to check this
Q25 = p25.s1.disAh[-1] + p25.s2.disAh[-1] - p25.s1.chgAh[-1] - p25.s2.chgAh[-1]
Q[k] = Q25

# discharge
indD = np.where(p25.s1.step == 2)[0]                            # slow discharge step
IR1Da = p25.s1.voltage[indD[0]-1] - p25.s1.voltage[indD[0]]     # the i*R voltage drop at beginning of discharge
IR2Da = p25.s1.voltage[indD[-1]+1] - p25.s1.voltage[indD[-1]]   # the i*R voltage drop at end of discharge

# charge
indC = np.where(p25.s3.step == 2)[0]                            # slow charge step
IR1Ca = p25.s3.voltage[indC[0]] - p25.s3.voltage[indC[0]-1]     # the i*R voltage rise at beginning of charge
IR2Ca = p25.s3.voltage[indC[-1]] - p25.s3.voltage[indC[-1]+1]   # the i*R voltage rise at end of charge

# put bounds on R
IR1D = min(IR1Da, 2*IR2Ca)
IR2D = min(IR2Da, 2*IR1Ca)
IR1C = min(IR1Ca, 2*IR2Da)
IR2C = min(IR2Ca, 2*IR1Da)

# discharge
blendD = np.linspace(0, 1, len(indD))   # linear blending from 0 to 1 for discharge
IRblendD = IR1D + (IR2D - IR1D)*blendD  # blend resistances for discharge
disV = p25.s1.voltage[indD] + IRblendD  # approximate discharge voltage at each point
disZ = 1 - p25.s1.disAh[indD]/Q25       # approximate SOC at each point
disZ = disZ + (1 - disZ[0])

# charge
blendC = np.linspace(0, 1, len(indC))   # linear blending from 0 to 1 for charge
IRblendC = IR1C + (IR2C - IR1C)*blendC  # blend resistances for charge
chgV = p25.s3.voltage[indC] - IRblendC  # approximate charge voltage at each point
chgZ = p25.s3.chgAh[indC]/Q25           # approximate SOC at each point
chgZ = chgZ - chgZ[0]

# compute voltage difference between charge and discharge at 50% SOC force i*R
# compensated curve to pass half-way between each charge and discharge at this
# point notice that vector chgZ and disZ must be increasing
deltaV50 = np.interp(0.5, chgZ, chgV) - np.interp(0.5, disZ[::-1], disV[::-1])
ind = np.where(chgZ < 0.5)[0]
vChg = chgV[ind] - chgZ[ind]*deltaV50
zChg = chgZ[ind]
ind = np.where(disZ > 0.5)[0]
vDis = disV[ind] + (1 - disZ[ind])*deltaV50
zDis = disZ[ind]

# rawocv now has our best guess of true ocv at this temperature
rawocv = np.interp(SOC, np.concatenate([zChg, zDis[::-1]]), np.concatenate([vChg, vDis[::-1]]))

# store calculated data into filedata object
filedata[k] = FileData(p25.s1.voltage[indD], disZ, p25.s3.voltage[indC], chgZ, rawocv, temps[k])

# Process Other Temperatures to Find Raw OCV Relationship and Eta
# Everything that follows is same as at 25 degC, except we need to compensate
# for different coulombic efficiencies eta at different temperatures.
# ------------------------------------------------------------------------------

not25, = np.where(temps != 25)

for k in not25:
    # adjust charge Ah  per eta25
    data[k].s2.chgAh = data[k].s2.chgAh * eta25
    data[k].s4.chgAh = data[k].s4.chgAh * eta25

    # coulombic efficiency
    eta[k] = ((data[k].s1.disAh[-1] + data[k].s2.disAh[-1] + data[k].s3.disAh[-1] 
             + data[k].s4.disAh[-1] - data[k].s2.chgAh[-1] - data[k].s4.chgAh[-1]) 
             / (data[k].s1.chgAh[-1] + data[k].s3.chgAh[-1]))

    # adjust charge Ah per eta at current temp
    data[k].s1.chgAh = data[k].s1.chgAh * eta[k]
    data[k].s3.chgAh = data[k].s3.chgAh * eta[k]

    # compute cell capacity
    Q[k] = data[k].s1.disAh[-1] + data[k].s2.disAh[-1] - data[k].s1.chgAh[-1] - data[k].s2.chgAh[-1]

    # discharge
    indD = np.where(data[k].s1.step == 2)[0]                                # slow discharge step
    IR1D = data[k].s1.voltage[indD[0]-1] - data[k].s1.voltage[indD[0]]      # the i*R voltage drop at beginning of discharge
    IR2D = data[k].s1.voltage[indD[-1]+1] - data[k].s1.voltage[indD[-1]]    # the i*R voltage drop at end of discharge

    # charge
    indC = np.where(data[k].s3.step == 2)[0]                             # slow charge step
    IR1C = data[k].s3.voltage[indC[0]] - data[k].s3.voltage[indC[0]-1]   # the i*R voltage rise at beginning of charge
    IR2C = data[k].s3.voltage[indC[-1]] - data[k].s3.voltage[indC[-1]+1] # the i*R voltage rise at end of charge

    # put bounds on R
    IR1D = min(IR1D, 2*IR2C)
    IR2D = min(IR2D, 2*IR1C)
    IR1C = min(IR1C, 2*IR2D)
    IR2C = min(IR2C, 2*IR1D)

    # discharge
    blend = np.linspace(0, 1, len(indD))        # linear blending from 0 to 1 for discharge
    IRblend = IR1D + (IR2D - IR1D)*blend        # blend resistances for discharge
    disV = data[k].s1.voltage[indD] + IRblend   # approximate discharge voltage at each point
    disZ = 1 - data[k].s1.disAh[indD]/Q25       # approximate SOC at each point
    disZ = disZ + (1 - disZ[0])

    # charge
    blend = np.linspace(0, 1, len(indC))        # linear blending from 0 to 1 for charge
    IRblend = IR1C + (IR2C - IR1C)*blend        # blend resistances for charge
    chgV = data[k].s3.voltage[indC] - IRblend   # approximate charge voltage at each point
    chgZ = data[k].s3.chgAh[indC]/Q25           # approximate SOC at each point
    chgZ = chgZ - chgZ[0]

    # compute voltage difference between charge and discharge at 50% SOC force i*R
    # compensated curve to pass half-way between each charge and discharge at this
    # point notice that vector chgZ and disZ must be increasing
    deltaV50 = np.interp(0.5, chgZ, chgV) - np.interp(0.5, disZ[::-1], disV[::-1])
    ind = np.where(chgZ < 0.5)[0]
    vChg = chgV[ind] - chgZ[ind]*deltaV50
    zChg = chgZ[ind]
    ind = np.where(disZ > 0.5)[0]
    vDis = disV[ind] + (1 - disZ[ind])*deltaV50
    zDis = disZ[ind]

    # rawocv now has our best guess of true ocv at this temperature
    rawocv = np.interp(SOC, np.concatenate([zChg, zDis[::-1]]), np.concatenate([vChg, vDis[::-1]]))

    # store calculated data into filedata object
    filedata[k] = FileData(data[k].s1.voltage[indD], disZ, data[k].s3.voltage[indC], chgZ, rawocv, temps[k])

# Use the SOC versus OCV data now available at each individual
# temperature to compute an OCV0 and OCVrel relationship
# ------------------------------------------------------------------------------

# compile the voltages and temperatures into single arrays rather than structures
postemps = temps[temps > 0]     # temps > 0
numtempskept = len(postemps)    # number of temps > 0

nocv = len(filedata[5].rawocv)          # number of rawocv values based on 25 degC results
Vraw = np.zeros([numtempskept, nocv])   # initialize rawocv array
idxpos = np.where(temps > 0)[0]         # indices of positive file temperatures

for k in range(numtempskept):
    Vraw[k] = filedata[idxpos[k]].rawocv

# use linear least squares to determine best guess for OCV at 0 degC
# and then the per-degree OCV change
OCV0 = np.zeros(len(SOC))
OCVrel = np.zeros(len(SOC))
H = np.ones([numtempskept, 2])
H[:, 1] = postemps

for k in range(len(SOC)):
    X = np.linalg.lstsq(H, Vraw[:, k])
    OCV0[k] = X[0][0]
    OCVrel[k] = X[0][1]

modelocv = ModelOcv(OCV0, OCVrel, SOC, 0, 0, 0, 0, 0)

# Make SOC0 and SOCrel
# Do same kind of analysis to find soc as a function of ocv
# ------------------------------------------------------------------------------

z = np.arange(-0.1, 1.1, 0.01)     # test soc vector
v = np.arange(minV-0.01, maxV+0.02, 0.01)
socs = np.zeros((len(temps), len(v)))

for k in range(len(temps)):
    T = temps[k]
    v1 = OCVfromSOCtemp(z, T, modelocv)
    socs[k, :] = np.interp(v, v1, z)

SOC0 = np.zeros(len(v))
SOCrel = SOC0
H = np.ones([len(temps), 2])
H[:, 1] = temps

for k in range(len(v)):
    X = np.linalg.lstsq(H, socs[:, k])  # fit SOC(v,T) = 1*SOC0(v) + T*SOCrel(v)
    SOC0[k] = X[0][0]
    SOCrel[k] = X[0][1]

# store ocv results in model object
modelocv = ModelOcv(OCV0, OCVrel, SOC, v, SOC0, SOCrel, eta, Q)

# Plot Results
# ------------------------------------------------------------------------------

plt.close('all')
plt.ion()

# k = 0   # index for -25 degC
# err = filedata[k].rawocv - OCVfromSOCtemp(SOC, filedata[k].temp, modelocv)
# rmserr = np.sqrt(np.mean(err**2))

# plt.figure(1)
# plt.plot(100*SOC, OCVfromSOCtemp(SOC, filedata[k].temp, modelocv), 'k')
# plt.plot(100*SOC, filedata[k].rawocv, 'r')
# plt.plot(100*filedata[k].disZ, filedata[k].disV, 'g--')
# plt.plot(100*filedata[k].chgZ, filedata[k].chgV, 'b--')
# plt.text(2, maxV-0.15, f'RMS error = {rmserr*1000:.01f} mV')
# plt.ylim(minV-0.2, maxV+0.2)
# plt.title('A123 OCV relationship at temp = -25')
# plt.xlabel('SOC (%)')
# plt.ylabel('OCV (V)')
# plt.grid()

# k = 5   # index for 25 degC
# err = filedata[k].rawocv - OCVfromSOCtemp(SOC, filedata[k].temp, modelocv)
# rmserr = np.sqrt(np.mean(err**2))

# plt.figure(5)
# plt.plot(100*SOC, OCVfromSOCtemp(SOC, filedata[k].temp, modelocv), 'k')
# plt.plot(100*SOC,filedata[k].rawocv, 'r')
# plt.plot(100*filedata[k].disZ, filedata[k].disV, 'g--')
# plt.plot(100*filedata[k].chgZ, filedata[k].chgV, 'b--')
# plt.text(2, maxV-0.15, f'RMS error = {rmserr*1000:.01f} mV')
# plt.ylim(minV-0.2, maxV+0.2)
# plt.title('A123 OCV relationship at temp = 25')
# plt.xlabel('SOC (%)')
# plt.ylabel('OCV (V)')
# plt.grid()

for k in range(len(temps)):
    err = filedata[k].rawocv - OCVfromSOCtemp(SOC, filedata[k].temp, modelocv)
    rmserr = np.sqrt(np.mean(err**2))

    plt.figure(k+1)
    plt.plot(100*SOC, OCVfromSOCtemp(SOC, filedata[k].temp, modelocv), 'k', label='model')
    plt.plot(100*SOC,filedata[k].rawocv, 'r', label='approx')
    plt.plot(100*filedata[k].disZ, filedata[k].disV, 'g--', label='dis')
    plt.plot(100*filedata[k].chgZ, filedata[k].chgV, 'b--', label='chg')
    plt.text(2, maxV-0.15, f'RMS error = {rmserr*1000:.01f} mV')
    plt.ylim(minV-0.2, maxV+0.2)
    plt.title(f'A123 OCV relationship at temp = {temps[k]}')
    plt.xlabel('SOC (%)')
    plt.ylabel('OCV (V)')
    plt.legend(numpoints=1, loc='lower right')
    plt.grid()

