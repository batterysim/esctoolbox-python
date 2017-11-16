"""
Plot the data from the four A123 battery tests conducted at a range of
temperatures. Change the tc variable to view plots from the other test data.
For example, change the tc string to N05 to create plots for the CSV files
named A123_OCV_N05_S1, A123_OCV_N05_S2, A123_OCV_N05_S3, and A123_OCV_N05_S4.
"""

import matplotlib.pyplot as plt
import pandas as pd

# Data
# ------------------------------------------------------------------------------

# string that represents temperature of battery test to determine which csv
# files to read, values can be N05, N15, N25, P05, P15, P25, P35, or P45
tc = 'N05'

df1 = pd.read_csv('../ocv_data/A123_OCV_' + tc + '_S1.csv')
test_time1 = df1['Test_Time(s)']
current1 = df1['Current(A)']
voltage1 = df1['Voltage(V)']

df2 = pd.read_csv('../ocv_data/A123_OCV_' + tc + '_S2.csv')
test_time2 = df2['Test_Time(s)']
current2 = df2['Current(A)']
voltage2 = df2['Voltage(V)']

df3 = pd.read_csv('../ocv_data/A123_OCV_' + tc + '_S3.csv')
test_time3 = df3['Test_Time(s)']
current3 = df3['Current(A)']
voltage3 = df3['Voltage(V)']

df4 = pd.read_csv('../ocv_data/A123_OCV_' + tc + '_S4.csv')
test_time4 = df4['Test_Time(s)']
current4 = df4['Current(A)']
voltage4 = df4['Voltage(V)']

# Compare Temperature Data
# ------------------------------------------------------------------------------

temps = ['N05', 'N15', 'N25', 'P05', 'P15', 'P25', 'P35', 'P45']
times = []
volts = []

for t in temps:
    df = pd.read_csv('../ocv_data/A123_OCV_' + t + '_S1.csv')
    time = df['Test_Time(s)'].values
    voltage = df['Voltage(V)'].values
    times.append(time)
    volts.append(voltage)

# Plot
# ------------------------------------------------------------------------------

plt.ion()
plt.close('all')

# Figure 1
plt.figure(1)
plt.plot(times[0], volts[0], label='-5$^{\circ}$C')
plt.plot(times[1], volts[1], label='-15$^{\circ}$C')
plt.plot(times[2], volts[2], label='-25$^{\circ}$C')
plt.plot(times[3], volts[3], label='5$^{\circ}$C')
plt.plot(times[4], volts[4], label='15$^{\circ}$C')
plt.plot(times[5], volts[5], label='25$^{\circ}$C')
plt.plot(times[6], volts[6], label='35$^{\circ}$C')
plt.plot(times[7], volts[7], label='45$^{\circ}$C')
plt.legend(loc='best')
plt.xlabel('Time (s)')
plt.ylabel('Voltage (V)')
plt.title('A123 battery cell')

# Figure 2
fig, ax1 = plt.subplots()
plt.title('A123_OCV_' + tc + '_S1')

ax1.plot(test_time1, current1, color='b', lw=2, label='current')
ax1.set_xlabel('Test Time (s)')
ax1.set_ylabel('Current (A)', color='b')
ax1.tick_params('y', colors='b')

ax2 = ax1.twinx()
ax2.plot(test_time1, voltage1, color='r', lw=2, label='voltage')
ax2.set_ylabel('Voltage (V)', color='r')
ax2.tick_params('y', colors='r')

# Figure 3
fig, ax1 = plt.subplots()
plt.title('A123_OCV_' + tc + '_S2')

ax1.plot(test_time2, current2, color='b', lw=2, label='current')
ax1.set_xlabel('Test Time (s)')
ax1.set_ylabel('Current (A)', color='b')
ax1.tick_params('y', colors='b')

ax2 = ax1.twinx()
ax2.plot(test_time2, voltage2, color='r', lw=2, label='voltage')
ax2.set_ylabel('Voltage (V)', color='r')
ax2.tick_params('y', colors='r')

# Figure 4
fig, ax1 = plt.subplots()
plt.title('A123_OCV_' + tc + '_S3')

ax1.plot(test_time3, current3, color='b', lw=2, label='current')
ax1.set_xlabel('Test Time (s)')
ax1.set_ylabel('Current (A)', color='b')
ax1.tick_params('y', colors='b')

ax2 = ax1.twinx()
ax2.plot(test_time3, voltage3, color='r', lw=2, label='voltage')
ax2.set_ylabel('Voltage (V)', color='r')
ax2.tick_params('y', colors='r')

# Figure 5
fig, ax1 = plt.subplots()
plt.title('A123_OCV_' + tc + '_S4')

ax1.plot(test_time4, current4, color='b', lw=2, label='current')
ax1.set_xlabel('Test Time (s)')
ax1.set_ylabel('Current (A)', color='b')
ax1.tick_params('y', colors='b')

ax2 = ax1.twinx()
ax2.plot(test_time4, voltage4, color='r', lw=2, label='voltage')
ax2.set_ylabel('Voltage (V)', color='r')
ax2.tick_params('y', colors='r')

