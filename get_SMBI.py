import MDSplus as mds
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

# shot = int(input("shot: "))
# signal = input("signal: ")
whichone = 1
shot = 100000
cn  = mds.Connection('mds.ipp.ac.cn')
cn.openTree('east_1', shot)

ip = cn.get('\ipg')
ip = ip.data()
time_ip = cn.get('dim_of(\ipg)')
idx_begin = 6000
idx_end = np.argwhere(ip > 10)[-1,0] # The return value of argwhere is row * 1

gauges = np.array([['\pjs105', '\pjs103'],['\pjs205', '\pjs203'],['\pas105', '\pas103']])

Pam2P = 4.82e20
volume = 2.0431e-4
which_gauge = 0
if idx_end - idx_begin >= 1+1e5:
    volume = 2.0431e-4 + 3.78e-3
    which_gauge = 1
gauge = gauges[whichone,which_gauge]
pressure = cn.get(gauges[whichone,which_gauge])

if max(pressure)<0.1 and which_gauge == 0:
        volume = 2.0431e-4+3.78e-3
        which_gauge = 1
        gauge = gauges[whichone, which_gauge]
        pressure = cn.get(gauge)

if which_gauge == 0 and max(pressure) > 9.8:
        which_gauge = 1
        gauge = gauges[whichone, which_gauge]
        pressure = cn.get(gauge)

background = np.mean(pressure[:idx_begin])
end_pressure = np.mean(pressure[idx_end:idx_end + 1000])
pressure = pressure[idx_begin:idx_end]

pressure = savgol_filter(pressure, 1000,3)

if which_gauge == 0:
    pressure[pressure < background] = background
    pressure[pressure > end_pressure] = end_pressure
    kp = 2e3
    pressure = pressure - background
else:
    pressure[pressure > background] = background
    pressure[pressure < end_pressure] = end_pressure
    kp = 4e5
    pressure = background - pressure

pressure = pressure * kp * volume * Pam2P
pressure[pressure < 1e18] = 0
data = pressure
plt.plot(data)
plt.show()