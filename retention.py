import MDSplus as mds
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d

cn = mds.Connection('mds.ipp.ac.cn')

def get_retention(shot,rate=False):
    """shot:int, ratio = True when get retention rate"""
    inject = get_inj(shot)
    exhaust = get_exh(shot)
    if not rate:
        exhaust = np.cumsum(exhaust) / 1000
        res = inject - exhaust
    else:
        inject = np.diff(inject) * 1000
        inject = np.append([0],inject)
        res =  savgol_filter(inject-exhaust, 1000, 3)
    return res
def get_inj(shot):

    puffing = get_puff(shot)

    smbi = get_SMBI(shot)

    pi = get_pi(shot)

    n_plasma = get_plasma(shot) # ??


    return puffing + smbi +  pi


def get_exh(shot):
    pumping = get_pump(shot)
    N_NBI = get_NBI(shot)
    return  pumping - N_NBI

def get_time(shot):
    """
    Get the time which will be used to uniform the time of all time of different signal
    return: idx_begin : the idx of time when discharge begins, idx_end: the idx of time when discharge ends
    time_ip: the time axis of ip
    """

    cn.openTree('east_1', shot)
    ip = cn.get('\ipg').data()
    time_ip = cn.get('dim_of(\ipg)').data()

    idx_begin = 6000
    idx_end = np.argwhere(ip > 10)[-1, 0]
    cn.closeTree('east_1', shot)
    return (idx_begin,idx_end),time_ip

def get_puff(shot):
    def get_puff_signal(shot, signal):

        (idx_begin, idx_end), time_ip = get_time(shot)  # There's closeTree() inside get_time func, so put it before openTree func

        cn.openTree('east_1', shot)
        pressure = cn.get(f'\\{signal}').data()
        cn.closeTree('east_1', shot)

        if len(time_ip) <= 7000 and max(pressure) - min(pressure) <= 0.3:
            return 0

        background = np.mean(pressure[:idx_begin])
        end_pressure = np.mean(pressure[idx_end:idx_end + 1000])
        pressure = pressure[idx_begin:idx_end]
        pressure = savgol_filter(pressure, 1000, 3)
        pressure[pressure > background] = background
        pressure[pressure < end_pressure] = end_pressure
        pressure = background - pressure

        Pam2P = 4.82e20  # 每pa立方米的粒子数
        kp = 2.5e4  # 充气规管的转换系数

        signal_volumes = {
            "JHG1": 3.118e-4, "JHG2": 2.91e-4, "JHG3": 2.922e-4, "JHG4": 2.99e-4, "JHG5": 2.997e-4,
            "JHG6": 2.949e-4,
            "OUG1": 1.4687e-3, "ODG1": 1.4601e-3, "CDG1": 6.792e-3,
            "HDG1": 6.75e-3, "KHG1": 2.938e-4, "DHG1": 2.919e-4}
        volume = signal_volumes[signal]
        pressure *= kp * volume * Pam2P
        pressure[pressure < 1e18] = 0
        return pressure

    signals = ["JHG1", "JHG4", "JHG5", "CDG1", "DHG1"]
    puffing = 0
    for signal in signals:
        puffing += get_puff_signal(shot, signal)
    return puffing

def get_SMBI(shot):
    def get_nbi(shot , which_one):
        (idx_begin,idx_end), _ = get_time(shot)

        cn.openTree('east_1',shot)
        gauges = np.array([['\pjs105', '\pjs103'], ['\pjs205', '\pjs203'], ['\pas105', '\pas103']])
        Pam2P = 4.82e20
        volume = 2.0431e-4
        which_gauge = 0
        if idx_end - idx_begin >= 1 + 1e5:
            volume = 2.0431e-4 + 3.78e-3
            which_gauge = 1
        gauge = gauges[which_one, which_gauge]
        pressure = cn.get(gauges[which_one, which_gauge]).data()

        if max(pressure) < 0.1 and which_gauge == 0:
            volume = 2.0431e-4 + 3.78e-3
            which_gauge = 1
            gauge = gauges[which_one, which_gauge]
            pressure = cn.get(gauge).data()

        if which_gauge == 0 and max(pressure) > 9.8:
            which_gauge = 1
            gauge = gauges[which_one, which_gauge]
            pressure = cn.get(gauge).data()

        background = np.mean(pressure[:idx_begin])
        end_pressure = np.mean(pressure[idx_end:idx_end + 1000])
        pressure = pressure[idx_begin:idx_end]

        pressure = savgol_filter(pressure, 1000, 3)

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
        cn.closeTree('east_1',shot)
        return pressure
    # compute smbi2 and smbi3
    return get_nbi(shot,1) + get_nbi(shot,2)

def get_NBI(shot):
    kp = 1e19 / 1.6

    (idx_begin, idx_end), _ = get_time(shot)
    cn.openTree('east_1', shot)
    NBI1LHI = 0.8 * (cn.get('\\NBI1LHI')).data()
    NBI1RHI = 0.8 * (cn.get('\\NBI1RHI')).data()
    NBI2LHI = 0.8 * (cn.get('\\NBI2LHI')).data()
    NBI2RHI = 0.8 * (cn.get('\\NBI2RHI')).data()
    # time_NBI=cn.get('dim_of(\\NBI2RHI)').data()

    background_NBI1LHI = np.mean(NBI1LHI[:idx_begin])
    NBI1LHI = NBI1LHI - background_NBI1LHI
    background_NBI1RHI = np.mean(NBI1RHI[:idx_begin])
    NBI1RHI = NBI1RHI - background_NBI1RHI
    background_NBI2RHI = np.mean(NBI2RHI[:idx_begin])
    NBI2RHI = NBI2RHI - background_NBI2RHI
    background_NBI2LHI = np.mean(NBI2LHI[:idx_begin])
    NBI2LHI = NBI2LHI - background_NBI2LHI
    NBI = NBI1LHI + NBI1RHI + NBI2LHI + NBI2RHI
    NBI = NBI[idx_begin:idx_end]
    NBI = NBI * kp

    NBI[NBI < 5e18] = 0
    NBI[1:1000] = 0
    cn.closeTree('east_1',shot)
    return NBI


def get_pi(shot):
    (idx_begin,idx_end),_ = get_time(shot)
    cn.openTree('east_1', shot)
    kp = 2e20
    PI = cn.get('\\vpi20').data()
    PI = PI[idx_begin:idx_end]
    PI = np.diff(PI)
    PI = np.append(PI, 0)
    PI[PI <= 4] = 0
    PI[PI > 4] = 1
    PI = np.cumsum(PI)
    cn.closeTree('east_1', shot)
    return kp * PI


# part two: pump
def get_pump(shot):
    (idx_begin,idx_end),_ = get_time(shot)
    cn.openTree('east_1', shot)
    KD2 = 2.9
    Pam2P = 4.82e20
    HD = 0.5
    main_speed = 77
    up_speed = 67
    low_speed = 60.5

    G107 = cn.get('\\G107').data()
    if len(G107) < 7000:
        G107 = cn.get('\\G101').data()

    G107 = G107[idx_begin:idx_end]
    G107 = 10 ** (1.667 * G107 - 9.333)

    G106 = cn.get('\\G106').data()

    if max(G106) < 1:
        G106 = G107
    else:
        G106 = G106[idx_begin:idx_end]
        G106 = 10 ** (1.667 * G106 - 9.333)

    G109 = cn.get('\\G109').data()

    if max(G109) < 1:
        G109 = G107
    else:
        G109 = G109[idx_begin:idx_end]
        G109 = 10 ** (1.667 * G109 - 9.333)

    data = up_speed * G109 + low_speed * G106 + main_speed * G107
    data = HD * KD2 * Pam2P * data

    data[data < 5e18] = 0
    data[1:1000] = 0

    # return savgol_filter(data,1000,3)
    cn.closeTree('east_1', shot)
    return data

def get_plasma(shot):
    kp = 2e13 / 3
    (idx_begin,idx_end),time_ip = get_time(shot)
    cn.openTree('pcs_east', shot)
    ne = cn.get('\\dfsdev2').data()
    t_ne = cn.get('dim_of(\\dfsdev2)').data()
    f_ne = interp1d(t_ne, ne, kind='cubic')  # Original freq is 500HZ, convert it to 1000HZ
    t_ne2 = np.arange(t_ne[0], t_ne[-1], 0.001)  # t_ne start from -4s, so t[3000] is -1s
    ne2 = f_ne(t_ne2)  # so idx_ne_begin = 3000, idx_ne_end = idx_end - idx_begin + idx_ne_begin
    idx_ne_begin = 3000
    idx_ne_end = idx_end - idx_begin + idx_ne_begin
    if len(ne2) < idx_ne_end:  # Uniform time
        ne_plasma = np.append(ne2[idx_ne_begin:], np.zeros(idx_ne_end - len(ne2)))
    else:
        ne_plasma = ne2[idx_ne_begin:idx_ne_end]
    cn.closeTree('pcs_east', shot)

    # Compute volume
    cn.openTree('efit_east', shot)
    volume = cn.get('\\volume').data()
    t_vol = cn.get('dim_of(\\volume)').data()
    f_v = interp1d(t_vol, volume, kind='cubic')
    old_s = int((t_vol[0] - time_ip[0]) * 1000)
    old_e = int((t_vol[-1] - time_ip[0]) * 1000)
    V_plasma = f_v(time_ip[old_s:old_e])
    V_plasma_long = np.zeros(idx_end - idx_begin)
    V_plasma_long[old_s - idx_begin:old_e - idx_begin] = V_plasma
    cn.closeTree('efit_east', shot)
    return kp * ne_plasma * V_plasma_long

shot  = 106912
# shot = 100000
(idx_begin,idx_end) , time_ip = get_time(shot)
# puffing = get_puff(shot)
# print(len(puffing))
# plt.plot(time_ip[idx_begin:idx_end],puffing)
# plt.show()

#SMBI
# smbi = get_SMBI(shot)
# print(len(smbi))
# plt.plot(time_ip[idx_begin:idx_end],smbi)
# plt.show()

#NBI
# NBI = get_NBI(shot)
# print(len(NBI))
# plt.plot(time_ip[idx_begin:idx_end],NBI)
# plt.show()

# PI
# PI = get_pi(shot)
# print(len(PI))
# plt.plot(time_ip[idx_begin:idx_end],PI)
# plt.show()

# pump
# pump = get_pump(shot)
# print(len(pump))
# plt.plot(time_ip[idx_begin:idx_end],pump)
# plt.show()

# Plasma
# plasma = get_plasma(shot)
# print(len(plasma))
# plt.plot(time_ip[idx_begin:idx_end],plasma)
# plt.show()


# Injection
# igt = get_inj(shot)
# print(len(igt))
# plt.plot(time_ip[idx_begin:idx_end],igt)
# plt.show()

# Exhaust
# eut = get_exh(shot)
# print(len(eut))
# plt.plot(time_ip[idx_begin:idx_end],eut)
# plt.show()

# Retention
ret = get_retention(shot,rate=True)
print(len(ret))
plt.plot(time_ip[idx_begin:idx_end],ret)
plt.show()

# ret = get_retention(shot,rate=True)
# print(len(ret))
# plt.plot(time_ip[idx_begin:idx_end],ret)
# plt.show()
# plt.plot(retention)
