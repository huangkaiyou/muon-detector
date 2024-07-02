#%%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scienceplots
from iminuit import Minuit
from iminuit.cost import LeastSquares
from scipy.signal import find_peaks
from matplotlib import ticker
def horizon(x, offset) :
    return 0 * x + offset
def fitting(t, a, b, tau): #charge 
    return a * np.exp(-t/ tau) + b

plt.style.use(['science', 'notebook'])
# from OuOb import pltSty2


# xs = np.loadtxt('xs_data.txt')
# ys = np.loadtxt('ys_data.txt')

plt.rcParams['figure.figsize'] = (10.0, 5.0) 
### initial parameter ###
threshold_first_ch1 = 2
threshold_first_ch2 = 1
stable_level = 0.03
prominence= 2
num_sigma = 2
ground_sigma_num = 2
time_interval = 1e-9
file_name = '20240418_01'
file_path = "./lifetimedatatxt/Lifetime_data(no B)/" + file_name + ".txt"
f= open(file= file_path, mode= 'r')

line = 6455
line = 4
start = 0.8e-6
end = 2.1e-6
# line = 2
count = 0
while True:
    count+=1
    #channel 1 data read
    raw1 = f.readline()
    raw1 = raw1.split('\t')
    raw1[-1] = raw1[-1].strip()

    #store the trimed data
    ch1 = []
    ch2 = []

    for i in raw1[:-1]:
        ch1.append(float(i))

    #channel 2 data read
    raw2 = f.readline()
    raw2 = raw2.split('\t')
    raw2[-1] = raw2[-1].strip()

    for i in raw2[:-1]:
        ch2.append(float(i))

    ch1 = pd.Series(ch1)
    ch2 = pd.Series(ch2)

    #-------------------------------#  
    # Smooth the data with a rolling window to reduce noise
    window_size = 9
    smoothed1 = ch1.rolling(window_size, center= True, min_periods=1).mean()
    smoothed2 = ch2.rolling(window_size, center= True, min_periods=1).mean()

    #first peak finding 
    ch1_peak_index, ch1_peak_info = find_peaks(smoothed1, height= threshold_first_ch1)
    ch2_first_peak_index, ch2_first_peak_info = find_peaks(smoothed2, height= threshold_first_ch2)

    #filter zero first peak cases
    if (len(ch1_peak_index) == 0 or len(ch2_first_peak_index) == 0):
        continue
    
    # -----check the first peak--------------#
    if abs(ch2_first_peak_index[0] - ch1_peak_index[0]) > 1:
        continue
    
    if ch1_peak_info['peak_heights'] < 8 :
        continue
    ys = smoothed2
    xs = [i*time_interval for i in range(10000)]
    # plt.scatter(xs, smoothed2, label= "smoothed", alpha= 0.5, s=5)
    start_index = int(start/time_interval)
    end_index = int(end/time_interval)

    xs = np.array(xs[start_index:end_index])
    xs -= xs[0]
    smoothed1 = smoothed1[start_index:end_index]
    smoothed2 = smoothed2[start_index:end_index]

    # Create a figure with subplots and set the figure size
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    #----more info-----------------#
    ax1.scatter(xs, smoothed1, label = 'CH1', color = 'black', s=5)
    ax2.scatter(xs, smoothed2, label = 'CH2', color = 'black', s=5)

    xtick = np.arange(0, int(end_index/1000), 1)
    ax1.set_xticks(xtick * 1e-6)
    ax1.set_xticklabels(xtick)

    ax2.set_xticks(xtick * 1e-6)
    ax2.set_xticklabels(xtick)

    ax1.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.2f}"))
    ax2.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.2f}"))

    # ax1.set_xlabel(r'Time $(\mu s)$')
    ax2.set_xlabel(r'Time $(\mu s)$', loc='right')
    ax1.set_ylabel('Voltage (v)', loc= 'top')
    # ax2.set_ylabel('Voltage (v)')
    
    # ax1.set_yscale('log')
    # ax2.set_yscale('log')
    ax1.legend(loc = 'best', fontsize = 15)
    ax2.legend(loc = 'best', fontsize = 15)
    print(count)
    plt.show()
# print(ch2[15]-ch2[14])
#%%
