#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import find_peaks

threshold_second = 0.15
threshold_first_ch1 = 2
threshold_first_ch2 = 1
trigger_level = 0.15
window_size = 9
stable_level = 0.03
num_sigma = 10
file_name = '20240429'
file_path = "lifetimedatatxt/Lifetime_data(B)/" + file_name + ".txt"
f= open(file= file_path, mode= 'r')

line = 2657
line = 5317-1
for i in range(line):
    f.readline()
    f.readline()
#channel 1 data read
raw1 = f.readline()
raw1 = raw1.split('\t')
raw1[-1] = raw1[-1].strip()

#store the trimed data
ch1 = []
ch2 = []

for i in raw1:
    ch1.append(float(i))

#channel 2 data read
raw2 = f.readline()
raw2 = raw2.split('\t')
raw2[-1] = raw2[-1].strip()

for i in raw2:
    ch2.append(float(i))

ch1 = pd.Series(ch1)
ch2 = pd.Series(ch2)
#-------------------------------#  
# Smooth the data with a rolling window to reduce noise
window_size = 9
smoothed1 = ch1.rolling(window_size, center= True, min_periods=1).mean()
smoothed2 = ch2.rolling(window_size, center= True, min_periods=1).mean()
            
# peak1, info1 = find_peaks(smoothed1, height= 2)
# peak2, info2 = find_peaks(smoothed2, height= 0.15, prominence=0.15)

for i in range(1):
    ch1_first_peak_index = []
    ch2_first_peak_index = []
    
    # Smooth the data with a rolling window to reduce noise
    window_size = 9
    smoothed1 = ch1.rolling(window_size, center= True, min_periods=1).mean()
    smoothed2 = ch2.rolling(window_size, center= True, min_periods=1).mean()
    
    #first peak finding 
    ch1_first_peak_index, ch1_first_peak_info = find_peaks(smoothed1, height= threshold_first_ch1)
    ch2_first_peak_index, ch2_first_peak_info = find_peaks(smoothed2, height= threshold_first_ch2)

    #filter two first peak cases
    if (len(ch1_first_peak_index) != 1 or len(ch2_first_peak_index) != 1):
        continue

    ch1_first_peak_height = ch1_first_peak_info['peak_heights'][0]
    ch2_first_peak_height = ch2_first_peak_info['peak_heights'][0]

    #look for smoothed
    for i in range(ch2_first_peak_index[0], 10000):
        if smoothed2[i] < stable_level:
            second_start = i
            break
    
    raw3 = smoothed1.iloc[second_start:]
    raw4 = smoothed2.iloc[second_start:]
    
    
    # plt.scatter(xs[second_start], ys[second_start], s=10)
    #second peak finding 
    ch1_second_peak_index, ch1_second_peak_info = find_peaks(raw3, height= np.mean(raw3) + np.std(raw3) * num_sigma, prominence= np.std(raw3) * num_sigma)
    ch2_second_peak_index, ch2_second_peak_info = find_peaks(raw4, height= np.mean(raw4) + np.std(raw4) * num_sigma, prominence= np.std(raw4) * num_sigma)

    ch1_second_peak_index +=second_start
    ch2_second_peak_index += second_start
    # print('a')
    if len(ch2_second_peak_index) != 1 :
        continue
    # print('b')
    ch1_second_peak_height = ch1_second_peak_info['peak_heights']
    ch2_second_peak_height = ch2_second_peak_info['peak_heights']
    
    #combine two segment
    peak_index_ch1 = np.append(ch1_first_peak_index, ch1_second_peak_index)
    peak_index_ch2 = np.append(ch2_first_peak_index, ch2_second_peak_index)
    peak_value_ch1 = np.append(ch1_first_peak_height, ch1_second_peak_height)
    peak_value_ch2 = np.append(ch2_first_peak_height, ch2_second_peak_height)
    
    # peak2 = np.append(ch2_first_peak_index, ch2_second_peak_index)
    print(peak_index_ch2)
    
    peak2 = peak_index_ch2
    # first_peak_index = peak1[0]
    # peak2, info = find_peaks(smoothed.iloc[first_peak_index:], height= threshold)
    # peak2 += first_peak_index
    xs = [0+ i*1e-9 for i in range(10000)]
    ys = smoothed2 
    # np.savetxt('xs_data.txt', xs)
    # np.savetxt('ys_data.txt', ys)
    plt.scatter(xs, ys, s=5, alpha= 0.5)
    plt.scatter([xs[i] for i in peak2], [ys[i] for i in peak2], color='red', label='CH2 : Peaks')
    plt.legend()
    plt.tight_layout()
    # plt.savefig('./analysis_result/lifetime/scipy/0.15v/frame' + '/test.jpg')
    plt.show()
    # np.savetxt('ys.txt', smoothed2)
# %%
