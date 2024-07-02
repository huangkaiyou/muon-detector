#%%
import pandas as pd
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import os
from iminuit import Minuit
from iminuit.cost import LeastSquares
import scienceplots


plt.rcParams['figure.figsize'] = (15.0, 10.0)
### initial parameter ###
threshold_first_ch1 = 2
threshold_first_ch2 = 1
stable_level = 0.03
prominence= 0
num_sigma = 7
ground_sigma_num = 5
time_interval = 1e-9

mode = '_stable' + ' (5tau)_fitting' + 'second_no_sigma'
prominice_value = '_' + 'prominence(sigma)'
ground_value = '_' + f'ground({ground_sigma_num})'
# file_names = ['20240424', '20240428', '20240429', '20240429_02', '20240430', '20240501', '20240502', '20240504', '20240505', '20240507', '20240510']
file_name = '20240510'
#----------------------------------------------------------------#
def fitting(t, a, b, tau): #charge 
    return a * np.exp(-t/ tau) + b


#-------------------------------------------------------#
file_path = "./new_lifetimedatatxt/Lifetime_data(B)/" + file_name + ".txt"
f= open(file= file_path, mode= 'r')
#%%
#----------------------------------------------------------------#
ys = []
left_list_ch1 = []
left_list_ch2 = []
count = 0
Times = 0
frames = 1
check = 0
try:
    # while True:
    while Times < frames :
        #process...
        check += 1
        if (count) % 1000 == 0:
                print("count: ", count)
        count+= 1
        #-------------------------------#
        #store the trimed data
        ch1 = []
        ch2 = []
        
        #channel 1 data read
        raw1 = f.readline()
        
        #check whether have data to read
        if raw1 == '':
            break
        raw1 = raw1.split('\t')
        raw1[-1] = raw1[-1].strip()
        
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

        # -----check the first peak--------------#
        #filter zero first peak cases
        if (len(ch1_peak_index) == 0 or len(ch2_first_peak_index) == 0):
            continue
        
        if abs(ch1_peak_index[0] - ch2_first_peak_index[0]) > 1:
            continue
        
        ch1_peak_index = ch1_peak_index
        ch2_first_peak_index = ch2_first_peak_index[0]
        ch1_peak_height = ch1_peak_info['peak_heights']
        ch2_first_peak_height = ch2_first_peak_info['peak_heights'][0]
        
        
        if (ch2_first_peak_height > 8):
            continue
        #look for smoothed
        for i in range(ch2_first_peak_index, 10000):
            if smoothed2[i] < stable_level:
                second_start = i
                break
        
        ys = smoothed2
        xs = [i*time_interval for i in range(10000)]
        #look actual smooth with fitting exponential
        #---------- -fitting part--------------------------                         
        yerror = 0.007
        start = int((ch2_first_peak_index * 4+second_start * 1)/ 5)
        # start = ch1_first_peak_index[0]
        xs_fit = np.array(xs[start : second_start])
        ys_fit = ys[start: second_start]

        #rescale x
        # print('i should see you')
        xs_temp = xs_fit - xs_fit[0]

        Tau = (second_start - start) * time_interval / 2 #decay about 80%
        least_square= LeastSquares(xs_temp, ys_fit, yerror= yerror, model= fitting)
        m=Minuit(least_square, a= ch2_first_peak_height, b= 0, tau= Tau)
        m.limits['a'] = (0, 10)
        m.limits['b'] = (-1, 1)
        m.limits['tau'] =  (1e-9, 1e-6)

        m.migrad() 
        m.hesse()

        #Judgement whether the fitting is great
        chi_square = m.fmin.reduced_chi2
        if chi_square > 10 :
            continue
        
        Tau = m.values['tau']
        cut_index = ch2_first_peak_index + 5 * int(Tau / time_interval)
        #-------------find the second peak----------------------
        raw3 = smoothed1.iloc[cut_index:]
        raw4 = smoothed2.iloc[cut_index:]
        # plt.hist(ys)
        # plt.show()
        left = raw3.to_list()
        left_list_ch1.extend(left)
        
        left = raw4.to_list()
        left_list_ch2.extend(left)
        
        Times +=1
        # if max(smoothed[:200])> 0.05:
        #     print(file_path)
finally:
    f.close()
# %%
import math

num_event = 0
#histogram
log_scale= True
counts, bins, patch = plt.hist(left_list_ch2, bins= np.linspace(-0.06, 0.1, 50), alpha=1, color='skyblue', edgecolor='navy', label=  'data', histtype= 'step', log= log_scale)
plt.clf()

bin_centers = (bins[:-1] + bins[1:]) / 2
plt.scatter(bin_centers, counts, alpha=0.5, label = 'no prominence')
plt.yscale('log')
num_event = len(left_list_ch2)
# print("the ratio: ", num_event/total_event*100, '(%)')
# Add labels and title
plt.xlabel('Voltage')
plt.ylabel('Entries')
plt.title('total data num: ' + str(num_event))
plt.legend()
# plt.xlim(0, 1.5)
# Show the plot
 

analysis_figure_dir = './analysis_result/noise/'
# plt.savefig(analysis_figure_dir + f'frame_num {Times} ' + '.jpg')
plt.show() 
plt.clf()
