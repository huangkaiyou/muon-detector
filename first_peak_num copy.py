#%%
import pandas as pd
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import os
from iminuit import Minuit
from iminuit.cost import LeastSquares
import scienceplots
# plt.ioff()

def horizon(x, offset) :
    return 0 * x + offset
plt.style.use('science')
plt.rcParams['figure.figsize'] = (15.0, 10.0) 

### initial parameter ###
threshold_first_ch1 = 2
threshold_first_ch2 = 1
stable_level = 0.03
prominence= 0
num_sigma = 7
ground_sigma_num = 5
time_interval = 1e-9
gate_interval = 10
file_names = ['20240424', '20240428', '20240429', '20240429_02', '20240430', '20240501', '20240502', '20240504', '20240505', '20240507', '20240510']
# file_name = '20240515'
#----------------------------------------------------------------#
def fitting(t, a, b, tau): #charge 
    return a * np.exp(-t/ tau) + b


'''
total_index_sigma_ch1 = [[] for i in range(15)]
total_index_sigma_ch2 = [[] for i in range(15)]
total_value_sigma_ch1 = [[] for i in range(15)]
total_value_sigma_ch2 = [[] for i in range(15)]

no_total_index_sigma_ch1 = [[] for i in range(15)]
no_total_index_sigma_ch2 = [[] for i in range(15)]
no_total_value_sigma_ch1 = [[] for i in range(15)]
no_total_value_sigma_ch2 = [[] for i in range(15)]
'''
lifetime_list = [[] for i in range(15)]

new_index_ch1 = []
new_index_ch2 = []

count = 0
time = 0
for file_name in file_names:
    while time < 10:
    # for count in range(1000):
        #process...
        file_path = "./new_lifetimedatatxt/Lifetime_data(B)/" + file_name + ".txt"
        f= open(file= file_path, mode= 'r')
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
        # smoothed1 = ch1.rolling(window_size, center= True, min_periods=1).mean()
        smoothed2 = ch2.rolling(window_size, center= True, min_periods=1).mean()
        
        #first peak finding 
        # ch1_peak_index, ch1_peak_info = find_peaks(smoothed1, height= threshold_first_ch1)
        ch2_first_peak_index, ch2_first_peak_info = find_peaks(smoothed2, height= threshold_first_ch2)

        # print(time)

        # ch1_peak_index = ch1_peak_index
        ch2_first_peak_index = ch2_first_peak_index[0]
        # ch1_peak_height = ch1_peak_info['peak_heights']
        ch2_first_peak_height = ch2_first_peak_info['peak_heights'][0]
        
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


        Tau = m.values['tau']
        cut_index = ch2_first_peak_index + 5 * int(Tau / time_interval)
        #-------------find the second peak----------------------
        # raw3 = smoothed1.iloc[cut_index:]
        raw4 = smoothed2.iloc[cut_index:]
        #second peak finding
        for i in range(1) :
            loop = i
        # for i in [9]:
            # ch1_second_peak_index, ch1_second_peak_info = find_peaks(raw3, height= np.mean(raw3) + np.std(raw3) * (i + 1), prominence= np.std(raw3) * ((i + 1)-3))
            # ch2_second_peak_index, ch2_second_peak_info = find_peaks(raw4, height= np.mean(raw4) + np.std(raw4) * (i + 1), prominence= np.std(raw4) * ((i + 1)-3))
            
            # ch1_second_peak_index, ch1_second_peak_info = find_peaks(raw3, height= np.mean(raw3) + np.std(raw3) * (i + 1), prominence= np.std(raw3) * (2))
            ch2_second_peak_index, ch2_second_peak_info = find_peaks(raw4, height= np.mean(raw4) + np.std(raw4) * (i + 1), prominence= np.std(raw4) * (2))
            
            # ch1_second_peak_height = ch1_second_peak_info['peak_heights']
            ch2_second_peak_height = ch2_second_peak_info['peak_heights']
            
            # for i in range(len(ch1_second_peak_index)):
            #     index = i
            #     if raw3[i + cut_index] > raw3[i + cut_index + 2] :
            #         new_index_ch1.append(index)
            '''
            for i in range(len(ch2_second_peak_index)):
                    index = i
                    gate_condition = True
                    for j in range(gate_interval):
                        if (raw4[i + cut_index] < raw4[i + cut_index + j+1]) and (raw4[i + cut_index] < raw4[i + cut_index - (j+1)]):
                            gate_condition = False
                            break
                    if gate_condition:
                        new_index_ch2.append(index)
                    
            # ch1_second_peak_index = new_index_ch1
            ch2_second_peak_index = new_index_ch2
            
            # ch1_second_peak_height = [ch1_second_peak_height[i] for i in ch1_second_peak_index]
            ch2_second_peak_height = [raw4[i + cut_index] for i in ch2_second_peak_index]
            '''
            
            if len(ch2_second_peak_height) == 0:
                continue
            
            second_peak_condition= False
            for peak in ch2_second_peak_height:
                if (ch2_first_peak_height < peak) :
                    second_peak_condition = True
                    break
            if second_peak_condition:
                continue
            
            
            # ch1_second_peak_index += cut_index
            ch2_second_peak_index += cut_index
            
            lifetime= (ch2_second_peak_index[0]-ch2_first_peak_index)*time_interval
            if lifetime < 0.5e-6:
                continue
            # total_index_sigma_ch1[loop].append((ch1_second_peak_index))
            # total_index_sigma_ch2[loop].append((ch2_second_peak_index))
            
            # total_value_sigma_ch1[loop].append(ch1_second_peak_height)
            # total_value_sigma_ch2[loop].append(ch2_second_peak_height[0])
            lifetime_list[loop].append(lifetime)
            # peak_index_ch1 = np.array(ch1_peak_index)
            # peak_index_ch2 = np.append(ch2_first_peak_index, ch2_second_peak_index)
            # peak_value_ch1 = np.array(ch1_peak_height)
            # peak_value_ch2 = np.append(ch2_first_peak_height, ch2_second_peak_height)
            
            # print("Ch1 peak index: ", peak_index_ch1)
            # print("Ch2 peak index: ", peak_index_ch2)
            # print("Ch1 peak height: ", peak_value_ch1)
            # print("Ch2 peak height: ", peak_value_ch2)
            
            #----more info-----------------#
            # ys = smoothed2
            # xs = [i*time_interval for i in range(10000)]
            # plt.vlines(xs[cut_index], -0.05, 3, label= 'cut')
            # horizon_index = np.array(xs[cut_index:])
            time += 1

'''
peak_num_list = [sum(i) for i in total_index_sigma_ch2]
peak_num_per = [num/len(total_index_sigma_ch2[0]) for num in peak_num_list]
plt.scatter([i+1 for i in range(15)], peak_num_list, label = 'Peak Num')
plt.xticks(ticks=[2*i for i in range(8)], labels=[2*i for i in range(8)])
plt.xlabel('Sigma')
plt.ylabel('Entries')
plt.legend()
plt.show()
'''
#%%
for i in range(1):
    # no_combined_list = no_total_value_sigma_ch2[i]
    # combined_list = total_value_sigma_ch2[i]
    
    lifetime_lists = lifetime_list[i]

    # no_hists, no_bins,_ = plt.hist(no_combined_list, bins=np.linspace(0, 0.6, 10))
    hists, bins,_ = plt.hist(lifetime_lists, bins=np.linspace(0, 0.5e-5, 100))
    plt.clf()

    # bin_centers = (no_bins[:-1] + no_bins[1:]) / 2
    # plt.scatter(bin_centers, no_hists,alpha=0.5, label = 'no prominence')

    bin_centers = (bins[:-1] + bins[1:]) / 2
    # plt.scatter(bin_centers, hists,alpha=0.5, label = 'prominence: 2 sigma')
    plt.plot(bin_centers, hists,alpha=0.5, label = 'prominence: 2 sigma')
    # plt.xlim(0, 1)
    plt.xlabel('time')
    plt.ylabel('Entries')
    plt.yscale('log')
    plt.legend()
    analysis_figure_dir = './analysis_result/all_peak_num/'
    # plt.savefig(analysis_figure_dir + f'{i}' + '.jpg')
    plt.show()
    plt.clf()
   

print('finish')
# %%
