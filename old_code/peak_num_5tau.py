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
# plt.rcParams['figure.figsize'] = (15.0, 10.0) 

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
file_name = '20240424'
#----------------------------------------------------------------#
def fitting(t, a, b, tau): #charge 
    return a * np.exp(-t/ tau) + b

#-------------------------------------------------------#
file_path = "./lifetimedatatxt/Lifetime_data(B)/" + file_name + ".txt"
f= open(file= file_path, mode= 'r')

final = {}
df_final = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in final.items() ]))

#-------------------------------------------------------#
#%%
count = 0
timeee= 0

total_index_sigma_ch1 = [[] for i in range(15)]
total_index_sigma_ch2 = [[] for i in range(15)]
total_value_sigma_ch1 = [[] for i in range(15)]
total_value_sigma_ch2 = [[] for i in range(15)]

fit_total_index_sigma_ch1 = [[] for i in range(15)]
fit_total_index_sigma_ch2 = [[] for i in range(15)]
fit_total_value_sigma_ch1 = [[] for i in range(15)]
fit_total_value_sigma_ch2 = [[] for i in range(15)]

lifetime_list = [[] for i in range(15)]
fit_lifetime_list = [[] for i in range(15)]
# while True:
for time in range(30000):
    #process...
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
    
    #first peak finding 
    ch1_peak_index, ch1_peak_info = find_peaks(smoothed1, height= threshold_first_ch1)
    ch2_first_peak_index, ch2_first_peak_info = find_peaks(smoothed2, height= threshold_first_ch2)

    #filter zero first peak cases
    if (len(ch1_peak_index) == 0 or len(ch2_first_peak_index) == 0):
        continue
    
    if abs(ch1_peak_index[0] - ch2_first_peak_index[0]) > 1 :
        continue
    # timeee += 1
    # print(time)

    ch1_peak_index = ch1_peak_index
    ch2_first_peak_index = ch2_first_peak_index[0]
    ch1_peak_height = ch1_peak_info['peak_heights']
    ch2_first_peak_height = ch2_first_peak_info['peak_heights'][0]
    
    #look for smoothed
    for i in range(ch2_first_peak_index, 10000):
        if smoothed2[i] < stable_level:
            second_start = i
            break
    
    ys = smoothed2
    xs = [i*time_interval for i in range(10000)]
    
    #look actual smooth with fitting exponential
    #-----------fitting part--------------------------                         
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
    timeee +=1
    Tau = m.values['tau']
    # fitting_end = start + 5 * int(tau/time_interval)
    #--------------------------------------#
    # raw3 = smoothed1.iloc[fitting_end:]
    # raw4 = smoothed2.iloc[fitting_end:]
    
    
    #-------------find the second peak
    # cut_index = ch2_first_peak_index + 5 * int(Tau / time_interval)
    cut_index = ch2_first_peak_index + 1
    raw3 = smoothed1.iloc[cut_index:]
    raw4 = smoothed2.iloc[cut_index:]
    #second peak finding
    for i in range(15) :
    # for i in [9]:
        # ch1_second_peak_index, ch1_second_peak_info = find_peaks(raw3, height= np.mean(raw3) + np.std(raw3) * (i + 1), prominence= np.std(raw3) * ((i + 1)-3))
        # ch2_second_peak_index, ch2_second_peak_info = find_peaks(raw4, height= np.mean(raw4) + np.std(raw4) * (i + 1), prominence= np.std(raw4) * ((i + 1)-3))
        
        ch1_second_peak_index, ch1_second_peak_info = find_peaks(raw3, height= np.mean(raw3) + np.std(raw3[1000:]) * (i + 1), prominence= np.std(raw3[1000:]) * (2))
        ch2_second_peak_index, ch2_second_peak_info = find_peaks(raw4, height= np.mean(raw4) + np.std(raw4[1000:]) * (i + 1), prominence= np.std(raw4[1000:]) * (2))
        
        ch1_second_peak_height = ch1_second_peak_info['peak_heights']
        ch2_second_peak_height = ch2_second_peak_info['peak_heights']
        
        
        ch1_second_peak_index += cut_index
        ch2_second_peak_index += cut_index
        
        if len(ch2_second_peak_index) < 1 :
            continue
        # total_index_sigma_ch1[i].append((ch1_second_peak_index[0]))
        total_index_sigma_ch2[i].append((ch2_second_peak_index[0]))
        
        # total_value_sigma_ch1[i].append(ch1_second_peak_height[0])
        total_value_sigma_ch2[i].append(ch2_second_peak_height[0])
        
        
        peak_index_ch1 = np.array(ch1_peak_index)
        peak_index_ch2 = np.append(ch2_first_peak_index, ch2_second_peak_index[0])
        peak_value_ch1 = np.array(ch1_peak_height)
        peak_value_ch2 = np.append(ch2_first_peak_height, ch2_second_peak_height[0])
        
        # print("Ch1 peak index: ", peak_index_ch1)
        # print("Ch2 peak index: ", peak_index_ch2)
        # print("Ch1 peak height: ", peak_value_ch1)
        # print("Ch2 peak height: ", peak_value_ch2)
        
        #----more info-----------------#
        ys = smoothed2
        xs = [i*time_interval for i in range(10000)]
        # plt.vlines(xs[cut_index], -0.05, 3, label= 'cut')
        horizon_index = np.array(xs[cut_index:])
        
        if len(peak_index_ch2) < 2 :
            continue       
        lifetime = time_interval * (peak_index_ch2[1] - peak_index_ch2[0])
        lifetime = float(lifetime)
        lifetime_list[i].append(lifetime)
        '''
        plt.plot(xs[cut_index:], horizon(horizon_index, offset= np.mean(raw4) +np.std(raw4[1000:]) * (i+1)), label= f"${(i+1)}\sigma$ : Second Height", linestyle= 'dashed', color = 'brown', lw= 3)
        plt.scatter(xs, smoothed2, label='CH2 : Smothed', alpha=0.5, s=5)
        peak_index_ch2 = list(peak_index_ch2)
        plt.scatter([xs[i] for i in peak_index_ch2], [ys[i] for i in peak_index_ch2], color='red', label='CH2 : Peaks')
        plt.xticks(ticks=[i*1e-6 for i in range(10)], labels=[i for i in range(10)])
        plt.xlabel('Time' + r'($\mu$s)')
        plt.ylabel('Voltage')
        plt.legend()
        plt.tight_layout()
        plt.xlim(850e-9, 1300e-9)
        plt.show()
        '''
        '''
        #show the plot  
        
        # print(count)
        plt.title(f"lifetime: {lifetime:.3e} (s)")
       
        # plt.scatter(xs, smoothed1, label='CH1 : Smothed', alpha=0.5, s=5)
        plt.plot(xs[cut_index:], horizon(horizon_index, offset= np.mean(raw4) +np.std(raw4[1000:]) * (i+1)), label= f"${(i+1)}\sigma$ : Second Height", linestyle= 'dashed', color = 'brown', lw= 3)
        plt.scatter(xs, smoothed2, label='CH2 : Smothed', alpha=0.5, s=5)
        peak_index_ch2 = list(peak_index_ch2)
        plt.scatter([xs[i] for i in peak_index_ch2], [ys[i] for i in peak_index_ch2], color='red', label='CH2 : Peaks')
        plt.xticks(ticks=[i*1e-6 for i in range(10)], labels=[i for i in range(10)])

        plt.xlabel('Time' + r'($\mu$s)')
        plt.ylabel('Voltage')
        plt.legend()
        plt.tight_layout()
        # plt.show()
        # plt.clf()
        
        fit_info = [
            f"$\\chi^2$ / $n_\\mathrm{{dof}}$ = {m.fval:.1f} / {len(xs_fit) - m.nfit}",
        ]
        for p, v, e in zip(m.parameters, m.values, m.errors):
            fit_info.append(f"{p} = ${v:.3e} \\pm {e:.3e}$")
        # plt.errorbar(xs, ys, yerr=yerror, fmt="o", label="data",ms= 3, alpha=0.6, color='brown')
        #---------rescale x back------------------
        ax1 = plt.subplot(1,1,1)
        plt.plot(xs_fit, fitting(xs_temp, *m.values), label="Fitting", linestyle= 'dashed', color = 'red')
        ax1.legend(loc = 'best', edgecolor = '#7e7474', fontsize = 12)
        plt.legend(title="\n".join(fit_info))
        plt.show()
        # print(i)
        # plt.savefig('test.png')
        '''
    
    cut_index = ch2_first_peak_index + 5 * int(Tau / time_interval)
    # cut_index = ch2_first_peak_index + 1
    raw3_fit = smoothed1.iloc[cut_index:]
    raw4_fit = smoothed2.iloc[cut_index:]
    for i in range(15) :
        
        fit_ch1_second_peak_index, fit_ch1_second_peak_info = find_peaks(raw3_fit, height= np.mean(raw3) + np.std(raw3[1000:]) * (i + 1), prominence= np.std(raw3[1000:]) * (2))
        fit_ch2_second_peak_index, fit_ch2_second_peak_info = find_peaks(raw4_fit, height= np.mean(raw4) + np.std(raw4[1000:]) * (i + 1), prominence= np.std(raw4[1000:]) * (2))
        
        fit_ch1_second_peak_height = fit_ch1_second_peak_info['peak_heights']
        fit_ch2_second_peak_height = fit_ch2_second_peak_info['peak_heights']
        
        
        fit_ch1_second_peak_index += cut_index
        fit_ch2_second_peak_index += cut_index
        
        if len(fit_ch2_second_peak_index) < 1 :
            continue
        
        # fit_total_index_sigma_ch1[i].append((fit_ch1_second_peak_index[0]))
        fit_total_index_sigma_ch2[i].append((fit_ch2_second_peak_index[0]))
        
        # fit_total_value_sigma_ch1[i].append(fit_ch1_second_peak_height[0])
        fit_total_value_sigma_ch2[i].append(fit_ch2_second_peak_height[0])
        '''
        if len(fit_ch1_second_peak_height) != 0 :
            fit_total_value_sigma_ch1[i].append(max(fit_ch1_second_peak_height))
        if len(fit_ch2_second_peak_height) != 0 :
            fit_total_value_sigma_ch2[i].append(max(fit_ch2_second_peak_height))
        '''
        
        peak_index_ch2 = np.append(ch2_first_peak_index, fit_ch2_second_peak_index[0])
        # peak_value_ch2 = np.append(ch2_first_peak_height, fit_ch2_second_peak_height[0])
        
        # print("Ch1 peak index: ", peak_index_ch1)
        # print("Ch2 peak index: ", peak_index_ch2)
        # print("Ch1 peak height: ", peak_value_ch1)
        # print("Ch2 peak height: ", peak_value_ch2)
        
        #----more info-----------------#
        ys = smoothed2
        xs = [i*time_interval for i in range(10000)]
        # plt.vlines(xs[cut_index], -0.05, 3, label= 'cut')
        horizon_index = np.array(xs[cut_index:])
        
        if len(peak_index_ch2) < 1 :
            continue       
        lifetime = time_interval * (peak_index_ch2[1] - peak_index_ch2[0])
        lifetime = float(lifetime)
        fit_lifetime_list[i].append(lifetime)
        
        # print('no fit: ', ch2_second_peak_index[:5])
        # print('fit: ', fit_ch2_second_peak_index[:5])
        # plt.scatter(xs, smoothed2)
        # plt.scatter([xs[i] for i in ch2_second_peak_index], [ys[i] for i in ch2_second_peak_index], label = 'no fit', alpha= .5)
        # plt.scatter([xs[i] for i in fit_ch2_second_peak_index], [ys[i] for i in fit_ch2_second_peak_index], label = 'fit', alpha= .5)
        # plt.legend()
        # plt.show()
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
no_combined_list = lifetime_list[0]
combined_list = fit_lifetime_list[0]

# no_hists, no_bins,_ = plt.hist(no_combined_list, bins=[1e-7*i for i in range(0)])
# hists, bins,_ = plt.hist(combined_list, bins=[1e-7*i for i in range(40)])
no_hists, no_bins,_ = plt.hist(no_combined_list, bins=20)
hists, bins,_ = plt.hist(combined_list, bins=20)
plt.clf()

bin_centers = (no_bins[:-1] + no_bins[1:]) / 2
plt.scatter(bin_centers, no_hists, label = 'no fit',alpha= 0.5)

bin_centers = (bins[:-1] + bins[1:]) / 2
plt.scatter(bin_centers, hists, label = 'fit', alpha=.5)
# plt.xlim(0, 1)
plt.xlabel('lifetime')
plt.ylabel('Entries')
# plt.yscale('log')
plt.legend()
plt.show()


print('finish')
# %%
