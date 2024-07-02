#%%
import pandas as pd
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import os
from iminuit import Minuit
from iminuit.cost import LeastSquares
import scienceplots

def horizon(x, offset) :
    return 0 * x + offset

# plt.rcdefaults
plt.style.use(['science', 'notebook'])
plt.rcParams['figure.figsize'] = (10.0, 10.0) 
### initial parameter ###
threshold_first_ch1 = 2
threshold_first_ch2 = 1
stable_level = 0.03
prominence= 0
time_interval = 1e-9

# file_names = ['20240416_01']
file_names = ['20240424']
# '20240416_01','20240417_01','20240417_02', '20240418_01','20240419_01','20240421_01','20240422_01','20240422_02','20240423_01'
# file_names = ['20240416_01']

#----------------------------------------------------------------#

#%%
def fitting(t, a, b, tau): #charge 
    return a * np.exp(-t/ tau) + b

tau_list = []
decay_list = []
peak_list = []
for file_name in file_names:
           
    count = 0
    file_path = "./lifetimedatatxt/Lifetime_data(B)/" + file_name + ".txt"
    #-------------------------------------------------------#
    f= open(file= file_path, mode= 'r')
    #-------------------------------------------------------#
    # while True:
    for i in range(50000):
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
        ch1_first_peak_index, ch1_first_peak_info = find_peaks(smoothed1, height= threshold_first_ch1)
        ch2_first_peak_index, ch2_first_peak_info = find_peaks(smoothed2, height= threshold_first_ch2)

        #filter two first peak cases
        if (len(ch1_first_peak_index) != 1 or len(ch2_first_peak_index) != 1):
            continue
        
        #filter first peaks are not happend at a time
        if abs(ch1_first_peak_index[0] - ch2_first_peak_index[0]) > 1 : 
            continue

        ch1_first_peak_height = ch1_first_peak_info['peak_heights'][0]
        ch2_first_peak_height = ch2_first_peak_info['peak_heights'][0]

        #look for smoothed
        for i in range(ch2_first_peak_index[0], 10000):
            if smoothed2[i] < stable_level:
                second_start = i
                break
        
        ys = smoothed2
        xs = [i*time_interval for i in range(10000)]
        
        #look actual smooth with fitting exponential
        #-----------fitting part--------------------------                         
        yerror = 0.007
        start = int((ch2_first_peak_index[0] * 4+second_start * 1)/ 5)
        # start = ch1_first_peak_index[0]
        xs_fit = np.array(xs[start : second_start])
        ys_fit = ys[start: second_start]
        
        #rescale x
        # print('i should see you')
        xs_temp = xs_fit - xs_fit[0]
        
        tau = (second_start - start) * time_interval / 2 #decay about 80%
        least_square= LeastSquares(xs_temp, ys_fit, yerror= yerror, model= fitting)
        m=Minuit(least_square, a= ch2_first_peak_height, b= 0, tau= 1e-7)
        m.limits['a'] = (0, 10)
        m.limits['b'] = (-1, 1)
        m.limits['tau'] =  (1e-9, 1e-7)
        # plt.errorbar(xs, ys, yerr= stdev(ys)/ , fmt='')
        m.migrad() 
        m.hesse()
        
        #Judgement whether the fitting is great
        chi_square = m.fmin.reduced_chi2
        if chi_square > 10 :
            continue
        # fitting_end = start + 5 * int(tau/time_interval)
        #--------------------------------------#
        # raw3 = smoothed1.iloc[fitting_end:]
        # raw4 = smoothed2.iloc[fitting_end:]
        if ch2_first_peak_height > 8 :
            continue
        cut_index = ch2_first_peak_index[0] + 5 * int(tau / time_interval)
        
        # tau vs. peak height
        peak_list.append(ch2_first_peak_height)
        tau_value =  m.values['tau']
        tau_list.append(tau_value)
        
        
        # plt.scatter(peak_list, tau_list)
        # plt.xlabel('Peak Height')
        # plt.ylabel('Tau Value')
        # plt.grid()
        # plt.show()
        
        
        # decay value vs. peak height
        decay_value = np.mean([smoothed2[cut_index]])
        decay_list.append(decay_value)
        
        
        
        
        '''
        #----more info-----------------#
        plt.vlines(xs[cut_index], -0.05, 3, label= 'cut')
        
        #show the plot            
        lifetime = time_interval * (peak_index_ch2[1] - peak_index_ch2[0])
        lifetime = float(lifetime)
        # print(count)
        plt.title(f"lifetime: {lifetime:.3e} (s)")
        plt.scatter(xs, ys, label='CH2 : Smothed', alpha=0.5, s=5)
        peak_index_ch2 = list(peak_index_ch2)
        plt.scatter([xs[i] for i in peak_index_ch2], [ys[i] for i in peak_index_ch2], color='red', label='CH2 : Peaks')
        plt.legend()
        
        
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
        
        
        plt.tight_layout()
        # plt.show()
        plt.savefig(frame_path + f'{count}' + '.jpg')
        plt.clf()
        '''
        
            # print('next')
    # except :      # works on python 3.x
        # print("interrupt")
# %%
# plt.scatter(peak_list, decay_list)
plt.hexbin(peak_list, decay_list, gridsize=20, cmap='Blues')
plt.colorbar(label='Counts')
# plt.title('Relation between First Peak Height and Decay Value')
plt.xlabel('First Peak Height (V)')
plt.ylabel('Decay Value (V)')
plt.grid()
plt.show()

plt.hexbin(peak_list, tau_list, gridsize=20, cmap='Blues')
plt.colorbar(label='Counts')
# plt.title('Relation between First Peak Height and Decay Value')
plt.xlabel('First Peak Height (V)')
plt.ylabel(r'Tau value ($\mu$)')
plt.yticks([i*0.25e-8 for i in range(0, 15, 2)], [i*0.25 for i in range(0, 15, 2)])
plt.grid()
plt.show()

# %%
