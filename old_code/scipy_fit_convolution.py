#%%
import pandas as pd
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import os
import scienceplots
from iminuit import Minuit
from iminuit.cost import LeastSquares
import sympy as sp
import scipy.special as sse

def horizon(x, offset) :
    return 0 * x + offset

def exp(t, a, b, tau) : 
    return a* sp.exp(-t/ tau) + b

def gaussian(t, sigma, mu) :
    return sp.exp(-1/2 * ((t-mu)/sigma)**2)

# u : mean value
# l : lamdad (1/tau)
# s : sigma
def fitting(t, l, sigma, mu, a, b, c) :
     return a * (l/2) * np.exp((l/2) * (2*mu + l*(sigma**2) - 2*(t - b))) * sse.erfc((mu + l*(sigma**2) - (t - b)) / (np.sqrt(2) * sigma)) + c
#the fucking trial of sympy
# t, a, b, tau, sigma, mu, x = sp.symbols('t, a, b, tau, sigma, mu, x')
# fitting = sp.integrate(gaussian(t - x, sigma, mu) * exp(x, a, b, tau), (x, -sp.oo, sp.oo))

# plt.rcParams.update(plt.rcParamsDefault)
# plt.rcParams.update({
#     'text.usetex': True,
#     'font.family': 'serif',})

plt.style.use(['science', 'no-latex'])
plt.rcParams['figure.figsize'] = (10.0, 4.0) 
threshold_first_ch1 = 2
threshold_first_ch2 = 1
stable_level = 0.03
prominence= 0
num_sigma = 10
relevant_sigma_num = 7
time_interval = 1e-9
file_name = '20240424'
file_path = "./lifetimedatatxt/Lifetime_data(B)/" + file_name + ".txt"
f= open(file= file_path, mode= 'r')

line = 2078-1

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
#%%
#-------------------------------#  
# Smooth the data with a rolling window to reduce noise
window_size = 9
smoothed1 = ch1.rolling(window_size, center= True, min_periods=1).mean()
smoothed2 = ch2.rolling(window_size, center= True, min_periods=1).mean()
            
# peak1, info1 = find_peaks(smoothed1, height= 2)
# peak2, info2 = find_peaks(smoothed2, height= 0.15, prominence=0.15)

for i in range(1):
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
    '''
    #filter two first peak cases
    if (len(ch1_first_peak_index) == 0 or len(ch2_first_peak_index) == 0):
        print('first peak is 0')
        continue
    
    #filter first peaks are not happend at a time
    if abs(ch1_first_peak_index[0] - ch2_first_peak_index[0]) > 20 : 
        print('first peak is not enough high')
        continue
    '''
    # ch1_first_peak_height = ch1_first_peak_info['peak_heights'][0]
    # ch2_first_peak_height = ch2_first_peak_info['peak_heights'][0]

    
    #look for smoothed
    for j in range(0, 10000):
        if smoothed2.iloc[j] > stable_level + 0.01:
            # print(j)
            first_start = j
            break
        
    for i in range(first_start, 10000):
        if smoothed2.iloc[i] < stable_level:
            first_end = i
            break
    '''
    #look for smoothed
    for i in range(ch2_first_peak_index[0], 10000):
        if smoothed2[i] < stable_level:
            second_start = i
            break
    '''
    ys = smoothed2
    xs = [i*time_interval for i in range(10000)]
    # plt.scatter(xs[first_start], smoothed2.iloc[first_end], label = 'start', s= 10)
    # plt.scatter(xs[first_end], smoothed2.iloc[first_end], label = 'end', s= 15)
    # plt.scatter(xs, smoothed2, s= 5, alpha= 0.5)
    # plt.legend()
    # plt.show()
    #look actual smooth with fitting exponential
    #-----------fitting part--------------------------                         
    yerror = 0.007
    # start = int((ch2_first_peak_index[0] * 4+second_start * 1)/ 5)
    start = first_start
    end = first_end
    # start = ch1_first_peak_index[0]
    xs_fit = np.array(xs[start : end])
    ys_fit = ys[start: end]
    
    #rescale x
    # print('i should see you')
    xs_temp = xs_fit - xs_fit[0]
    least_square= LeastSquares(xs_temp, ys_fit, yerror= yerror, model= fitting)
    m=Minuit(least_square, a= 15e-7, b= 1e-6, c= 0, l= 1/2e-7, sigma= 1e-6, mu= 2e-7)
    # m.limits['a'] = (0, 10)
    # m.limits['b'] = (-1, 1)
    # m.limits['l'] =  (1/1e-7, 1/1e-9)
    # plt.errorbar(xs, ys, yerr= stdev(ys)/ , fmt='')
    m.migrad() 
    m.hesse()
    
    #Judgement whether the fitting is great
    chi_square = m.fmin.reduced_chi2
    # if chi_square > 10 :
    #     continue
    '''
    tau = 1/m.values['l']
    # fitting_end = start + 5 * int(tau/time_interval)
    #--------------------------------------#
    # raw3 = smoothed1.iloc[fitting_end:]
    # raw4 = smoothed2.iloc[fitting_end:]
    cut_index = ch2_first_peak_index[0] + 5 * int(tau / time_interval)
    raw3 = smoothed1.iloc[cut_index:]
    raw4 = smoothed2.iloc[cut_index:]
    #second peak finding 
    ch1_second_peak_index, ch1_second_peak_info = find_peaks(raw3, height= np.mean(raw3) + np.std(raw3) * num_sigma)
    ch2_second_peak_index, ch2_second_peak_info = find_peaks(raw4, height= np.mean(raw4) + np.std(raw4) * num_sigma)

    # ch1_second_peak_index +=fitting_end
    # ch2_second_peak_index += fitting_end
    
    ch1_second_peak_height = ch1_second_peak_info['peak_heights']
    ch2_second_peak_height = ch2_second_peak_info['peak_heights']
    
    if len(ch2_second_peak_index) != 1 :
        print('ch2_second_peak_index is not 1')
        # print(ch2_second_peak_index)
        continue
    
    #---------relevant fake peak------------
    # ch1_Fake_Height = ch1_second_peak_height[0] - relevant_sigma_num * np.std(raw3)
    ch2_Fake_Height = ch2_second_peak_height[0] - relevant_sigma_num * np.std(raw4)
    
    # ch1_second_peak_relevant_index, ch1_second_peak_relevant_info = find_peaks(raw3, height= ch1_Fake_Height, prominence= np.std(raw3) * (num_sigma - relevant_sigma_num))
    ch2_second_peak_relevant_index, ch2_second_peak_relevant_info = find_peaks(raw4, height= ch2_Fake_Height)

    if len(ch2_second_peak_relevant_index) != 1 :
        print('ch2_second_peak_relevant_index is not 1')
        continue
    
    # ch1_second_peak_relevant_index += fitting_end
    ch2_second_peak_relevant_index += cut_index
    
    # ch1_second_peak_height = ch1_second_peak_relevant_info['peak_heights']
    ch2_second_peak_relevant_height = ch2_second_peak_relevant_info['peak_heights']
    #----more info-----------------#
    # plt.scatter(xs, ys, label='CH2 : Smothed', alpha=0.5, s=5)
    peak_index_ch2 = list(peak_index_ch2)
    # plt.scatter([xs[i] for i in peak_index_ch2], [ys[i] for i in peak_index_ch2], color='red', label='CH2 : Peaks')
    plt.legend()
'''
    
    fit_info = [
        f"$\\chi^2$ / $n_\\mathrm{{dof}}$ = {m.fval:.1f} / {len(xs_fit) - m.nfit}",
    ]
    for p, v, e in zip(m.parameters, m.values, m.errors):
        fit_info.append(f"{p} = ${v:.3e} \\pm {e:.3e}$")
    # plt.errorbar(xs, ys, yerr=yerror, fmt="o", label="data",ms= 3, alpha=0.6, color='brown')
    #---------rescale x back------------------
    # ax1 = plt.subplot(1,1,1)
    plt.scatter(xs_fit, ys_fit, label = 'raw', s= 5, alpha= 0.5)
    plt.plot(xs_fit, fitting(xs_fit, a= 15e-7, b= 1e-6, c= 0, l= 1/2e-6, sigma= 1e-6, mu= 2e-7), label = 'theory')
    plt.plot(xs_fit, fitting(xs_temp, *m.values), label="Fitting", linestyle= 'dashed', color = 'red')
    # ax1.legend(loc = 'best', edgecolor = '#7e7474', fontsize = 12)
    plt.legend(title="\n".join(fit_info))
    plt.tight_layout()
    plt.show()
    # plt.savefig('./analysis_result/lifetime/scipy/0.15v/frame' + '/test.jpg')
# %%
