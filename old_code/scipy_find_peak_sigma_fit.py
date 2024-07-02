#%%
import pandas as pd
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import os
import scienceplots
from iminuit import Minuit
from iminuit.cost import LeastSquares
from matplotlib import ticker
def horizon(x, offset) :
    return 0 * x + offset

def fitting(t, a, b, tau): #charge 
    return a * np.exp(-t/ tau) + b

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
relevant_sigma_num = 0
time_interval = 1e-9
file_name = '20240424'
file_path = "new_lifetimedatatxt/Lifetime_data(B)/" + file_name + ".txt"
f= open(file= file_path, mode= 'r')


ax_start = 0.8e-6
ax_end = 3.1e-6

line = 32596-1
line = 5
for i in range(line):
    f.readline()
    f.readline()
    
#channel 1 data read
raw1 = f.readline()
raw1 = raw1.split('\t')
raw1[-1] = raw1[-1].strip()

#%%
            
# peak1, info1 = find_peaks(smoothed1, height= 2)
# peak2, info2 = find_peaks(smoothed2, height= 0.15, prominence=0.15)

for i in range(1000):
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
    
    
    
    
    
    
    ys = smoothed2
    xs = [i*time_interval for i in range(10000)]
    # start_index = int(ax_start/time_interval)
    # end_index = int(ax_end/time_interval)

    # xs = np.array(xs[start_index:end_index])
    # xs -= xs[0]
    # smoothed1 = smoothed1[start_index:end_index]
    # smoothed2 = smoothed2[start_index:end_index]
    
    
    
    
    
    
    
    #first peak finding 
    ch1_peak_index, ch1_peak_info = find_peaks(smoothed1, height= threshold_first_ch1)
    ch2_peak_index, ch2_peak_info = find_peaks(smoothed2, height= threshold_first_ch2)
    # print(ch2_first_peak_index)

    #filter zero first peak cases
    if (len(ch1_peak_index) == 0 or len(ch2_peak_index) == 0):
        # print('not you')
        continue
    
    # -----check the first peak--------------#
    if abs(ch2_peak_index[0] - ch1_peak_index[0]) > 1:
        # print('also not you')
        continue
    
    # #filter first peaks are not happend at a time
    # if abs(ch1_first_peak_index[0] - ch2_first_peak_index[0]) > 20 : 
    #     continue

    ch1_peak_index = ch1_peak_index
    ch2_first_peak_index = ch2_peak_index[0]
    ch1_peak_height = ch1_peak_info['peak_heights']
    ch2_first_peak_height = ch2_peak_info['peak_heights'][0]
    
    if (ch2_first_peak_height > 8): #no too high signal
        continue
    
    #look for smoothed
    for i in range(ch2_first_peak_index, 10000):
        if smoothed2.iloc[i] < stable_level:
            second_start = i
            break    
    
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
        # print('heheehehe')
        continue
    
    Tau = m.values['tau']
    # fitting_end = start + 5 * int(tau/time_interval)
    #--------------------------------------#
    # raw3 = smoothed1.iloc[fitting_end:]
    # raw4 = smoothed2.iloc[fitting_end:]
    cut_index = ch2_first_peak_index + 5 * int(Tau / time_interval)
    raw3 = smoothed1.iloc[cut_index:]
    raw4 = smoothed2.iloc[cut_index:]
    #second peak finding
    ch1_second_peak_index, ch1_second_peak_info = find_peaks(raw3, height= np.mean(raw3) + np.std(raw3) * num_sigma, prominence= np.std(raw3) * (num_sigma-3))
    ch2_second_peak_index, ch2_second_peak_info = find_peaks(raw4, height= np.mean(raw4) + np.std(raw4) * num_sigma, prominence= np.std(raw4) * (num_sigma-3))
    
    ch1_second_peak_height = ch1_second_peak_info['peak_heights']
    ch2_second_peak_height = ch2_second_peak_info['peak_heights']
    
    if len(ch2_second_peak_index) != 1 :
        # print(ch2_second_peak_index)
        continue
    
    #---------relevant fake peak------------
    # ch1_Fake_Height = ch1_second_peak_height[0] - relevant_sigma_num * np.std(raw3)
    ch2_Fake_Height = ch2_second_peak_height[0] - relevant_sigma_num * np.std(raw4)
    
    # ch1_second_peak_relevant_index, ch1_second_peak_relevant_info = find_peaks(raw3, height= ch1_Fake_Height, prominence= np.std(raw3) * (num_sigma - relevant_sigma_num))
    ch2_second_peak_relevant_index, ch2_second_peak_relevant_info = find_peaks(raw4, height= ch2_Fake_Height, prominence= np.std(raw4) * 1)

    if len(ch2_second_peak_relevant_index) != 1 :
        # print('there is fake peak', 'count: ', count)
        continue
    
    ch1_second_peak_index += cut_index
    ch2_second_peak_index += cut_index
    
    # ch1_second_peak_height = ch1_second_peak_relevant_info['peak_heights']
    ch2_second_peak_height = ch2_second_peak_info['peak_heights']
    
    peak_index_ch1 = np.append(ch1_peak_index, ch1_second_peak_index)
    peak_index_ch2 = np.append(ch2_first_peak_index, ch2_second_peak_index)
    peak_value_ch1 = np.append(ch1_peak_height, ch1_second_peak_height)
    peak_value_ch2 = np.append(ch2_first_peak_height, ch2_second_peak_height)
    
    fit_info = [
        f"$\\chi^2$ / $n_\\mathrm{{dof}}$ = {m.fval:.1f} / {len(xs_fit) - m.nfit}",
    ]
    for p, v, e in zip(m.parameters, m.values, m.errors):
        fit_info.append(f"{p} = ${v:.3e} \\pm {e:.3e}$")
    # plt.errorbar(xs, ys, yerr=yerror, fmt="o", label="data",ms= 3, alpha=0.6, color='brown')
    #---------rescale x back------------------
    # Create a figure with subplots and set the figure size
    
    fig, ax = plt.subplots(1, 1, figsize=(10, 5)) 
    #----more info-----------------#
    
    filter = ch2_second_peak_height[0] - 5 * np.std(raw4)      
    lifetime = time_interval * (peak_index_ch2[1] - peak_index_ch2[0])
    lifetime = float(lifetime)
    # print(count)
    
    plt.title(f"lifetime: {lifetime:.3e} (s)")
    # plt.scatter(xs, ys, label='CH2 : Smothed', alpha=0.5, s=5)
    peak_index_ch2 = list(peak_index_ch2)
    ax.scatter([xs[i] for i in peak_index_ch2], [ys[i] for i in peak_index_ch2], color='red', label='CH2 : Peaks', s=80)
    #----more info-----------------#
    # ax.scatter(xs, smoothed1, label = 'CH1', color = 'black', s=5)    
    ax.errorbar(xs, smoothed2, yerr= yerror, fmt= '.', color= 'black')
    ax.plot(xs_fit, fitting(xs_temp, *m.values), label= 'fit', color = 'red', lw= 4)
    # ax.scatter(xs[start], smoothed2.iloc[start], label = 'Fit_start', s= [100], alpha= 0.5)
    # ax.scatter(xs[second_start], smoothed2[second_start], label = 'stable', s= 80, alpha= 0.5)
    # ax.scatter(xs, smoothed2, label = 'CH2', color = 'black', s=5)
    # ax.vlines(xs[cut_index], -0.05, peak_value_ch2[0], label= 'cut')
    
    ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.2f}"))
    ax.set_xlabel(r'Time $(\mu s)$', loc='right')
    ax.set_ylabel('Voltage (v)', loc= 'top')
    
    # ax2.set_ylabel('Voltage (v)')
    ax.legend(loc = 'best', fontsize = 15)
    ax.legend(title="\n".join(fit_info))
    
    xtick = np.arange(int(ax_start * 1e6), int(ax_end * 1e6), 1)
    ax.set_xticks(xtick * 1e-6)
    ax.set_xticklabels(xtick)
    ax.set_xlim(ax_start, ax_end)
    
    plt.show()
    # plt.savefig('./analysis_result/lifetime/scipy/0.15v/frame' + '/test.jpg')
# %%
