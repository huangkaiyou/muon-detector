#%%
import pandas as pd
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import os
from iminuit import Minuit
from iminuit.cost import LeastSquares
import scienceplots

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
file_path = "./lifetimedatatxt/Lifetime_data(B)/" + file_name + ".txt"
f= open(file= file_path, mode= 'r')
#%%
#----------------------------------------------------------------#
ys = []
left_list_ch1 = []
left_list_ch2 = []
count = 0
Times = 0
frames = 1
# while True:
while Times < frames :
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
    # print(time)
    # #filter first peaks are not happend at a time
    # if abs(ch1_first_peak_index[0] - ch2_first_peak_index[0]) > 20 : 
    #     continue

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
    # cut_index = ch2_first_peak_index + 5 * int(Tau / time_interval)
    cut_index = ch2_first_peak_index + 1
    #-------------find the second peak----------------------
    raw3 = smoothed1.iloc[cut_index:]
    raw4 = smoothed2.iloc[cut_index:]
    # plt.hist(ys)
    # plt.show()
    raw3 = (raw3 - raw3.mean()) / raw3.std()
    raw4 = (raw4 - raw4.mean()) / raw4.std()
    
    left = raw3.to_list()
    left_list_ch1.append(left)
    
    left = raw4.to_list()
    left_list_ch2.append(left)
    
    Times +=1
    # if max(smoothed[:200])> 0.05:
    #     print(file_path)

# %%
# left_list_ch2 = np.array(left_list_ch2)
# left_list_ch2 = left_list_ch2.flatten()

num_event = 0
left_list_ch2 = np.array(left_list_ch2)
left_list_ch2 = left_list_ch2[:200].reshape(-1)
#histogram
log_scale= True
counts, bins, patch = plt.hist(left_list_ch2, bins= 30, alpha=1, color='skyblue', edgecolor='navy', label=  'data', histtype= 'step', log= log_scale)
# plt.ylim((0, 299))
'''
for i in range(1,15):
    if bins[i]<= threshold:
        num_event += counts[i-1]
    else:
        # print("bins:", bins[i])
        break
'''
total_event = len(left_list_ch2)
# print("the ratio: ", num_event/total_event*100, '(%)')
# Add labels and title
plt.xlabel('Sigma')
plt.ylabel('Entries')
plt.title('total data num: ' + str(total_event))
plt.legend()
# Show the plot
# plt.show()  

analysis_figure_dir = './analysis_result/noise/'
plt.savefig(analysis_figure_dir + f'frame_num {Times} (sigma) ' + '.jpg')
plt.clf()



#%%
'''
start = 0
frame = 20
xs = bins[1+start:frame]
ys = counts[start:frame-1]


plt.scatter(xs, ys, label="data",s= 5)
ax1 = plt.subplot(1,1,1)
pltSty2(xName = 'Voltage v', yName = 'Entities')
# ax1.scatter(xs, ys, color = '#070d58', label = 'x(actual)', s = 5)
ax1.legend(loc = 'best', edgecolor = '#7e7474', fontsize = 12)

plt.tight_layout()
plt.draw()
plt.grid(color= 'gray', visible=True)
# plt.xlabel('t (s)', fontsize= 15)
# plt.ylabel('V (volt)', fontsize= 15)
plt.legend()
# plt.savefig('test.png',dpi= 500)
plt.show()


# %%
#-----fitting gaussain distribution-----#
from iminuit import Minuit
from iminuit.cost import LeastSquares
# from OuOb import pltSty2

start = 0
frame = 20
xs = bins[1+start:frame]
ys = counts[start:frame-1]
print()
# yerror= [np.std(i) for i in data]
yerror= 1


def gaussian(x, mu, sigma, N0): #charge 
    return N0 * np.exp(-((x-mu)/sigma)**2 / 2)

least_square= LeastSquares(xs, ys, yerror= yerror, model= gaussian)
m=Minuit(least_square, mu= 0.05, sigma= 0.001, N0= 12e4, )
# plt.errorbar(xs, ys, yerr= stdev(ys)/ , fmt='')
m.migrad() 
m.hesse()
# plt.errorbar(xs, ys, yerr=yerror, fmt="o", label="data",ms= 3, alpha=0.6, color='brown')

plt.errorbar(xs, ys, yerr=yerror, fmt="o", label="data",ms= 2)
plt.plot(xs, gaussian(xs, *m.values), label="gaussian", color= 'blue', linestyle= 'dashed')
#for test
# plt.scatter(xs, ys, label= 'mean',s=2)


fit_info = [
    f"$\\chi^2$ / $n_\\mathrm{{dof}}$ = {m.fval:.1f} / {len(xs) - m.nfit}",
]
for p, v, e in zip(m.parameters, m.values, m.errors):
    fit_info.append(f"{p} = ${v:.3e} \\pm {e:.3e}$")


ax1 = plt.subplot(1,1,1)
pltSty2(xName = 'Voltage v', yName = 'Entities')
# ax1.scatter(xs, ys, color = '#070d58', label = 'x(actual)', s = 5)
ax1.legend(loc = 'best', edgecolor = '#7e7474', fontsize = 12)

plt.tight_layout()
plt.draw()
plt.grid(color= 'gray', visible=True)
# plt.xlabel('t (s)', fontsize= 15)
# plt.ylabel('V (volt)', fontsize= 15)
plt.legend(title="\n".join(fit_info))
# plt.savefig('test.png',dpi= 500)
plt.show()
'''
# %%
