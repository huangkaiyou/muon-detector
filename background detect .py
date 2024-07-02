#%%
import pandas as pd
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import os
from iminuit import Minuit
from iminuit.cost import LeastSquares
import scienceplots
import glob

plt.rcParams['figure.figsize'] = (15.0, 10.0) 

### initial parameter ###
threshold_first_ch1 = 2
threshold_first_ch2 = 1
stable_level = 0.03
prominence= 0
num_sigma = 7
time_interval = 1e-9
plt_condition = True
# file_names = ['20240424', '20240428', '20240429', '20240429_02', '20240430', '20240501', '20240502', '20240504', '20240505', '20240507', '20240510']
# file_names = ['20240415' , '20240417_01', '20240417_02','20240418_01','20240419_01', '20240421_01']
# path = r'C:\Users\user\code\python 練習\實驗物理功課\二下小年會\lifetimedatatxt\Lifetime_data(B)'
# file_names = glob.glob( path+ r'\*.txt')
#'20240415' , '20240417_01', '20240417_02','20240418_01','20240419_01', '20240421_01'
times = 10 # how many data we want
mode = 'no B'
#-------ensure the save directorys are exist

if mode == 'with B':
    file_names = ['20240514','20240515','20240516','20240518','20240520',]
    
    analysis_result_dir = './analysis_result/background/with B/no ch1 first peak/' + f'num_sigma({num_sigma})/'
elif mode == 'no B':
    file_names =  ['20240415' , '20240417_01', '20240417_02','20240418_01','20240419_01', '20240421_01']

    analysis_result_dir = './analysis_result/background/no B/no ch1 first peak/' + f'num_sigma({num_sigma})/'

# Check whether the specified path exists or not
isExist = os.path.exists(analysis_result_dir)
if not isExist:
    # Create a new directory because it does not exist
    os.makedirs(analysis_result_dir)
    print(f"The new directory({analysis_result_dir}) is created!")
    
#----------------------------------------------------------------#
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

ch2_first_list= []
ch1_lifetime = []
Times = 0
def fitting(t, a, b, tau): #charge 
    return a * np.exp(-t/ tau) + b

for file_name in file_names:
    
    count = 0
    lifetime_list_temp = []
    if mode == 'with B':
        file_path = "./lifetimedatatxt/Lifetime_data(B)/" + file_name + ".txt"
    elif mode == 'no B':
        file_path = "./lifetimedatatxt/Lifetime_data(no B)/" + file_name + ".txt"
    #-------------------------------------------------------#
    f= open(file= file_path, mode= 'r')
    
    final = {}
    df_final = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in final.items() ]))

    #-------------------------------------------------------#
    try:
        while True:
        # while Times < times :
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
            
            # #filter first peaks are not happend at a time
            # if abs(ch1_first_peak_index[0] - ch2_first_peak_index[0]) > 20 : 
            #     continue

            if len(ch1_first_peak_index) > 0 :
                continue
            
            if len(ch2_first_peak_index) == 0 :
                continue
            
            ch2_first_peak_height = ch2_first_peak_info['peak_heights']
            
            # if (ch1_peak_value_list[0] > 8) or (ch2_peak_value_list[0] > 8):
            
            
            if len(ch2_first_peak_index) == 2:
                if ch2_first_peak_height[0] < ch2_first_peak_height[1]:
                    continue
            elif len(ch2_first_peak_index) == 3:
                if ch2_first_peak_height[0] < ch2_first_peak_height[1] or ch2_first_peak_height[0] < ch2_first_peak_height[2]:
                    print('3 first peak')
                    continue
            elif len(ch2_first_peak_index) == 4:
                if ch2_first_peak_height[0] < ch2_first_peak_height[1] or ch2_first_peak_height[0] < ch2_first_peak_height[2] or ch2_first_peak_height[0] < ch2_first_peak_height[3]:
                    print('4 first peak')
                    continue
                
            ch2_first_peak_index = ch2_first_peak_index[0]
            ch2_first_peak_height = ch2_first_peak_info['peak_heights'][0]
            
            
            #look for smoothed
            for i in range(ch2_first_peak_index, 10000):
                if smoothed2[i] < stable_level:
                    second_start = i
                    break
            
            if i == 9999:
                continue
            
            ys = smoothed2
            xs = [i*time_interval for i in range(10000)]
            
            #look actual smooth with fitting exponential
            #-----------fitting part--------------------------                         
            yerror = 0.007
            start = int((ch2_first_peak_index * 4+second_start * 1)/ 5)
            # start = ch2_first_peak_index[0]
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
            # fitting_end = start + 5 * int(tau/time_interval)
            #--------------------------------------#
            # raw3 = smoothed1.iloc[fitting_end:]
            # raw4 = smoothed1.iloc[fitting_end:]
            cut_index = ch2_first_peak_index + 5 * int(Tau / time_interval)
            raw3 = smoothed1.iloc[cut_index:]
            raw4 = smoothed2.iloc[cut_index:]
            #second peak finding
            ch1_second_peak_index, ch1_second_peak_info = find_peaks(raw3, height= np.mean(raw3) + np.std(raw3) * num_sigma, prominence= np.std(raw3) * 2)
            ch2_second_peak_index, ch2_second_peak_info = find_peaks(raw4, height= np.mean(raw4) + np.std(raw4) * num_sigma, prominence= np.std(raw4) * 2)
           
            ch1_second_peak_height = ch1_second_peak_info['peak_heights']
            ch2_second_peak_height = ch2_first_peak_info['peak_heights']
            if len(ch2_second_peak_index) != 1 :
                # print(ch2_second_peak_index)
                continue
            
            ch1_second_peak_index += cut_index
            ch2_second_peak_index += cut_index
            
            
            peak_index_ch2 = np.append([ch2_first_peak_index], ch2_second_peak_index)
            peak_value_ch2 = np.append([ch2_first_peak_height], ch2_second_peak_height)
            
            
            ch2_first_list.append(ch2_first_peak_height)
            lifetime = (peak_index_ch2[1] - peak_index_ch2[0]) * time_interval
            lifetime_list_temp.append(lifetime)
            Times += 1
            
            #---------rescale x back------------------
            # if plt_condition :
            if False:
                    
                if lifetime > 3e-6:
                    continue
                fit_info = [
                    f"$\\chi^2$ / $n_\\mathrm{{dof}}$ = {m.fval:.1f} / {len(xs_fit) - m.nfit}",
                ]
                for p, v, e in zip(m.parameters, m.values, m.errors):
                    fit_info.append(f"{p} = ${v:.3e} \\pm {e:.3e}$")
                
                
                ax1.scatter(xs, ys, label= 'CH1', s= 8)
                ax1.scatter([xs[i] for i in peak_index_ch1], [ys[i] for i in peak_index_ch1], label= 'Peak', s= 50, alpha= 0.5)
                ax1.plot(xs_fit, fitting(xs_temp, *m.values), label="Fitting", linestyle= 'dashed', color = 'red')
                
                ax2.scatter(xs, smoothed2, label= 'CH2', s= 8)
                ax2.scatter([xs[i] for i in peak_index_ch2], [smoothed2[i] for i in peak_index_ch2], label= 'Peak', s= 50, alpha= 0.5)
                ax2.legend()
                
                ax1.legend(loc = 'best', edgecolor = '#7e7474', fontsize = 12)
                ax1.legend(title="\n".join(fit_info))
                ax1.set_title(f"lifetime: {lifetime:.3e} (s)")
                
                # ax1.tight_layout()
                # plt.show()
                plt.savefig(analysis_fig_dir +'before 3 us/' + f'{count}' + '.jpg')
                ax1.clear()
                ax2.clear()
                
            
            
            
        # print('next')
    # except :     
        # print("interrupt")
        
    finally:
        print("finish")
        f.close()
        lifetime_list_temp = pd.Series(lifetime_list_temp)
        lifetime_list_temp.transpose()
        pd.DataFrame.to_csv(lifetime_list_temp, analysis_result_dir + f'{file_name}' + '.csv', header= True)  

        # pd.DataFrame.to_csv(df_final, analysis_result_dir + f'/{file_name}.csv')


# yoyo  真的好頂

# %%
