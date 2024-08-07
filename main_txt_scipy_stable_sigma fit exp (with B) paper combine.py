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

plt.rcParams['figure.figsize'] = (15.0, 10.0) 

### initial parameter ###
threshold_first_ch1 = 2
threshold_first_ch2 = 1
stable_level = 0.03
prominence= 0
num_sigma = 7
ground_sigma_num = 0
time_interval = 1e-9

mode = '_stable' + ' (5tau)_fitting' + 'second_no_sigma_paper'
prominice_value = '_' + 'prominence(sigma)'
ground_value = '_' + f'ground({ground_sigma_num})'
file_names = ['20240424', '20240428', '20240429', '20240429_02', '20240430', '20240501', '20240502', '20240504', '20240505', '20240507', '20240510'
              ,'20240514', '20240515', '20240516', '20240518', '20240520', '20240524']
# file_names = ['20240524']
#-------ensure the save directorys are exist
analysis_result = './analysis_result/lifetime/with B/scipy/' + f'threshold_second (mean, num_sigma({num_sigma}))' + mode + prominice_value + ground_value 
analysis_result_dir = analysis_result + '/analysis'
analysis_figure_dir = analysis_result + '/frame'
# Check whether the specified path exists or not
isExist = os.path.exists(analysis_result_dir)
if not isExist:
    # Create a new directory because it does not exist
    os.makedirs(analysis_result_dir)
    print(f"The new directory({analysis_result_dir}) is created!")
    
#----------------------------------------------------------------#
lifetime_list = [] # all files lifetimedata
first_peak_interval = []

def fitting(t, a, b, tau): #charge 
    return a * np.exp(-t/ tau) + b

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))


for file_name in file_names:    
    #directory for frame
    # frame2_path = analysis_figure_dir + '/CH2/' + f'/{file_name}/' 
    # for frame_path in [frame2_path]:
    #     isExist = os.path.exists(frame_path)
    #     if not isExist:
    #         # Create a new directory because it does not exist
    #         os.makedirs(frame_path)
    #         print(f"The new directory({frame_path}) is created!")
            
    count = 0
    file_path = "./new_lifetimedatatxt/Lifetime_data(B)/" + file_name + ".txt"
    #-------------------------------------------------------#
    f= open(file= file_path, mode= 'r')
    
    final = {}
    df_final = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in final.items() ]))
    lifetime_list_temp = [] #clear up for next file lifetime
    #-------------------------------------------------------#
    try:
        while True:
        # for i in range(10000):
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
            ch1_second_peak_index, ch1_second_peak_info = find_peaks(raw3, height= np.mean(raw3) + np.std(raw3) * num_sigma, prominence= np.std(raw3) * 2)
            ch2_second_peak_index, ch2_second_peak_info = find_peaks(raw4, height= np.mean(raw4) + np.std(raw4) * num_sigma, prominence= np.std(raw4) * 2)
            
            ch1_second_peak_height = ch1_second_peak_info['peak_heights']
            ch2_second_peak_height = ch2_second_peak_info['peak_heights']
            # print(ch2_second_peak_index)
            if len(ch2_second_peak_index) == 0 :
                # print(ch2_second_peak_index)
                continue
            
            ch1_second_peak_index += cut_index
            ch2_second_peak_index += cut_index
            
            # ch1_second_peak_height = ch1_second_peak_ground_info['peak_heights']
            ch2_second_peak_height = ch2_second_peak_info['peak_heights']
            
            peak_index_ch1 = np.array(ch1_peak_index)
            peak_index_ch2 = np.append(ch2_first_peak_index, ch2_second_peak_index)
            peak_value_ch1 = np.array(ch1_peak_height)
            peak_value_ch2 = np.append(ch2_first_peak_height, ch2_second_peak_height)
            '''
            print("Ch1 peak index: ", peak_index_ch1)
            print("Ch2 peak index: ", peak_index_ch2)
            print("Ch1 peak height: ", peak_value_ch1)
            print("Ch2 peak height: ", peak_value_ch2)
            '''            
            
            #-----check the second peak--------------#
            second_peak_condition= False
            for peak_i in peak_index_ch1:
                if abs(peak_index_ch2[1] - peak_i) <= 2:
                    second_peak_condition = True
                    break
                
            if second_peak_condition:
                continue
            
            
            if (peak_value_ch2[0] <= peak_value_ch2[1]) : #no too high second signal
                continue
            
            lifetime= (peak_index_ch2[1] - peak_index_ch2[0]) * time_interval
            lifetime_list_temp.append(lifetime)
                
            # Store the peak info
            temp_dic = {
                        f'ch1_result_index({count})' : peak_index_ch1
                        , f'ch2_result_index({count})' : peak_index_ch2
                        , f'ch1_result_value({count})' : peak_value_ch1
                        , f'ch2_result_value({count})' : peak_value_ch2              
                    }

            df_temp = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in temp_dic.items() ]))

            df_final = pd.concat([df_final, df_temp], axis= 1, join='outer')
            
            '''
            # Create a figure with subplots and set the figure size
            # fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
            #----more info-----------------#
            ax2.vlines(xs[cut_index], -0.05, peak_value_ch2[0], label= 'cut')
            horizon_index = np.array(xs[cut_index:])
            # ax2.plot(xs[cut_index:], horizon(horizon_index, offset= ch2_Fake_Height), label=f"${ground_sigma_num}\sigma$ : Fake Height", linestyle= 'dashed', color = 'brown', lw= 3)
            
            # ax2.plot(xs[cut_index:], horizon(horizon_index, offset= np.mean(raw4) + np.std(raw4) * 3), label="cut line", linestyle= 'dashed', color = 'brown', lw= 3)
            # ax2.plot(xs[cut_index:], horizon(horizon_index, offset= np.mean(raw4) +np.std(raw4) * 5), label=" $5sigma$", linestyle= 'dashed', color = 'brown', lw= 3)
            # ax2.plot(xs[cut_index:], horizon(horizon_index, offset= np.mean(raw4) +np.std(raw4) * 7), label="cut line", linestyle= 'dashed', color = 'brown', lw= 3)
            ax2.plot(xs[cut_index:], horizon(horizon_index, offset= np.mean(raw4) +np.std(raw4) * num_sigma), label= f"${num_sigma}\sigma$ : Second Height", linestyle= 'dashed', color = 'brown', lw= 3)
            '''
            # ax2.vlines(xs[ch2_first_peak_index[0] + 6 * int(Tau / time_interval)], -0.05, peak_value_ch2[0], label= 'cut')
            # ax2.vlines(xs[ch2_first_peak_index[0] + 8 * int(Tau / time_interval)], -0.05, peak_value_ch2[0], label= 'cut')
            # ax2.vlines(xs[ch2_first_peak_index[0] + 10 * int(Tau / time_interval)], -0.05, peak_value_ch2[0], label= 'cut')
            # ax2.vlines(xs[ch2_first_peak_index[0] + 15 * int(Tau / time_interval)], -0.05, peak_value_ch2[0], label= 'cut')
            '''
            #show the plot            
            lifetime = time_interval * (peak_index_ch2[1] - peak_index_ch2[0])
            lifetime = float(lifetime)
            # print(count)
            ax2.set_title(f"lifetime: {lifetime:.3e} (s)")
            # ax2.scatter(xs, smoothed1, label='CH1 : Smothed', alpha=0.5, s=5)
            ax2.scatter(xs, smoothed2, label='CH2 : Smothed', alpha=0.5, s=5)
            peak_index_ch2 = list(peak_index_ch2)
            ax2.scatter([xs[i] for i in peak_index_ch2], [ys[i] for i in peak_index_ch2], color='red', label='CH2 : Peaks')
            ax2.legend()
            
            
            fit_info = [
                f"$\\chi^2$ / $n_\\mathrm{{dof}}$ = {m.fval:.1f} / {len(xs_fit) - m.nfit}",
            ]
            for p, v, e in zip(m.parameters, m.values, m.errors):
                fit_info.append(f"{p} = ${v:.3e} \\pm {e:.3e}$")
            # ax2.errorbar(xs, ys, yerr=yerror, fmt="o", label="data",ms= 3, alpha=0.6, color='brown')
            #---------rescale x back------------------
            ax2.plot(xs_fit, fitting(xs_temp, *m.values), label="Fitting", linestyle= 'dashed', color = 'red')
            ax2.legend(loc = 'best', edgecolor = '#7e7474', fontsize = 12)
            ax2.legend(title="\n".join(fit_info))
            # ax2.tight_layout()
            
            ax1.scatter(xs, smoothed1, label='CH1 : Smothed', alpha=0.5, s=5)
            ax1.legend()
            # ax1.tight_layout()
            # plt.show()
            
            plt.savefig(frame2_path + f'{count}' + '.jpg')
            ax1.clear()
            ax2.clear()
            '''
            # print('next')
    # except :     
        # print("interrupt")
        
    finally:
        print(f'{file_name} ' + 'finish!')
        f.close()
        pd.DataFrame.to_csv(df_final, analysis_result_dir + f'/{file_name}.csv')
        lifetime_list.extend(lifetime_list_temp)
        
        savecsv_path = analysis_result + '/lifetimedata/'
        isExist = os.path.exists(savecsv_path)
        if not isExist:
            # Create a new directory because it does not exist
            os.makedirs(savecsv_path)
            print(f"The new directory({savecsv_path}) is created!")
        lifetime_list_temp = pd.Series(lifetime_list_temp)
        lifetime_list_temp.transpose()
        pd.DataFrame.to_csv(lifetime_list_temp, savecsv_path + f'{file_name}' + '.csv', header= True)  

# yoyo  真的好頂

# %%
### plot the histogram ###
log_scale = False
(n, bin, patch) = plt.hist(lifetime_list, bins=50, alpha=0.5, color='blue', edgecolor='black', histtype= "step", log= log_scale)
total_event = len(lifetime_list)
# Add labels and title
plt.xlabel('Time')
plt.ylabel('Entities')
plt.title('total event: ' + str(total_event))
# plt.figure(figsize= (15, 10))
# Show the plot
plt.show()
# %%
