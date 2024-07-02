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
num_sigma = 10
relevant_sigma_num = 3
time_interval = 1e-9

mode = '_stable' + ' (5tau)_fitting' + 'second_no_sigma'
prominice_value = '_' + 'prominence(sigma)'
relevant_value = '_' + f'relevant({relevant_sigma_num})'
# file_names = ['20240416_01']
file_names = ['20240416_01']
# '20240416_01','20240417_01','20240417_02', '20240418_01','20240419_01','20240421_01','20240422_01','20240422_02','20240423_01'

analysis_result_dir = './analysis_result/lifetime/no B/scipy/' + f'threshold_second (mean, num_sigma({num_sigma}))' + mode + prominice_value + relevant_value + '/analysis'
analysis_figure_dir = './analysis_result/lifetime/no B/scipy/' + f'threshold_second (mean, num_sigma({num_sigma}))' + mode + prominice_value + relevant_value + '/frame'
# Check whether the specified path exists or not
isExist = os.path.exists(analysis_result_dir)
if not isExist:
    # Create a new directory because it does not exist
    os.makedirs(analysis_result_dir)
    print(f"The new directory({analysis_result_dir}) is created!")
    
#----------------------------------------------------------------#

#%%
def fitting(t, a, b, tau): #charge 
    return a * np.exp(-t/ tau) + b

for file_name in file_names:
    
    #directory for frame
    frame_path = analysis_figure_dir + f'/{file_name}/' 
    isExist = os.path.exists(frame_path)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(frame_path)
        print(f"The new directory({frame_path}) is created!")
        
    count = 0
    file_path = "./lifetimedatatxt/Lifetime_data(no B)/" + file_name + ".txt"
    #-------------------------------------------------------#
    f= open(file= file_path, mode= 'r')
    
    final = {}
    df_final = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in final.items() ]))

    #-------------------------------------------------------#
    try:
        # while True:
        for i in range(1000):
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

            #filter zero first peak cases
            if (len(ch1_first_peak_index) == 0 or len(ch2_first_peak_index) == 0):
                continue
            
            # #filter first peaks are not happend at a time
            # if abs(ch1_first_peak_index[0] - ch2_first_peak_index[0]) > 20 : 
            #     continue

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
            # raw4 = smoothed2.iloc[fitting_end:]
            cut_index = ch2_first_peak_index[0] + 5 * int(Tau / time_interval)
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
            ch2_second_peak_relevant_index, ch2_second_peak_relevant_info = find_peaks(raw4, height= ch2_Fake_Height)

            if len(ch2_second_peak_relevant_index) != 1 :
                # print(ch2_second_peak_index)
                continue
            
            ch1_second_peak_index += cut_index
            ch2_second_peak_index += cut_index
            
            # ch1_second_peak_height = ch1_second_peak_relevant_info['peak_heights']
            ch2_second_peak_height = ch2_second_peak_info['peak_heights']
            
            peak_index_ch1 = np.append(ch1_first_peak_index, ch1_second_peak_index)
            peak_index_ch2 = np.append(ch2_first_peak_index, ch2_second_peak_index)
            peak_value_ch1 = np.append(ch1_first_peak_height, ch1_second_peak_height)
            peak_value_ch2 = np.append(ch2_first_peak_height, ch2_second_peak_height)
            
            # Store the peak info
            temp_dic = {
                        f'ch1_result_index({count})' : peak_index_ch1
                        , f'ch2_result_index({count})' : peak_index_ch2
                        , f'ch1_result_value({count})' : peak_value_ch1
                        , f'ch2_result_value({count})' : peak_value_ch2                 
                    }

            df_temp = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in temp_dic.items() ]))

            df_final = pd.concat([df_final, df_temp], axis= 1, join='outer')
            
            #----more info-----------------#
            plt.vlines(xs[cut_index], -0.05, peak_value_ch2[0], label= 'cut')
            filter = ch2_second_peak_height[0] - 5 * np.std(raw4)
            horizon_index = np.array(xs[cut_index:])
            plt.plot(xs[cut_index:], horizon(horizon_index, offset= filter), label="cut line", linestyle= 'dashed', color = 'brown', lw= 5)
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
            # print('next')
    # except :      # works on python 3.x
        # print("interrupt")
        
    finally:
        print("finish")
        f.close()
        pd.DataFrame.to_csv(df_final, analysis_result_dir + f'/{file_name}.csv')


# yoyo  真的好頂

# %%
