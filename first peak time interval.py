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
time_interval = 1e-9

mode = '_stable' + ' (5tau)_fitting' + 'second_no_sigma_paper'

file_names = ['20240424', '20240428', '20240429', '20240429_02', '20240430', '20240501', '20240502', '20240504', '20240505', '20240507', '20240510']
# file_names = []
#-------ensure the save directorys are exist

analysis_result_dir = './analysis_result/first_peak_interval/'
# Check whether the specified path exists or not
isExist = os.path.exists(analysis_result_dir)
if not isExist:
    # Create a new directory because it does not exist
    os.makedirs(analysis_result_dir)
    print(f"The new directory({analysis_result_dir}) is created!")
    
#----------------------------------------------------------------#

#%%
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))


for file_name in file_names:            
    count = 0
    file_path = "./lifetimedatatxt/Lifetime_data(B)/" + file_name + ".txt"
    #-------------------------------------------------------#
    f= open(file= file_path, mode= 'r')
    
    final = {}
    df_final = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in final.items() ]))
    peak_interval_temp = []
    #-------------------------------------------------------#
    try:
        while True:
        # for i in range(5000):
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
            ch2_peak_index, ch2_first_peak_info = find_peaks(smoothed2, height= threshold_first_ch2)
            
            #filter zero first peak cases
            if (len(ch1_peak_index) == 0 or len(ch2_peak_index) == 0):
                # print('not you')
                continue
            
            ch1_first_peak_height = ch1_peak_info['peak_heights'][0]
            ch2_first_peak_height = ch2_first_peak_info['peak_heights'][0]

            if (ch1_first_peak_height > 8) or (ch2_first_peak_height > 8): #no too high signal
                continue
            '''
            # -----check the first peak--------------#
            if abs(ch2_first_peak_index[0] - ch1_peak_index[0]) > 1:
                # print('also not you')
                continue
            '''
            
            # #filter first peaks are not happend at a time
            # if abs(ch1_first_peak_index[0] - ch2_first_peak_index[0]) > 20 : 
            #     continue

            peak_interval = ch2_peak_index[0] - ch1_peak_index[0]
            peak_interval_temp.append(peak_interval)
            
            
            
            # print('next')
    # except :     
        # print("interrupt")
        
    finally:

        f.close()
            
        peak_interval_temp = pd.Series(peak_interval_temp)
        peak_interval_temp.transpose()
        pd.DataFrame.to_csv(peak_interval_temp, analysis_result_dir + f'{file_name}' + '.csv', header= True)  
        print(f'{file_name} ' + 'finish!')
# yoyo  真的好頂

# %%
