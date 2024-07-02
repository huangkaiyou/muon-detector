#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import scienceplots
import pandas as pd

### initial parameter ###
threshold_first_ch1 = 2
threshold_first_ch2 = 1
stable_level = 0.03
prominence= 0
num_sigma = 10
ground_sigma_num = 5
time_interval = 1e-9

mode = '_stable' + ' (5tau)_fitting' + 'second_no_sigma'
prominice_value = '_' + 'prominence(sigma)'
ground_value = '_' + f'ground({ground_sigma_num})'

plt.rcdefaults
# mode = ''
lifetime_list = []
first_peak_interval = []
#the folder data used
# folder_paths = ["lifetimedata1(V=0.1v)","lifetimedata2(V=0.1v)","lifetimedata3(V=0.1v)"]
# folder_paths = ['20240416_01']
folder_paths = ['20240416_01','20240417_01','20240418_01','20240419_01','20240421_01','20240422_01','20240422_02','20240423_01']
folder_paths = ['20240424', '20240428', '20240429', '20240429_02', '20240430', '20240501', '20240502', '20240504', '20240505', '20240507', '20240510']
ch2_peak_value_lists = [] #store ch2 peak value

analysis_dir = './analysis_result/lifetime/with B/scipy/' + f'threshold_second (mean, num_sigma({num_sigma}))' + mode + prominice_value + ground_value + '/analysis'
for file_name in folder_paths:
    raw = pd.read_csv(analysis_dir + f'/{file_name}' +'.csv')
    try :
        Times = int((raw.shape[1] - 1) / 4)
        # while True:
        for i in range(Times) :            
            #read a new frame info    
            ch1_peak_index_list = raw.iloc[:,4* i + 1].dropna()
            ch2_peak_index_list = raw.iloc[:,4* i + 2].dropna()
            ch1_peak_value_list = raw.iloc[:,4* i + 3].dropna()
            ch2_peak_value_list = raw.iloc[:,4* i + 4].dropna()
            # -----check the first peak--------------#
            
            if abs(ch1_peak_index_list[0] - ch2_peak_index_list[0]) > 5 : 
                continue
            
            if (ch1_peak_value_list[0] > 8) or (ch2_peak_value_list[0] > 8 ) : 
                continue
            
            first_peak_interval.append(ch1_peak_index_list[0] - ch2_peak_index_list[0])
            
            #-----check the second peak--------------#
            second_peak_condition= True
            for i in ch1_peak_index_list:
                if abs(ch2_peak_index_list[1] - i) < 5:
                    second_peak_condition = False
                    break
            if second_peak_condition:
                lifetime= (ch2_peak_index_list[1] - ch2_peak_index_list[0]) * time_interval
                lifetime_list.append(lifetime)
                ch2_peak_value_lists.append(ch2_peak_value_list)
    finally :   
        print(f'{file_name} ' + 'finish!')
        ch2_second_peak_list = [i[1] for i in ch2_peak_value_lists]
        savecsv_path = './analysis_result/lifetime/with B/scipy/' + f'threshold_second (mean, num_sigma({num_sigma}))' + mode + prominice_value + ground_value + '/lifetimedata/'
        isExist = os.path.exists(savecsv_path)
        if not isExist:
            # Create a new directory because it does not exist
            os.makedirs(savecsv_path)
            print(f"The new directory({savecsv_path}) is created!")
        lifetime_list_temp = pd.Series(lifetime_list)
        lifetime_list_temp.transpose()
        pd.DataFrame.to_csv(lifetime_list_temp, savecsv_path + f'{file_name}' + '.csv', header= True)
#%%       
#-----------------------------------------------------------------#
### plot the histogram ###
log_scale = True
(n, bin, patch) = plt.hist(lifetime_list, bins=20, alpha=0.5, color='blue', edgecolor='black', histtype= "step", log= log_scale)
total_event = len(lifetime_list)
# Add labels and title
plt.xlabel('Time')
plt.ylabel('Entities')
plt.title('total event: ' + str(total_event))
# plt.figure(figsize= (15, 10))
# Show the plot
plt.show()



# print((lifetime_list))
# %%
