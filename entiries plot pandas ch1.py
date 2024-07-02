#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import scienceplots
import pandas as pd
from iminuit import Minuit
from iminuit.cost import LeastSquares
from matplotlib.colors import LogNorm
### initial parameter ###
threshold_first_ch1 = 2
threshold_first_ch2 = 1
stable_level = 0.03
threshold_second = 0.1
prominence= 0
num_sigma = 7
ground_sigma_num = 0
time_interval = 1e-9

lifetime_start = 1e-6
lifetime_end = 10e-6
Count = 0


mode = '_stable' + ' (5tau)_fitting' + 'second_no_sigma_paper'
prominice_value = '_' + 'prominence(sigma)'
ground_value = '_' + f'ground({ground_sigma_num})'

plt.rcdefaults
# mode = ''
lifetime_list = [] # all files lifetimedata
value_list = []
first_peak_interval = []
#the folder data used
# folder_paths = ["lifetimedata1(V=0.1v)","lifetimedata2(V=0.1v)","lifetimedata3(V=0.1v)"]
# folder_paths = ['20240416_01']
# folder_paths = ['20240415','20240416_01','20240417_01','20240417_02','20240418_01','20240419_01','20240421_01','20240422_01','20240422_02','20240423_01'
#                 ,'20240511','20240513','20240521_01','20240522']
folder_paths = [ '20240424', '20240428', '20240429', '20240429_02', '20240430', '20240501', '20240502', '20240504', '20240505'
                , '20240507', '20240510', '20240514', '20240515', '20240516', '20240518', '20240520', '20240524']
# folder_paths = ['20240505']
ch2_peak_value_lists = [] #store ch2 peak value
analysis_result = './analysis_result/lifetime/with B/scipy/' + f'threshold_second (mean, num_sigma({num_sigma}))' + mode + prominice_value + ground_value 

# analysis_result = './analysis_result/lifetime/with B/scipy/' + f'threshold_second ({threshold_second})' + mode + prominice_value + ground_value

analysis_result_dir = analysis_result + '/backup/analysis'

for file_name in folder_paths:
    
    raw = pd.read_csv(analysis_result_dir + f'/{file_name}' +'.csv')
    Times = int((raw.shape[1] - 1) / 4)
    lifetime_list_temp = [] #clear up for next file lifetime
    # while True:
    for time in range(Times) :            
        #read a new frame info    
        ch1_peak_index_list = raw.iloc[:,4* time + 1].dropna()
        ch2_peak_index_list = raw.iloc[:,4* time + 2].dropna()
        ch1_peak_value_list = raw.iloc[:,4* time + 3].dropna()
        ch2_peak_value_list = raw.iloc[:,4* time + 4].dropna()
        
        # if (len(ch2_peak_index_list) != 2) and (ch2_peak_index_list[1] != ch2_peak_index_list[2]):
        #     print((time))
        
            
        # # -----check the first peak--------------#
        # if abs(ch1_peak_index_list[0] - ch2_peak_index_list[0]) > 1:
        #     continue
        
        # if (ch1_peak_value_list[0] > 8) or (ch2_peak_value_list[0] > 8):
        if (ch2_peak_value_list[0] > 8):
            Count +=1
            continue
        
        if (ch2_peak_value_list[0] <= ch2_peak_value_list[1]) :
            continue
        
        
        first_peak_interval.append(ch1_peak_index_list[0] - ch2_peak_index_list[0])
        
        
        #-----check the second peak--------------#
        second_peak_condition= False
        for i in ch1_peak_index_list:
            if abs(ch2_peak_index_list[1] - i) <= 2:
                second_peak_condition = True
                # print('noooooooooo')
                break
        if second_peak_condition:
            continue
        
        if len(ch1_peak_index_list) <= 2:
            continue
        lifetime = (ch1_peak_index_list[1]- ch1_peak_index_list[0]) * time_interval
        if (lifetime < lifetime_start)  or  (lifetime >lifetime_end) :
            continue

        lifetime_list_temp.append(lifetime)
        # print(ch1_peak_index_list)
        ch2_peak_value_lists.append(ch2_peak_value_list)
        # value_list.append(ch2_peak_value_list[1])
    print(f'{file_name} ' + 'finish!')
    lifetime_list.extend(lifetime_list_temp)
    # print(len(lifetime_list_temp))
    ch2_second_peak_list = [i[1] for i in ch2_peak_value_lists]
    
    savecsv_path = './analysis_result/background/no B/ch1' + '/lifetimedata/'
    isExist = os.path.exists(savecsv_path)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(savecsv_path)
        print(f"The new directory({savecsv_path}) is created!")
    lifetime_list_temp = pd.Series(lifetime_list_temp)
    lifetime_list_temp.transpose()
    # pd.DataFrame.to_csv(lifetime_list_temp, savecsv_path + f'{file_name}' + '.csv', header= True)   
    
#-----------------------------------------------------------------#
#%%
### plot the histogram ###

log_scale = True
(n, bin, patch) = plt.hist(lifetime_list, bins=np.linspace(0, 10e-6, 50, endpoint= False), alpha=0.5, color='blue', edgecolor='black', histtype= "step", log= log_scale)
# total_event = len(lifetime_list)
# plt.hist2d(lifetime_list, value_list, bins=50, alpha=1, cmap= 'Blues', norm=LogNorm(), range= [[0., 10e-6], [0, 1]])
total_event = len(lifetime_list)
# total_event = len(value_list)
print(n)
plt.clf()
bin_centers = (bin[:-1] + bin[1:]) / 2
# plt.plot(bin_centers, n,alpha=0.5, label = 'prominence: 2 sigma')
plt.scatter(bin_centers, n, alpha=0.5, label = 'prominence: 2 sigma')
# plt.colorbar()
# Add labels and title
plt.xlabel('time')
plt.ylabel('Entries')
# plt.semilogy()
plt.title('total event: ' + str(total_event) + "/"
        + 'lifetime_range: ' + f'{lifetime_start}' + '--' + f'{lifetime_end}'
        ,fontsize = 20)
# plt.figure(figsize= (15, 10))
# Show the plot
plt.yscale('log')
plt.show()
# plt.savefig(fname = r'C:\Users\user\Downloads\muon\second_peak_dis_new\0.01v' 
            # + f'/{lifetime_start}' + '--' + f'{lifetime_end}'
            # + '.png')
# print((lifetime_list))

# %%
