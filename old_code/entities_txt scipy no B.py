#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import scienceplots

# plt.style.use(['science', 'nature'])
### initial parameter ###
threshold_first_ch1 = 2
threshold_first_ch2 = 1
threshold_second = 0.1
stable_level = 0.05
prominence= 0.05
time_interval = 1e-9
mode = '_stable' + f' ({stable_level})'
prominice_value = '_' + f'prominence({prominence})'

# mode = ''
lifetime_list = []
first_peak_interval = []
#the folder data used
# folder_paths = ["lifetimedata1(V=0.1v)","lifetimedata2(V=0.1v)","lifetimedata3(V=0.1v)"]
folder_paths = ['20240416_01','20240417_01','20240417_02', '20240418_01','20240419_01','20240421_01','20240422_01','20240422_02','20240423_01']
# folder_paths = [folder_paths_list[0]]
ch2_peak_value_lists = [] #store ch2 peak value
for file_name in folder_paths:
    count=0
    ch1_peak_index_file = open('./analysis_result/lifetime/no B/scipy/' + f'threshold_second ({threshold_second}v)' + mode + prominice_value + '/peak_index/ch1/' + f'{file_name}' + '.txt', 'r')
    ch2_peak_index_file = open('./analysis_result/lifetime/no B/scipy/' + f'threshold_second ({threshold_second}v)' + mode + prominice_value + '/peak_index/ch2/' + f'{file_name}' + '.txt', 'r')
    ch1_peak_value_file = open('./analysis_result/lifetime/no B/scipy/' + f'threshold_second ({threshold_second}v)' + mode + prominice_value + '/peak_value/ch1/' + f'{file_name}' + '.txt', 'r')
    ch2_peak_value_file = open('./analysis_result/lifetime/no B/scipy/' + f'threshold_second ({threshold_second}v)' + mode + prominice_value + '/peak_value/ch2/' + f'{file_name}' + '.txt', 'r')
    
    # file_path_list_temp = open(f'./analysis_result/lifetime/file_path/{file_name}.txt', 'r')
    try:
        while True:
            count+=1
            #read a new frame info
            ch1_peak_index_list = []
            ch2_peak_index_list = []
            ch1_peak_value_list = []
            ch2_peak_value_list = []
            #index_read
            raw1_index = ch1_peak_index_file.readline()
            
            if raw1_index == '':
                break
            
            raw1 = raw1_index.split('\t')
            for i in raw1[:-1]:
                ch1_peak_index_list.append(float(i))
                
            raw2_index = ch2_peak_index_file.readline()
            raw2 = raw2_index.split('\t')
            for i in raw2[:-1]:
                ch2_peak_index_list.append(float(i))
            
            #value_reade
            raw1_value = ch1_peak_value_file.readline()
            raw1 = raw1_value.split('\t')
            for i in raw1[:-1]:
                ch1_peak_value_list.append(float(i))
                
            raw2_value = ch2_peak_value_file.readline()
            raw2 = raw2_value.split('\t')
            for i in raw2[:-1]:
                ch2_peak_value_list.append(float(i))
            # -----check the first peak--------------#
            if abs(ch1_peak_index_list[0] - ch2_peak_index_list[0]) >5:
                continue
            
            if (ch1_peak_value_list[0] > 8) and (ch2_peak_value_list[0] > 8):
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
           
            
    # except:
    #     print('interrupted')
    #     ch1_peak_index_file.close()
    #     ch2_peak_index_file.close()
    #     ch1_peak_value_file.close()
    #     ch2_peak_value_file.close()
    finally:
        print('finish!')
        ch1_peak_index_file.close()
        ch2_peak_index_file.close()
        ch1_peak_value_file.close()
        ch2_peak_value_file.close()
        
        ch2_second_peak_list = [i[1] for i in ch2_peak_value_lists]
        savecsv_path = './analysis_result/lifetime/scipy/' + f'threshold_second ({threshold_second}v)' + mode + prominice_value + '/lifetimedata/'
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
(n, bin, patch) = plt.hist(lifetime_list, bins=30, alpha=0.5, color='blue', edgecolor='black', histtype= "step", log= log_scale)
total_event = len(lifetime_list)
# Add labels and title
plt.xlabel('Time')
plt.ylabel('Entities')
plt.title('total event: ' + str(total_event))
# Show the plot
plt.show()



# print((lifetime_list))
# %%
