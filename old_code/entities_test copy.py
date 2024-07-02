#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#the folder data used
# folder_paths = ["lifetimedata1(V=0.1v)","lifetimedata2(V=0.1v)","lifetimedata3(V=0.1v)"]
folder_paths = ['lifetimedata1(V=0.1v)']
peak_index_list = np.zeros((9))
peak_value_list = np.zeros((9))
file_path_list = np.zeros(1)
for file_name in folder_paths:
    peak_index_list_temp = np.loadtxt(f'./analysis_result/lifetime/peak_index/{file_name}.txt')
    peak_value_list_temp = np.loadtxt(f'./analysis_result/lifetime/peak_value/{file_name}.txt')
    file_path_list_temp = np.loadtxt(f'./analysis_result/lifetime/file_path/{file_name}.txt', dtype= "str")
    
    peak_index_list = np.append(peak_index_list, peak_index_list_temp)
    peak_value_list = np.append(peak_value_list, peak_value_list_temp)
    file_path_list = np.append(file_path_list, file_path_list_temp)

### initial parameter ###
threshold_second = 0.15
time_interval = 1e-9
lifetime_list = []
#%%
#----------------------------------------------------------------#
for i in range(1,int(len(peak_index_list)/9)):
    peak_index = peak_index_list[i*9:(i+1)*9]
    peak_value = peak_value_list[i*9:(i+1)*9]
    file_path = file_path_list[i]
    
    peak_count = len([x for x in peak_value if x > threshold_second])
    if peak_count == 2:
        for i in range(1, 9):
            if peak_value[i] > threshold_second:
                # print("which is: ", num)
                # print("peak_indexs: ", peak_index)
                # print("peak_values: ", peak_value)
                # print("file path: ",file_path)
                lifetime= (peak_index[i] - peak_index[0]) * time_interval
                lifetime_list.append(lifetime)
                
                if lifetime > 2e-6:
                    print(lifetime)
                    #--------print the waveform------------#
                    ### read csv ###
                    df = pd.read_csv(file_path, header= 25) 
                    df = df.iloc[:, :-1]
                    df.columns = ["Time1", "Amplitude1", "Time2", "Amplitude2"]

                    # Smooth the data with a rolling window to reduce noise
                    window_size = 9
                    smoothed = df['Amplitude2'].rolling(window_size, center= True, min_periods=1).mean()
                    '''
                    plt.title(f"lifetime: {lifetime:.3e} (s)")
                    plt.scatter(df['Time2'], df['Amplitude2'], label='Original Waveform', alpha=0.5)
                    plt.scatter(df['Time2'], smoothed, label='Smothed', alpha=0.5, s=5)
                    plt.scatter(df['Time2'].loc[peak_index], df['Amplitude2'].loc[peak_index], color='red', label='Peaks')
                    plt.legend()
                    plt.show()
                    '''
#%%       
#-----------------------------------------------------------------#
### plot the histogram ###
(n, bin, patch) = plt.hist(lifetime_list, bins=30, alpha=0.5, color='blue', edgecolor='black', histtype= "step", log= True)
total_event = len(lifetime_list)
# Add labels and title
plt.xlabel('Time')
plt.ylabel('Entities')
plt.title('total event: ' + str(total_event))
# Show the plot
plt.show()  
# %%
# print((lifetime_list))
# %%
