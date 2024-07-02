#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob

### initial parameter ###
threshold_first = 1
threshold_second = 0.20
time_interval = 1e-9
lifetime_list = []
peak_list = []
file_path_list = []
folder_paths = ["lifetimedata2(V=0.1v)", 'lifetimedata(V=1v)']
#----------------------------------------------------------------#
for folder_path in folder_paths :
    ### read csv ###
    # Find all CSV files in the folder
    file_paths = glob.glob(folder_path + "/*.CSV")
    num = len(file_paths)
    count = 0
    # Iterate over each file and read it
    for file_path in file_paths:
        is_peak = []
        if count % 500 == 0: 
            print("percent: ", count/num * 100)
        count+=1
        ### read csv ###
        df = pd.read_csv(file_path, header= 25) 
        df = df.iloc[:, :-1]
        df.columns = ["Time1", "Amplitude1", "Time2", "Amplitude2"]

        x = df['Time1']
        y1 = df['Amplitude1']
        y2 = df['Amplitude2']
        #exam whether both channel have first signal
        if ((df['Amplitude1'] > threshold_first).any() and (df['Amplitude2'] > threshold_first).any()) : 
            
            # Smooth the data with a rolling window to reduce noise
            window_size = 9
            smoothed = df['Amplitude2'].rolling(window_size, center= True, min_periods=1).mean()
            data = smoothed.iloc[500:1200]
            if data.max() > threshold_first:
                is_peak.append(data.idxmax())
                first_max = data.max()
            
                data = smoothed.iloc[1200:]
                if (data.max() > threshold_second) and (data.max() < first_max):
                    if (data.max() > data.iloc[data.idxmax()-1]) and (data.max() > data.iloc[data.idxmax()+1]):
                        is_peak.append(data.idxmax())
            # Get the time and amplitude values of the peaks
            peaks = df.iloc[is_peak]
            # print(is_peak)
            peak_index = np.array(is_peak)
            if len(peak_index) >1:
                lifetime= (peak_index[1] - peak_index[0]) * time_interval
                
                lifetime_list.append(lifetime)
                peak_list.append(is_peak)
                file_path_list.append(file_path)
                print("file path: ",file_path)
                plt.scatter(x, df['Amplitude2'], label='Original Waveform', alpha=0.5)
                plt.scatter(x, smoothed, label='Smothed', alpha=0.5, s=5)
                plt.scatter(peaks['Time2'], peaks['Amplitude2'], color='red', label='Peaks')
                plt.legend()
                plt.show()
                
        else:
            # print('no peak')
            pass
            

#-----------------------------------------------------------------#

file_name= folder_path
np.savetxt(f'./analysis_result/lifetime/{file_name}.txt', lifetime_list)
# np.savetxt(f'./peak_data/{file_name}.txt', peak_list)
#------------------------------------------------------------------#

### plot the histogram ###
plt.hist(lifetime_list, bins=30, alpha=0.5, color='blue', edgecolor='black')
total_event = len(lifetime_list)
# Add labels and title
plt.xlabel('Time')
plt.ylabel('Count')
plt.title('total event: ' + str(total_event))

# Show the plot
plt.show()  


# %%
