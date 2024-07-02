#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob

### initial parameter ###
threshold = 1
# Define the noise range to set to 0
lower_bound = -0.02
upper_bound = 0.2
time_interval = 10e-9
lifetime = -1
lifetime_list = []
peak_list = np.array([])
folder_path = "lifetimedata"
Condition_two_peak = True #consider the logic of two scintilletor

dead_time = int(100e-9/time_interval)
#----------------------------------------------------------------#
### read csv ###
# Find all CSV files in the folder
file_paths = glob.glob(folder_path + "/*.CSV")
is_peak= []
#%%
try:
    # Iterate over each file and read it
    for file_path in file_paths:
        peak_index = []
        ### read csv ###
        df = pd.read_csv(file_path, header= 25) 
        df = df.iloc[:, :-1]
        df.columns = ["Time1", "Amplitude1", "Time2", "Amplitude2"]

        x = df['Time1']
        y1 = df['Amplitude1']
        y2 = df['Amplitude2']
        
        if ((df['Amplitude1'] > threshold).any() and (df['Amplitude2'] > threshold).any()) and Condition_two_peak : 
            ## find the second peak ###
            # plt.scatter(x, df['Amplitude2'], label= "raw", s=20, alpha=1)

            
            # Smooth the data with a rolling window to reduce noise
            window_size = 9
            smoothed = df['Amplitude2'].rolling(window_size, center= True, min_periods=1).mean()

            for i in range(10, len(smoothed)-10):
                if (smoothed.iloc[i] > smoothed.iloc[i-5]) and (smoothed.iloc[i] > smoothed.iloc[i+5]):
                    is_peak.append(i)
            # Get the time and amplitude values of the peaks
            peaks = df.iloc[is_peak]
            print(is_peak)
            peak_index = np.array(is_peak)
            
            plt.scatter(x, df['Amplitude2'], label='Original Waveform', alpha=0.5)
            # print("file path: ", file_path, " | ", "peak_index: ", peak_index)
            plt.scatter(peaks['Time2'], peaks['Amplitude2'], color='red', label='Peaks')
            plt.legend()
            plt.show()
            if len(peak_index) > 1:
                lifetime= (peak_index[1] - peak_index[0]) * time_interval
                
                # Plot the peaks
                # print(lifetime)
                # print(peak_index)
                # plt.scatter(x, df['Amplitude2'], label='Original Waveform')
                # plt.scatter(peaks['Time2'], peaks['Amplitude2'], color='red', label='Peaks')
                # plt.legend()
                # plt.show()
        else:
            print('no peak')
        
        lifetime_list.append(lifetime)
        peak_list = np.append(peak_list, peak_index)
finally:        
    #-----------------------------------------------------------------#
    file_name= folder_path
    # np.savetxt(f'.\lifetime_data\{file_name}.txt', lifetime_list)
    np.savetxt(f'./peak_data/{file_name}.txt', peak_list)
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

