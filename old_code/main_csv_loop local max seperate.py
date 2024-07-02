#%%
import pandas as pd
import numpy as np
import glob
import os

### initial parameter ###
threshold_first = 1
time_interval = 1e-9
folder_paths = ["lifetimedata1(V=0.1v)"]

directorys = ['./analysis_result/lifetime/peak_index/', './analysis_result/lifetime/peak_value/', './analysis_result/lifetime/file_path/']
for dir in directorys:
    # Check whether the specified path exists or not
    isExist = os.path.exists(dir)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(dir)
        print(f"The new directory({dir}) is created!")
#%%
#----------------------------------------------------------------#
for folder_path in folder_paths :
    result_index=open(f'./analysis_result/lifetime/peak_index/{folder_path}.txt','w')
    result_value=open(f'./analysis_result/lifetime/peak_value/{folder_path}.txt','w')
    result_file_path=open(f'./analysis_result/lifetime/file_path/{folder_path}.txt','w')

    ### read csv ###
    # Find all CSV files in the folder
    file_paths = glob.glob(folder_path + "/*.CSV")
    num = len(file_paths)
    count = 0
    # Iterate over each file and read it
    for file_path in file_paths[:100]:
        peak_index = [0] * 9
        peak_value = [0] * 9
        if count % 500 == 0: 
            print("percent: ", count/num * 100)
        count+=1
        ### read csv ###
        df = pd.read_csv(file_path, header= 25) 
        df = df.iloc[:, :-1]
        df.columns = ["Time1", "Amplitude1", "Time2", "Amplitude2"]

        #exam whether both channel have first signal
        if ((df['Amplitude1'] > threshold_first).any() and (df['Amplitude2'] > threshold_first).any()) : 
            
            # Smooth the data with a rolling window to reduce noise
            window_size = 9
            smoothed = df['Amplitude2'].rolling(window_size, center= True, min_periods=1).mean()
            
            #first peak finding 
            data = smoothed.iloc[0:1200]
            peak_index[0] = (data.idxmax())
            first_max = data.max()
            peak_value[0] = (first_max)
            
            #second peak judgement 
            section_num = 1100
            section = int((10000-1200)/section_num) #slice multiply section
            for i in range(section):
                data = smoothed.iloc[1200+i*section_num:1200+(i+1)*section_num]
                if (data.max() < first_max):
                    if data.idxmax() != 1200+i*section_num:
                        peak_index[i+1] = (data.idxmax())
                        peak_value[i+1] = (data.max())
            # # Get the time and amplitude values of the peaks
            # peaks = df.loc[peak_index]
            # print(peak_index)
            print("file path: ",file_path)
            
            file_name = folder_path
            np.savetxt(result_index, peak_index)
            np.savetxt(result_value, peak_value)
            np.savetxt(result_file_path, [file_path], fmt= "%s")
    
    result_file_path.close()
    result_index.close()
    result_value.close()  
    
print("finish!")

# %%
