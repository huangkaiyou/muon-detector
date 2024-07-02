#%%
import pandas as pd
import numpy as np
import os

### initial parameter ###
threshold_first_ch1 = 2
threshold_first_ch2 = 1
time_interval = 1e-9
file_names = ['20240416_01','20240417_01','20240417_02','20240418_01','20240419_01','20240421_01','20240422_01','20240422_02','20240423_01']
#-------ensure the save directorys are exist
directorys = ['./analysis_result/lifetime/peak_index/', './analysis_result/lifetime/peak_value/', './analysis_result/lifetime/file_path/']
for dir in directorys:
    # Check whether the specified path exists or not
    isExist = os.path.exists(dir)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(dir)
        print(f"The new directory({dir}) is created!")
#----------------------------------------------------------------#

#%%
for file_name in file_names:
    count = 0
    
    file_path = "lifetimedatatxt/" + file_name + ".txt"
    result_index=open(f'./analysis_result/lifetime/peak_index/{file_name}.txt','w')
    result_value=open(f'./analysis_result/lifetime/peak_value/{file_name}.txt','w')
    result_frame_index=open(f'./analysis_result/lifetime/file_path/{file_name}.txt','w')
    f= open(file= file_path, mode= 'r')
    try:
        while True:
            #process...
            if (count) % 1000 == 0:
                    print("count: ", count)
            count+= 1

            #channel one data read
            raw1 = f.readline()
            raw1 = raw1.split('\t')
            raw1[-1] = raw1[-1].strip()
            
            ch1 = []
            ch2 = []
            
            for i in raw1:
                ch1.append(float(i))
            #channel two data read
            raw2 = f.readline()
            raw2 = raw2.split('\t')
            raw2[-1] = raw2[-1].strip()

            for i in raw2:
                ch2.append(float(i))

            peak_index = [0]*9
            peak_value = [0]*9
                
            ch1 = pd.Series(ch1)
            ch2 = pd.Series(ch2)
                    
            #exam whether both channel have first signal
            if ((ch1.max() > threshold_first_ch1) and (ch2.max() > threshold_first_ch2)) : 
                
                # Smooth the data with a rolling window to reduce noise
                window_size = 9
                smoothed = ch2.rolling(window_size, center= True, min_periods=1).mean()
                
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
                        if (data.idxmax() != 1200+i*section_num) and (data.idxmax() != 1200+(i+1)*section_num-1):
                            peak_index[i+1] = (data.idxmax())
                            peak_value[i+1] = (data.max())
                
                # Store the peak info
                # print("file path: ",file_path)
                np.savetxt(result_index, peak_index)
                np.savetxt(result_value, peak_value)
                np.savetxt(result_frame_index, [count], fmt= '%i')
    except:
        result_index.close()
        result_frame_index.close()
        result_value.close()
        f.close()
    finally:
        print("finish")
        result_index.close()
        result_frame_index.close()
        result_value.close()
        f.close()
# %%
