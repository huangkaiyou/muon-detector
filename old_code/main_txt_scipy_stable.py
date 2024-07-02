#%%
import pandas as pd
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import os


### initial parameter ###
threshold_first_ch1 = 2
threshold_first_ch2 = 1
threshold_second = 0.1
stable_level = 0.03
prominence= 0.0075
time_interval = 1e-9
mode = '_stable' + f' ({stable_level})'
# prominice_value = '_' + f'prominence({prominence})'

prominice_value = '_' + 'prominence(sigma)'

file_names = ['20240416_01','20240417_01','20240417_02']
# file_names = ['20240416_01']
#-------ensure the save directorys are exist
'''
ch1_peak_index_file = './analysis_result/lifetime/no B/scipy/' + f'threshold_second ({threshold_second}v)' + mode + prominice_value + '/peak_index/ch1/'
ch1_peak_value_file = './analysis_result/lifetime/no Bscipy/' + f'threshold_second ({threshold_second}v)' + mode + prominice_value + '/peak_value/ch1/'

ch2_peak_index_file = './analysis_result/lifetime/no B/scipy/' + f'threshold_second ({threshold_second}v)' + mode + prominice_value + '/peak_index/ch2/'
ch2_peak_value_file = './analysis_result/lifetime/no B/scipy/' + f'threshold_second ({threshold_second}v)' + mode + prominice_value + '/peak_value/ch2/'

peak_path_file = './analysis_result/lifetime/no Bscipy/' + f'threshold_second ({threshold_second}v)' + mode + prominice_value + '/file_path/'
'''
ch1_peak_index_file = './analysis_result/lifetime/no B/scipy/' + 'threshold_second (mean_all)' + mode + prominice_value + '/peak_index/ch1/'
ch1_peak_value_file = './analysis_result/lifetime/no B/scipy/' + 'threshold_second (mean_all)' + mode + prominice_value + '/peak_value/ch1/'

ch2_peak_index_file = './analysis_result/lifetime/no B/scipy/' + 'threshold_second (mean_all)' + mode + prominice_value + '/peak_index/ch2/'
ch2_peak_value_file = './analysis_result/lifetime/no B/scipy/' + 'threshold_second (mean_all)' + mode + prominice_value + '/peak_value/ch2/'

peak_path_file = './analysis_result/lifetime/no B/scipy/' + 'threshold_second (mean_all)' + mode + prominice_value + '/file_path/'
directorys = [ch1_peak_index_file
              , ch1_peak_value_file
              ,ch2_peak_index_file
              , ch2_peak_value_file
              , peak_path_file
              ]

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
    #directory for frame
    frame_path = './analysis_result/lifetime/no B/scipy/' + 'threshold_second (mean_all)' + mode + prominice_value + f'/frame/{file_name}/' 
    isExist = os.path.exists(frame_path)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(frame_path)
        print(f"The new directory({frame_path}) is created!")
        
    count = 0
    file_path = "lifetimedatatxt/Lifetime_data(no B)" + file_name + ".txt"
    #-------------------------------------------------------#
    ch1_result_index=open(ch1_peak_index_file + f'{file_name}.txt','w')
    ch1_result_value=open(ch1_peak_value_file + f'{file_name}.txt','w')
    
    ch2_result_index=open(ch2_peak_index_file + f'{file_name}.txt','w')
    ch2_result_value=open(ch2_peak_value_file + f'{file_name}.txt','w')
    
    result_path=open(peak_path_file + f'{file_name}.txt','w')
    f= open(file= file_path, mode= 'r')
    #-------------------------------------------------------#
    try:
        while True:
        # for i in range(10000):
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

            #filter two first peak cases
            if (len(ch1_first_peak_index) != 1 or len(ch2_first_peak_index) != 1):
                continue

            ch1_first_peak_height = ch1_first_peak_info['peak_heights'][0]
            ch2_first_peak_height = ch2_first_peak_info['peak_heights'][0]

            #look for smoothed
            for i in range(ch2_first_peak_index[0], 10000):
                if smoothed2[i] < stable_level:
                    second_start = i
                    break
            
            #second peak finding 
            ch1_second_peak_index, ch1_second_peak_info = find_peaks(smoothed1.iloc[second_start:], height= threshold_second, prominence= prominence)
            ch2_second_peak_index, ch2_second_peak_info = find_peaks(smoothed2.iloc[second_start:], height= threshold_second, prominence= prominence)

            ch1_second_peak_index +=second_start
            ch2_second_peak_index += second_start
            
            if len(ch2_second_peak_index) != 1 :
                continue
            
            ch1_second_peak_height = ch1_second_peak_info['peak_heights']
            ch2_second_peak_height = ch2_second_peak_info['peak_heights']
            
            peak_index_ch1 = np.append(ch1_first_peak_index, ch1_second_peak_index)
            peak_index_ch2 = np.append(ch2_first_peak_index, ch2_second_peak_index)
            peak_value_ch1 = np.append(ch1_first_peak_height, ch1_second_peak_height)
            peak_value_ch2 = np.append(ch2_first_peak_height, ch2_second_peak_height)
            
            # Store the peak info
            for i in peak_index_ch1:
                ch1_result_index.write(f'{i}' + '\t')
            ch1_result_index.write("\n")
            
            for i in peak_value_ch1:
                ch1_result_value.write(f'{i}' + '\t')
            ch1_result_value.write("\n")
            
            for i in peak_index_ch2:
                ch2_result_index.write(f'{i}' + '\t')
            ch2_result_index.write("\n")
            
            for i in peak_value_ch2:
                ch2_result_value.write(f'{i}' + '\t')
            ch2_result_value.write("\n")
            
            result_path.write(f'{count}' + '\n')
            
            #show the plot
            ys = smoothed2
            xs = [i*time_interval for i in range(10000)]
            
            lifetime = time_interval * (peak_index_ch2[1] - peak_index_ch2[0])
            lifetime = float(lifetime)
            # print(count)
            plt.title(f"lifetime: {lifetime:.3e} (s)")
            plt.scatter(xs, ys, label='CH2 : Smothed', alpha=0.5, s=5)
            peak_index_ch2 = list(peak_index_ch2)
            plt.scatter([xs[i] for i in peak_index_ch2], [ys[i] for i in peak_index_ch2], color='red', label='CH2 : Peaks')
            plt.legend()
            # plt.show()
            
            plt.tight_layout()
            plt.savefig(frame_path + f'{count}' + '.jpg')
            plt.clf()
    except :      # works on python 3.x
        print("interrupt")
        
    finally:
        print("finish")
        ch1_result_index.close()
        ch2_result_index.close()
        ch1_result_value.close()
        ch2_result_value.close()
        result_path.close()
        f.close()


# yoyo  真的好頂

# %%
