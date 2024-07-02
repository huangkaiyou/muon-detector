import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob

### initial parameter ###
folder_name = 'testdata4 (300mv)'

threshold = 250e-3
Peak2_voltage = 250e-3
time_interval = 1e-9
lifetime = -1 #default lifetime
lifetime1 = lifetime
lifetime2 = lifetime

Condition_two_peak = True #consider the logic of two scintilletor

lifetime_Col1 = []
lifetime_Col2 = []

# Find all CSV files in the folder
file_paths = glob.glob(folder_name + "/*.csv")

# Iterate over each file and read it
for file_path in file_paths:
    # print(file_path)
    try:
        df = pd.read_csv(file_path, header= 25) 
    except: 
        continue
    df = df.iloc[:, :-1]
    df.columns = ["Time1", "Waveform1", "Time2", "Waveform2"]

    x = df['Time1'].to_numpy()
    y1 = df['Waveform1'].to_numpy()
    y2 = df['Waveform2'].to_numpy()

    
    raw = [y1, y2]
    # print(raw)            
    if ((raw[0] > threshold).any() and (raw[1] > threshold).any()) and Condition_two_peak : 
        ### see which file has peak at both channels ###
        # print(file_path) 
        
        ## find the second peak ###
        #Channel 1
        i = 0
        while i < len(raw[0]):
            if raw[0][i] <= Peak2_voltage:
                i += 1
                continue
        
            j = i
            while j < len(raw[0]) and raw[0][j] >= 0.1:
                j += 1
            
            if j == len(raw[0]):
                break
            
            k = j
            while k < len(raw[0]) and raw[0][k] <= Peak2_voltage:
                k += 1
                
            lifetime1 = (k - i) * time_interval
            # print("k:", k, "i: ", i)
            # print('ok')
            break
        
        #Channel 2
        i = 0
        while i < len(raw[1]):
            if raw[1][i] <= Peak2_voltage:
                i += 1
                continue
        
            j = i
            while j < len(raw[1]) and raw[1][j] >= 0.1:
                j += 1
            
            if j == len(raw[1]):
                break
            
            k = j
            while k < len(raw[1]) and raw[1][k] <= Peak2_voltage:
                k += 1

            lifetime1 = (k - i) * time_interval
            # print('ok')
            break
    
    ### Store the lifetime to a list ###
    lifetime_Col1.append(lifetime1)
    lifetime_Col2.append(lifetime2)
    ### initialize the lifetime
    lifetime1 = lifetime
    lifetime2 = lifetime   
    
#-----------------------------------------------------------------#
### save the lifetime list to csv file ###
# Convert the list to a DataFrame
data= {
    "CH1":lifetime_Col1,
    "CH2":lifetime_Col2
    }
df = pd.DataFrame(data)

# Save the DataFrame to a CSV file
df.to_csv(f'{folder_name}.csv', index=False)
#------------------------------------------------------------------#

### plot the histogram ###
plt.hist(lifetime_Col2, bins=30, alpha=0.5, color='blue', edgecolor='black')
total_event = len(lifetime_Col1)
# Add labels and title
plt.xlabel('Time')
plt.ylabel('Count')
plt.title('total event: ' + str(total_event))

# Show the plot
plt.show()  

