#%%
import pandas as pd
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import os
from iminuit import Minuit
from iminuit.cost import LeastSquares
import scienceplots

### initial parameter ###
threshold_first_ch1 = 2
threshold_first_ch2 = 1
stable_level = 0.03
prominence= 0
num_sigma = 7
ground_sigma_num = 5
time_interval = 1e-9

mode = '_stable' + ' (5tau)_fitting' + 'second_no_sigma'
prominice_value = '_' + 'prominence(sigma)'
ground_value = '_' + f'ground({ground_sigma_num})'
# file_names = ['20240415', '20240416_01','20240417_01','20240417_02', '20240418_01','20240419_01','20240421_01','20240422_01','20240422_02','20240423_01','20240511','20240513']
# file_names = ['20240424', '20240428', '20240429', '20240429_02', '20240430', '20240501', '20240502', '20240504', '20240505', '20240507', '20240510']
file_names = ['20240524']
#----------------------------------------------------------------#
def fitting(t, a, b, tau): #charge 
    return a * np.exp(-t/ tau) + b

#%%
#-------------------------------------------------------#
for file_name in file_names:
    file_path = "./lifetimedatatxt/Lifetime_data(B)/" + file_name + ".txt"
    new_file_path = "./new_lifetimedatatxt/Lifetime_data(B)/" + file_name + ".txt"

    f= open(file= file_path, mode= 'r')
    new_f = open(file= new_file_path, mode= 'w')

    #----------------------------------------------------------------#
    count = 0
    Times = 0 #filtered length 
    frames = 1
    print('start run file: ', file_name)
    try:
        while True:
        # while Times < frames :
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
            ch1_peak_index, ch1_peak_info = find_peaks(smoothed1, height= threshold_first_ch1)
            ch2_first_peak_index, ch2_first_peak_info = find_peaks(smoothed2, height= threshold_first_ch2)

            # -----check the first peak--------------#
            #filter zero first peak cases
            if (len(ch1_peak_index) == 0 or len(ch2_first_peak_index) == 0):
                continue
            
            if abs(ch1_peak_index[0] - ch2_first_peak_index[0]) > 1:
                continue
            '''
            ch1_peak_index = ch1_peak_index
            ch2_first_peak_index = ch2_first_peak_index[0]
            ch1_peak_height = ch1_peak_info['peak_heights']
            ch2_first_peak_height = ch2_first_peak_info['peak_heights'][0]
            '''
            #store the waveform
            for i in ch1:
                new_f.write(f'{i}' + '\t')
            new_f.write('\n')
            
            for i in ch2:
                new_f.write(f'{i}' + '\t')
            new_f.write('\n')
            Times += 1
        
    finally:
        f.close()
        new_f.close()
        print(file_name, ' original length: ', count)
        print(file_name, 'length after filter: ', Times)




    # %%
    '''
    #-----fitting gaussain distribution-----#
    from iminuit import Minuit
    from iminuit.cost import LeastSquares
    # from OuOb import pltSty2

    start = 0
    frame = 20
    xs = bins[1+start:frame]
    ys = counts[start:frame-1]
    print()
    # yerror= [np.std(i) for i in data]
    yerror= 1


    def gaussian(x, mu, sigma, N0): #charge 
        return N0 * np.exp(-((x-mu)/sigma)**2 / 2)

    least_square= LeastSquares(xs, ys, yerror= yerror, model= gaussian)
    m=Minuit(least_square, mu= 0.05, sigma= 0.001, N0= 12e4, )
    # plt.errorbar(xs, ys, yerr= stdev(ys)/ , fmt='')
    m.migrad() 
    m.hesse()
    # plt.errorbar(xs, ys, yerr=yerror, fmt="o", label="data",ms= 3, alpha=0.6, color='brown')

    plt.errorbar(xs, ys, yerr=yerror, fmt="o", label="data",ms= 2)
    plt.plot(xs, gaussian(xs, *m.values), label="gaussian", color= 'blue', linestyle= 'dashed')
    #for test
    # plt.scatter(xs, ys, label= 'mean',s=2)


    fit_info = [
        f"$\\chi^2$ / $n_\\mathrm{{dof}}$ = {m.fval:.1f} / {len(xs) - m.nfit}",
    ]
    for p, v, e in zip(m.parameters, m.values, m.errors):
        fit_info.append(f"{p} = ${v:.3e} \\pm {e:.3e}$")


    ax1 = plt.subplot(1,1,1)
    pltSty2(xName = 'Voltage v', yName = 'Entities')
    # ax1.scatter(xs, ys, color = '#070d58', label = 'x(actual)', s = 5)
    ax1.legend(loc = 'best', edgecolor = '#7e7474', fontsize = 12)

    plt.tight_layout()
    plt.draw()
    plt.grid(color= 'gray', visible=True)
    # plt.xlabel('t (s)', fontsize= 15)
    # plt.ylabel('V (volt)', fontsize= 15)
    plt.legend(title="\n".join(fit_info))
    # plt.savefig('test.png',dpi= 500)
    plt.show()
    '''
# %%
