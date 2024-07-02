#%%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob
# from OuOb import pltSty2

# %matplotlib widget
folder_path = 'lifetimedata'
file_paths = glob.glob(folder_path + "/*.CSV")
threshold = 1
#%% ### read csv ###
for file_path in file_paths:
    df = pd.read_csv(file_path, header= 25)
    df = df.iloc[:, :-1]
    df.columns = ["Time1", "Amplitude1", "Time2", "Amplitude2"]


    # %%
    # df.plot(x='Time1', y=['Amplitude1', 'Amplitude2'], label=['CH 1', 'CH 2'], kind= 'scatter')
    # plt.show()
    if ((df['Amplitude1'] > threshold).any() and (df['Amplitude2'] > threshold).any()) :
        x = df['Time1'].to_numpy()
        y1 = df['Amplitude1'].to_numpy()
        y2 = df['Amplitude2'].to_numpy()
        print(file_path)


        # plt.scatter(x, y, label= "CH1")
        plt.scatter(x, y1, label= "CH1")
        plt.scatter(x, y2, label= "CH2")
        plt.xlabel('Time')
        plt.ylabel('Voltage')
        plt.legend()
        plt.show()

'''
from OuOb import pltSty2
import numpy as np
import matplotlib.pyplot as plt
from iminuit import Minuit
from iminuit.cost import LeastSquares

def f(x, tau, V, c): #charge 
    tau= 0.0051
    return V*(-np.exp(-((x-c)/(tau))))

index1 = np.argmax(y1)
num = 50
xs=x[index1:index1+num]
ys= y1[index1:index1+num]
yerror=0.01/(2*3**(1/2))

least_square= LeastSquares(xs, ys, yerror= yerror, model= f)
m=Minuit(least_square, tau=0.0051, V=0.2,c=-0.)
# plt.errorbar(xs, ys, yerr= stdev(ys)/ , fmt='')
m.migrad() 
m.hesse()
plt.errorbar(xs, ys, yerr=yerror, fmt="o", label="data",ms= 5, color= 'brown')
plt.plot(xs, f(xs, *m.values), label="fit")

fit_info = [
    f"$\\chi^2$ / $n_\\mathrm{{dof}}$ = {m.fval:.1f} / {len(xs) - m.nfit}",
]
for p, v, e in zip(m.parameters, m.values, m.errors):
    fit_info.append(f"{p} = ${v:.3e} \\pm {e:.3e}$")

ax1 = plt.subplot(1,1,1)
pltSty2(xName = 'Time t', yName = 'Voltage v')
# ax1.scatter(xs, ys, color = '#070d58', label = 'x(actual)', s = 5)
ax1.legend(loc = 'best', edgecolor = '#7e7474', fontsize = 12)
plt.legend()
plt.tight_layout()
plt.draw()
# plt.xlabel('t (s)', fontsize= 15)
# plt.ylabel('V (volt)', fontsize= 15)
# plt.legend(title="\n".join(fit_info))
# plt.savefig('test.png',dpi= 500)
plt.show()
'''
# %%
