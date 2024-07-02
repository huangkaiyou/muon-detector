#%%
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os 
import scienceplots
import numpy as np
from iminuit import Minuit
from iminuit.cost import LeastSquares


path1 = r"C:\Users\user\code\python 練習\實驗物理功課\二下小年會\analysis_result\background\no B\ch1\lifetimedata"
path2 = r"C:\Users\user\code\python 練習\實驗物理功課\二下小年會\analysis_result\lifetime\no B\scipy\threshold_second (mean, num_sigma(7))_stable (5tau)_fittingsecond_no_sigma_paper_prominence(sigma)_ground(0)\lifetimedata"
paths = [path1, path2]
name = 'start_stop(B)(01v)'


start_time = 1.e-6
end_time = 10e-6
window = 10e-6

labels= ['ch1: no B', 'ch2: no B']
colors= ['black', 'navy']
plt.style.use(['science', 'notebook'])
fig, ax = plt.subplots(1, 1, figsize=(8,6))
for cycle in range(2):
    raw_time = []
    time = []
    
    path = paths[cycle]
    files = glob.glob(path + r'\*.csv')
    for filename in files:
        df = pd.read_csv(filename)
        tem = df['0'].tolist()
        for i in tem:
            raw_time.append(i)



    for i in raw_time:
        if i>start_time and i <= end_time:
            time.append(i)

    ybin = pd.Series(time)
    q1 = ybin.quantile(0.25)
    q3 = ybin.quantile(0.75)
    iqr = q3 - q1

    bin_width = (2 * iqr) / (len(time) ** (1 / 3))
    # print('bin width : ' + str(bin_width))

    bin_width = 0.2e-6
    bin_number = int(window // bin_width)





    a = np.histogram(time, bins=bin_number, range=(0, 10e-6))
    entries = a[0]
    edges = a[1]


    x = [(edges[i] + edges[i+1]) / 2 for i in range(len(edges) - 1)]


    sigma = []
    for j in range(len(edges)-1):
        count = 0
        for i in time:
            if i>=edges[j] and i<=edges[j+1]:
                count+=1
        s = count * (1-(count/len(time)))
        m = np.sqrt(s)
        sigma.append(m)


    f = {
        'x' : x,
        'y' : entries,
        'sigma' : sigma   
    }

    df = pd.DataFrame(f)
    df = df[(df.x != 0) & (df.y != 0)  & (df.sigma != 0)]
    df.reset_index(inplace=True)

    x = df.x.tolist()
    y = df.y.tolist()
    sigma = df.sigma.tolist()

    print('total events : ' + str(sum(y)))



    def fit(t, N0, tau, c):
        I = N0*np.exp((-t)/tau) + c
        return I


    # def fit_seperate(t, N1, tau1, N2, tau2, c):
    #     N = N1*np.exp((-t)/tau1) + N2*np.exp((-t)/tau2) + c
    #     return N


    def fit_B(t, N0, tau, alpha, omega, delta, c):
        I = N0*np.exp((-t)/tau)*(1 + alpha*(omega*t + delta)) + c
        return I

    mode = 2
    if mode == 1:
        fit_formula = fit
    else:
        fit_formula = fit_B


    least_squares = LeastSquares(x, y, sigma, fit_formula) # fit
    # m = Minuit(least_squares, N0=np.amax(y), tau=2.2e-6, c=0)
    # m = Minuit(least_squares, N1=np.amax(y), tau1=2.2e-6, N2=np.amax(y), tau2=2.2e-6,c=0)
    m = Minuit(least_squares, N0=np.amax(y), tau=2.2e-6, alpha=0.06, omega=3.76*1e6, delta=-0.6, c=0)

    # m.limits["I0"] = (6000, 10000)
    m.limits["tau"] = (1.5e-6, 2.5e-6)
    # m.fixed['N0'] = True
    m.limits['omega'] = (1e5, 1e7)

    m.migrad()
    m.hesse()
    print(m.hesse())
    print(m.covariance.correlation())

    fitx = np.linspace(np.amin(x), np.amax(x), 1000)

    # parameter = (r'$I_0$', r'$\tau$', r'c')

    # # = {m.fval:.1e} / {m.ndof:.0e} 
    fit_info = [
        f"$\\chi^2$/$n_\\mathrm{{dof}}$ = {m.fmin.reduced_chi2:.1e}",
    ]
    for p, v, e in zip(m.parameters, m.values, m.errors):
        fit_info.append(f"{p} = ${v:.2e} \\pm {e:.2e}$")


    

    y /= np.max(y)
    s = ax.scatter(x, y, label= labels[cycle])
    # s = ax.errorbar(x, y, yerr=sigma, fmt='o', capsize=5, color= colors[cycle], label= labels[cycle])
    # s = ax.semilogy(x, y, marker='o', color='black', label='data', lw=0)
    # ax.plot(fitx, fit_formula(fitx, *m.values), color='red', label="fit")
   

    # ax.legend(title="\n".join(fit_info),fontsize=10, fancybox=False, edgecolor='black', frameon=True) 
    ax.legend(fontsize=15, fancybox=False, edgecolor='black', frameon=True) 
    
    ax.set_ylabel('Entries', fontsize=15, loc='top')
    ax.set_xlabel(r'Lifetime ($\mu$s)', fontsize=15, loc='right')
    ax.set_xticks(ticks=[i*1e-6 for i in range(0, 11, 2)], labels=(i for i in range(0, 11, 2)))
    ax.set_title('Muon lifetime')


# plt.savefig(fr"D:\University\Muon_Detection\Result_figure\{name}")
ax.set_yscale('log')
plt.show()
# %%
