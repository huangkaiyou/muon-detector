#%%
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os 
import scienceplots
import numpy as np
from iminuit import Minuit
from iminuit.cost import LeastSquares

def fit(t, N0, tau, c):
    I = N0*np.exp((-t)/tau) + c
    return I

def fit_B(t, N0, tau, alpha, omega, delta, c):
    I = N0*np.exp((-t)/tau)*(1 + alpha*np.cos(omega*t + delta)) + c
    return I

plt.style.use(['science', 'notebook'])
fig, ax = plt.subplots(1, 1, figsize=(8,6))

path = r"C:\Users\user\code\python 練習\實驗物理功課\二下小年會\analysis_result\lifetime\no B\scipy\threshold_second (mean, num_sigma(7))_stable (5tau)_fittingsecond_no_sigma_paper_prominence(sigma)_ground(0)\lifetimedata"
# path = r"C:\Users\user\code\python 練習\實驗物理功課\二下小年會\analysis_result\lifetime\with B\scipy\threshold_second (mean, num_sigma(7))_stable (5tau)_fittingsecond_no_sigma_paper_prominence(sigma)_ground(0)\lifetimedata"
files = glob.glob(path + r'\*.csv')

name = 'start_stop(B)(01v)'


start_time = 1e-6
end_time = 10e-6
window = 10e-6
mode = 1 #mode 1: no B ; mode 2: With B
cut_index = 10

raw_time = []
time = []

for filename in files:
    df = pd.read_csv(filename)
    tem = df['0'].tolist()
    for i in tem:
        raw_time.append(i)


#select data segment 
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
        if i>edges[j] and i<=edges[j+1]:
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

xs1 = x[:cut_index]
xs2 = x[cut_index:]
xs = [xs1, xs2]

ys1 = y[:cut_index]
ys2 = y[cut_index:]
ys = [ys1, ys2]

sigma1 = sigma[:cut_index]
sigma2 = sigma[cut_index:]
sigmas = [sigma1, sigma2]
#-----------------fit-----------------#
for i in range(2):
    x = xs[i]
    y = ys[i]
    sigma = sigmas[i]
    if mode == 1:
        fit_formula = fit
        least_squares = LeastSquares(x, y, sigma, fit_formula) # fit
        m = Minuit(least_squares, N0=np.amax(y), tau=2.2e-6, c=0)
        fig_title = 'Muon lifetime (No B)'
    elif mode == 2:
        fit_formula = fit_B
        least_squares = LeastSquares(x, y, sigma, fit_formula) # fit
        m = Minuit(least_squares, N0=np.amax(y), tau=2.2e-6, alpha=0.06, omega=3.76*1e6, delta=-0.6, c=0)
        fig_title = 'Muon lifetime (With B)'



    # m = Minuit(least_squares, N1=np.amax(y), tau1=2.2e-6, N2=np.amax(y), tau2=2.2e-6,c=0)

    # m.limits["I0"] = (6000, 10000)
    m.limits["tau"] = (1.e-7, 1.0e-5)
    # m.fixed['N0'] = True
    # m.limits['omega'] = (1e5, 1e7)

    m.migrad()
    m.hesse()
    print(m.hesse())
    print(m.covariance.correlation())

    fitx = np.linspace(np.amin(x), np.amax(x), 1000)

    # parameter = (r'$I_0$', r'$\tau$', r'c')

    # # = {m.fval:.1e} / {m.ndof:.0e} 

    parameter = [r'$N_0$', r'$\tau$', r'$\alpha$', r'$\omega$', r'$\delta$', r'c']

    fit_info = [
        f"$\\chi^2$/$n_\\mathrm{{dof}}$ = {m.fmin.reduced_chi2:.1e}",
    ]
    for p, v, e in zip(parameter, m.values, m.errors):
        fit_info.append(f"{p} = ${v:.2e} \\pm {e:.2e}$")


    s = ax.errorbar(x, y, yerr=sigma, fmt='o', capsize=5, color= 'black', label=['data', ''][i])
    # s = ax.semilogy(x, y, marker='o', color='black', label='data', lw=0)
    ax.plot(fitx, fit_formula(fitx, *m.values), color=['red', 'blue'][i], label=["", ''][i])
    # ax.plot(fitx, fit_formula(fitx, N0=np.amax(y), tau=2.2e-6, alpha=-2.43, omega=3.76*1e6, delta=-0.6, c=0), label = 'theory')

ax.legend(title_fontsize=12, fancybox=False, edgecolor='black', frameon=True, fontsize = 15) 
ax.set_ylabel('Entries', fontsize=15, loc='top')
ax.set_xlabel(r'Lifetime ($\mu$s)', fontsize=15, loc='right')
ax.set_xticks(ticks=[i*1e-6 for i in range(0, 11, 2)], labels=(i for i in range(0, 11, 2)))
ax.set_title(fig_title)

ax.set_yscale('log')
# plt.savefig(fr"D:\University\Muon_Detection\Result_figure\{name}")
plt.show()
# %%
