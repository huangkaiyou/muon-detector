#%%
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os 
import scienceplots
import numpy as np
from iminuit import Minuit
from iminuit.cost import LeastSquares

def channel(version, mode):
    global title_name
    if version == 'first peak stagger':
        title_name = 'first peak stagger: '
        if mode == 'no B':
            path = './analysis_result/background/no B/first peak stagger/' + 'num_sigma(7)/'
            
        elif mode == 'with B':
            path = './analysis_result/background/with B/first peak stagger/' + 'num_sigma(7))/'

        
    elif version == 'no ch1 first peak':
        title_name = 'no ch1 first peak: '
        if mode == 'no B':
            path = './analysis_result/background/no B/no ch1 first peak/' + 'num_sigma(7)/'
            
        elif mode == 'with B':
            path = './analysis_result/background/with B/no ch1 first peak/' + 'num_sigma(7)/'
    return path
# path = r"C:\Users\user\code\python 練習\實驗物理功課\二下小年會\analysis_result\lifetime\with B\scipy\threshold_second (mean, num_sigma(7))_stable (5tau)_fittingsecond_no_sigma_paper_prominence(sigma)_ground(0)\lifetimedata"
# path = r"C:\Users\user\code\python 練習\實驗物理功課\二下小年會\analysis_result\background\with B\num_sigma(7)"



'''========== this is channel ============='''
version = 'first peak stagger'
mode = 'no B'
'''========================================'''
files = glob.glob(channel(version, mode) + r'\*.csv')
name = 'start_stop(B)(01v)'


start_time = 1e-6
end_time = 10e-6
window = 10e-6
mode =  1   #mode 1: no B ; mode 2: With B ; 
            # mode 6: no B, seperate
            # mode 7: with B, seperate exp
            # mode 8: with B, seperate 


# Tau1 = 7.82e-7
# Alpha1 = 0.06
# Omega1 = 4.53e6
# Delta1 = -5.14e-01

# Tau1 = 2.2e-7

raw_time = []
time = []

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



def fit(t, N0, tau, c):
    I = N0*np.exp(-t/tau) + c
    return I

def fit_B(t, N0, tau, alpha, omega, delta, c):
    I = N0*np.exp((-t)/tau)*(1 + alpha*np.cos(omega*t + delta)) + c
    return I

def fit_seperate(t, N1, tau1, N2, tau2, c):
    N = N1*np.exp((-t)/tau1) + N2*np.exp(-t/tau2) + c
    return N

def fit_seperate_B(t, N1, tau1, alpha1, omega1, delta1, N2, tau2, alpha2, omega2, delta2, c):
    N = N1*np.exp((-t)/tau1)*(1 + alpha1*np.cos(omega1*t + delta1)) + N2*np.exp((-t)/tau2)*(1 + alpha2*np.cos(omega2*t + delta2)) + c
    return N
def fit_seperate_B_exp(t, N1, tau1, N2, tau2,  alpha2, omega2, delta2, c):
    N = N1*np.exp((-t)/tau1) + N2*np.exp((-t)/tau2)*(1 + alpha2*np.cos(omega2*t + delta2)) + c
    return N
def fit_seperate_B_sum(t, N1, tau1, N2, tau2, alpha1, omega1, delta1, c) :
    N = (N1*np.exp((-t)/tau1) + N2*np.exp((-t)/tau2)) *(1 + alpha1*np.cos(omega1*t + delta1))+ c
    return N

if mode == 1:
    fit_formula = fit
    least_squares = LeastSquares(x, y, sigma, fit_formula) # fit
    m = Minuit(least_squares, N0=np.amax(y), tau=2.2e-6, c=0)
    fig_title = 'Muon lifetime (No B)'
    m.limits["tau"] = (1e-8, 1.0e-5)
    parameter = [r'$N_0$', r'$\tau$', r'c']
elif mode == 2:
    fit_formula = fit_B
    least_squares = LeastSquares(x, y, sigma, fit_formula) # fit
    m = Minuit(least_squares, N0=np.amax(y), tau=2.2e-6, alpha=0.06, omega=3.76*1e6, delta=-0.6, c=0)
    fig_title = 'Muon lifetime (B)'
    m.limits["tau"] = (1e-8, 1.0e-5)
    parameter = [r'$N_0$', r'$\tau$', r'$\alpha$', r'$\omega$', r'$\delta$', r'c']
    
elif mode == 6:
    fit_formula = fit_seperate
    least_squares = LeastSquares(x, y, sigma, fit_formula) # fit
    m = Minuit(least_squares, N1= np.amax(y), tau1= 5.38e-7, N2= np.amax(y), tau2= 2.2e-6, c=0)
    fig_title = 'Muon lifetime (No B)'
    m.fixed["tau1"] = True
    # m.fixed['n'] = True
    m.limits["tau2"] = (1e-7, 1.0e-5)
    m.limits['c'] = (0, 1e4)
    parameter = [r'$N_1$', r'$\tau1$', r'$N_2$', r'$\tau2$', r'c']

elif mode == 7:
    fit_formula = fit_seperate_B_exp
    least_squares = LeastSquares(x, y, sigma, fit_formula) # fit
    m = Minuit(least_squares, N1=np.amax(y), tau1= 6.34e-7, N2= np.amax(y), tau2= 2.2e-6, alpha2= 0.1, omega2= 4.53e6, delta2= -0.6, c=0)
    fig_title = 'Muon lifetime (B)'
    # m.fixed["tau1"] = True
    m.limits["tau2"] = (1e-8, 1.0e-5)
    parameter = [r'$N_1$', r'$\tau1$'
                 , r'$N_2$', r'$\tau2$', r'$\alpha_2$', r'$\omega_2$', r'$\delta_2$', r'c']

elif mode == 8:
    fit_formula = fit_seperate_B
    least_squares = LeastSquares(x, y, sigma, fit_formula) # fit
    m = Minuit(least_squares, N1=np.amax(y), tau1= 7.65, alpha1= 0.06, omega1= 4.53e6, delta1= -0.6, N2= np.amax(y), tau2= 2.2e-6, alpha2= 0.06, omega2= 4.53e6, delta2= -0.6, c=0)
    fig_title = ' Muon lifetime (B)'
    m.fixed["tau1"] = True
    m.fixed['omega1'] =True
    # m.fixed['delta1'] = True
    m.limits["tau2"] = (1e-8, 1.0e-5)
    parameter = [r'$N_1$', r'$\tau1$', r'$\alpha_1$', r'$\omega_1$', r'$\delta_1$'
                 , r'$N_2$', r'$\tau2$', r'$\alpha_2$', r'$\omega_2$', r'$\delta_2$', r'c']
elif mode == 9:
    fit_formula = fit_seperate_B_sum
    least_squares = LeastSquares(x, y, sigma, fit_formula) # fit
    m = Minuit(least_squares, N1=np.amax(y), tau1= 2.2e-7, N2= np.amax(y), tau2= 2.2e-6, alpha1= 0.06, omega1= 4.53e6, delta1= -0.6, c=0)
    fig_title = 'Muon lifetime (B)'
    m.limits["tau2"] = (1e-8, 1.0e-5)
    parameter = [r'$N_1$', r'$\tau1$', r'$N_2$', r'$\tau2$'
                 , r'$\alpha_1$', r'$\omega_1$', r'$\delta_1$', r'c']

# m = Minuit(least_squares, N1=np.amax(y), tau1=2.2e-6, N2=np.amax(y), tau2=2.2e-6,c=0)

# m.limits["I0"] = (6000, 10000)

# m.fixed['N0'] = True
# m.limits['omega'] = (1e5, 1e7)

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
for p, v, e in zip(parameter, m.values, m.errors):
    fit_info.append(f"{p} = ${v:.2e} \\pm {e:.2e}$")


plt.style.use(['science', 'notebook'])
fig, ax = plt.subplots(1, 1, figsize=(8,6))


s = ax.errorbar(x, y, yerr=sigma, fmt='o', capsize=5, color='black', label='data')
# s = ax.semilogy(x, y, marker='o', color='black', label='data', lw=0)
ax.plot(fitx, fit_formula(fitx, *m.values), color='red', label="fit")
# ax.plot(fitx, fit(fitx, N0= 1.26e4, tau= 4.79e-7, c=0), label = 'theory_bgd')
# ax.plot(fitx, fit(fitx, N0= 1.95e3, tau= 1.78e-6, c=0), label = 'theory_sig')
# ax.plot(fitx, horizon(fitx, offset= 33.4), label= 'shift')


ax.legend(title="\n".join(fit_info),title_fontsize=10, fancybox=False, edgecolor='black', frameon=True, fontsize = 15) 
ax.set_ylabel('Entries', fontsize=15, loc='top')
ax.set_xlabel(r'Lifetime ($\mu$s)', fontsize=15, loc='right')
ax.set_xticks(ticks=[i*1e-6 for i in range(0, 11, 2)], labels=(i for i in range(0, 11, 2)))
ax.set_title(title_name + fig_title)

ax.set_yscale('log')
# ax.plot(x, fit(np.array(x), N0= 3.34e3, tau= 7.66e-7, c= -1))
# plt.savefig(fr"D:\University\Muon_Detection\Result_figure\{name}")
plt.show()