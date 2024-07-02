import numpy as np
import matplotlib.pyplot as plt
from iminuit import Minuit
from iminuit.cost import LeastSquares


format 
def f(t, a, b): #charge 
    return a + (-b*t)

least_square= LeastSquares(xs, ys, yerror= yerror, model= f)
m=Minuit(least_square, tau=0.0051, V=3,c=-0.)
# plt.errorbar(xs, ys, yerr= stdev(ys)/ , fmt='')
m.migrad() 
m.hesse()
plt.errorbar(xs, ys, yerr=yerror, fmt="o", label="data",ms= 3, alpha=0.6, color='brown')

# plt.errorbar(xs1, ys2, yerr=yerror, fmt="o", label="raspberry pi",ms= 2)
# plt.plot(xs, f(xs, *m.values), label="fit", color= 'green')
plt.plot(xs, g(xs, *m.values), label="theory", color= 'blue', linestyle= 'dashed')
#for test
# plt.scatter(xs, ys, label= 'mean',s=2)


fit_info = [
    f"$\\chi^2$ / $n_\\mathrm{{dof}}$ = {m.fval:.1f} / {len(xs) - m.nfit}",
]
for p, v, e in zip(m.parameters, m.values, m.errors):
    fit_info.append(f"{p} = ${v:.3e} \\pm {e:.3e}$")

print(m.values['V'])
ax1 = plt.subplot(1,1,1)
pltSty2(xName = 'Time t', yName = 'Voltage v')
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
