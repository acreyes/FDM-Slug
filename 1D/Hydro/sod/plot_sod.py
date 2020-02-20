import numpy as np
import matplotlib.pyplot as plt

def add_plot(ax, fname):
    analytic = np.loadtxt('sod_analytic.dat')
    data = np.loadtxt(fname)
    fvars = [1,2,3]
    labels = [r'$\rho$',r'$u$',r'$P$']
    colors = ['r','g','b']
    markers = ['o','D','s']
    for var, label, color, marker in zip(fvars, labels, colors, markers):
        #ax.plot(analytic[:,0],analytic[:,var],color=color)
        ax.plot(data[:,0],data[:,var],label=label,color=color,marker=marker,markeredgecolor='k',markersize=4,ls='')
        


fig, ax = plt.subplots()        
fname = "slug_sod_10_10001.dat"
plt.title('GP-R3')
add_plot(ax, fname)
plt.legend()
#plt.ylim(-0.02,1.02)
plt.show()
