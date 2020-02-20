import numpy as np
import matplotlib.pyplot as plt

def add_plot(ax, fname):
    #analytic = np.loadtxt('sod_analytic.dat')
    data = np.loadtxt(fname)
    fvars = [1,2,3]
    labels = [r'$\rho$',r'$u$',r'$P$']
    colors = ['r','g','b']
    markers = ['o','D','s']
    for var, label, color, marker in zip(fvars, labels, colors, markers):
        #ax.plot(analytic[:,0],analytic[:,var],color=color)
        ax.plot(data[:,0],data[:,var],label=label,color=color,marker=marker,markeredgecolor='k',markersize=4,ls='')
        


fig, ax = plt.subplots()        
fname = "slug__10001.dat"
gpr2 = np.loadtxt('slug_blast_gp_10001.dat')
weno = np.loadtxt('slug_blast_weno_10001.dat')
ref = np.loadtxt('slug_blast_ref_10001.dat')
plt.title('Interacting Blast Wave')
plt.plot(gpr2[:,0],gpr2[:,1],marker='s',label='GP-R2',mec='k',c='r')
plt.plot(weno[:,0],weno[:,1],marker='d',label='WENO-JS',mec='k',c='cyan')
plt.plot(ref[:,0],ref[:,1],'k',label='ref')
plt.xlim(0.4,1.)
plt.ylim(2.,6.5)
plt.legend()
#plt.ylim(-0.02,1.02)
plt.show()
