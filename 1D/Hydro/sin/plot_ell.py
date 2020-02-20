import numpy as np
import matplotlib.pyplot as plt

#weno = 1.405990e-02
#weno = np.loadtxt('WENO-RK3_Nred.dat')
NN = [25,50,100,200]
N25 = np.loadtxt("ell_RK3_25.dat")
'''
for N in NN:
    fname = "ell_RK4_%d.dat" % N
    data = np.loadtxt(fname)
    label = 'N=%d' % N
    plt.plot(data[:,0],data[:,1],label=label)
'''
colors = ['b','g','r','cyan']
for i in range(len(NN)):
    N = NN[i]
    color = colors[i]
    fname = "ell_RK3_%d.dat" % N
    data = np.loadtxt(fname)
    label = 'N=%d' % N
    plt.plot(data[:,0],data[:,1],label=label,color=color)
    #print N, weno[i,0]
    #plt.axhline(y=weno[i,1],color=color,ls='--')


#plt.axhline(y=weno,xmin=0,xmax=10.,c='r',ls='--',label='WENO-JS')
plt.yscale('log')
plt.xlabel(r'$\ell$')
plt.ylabel(r'$L_1$ Error')
plt.title('1D Smooth Advection')
plt.legend()
plt.show()
