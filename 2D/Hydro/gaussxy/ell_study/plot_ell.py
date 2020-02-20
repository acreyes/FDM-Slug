import numpy as np
import matplotlib.pyplot as plt

#weno = 1.405990e-02

NN = [25,50,100,200]

for N in NN:
    fname = "ell_RK3_gauss_%d.dat" % N
    data = np.loadtxt(fname)
    label = 'N=%d' % N
    plt.plot(data[:,0],data[:,1],label=label)


#plt.axhline(y=weno,xmin=0,xmax=10.,c='r',ls='--',label='WENO-JS')
plt.yscale('log')
plt.xlabel(r'$\ell$')
plt.ylabel(r'$L_1$ Error')
plt.title('Gaussian Advection')
plt.legend()
plt.show()
