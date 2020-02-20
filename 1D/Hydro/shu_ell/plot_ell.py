import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as intp

#weno = 1.405990e-02
#weno = np.loadtxt('WENO-RK3_Nred.dat')
NN = [64,128,256]
NN = [32,64,128,256,512]
L0 = 3.*9./32.
NNN = np.array([1./32.,1./64.,1./128.,1./256.,1./512.])
deltas = 9.*NNN
eldels = L0/np.array(deltas)

#N25 = np.loadtxt("ell_RK3_25.dat")
'''
for N in NN:
    fname = "ell_RK4_%d.dat" % N
    data = np.loadtxt(fname)
    label = 'N=%d' % N
    plt.plot(data[:,0],data[:,1],label=label)
'''
L1 = np.zeros(len(NN))
colors = ['b','g','r','cyan','k']
for i in range(len(NN)):
    N = NN[i]
    color = colors[i]
    fname = "ell_RK3_%d.dat" % N
    data = np.loadtxt(fname)
    label = 'N=%d' % N
    plt.plot(data[:,0],data[:,1],label=label,color=color)
    f = intp.interp1d(data[:,0],data[:,1],fill_value='extrapolate')
    L1[i] = f(eldels[i])
    #print N, weno[i,0]
    #plt.axhline(y=weno[i,1],color=color,ls='--')

plt.plot(eldels,L1,marker='*',color='orange',ms=10)


#plt.axhline(y=weno,xmin=0,xmax=10.,c='r',ls='--',label='WENO-JS')
#plt.yscale('log')
plt.xlim(1,10)
plt.xlabel(r'$\ell/\Delta$')
plt.ylabel(r'$L_1$ Error')
plt.title('ShuOsher Problem')
plt.legend()
plt.show()
