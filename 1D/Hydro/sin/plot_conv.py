import numpy as np
import matplotlib.pyplot as plt

def logline(x, x0, y0, m):
    return(y0*(x/x0)**-m)
def make_label(order):
    return(r'$\Delta^%d$' % order)

fname = "GP-RK3_red8.dat"
gpr2 = np.loadtxt(fname)

fname = "WENO-RK3_nored.dat"
weno = np.loadtxt(fname)

Nx = np.linspace(gpr2[0,0],gpr2[len(gpr2[:,0])-1,0],100)

plt.plot(gpr2[:,0], gpr2[:,1], marker='o',linestyle='',label='GP-R2')
order = 5
plt.plot(Nx,logline(Nx,gpr2[3,0],gpr2[3,1],order),'k--')#,label=make_label(order))

plt.plot(weno[:,0], weno[:,1], marker='*',linestyle='',label='WENO-JS')
order = 5
plt.plot(Nx,logline(Nx,100,weno[2,1],order),'k--',label=make_label(order))

plt.yscale('log')
plt.xscale('log')
plt.legend()


plt.xlabel('Resolution')
plt.ylabel(r'$L_1$ Error')
plt.show()


