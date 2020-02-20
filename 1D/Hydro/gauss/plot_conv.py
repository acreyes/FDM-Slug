import numpy as np
import matplotlib.pyplot as plt

def logline(x, x0, y0, m):
    return(y0*(x/x0)**-m)
def make_label(order):
    return(r'$\Delta^%d$' % order)

fname = "GPR1-RK4_red.dat"
gpr1 = np.loadtxt(fname)

fname = "GPR2-RK4_red.dat"
gpr2 = np.loadtxt(fname)

fname = "GPR3-RK4_red.dat"
gpr3 = np.loadtxt(fname)


fname = "WENO-RK4_red.dat"
weno = np.loadtxt(fname)

Nx = np.linspace(gpr2[0,0],gpr2[len(gpr2[:,0])-1,0],100)

plt.plot(gpr1[:,0], gpr1[:,1],c='b', marker='o',linestyle='',label='GP-R1',mec='k',ms=7)
order = 3
plt.plot(Nx,logline(Nx,gpr1[2,0],gpr1[2,1],order),'--',c='b')



plt.plot(gpr2[:,0], gpr2[:,1], marker='^',linestyle='',label='GP-R2',c='g',mec='k',ms=7)
order = 5
plt.plot(Nx,logline(Nx,gpr2[2,0],gpr2[2,1],order),'--',c='g')

plt.plot(gpr3[:,0], gpr3[:,1], marker='D',linestyle='',label='GP-R3',c='r',mec='k',ms=7)
order = 7
plt.plot(Nx,logline(Nx,gpr3[1,0],gpr3[1,1],order),'--',c='r')
'''
plt.plot(weno[:,0], weno[:,1], marker='*',linestyle='',label='WENO-JS')
order = 4
plt.plot(Nx,logline(Nx,100,weno[2,1],order),'--',label=make_label(order))
'''
plt.title('Gaussian Advection')
plt.yscale('log')
plt.xscale('log')
plt.legend()


plt.xlabel('Resolution')
plt.ylabel(r'$L_1$ Error')
plt.show()


