import numpy as np
import matplotlib.pyplot as plt

def get_data(fbase):
    f0 = fbase + '_10000.dat'
    data = np.loadtxt(f0)

    x = data[:,0]
    d0 = data[:,1] 

    f1 = fbase + '_10001.dat'
    data = np.loadtxt(f1)
    d1 = data[:,1] 
    return(x, d0, d1)

def add_plot(fbase, ax, **kwargs):
    x, d0, d1 = get_data(fbase)
    ax.plot(x, d0, '-k')
    ax.plot(x, d1, linestyle='', **kwargs)


fig0, ax0 = plt.subplots()
fig, ax = plt.subplots()
fbase = 'slug_weno'

add_plot(fbase, ax0, marker='o', label='WENO-JS')
add_plot(fbase, ax , marker='o', label='WENO-JS')

#ax0.set_title('WENO-JS')
ax0.set_ylim(0,1.2)
#plt.savefig('WENO-JS', dpi=200)
R = [1,2,3]
R = [2]

fbase = 'slug_gpr2_zero'
add_plot(fbase, ax0, marker='^', label='GP-R2: zero mean')
ax0.legend(loc='upper right')
fig0.savefig('zero-mean.png',dpi=200)
fbase = 'slug_gpr2_nonzero'
add_plot(fbase, ax, marker='^', label='GP-R2: non-zero mean')
ax.legend(loc='upper right')
ax.set_ylim(0.,1.2)
plt.savefig('nonzero-mean.png',dpi=200)
'''
for r in R:
    fig[r], ax[r] = plt.subplots()
    fbase = 'slug_scal_advE_cfl8_gpr%d' % r
    add_plot(fbase, ax[r], marker='o')
    
    title = 'GP-R%d: zero mean + Euler' % r
    ax[r].set_title(title)
    ax[r].set_ylim(1,2.2)
    title = 'GP-R%d' % r
    plt.savefig(title, dpi=200)
'''
plt.show()
    
