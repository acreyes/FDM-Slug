import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.collections import LineCollection
import numpy as np
import pickle

def plot_dat(fname, title):
    fid = open(fname, 'rb')
    freq, sigs_all, real, imag = pickle.load(fid)
    fid.close()

    fig, (ax1, ax2) = plt.subplots(2,sharex=True)
    fig.suptitle(title)
    
    sigs = [1., 1.5, 3., 5.5,4., 8.]

    ax1.set_xlim(np.min(freq),np.max(freq))
    ax1.set_ylim(-0.7,np.max(freq))
    segs = [np.column_stack([freq,real[sig]]) for sig in sigs[::-1]]
    lc = LineCollection(segs,cmap='rainbow')
    lc.set_array(np.asarray(sigs[::-1]))
    ax1.add_collection(lc)
    ax1.plot(freq,freq,'k--',linewidth=3,zorder=1)
    #ax1.axis('equal')
    ax1.set_title('Dispersion')
    ax2.set_title('Dissipation')
    ax2.set_xlim(np.min(freq),np.max(freq))
    ax2.set_ylim(-1.8,0.1)
    ax2.set_xlabel('wavenumber')
    segs = [np.column_stack([freq,imag[sig]]) for sig in sigs[::-1]]
    lc = LineCollection(segs,cmap='rainbow')
    lc.set_array(np.asarray(sigs[::-1]))
    ax2.add_collection(lc)

    cbar_ax = fig.add_axes([0.91,0.15,0.01,0.7])
    axcb = fig.colorbar(lc,cax=cbar_ax)
    axcb.set_label(r'$\sigma/\Delta$')
    [axcb.ax.hlines(axcb.norm(sig),0,1,color='white') for sig in sigs]
'''    
fid = open('gpr2_sigdel.dat', 'rb')
freq, sigs, real, imag = pickle.load(fid)
fid.close()

fig, ax = plt.subplots()
ax.set_xlim(np.min(freq),np.max(freq))
ax.set_ylim(np.min(freq),np.max(freq))
segs = [np.column_stack([freq,real[sig]]) for sig in sigs]
lc = LineCollection(segs,cmap='rainbow')
lc.set_array(np.asarray(sigs))
ax.add_collection(lc)

ax.plot(freq,freq, 'k--', linewidth = 3)

axcb = fig.colorbar(lc)
'''
for r in range(1,4):
    fname = 'gpr%d_sigdel.dat' % r
    leg = 'GP-R%d' % r
    plot_dat(fname,leg)
    plt.savefig('MwaveNum_gpr%d.pdf'%r, format='pdf')
#plot_dat('gpr4_sigdel.dat', 'GP-R4')

plt.show()
