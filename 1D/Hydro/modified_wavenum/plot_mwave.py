import pickle
import numpy
import matplotlib.pyplot as plt

fname = 'mod_wavenum.dat'
fid = open(fname, 'rb')
freq, real, imag = pickle.load(fid)

fig, (ax1, ax2) = plt.subplots(2,sharex=True)
ax1.set_title('Dispersion')
ax2.set_title('Diffusion')
ax2.set_xlabel('wavenumber')
ax1.plot(freq,freq,'k--',linewidth=3)

keys = ['weno', 'gpr1', 'gpr2', 'gpr3']
legs = ['WENO-JS', 'GP-R1', 'GP-R2', 'GP-R3']

for key, leg in zip(keys,legs):
    ax1.plot(freq, real[key], label=leg)
    ax2.plot(freq, imag[key], label=leg)
#ax1.legend()
ax2.legend()
#plt.show()
plt.savefig('MwaveNum_all.pdf',format='pdf')
