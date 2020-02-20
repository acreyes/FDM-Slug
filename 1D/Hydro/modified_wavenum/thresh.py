import pickle
import numpy as np
import matplotlib.pyplot as plt

class spec_prop():
    def __init__(self, freq, real, imag):
        self.freq = freq
        self.disp = real
        self.diff = imag

    def get_kc(self, eps):
        delta = np.abs(self.disp - self.freq)/freq
        w = freq[delta>eps]
        k_disp = w[0]/np.pi
        
        delta = np.abs(self.diff)
        w = freq[delta>eps]
        k_diff = w[0]
        return(k_disp, k_diff)

class gp_data():
    def __init__(self, fname, eps):
        fid = open(fname, 'rb')
        freq, sigs, real, imag = pickle.load(fid)
        fid.close()

        self.sigs = sigs
        self.eps = eps
        self.k_crit = np.zeros((len(sigs),2))
        for i in range(len(self.sigs)):
            sig = sigs[i]
            data = spec_prop(freq, real[sig], imag[sig])
            self.k_crit[i,:] = np.array(data.get_kc(eps))
           
        print fname, sigs[2], self.k_crit[2,0]

#get WENO data
eps = 0.1

            
fname = 'mod_wavenum.dat'
fid = open(fname, 'rb')
freq, real, imag = pickle.load(fid)
fid.close()

weno = spec_prop(freq, real['weno'], imag['weno'])
weno_disp, weno_diff = weno.get_kc(eps)
print 'WENO: ', weno_disp
gpr1 = gp_data('gpr1_sigdel.dat', eps)
gpr2 = gp_data('gpr2_sigdel.dat', eps)
gpr3 = gp_data('gpr3_sigdel.dat', eps)
gpr4 = gp_data('gpr4_sigdel.dat', eps)
plt.figure()
plt.axhline(weno_disp, color='k', linestyle = '--', label='WENO-JS')
plt.plot(gpr1.sigs, gpr1.k_crit[:,0], label='GP-R1', ls='-', marker='o')
plt.plot(gpr1.sigs, gpr2.k_crit[:,0], label='GP-R2', ls='-', marker='^')
plt.plot(gpr1.sigs, gpr3.k_crit[:,0], label='GP-R3', ls='-', marker='x')
#plt.plot(gpr1.sigs, gpr4.k_crit[:,0], label='GP-R4', ls='-', marker='+')
plt.legend(loc=1)

plt.xlabel(r'$\sigma/\Delta$')
plt.ylabel('Resolving Efficiency')
#plt.title(r'$\epsilon=%.02e$'%eps)
plt.title('')

plt.savefig('resolving_eff_eps-01.pdf',format='pdf')

