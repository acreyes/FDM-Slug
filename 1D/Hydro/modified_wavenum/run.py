import numpy as np
import os
import matplotlib.pyplot as plt
import pickle

class data_set():
    def __init__(self, pars):
        self.pars = pars
        self.name = pars.name
        #self.x, self.d, self.dp = self.get_data()
        #self.k, self.freq, self.f, self.fp = self.fft()
        self.k = []
        self.w = []

    def get_vhat(self, k, d):
        pars = self.pars
        vhat = 0.
        for i in range(pars.nx):
            vhat += d[i]*np.exp(-1j*2.*np.pi*k*i/pars.nx)
        return(vhat/pars.nx)

    def get_data(self):
        pars = self.pars
        kk = np.arange(1,pars.nx/2 )
        w = []
        ks = []
        dx = 1./pars.nx
        dt = 0.0001
        dtdx = dt/dx
        print dtdx
        #plt.figure()
        for k in kk:
            #if k > 5: break
            pars.k = k
            pars.run()
            freq = np.pi*k/kk[-1]
            x, d0, dn = pars.get_data()
            '''
            vhat  = np.fft.rfft(dn)
            vhat0 = np.fft.rfft(d0)
            #plt.plot(vhat)
            wn = (1j/dtdx)*np.log(vhat[k]/vhat0[k])
            '''
            wn = (1j/dtdx)*np.log(self.get_vhat(k,dn)/self.get_vhat(k,d0))
            w.append(wn)
            ks.append(freq)
        '''
        pars.k = 1
        pars.run()
        x, d0, dn = pars.get_data()
        vhat = np.fft.rfft(dn)
        vhat0 = np.fft.rfft(d0)
        #ks = np.fft.rfftfreq(dn.shape[-1])
        ks = np.pi*kk/kk[-1]
        w = 1j/dtdx * np.log(vhat/vhat0)
        '''
        #print (w[1].real-w[0].real)/(ks[1]-ks[0])
        #return(np.array(ks), np.array(w))
        self.k = np.array(ks)
        self.w = np.array(w)

    def filter_dat(self, d):
        eps = 0.00
        N = len(d)
        dn = np.copy(d)
        x = np.zeros(N)

        for i in range(2,N-2):
            mpp = d[i+2]-  d[i+1]
            mp  = d[i+1]-  d[i]
            m   = d[i+1] - d[i-1]
            mm  = d[i]   - d[i-1]
            mmm = d[i-1] - d[i-2]
            if  (mp*mm < 0.) and (mmm*m > -eps) :
                d[i] = 0.5*(d[i+1]+d[i-1])
        for i in range(2,N-2):
            mpp = d[i+2]-  d[i+1]
            mp  = d[i+1]-  d[i]
            m   = d[i+1] - d[i-1]
            mm  = d[i]   - d[i-1]
            mmm = d[i-1] - d[i-2]
            if  (mp*mm < 0.) and (mmm*m > -eps) :
                d[i] = 0.5*(d[i+1]+d[i-1])
    
        return(d)
    def plot_disp(self, ax):
        k = self.k
        wn = self.w
        label = self.name.split('_')
        print label
        ax.set_title('Dispersion')
        ax.plot(k, self.filter_dat(wn.real), label=label[0])

    def plot_diff(self, ax):
        k = self.k
        wn = self.w
        label = self.name.split('_')
        ax.set_title('Diffusion')
        ax.plot(k, self.filter_dat(wn.imag), label=label[0])
        
    
class run_pars():
    def __init__(self):
        self.order = 10
        self.R = 2
        self.ell = 0.08
        self.sigdel = 3.
        self.nx = 400
        self.name = self.get_name()
        self.k = 1

    def get_data(self):
        fname = 'slug_%s_10000.dat' % self.name
        data = np.loadtxt(fname,usecols=(0,1))
        x = data[:,0]
        d0 = data[:,1]
        fname = 'slug_%s_10001.dat' % self.name
        data = np.loadtxt(fname,usecols=(0,1))
        dn = data[:,1]
        return(x, d0, dn)
        
    def get_name(self):
        if self.order == 5:
            meth = 'WENO-JS'
        else :
            meth = 'GP-R%d' % self.R
        name = "%s_%d" % (meth, self.nx)
        return(name)

    def make_init(self):
        dt = 0.8/pars.nx
        
        self.name = self.get_name()
        fid = open('base.init', 'r')
        fwr = open('slug.init', 'w')

        for line in fid :
            fwr.write(line)

        nline = 'sim_order  %d\n' % self.order
        fwr.write(nline)
        nline = 'gp_radius  %d\n' % self.R
        fwr.write(nline)
        nline = 'gp_ell  %f\n' % self.ell
        fwr.write(nline)
        nline = 'gp_sigdel  %f\n' % self.sigdel
        fwr.write(nline)
        nline = 'gr_nx  %d\n' % self.nx
        fwr.write(nline)
        nline = "sim_name  '%s'\n\n" % self.name
        fwr.write(nline)
        nline = "sim_kwave  %f\n\n" % self.k
        fwr.write(nline)
        nline = "sim_dt %f\n\n" % dt
        fwr.write(nline)
        fid.close()
        fwr.close()

    def run(self):
        self.make_init()
        os.system('../slugEuler1d > dump.dat')

    def add_plot(self, ax):
        data = np.loadtxt('slug_%s_10001.dat'%self.name)
        ax.plot(data[:,0],data[:,1],label=self.name)

'''
fig, (ax, ax2) = plt.subplots(2, sharex=True)
#fig2, ax2 = plt.subplots()
#ax.set_xlim(0,100)
data_real = {}
data_imag = {}
pars = run_pars()
pars.nx = 512
pars.k = 10
pars.order = 5
pars.run()
d1 = np.loadtxt('slug_%s_10001.dat'%pars.name)
#pars.add_plot(ax)
weno = data_set(pars)
weno.get_data()
data_real['weno'] = weno.filter_dat(weno.w.real)
data_imag['weno'] = weno.filter_dat(weno.w.imag)
weno.plot_disp(ax)
weno.plot_diff(ax2)


kk = np.arange(1,pars.nx/2+1)
kk = kk*np.pi/kk[-1]
ax.plot(kk,kk,'r--')
pars.order = 10
freq = weno.k
limit = 10
#limit = weno.freq[-1]
#ax.set_xlim(0,limit)
#ax.set_ylim(0,limit)

for r in range(1,4):
    pars.R = r
    key = 'gpr%d' % r
    pars.run()
    d2 = np.loadtxt('slug_%s_10001.dat'%pars.name)
    #pars.add_plot(ax)
    gp = data_set(pars)
    gp.get_data()
    data_real[key] = gp.filter_dat(gp.w.real)
    data_imag[key] = gp.filter_dat(gp.w.imag)
    gp.plot_disp(ax)
    gp.plot_diff(ax2)
'''
pars = run_pars()
pars.order = 10
pars.nx = 512
pars.R = 4
sigs = np.arange(1., 12.5, .5)
gpr2_real = {}
gpr2_imag = {}
print sigs
for sig in sigs:
    print sig, 'of', sigs[-1]
    pars.sigdel = sig
    gp = data_set(pars)
    gp.get_data()
    gpr2_real[sig] = gp.filter_dat(gp.w.real)
    gpr2_imag[sig] = gp.filter_dat(gp.w.imag)


print sigs
freq = gp.k
fid = open('gpr4_sigdel.dat', 'wb')
pickle.dump((freq, sigs, gpr2_real, gpr2_imag), fid)
fid.close()

'''
fid = open('mod_wavenum.dat','wb')
pickle.dump((freq, data_real, data_imag),fid)
fid.close()
ax.legend()
ax2.legend()
plt.show()
'''

