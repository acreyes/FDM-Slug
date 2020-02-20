import numpy as np
import matplotlib.pyplot as plt

#ffmpeg -r 18 -i isen_%04d.png -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -b 8000k isen_pulse.mp4

class pulse():
    def __init__(self, N):
        self.N = N
        self.name = 'slug_gauss_1%04d.dat' % N
        self.t = 0.02*N
        self.data = np.loadtxt(self.name)
        self.Nx = len(self.data[:,0])

    def ent_err(self):
        total_S = 1.*self.Nx
        gamm = 1.66666666667
        gamm = 1.4
        S = 0.
        for i in range(self.Nx):
            S += self.data[i,3]/(self.data[i,1]**gamm)
        return (self.t, np.abs(S))

    def make_plot(self, t, err):
        tn, errn = self.ent_err()
        t.append(tn)
        err.append(errn)
        fig, axes = plt.subplots(2)
        axes[0].plot(self.data[:,0],self.data[:,1],color='k',marker='x')
        axes[1].plot(np.array(t)+1.e-12,err,color='r')
        axes[1].set_xlim(1.e-12,2.5)
        S = 1.*self.Nx
        axes[1].set_ylim(S-.01,S+0.05)
        axes[1].set_yscale('log')
        #axes[1].set_xscale('log')
        yloc = abs(errn-S)/(0.05)
        #axes[1].axvline(tn,ymin=yloc,color='k',ls='--')
        frame = 'isen_%04d.png' % self.N
        plt.savefig(frame,dpi=200)
#        if self.N == 50: plt.show()

        plt.close(fig)
        return(t,err)
t = []
err = []
Nmax = 125

for N in range(Nmax+1):
    print N
    dat = pulse(N)
    t, err = dat.make_plot(t,err)


