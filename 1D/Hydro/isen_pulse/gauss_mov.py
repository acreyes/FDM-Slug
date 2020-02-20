import numpy as np
import matplotlib.pyplot as plt

#ffmpeg -r 18 -i gauss_%04d.png -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -b 8000k gauss.mp4

class pulse():
    def __init__(self, N):
        self.N = N
        self.name5 = 'slug_gauss5_1%04d.dat' % N
        self.name1 = 'slug_gauss1_1%04d.dat' % N
        self.t = 0.02*N
        self.data1 = np.loadtxt(self.name1)
        self.data5 = np.loadtxt(self.name5)
        self.Nx = len(self.data1[:,0])

    def ent_err(self):
        total_S = 1.*self.Nx
        total_S = 101.359976486
        gamm = 1.66666666667
        S1 = 0.
        S5 = 0.
        for i in range(self.Nx):
            S1 += self.data1[i,3]/(self.data1[i,1]**gamm)
            S5 += self.data5[i,3]/(self.data5[i,1]**gamm)
        return (self.t, abs(S1-total_S), abs(S5-total_S))

    def make_frame(self, t, err1, err5):
        tn, errn1, errn5 = self.ent_err()
        print errn1, errn5
        t.append(tn)
        err1.append(errn1)
        err5.append(errn5)

        fig, axes = plt.subplots(2,2, sharex='row', sharey='row')

        axes[1,1].set_xlim(0.0001,5.)
        axes[1,1].set_ylim(1e-10,5.)
        axes[1,1].set_yscale('log')
        axes[1,1].set_xscale('log')
        axes[0,0].set_ylim(1.,2.)

        S = 1.*self.Nx

        axes[0,0].plot(self.data1[:,0],self.data1[:,1],color='k',marker='x')
        axes[0,0].set_title('1st Order')
        axes[0,1].plot(self.data5[:,0],self.data5[:,1],color='k',marker='x')
        axes[0,1].set_title('5th Order')

        axes[1,0].plot(t,err1,color='r')
        axes[1,1].plot(t,err5,color='r')
        frame = 'gauss_%04d.png' % self.N
        plt.savefig(frame,dpi=200)
        plt.close(fig)
        return(t, err1, err5)
    '''
    def make_plot(self, axD, axS):
        
        axes[0].plot(self.data[:,0],self.data[:,1],color='k',marker='x')
        axes[1].plot(t,err,color='r')
        axes[1].set_xlim(0.,2.5)
        S = 1.*self.Nx
        axes[1].set_ylim(S,S+0.05)
        axes[1].set_yscale('log')
        yloc = abs(errn-S)/(0.05)
        axes[1].axvline(tn,ymin=yloc,color='k',ls='--')
        frame = 'gauss_%04d.png' % self.N
        plt.savefig(frame,dpi=200)
        plt.close(fig)
        return(t,err)
'''
t = []
err1 = []
err5 = []
Nmax = 100
Nmax = 250
dat = pulse(Nmax)
print dat.ent_err()
Nmax = 250
for N in range(Nmax+1):
    dat = pulse(N)
    print N
    t, err1, err5 = dat.make_frame(t,err1,err5)

