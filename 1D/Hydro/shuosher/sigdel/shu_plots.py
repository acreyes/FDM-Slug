import matplotlib.pyplot as plt
import numpy as np


class shu_pars():
    def __init__(self, eldel, sigdel, nx, R):
        self.name = 'shu_default'
        self.eldel = eldel
        self.sigdel = sigdel
        self.nx = nx
        self.R = R
        self.get_name()

    def get_name(self):
        self.name = 'shu_%d_R%d_el%.1f_sig%.1f' % (self.nx, self.R, self.eldel, self.sigdel)

    def add_plot(self, ax, xlim, ylim):
        fname = 'slug_%s_10001.dat' % self.name
        data = np.loadtxt(fname[:35])
        ref_dat = np.loadtxt('slug_shu_ref_10001.dat')
        label = r'$\sigma/\Delta=%.1f$' % self.sigdel
        ax.plot(data[:,0], data[:,1], ls='-', marker='x',label=label)
        ax.plot(ref_dat[:,0], ref_dat[:,1],'r')
        ax.legend()
        ax.set_xlim(xlim[0],xlim[1])
        ax.set_ylim(ylim[0],ylim[1])

    def add_plot_sig(self, ax,xlim,ylim):
        fname = 'slug_%s_10001.dat' % self.name
        data = np.loadtxt(fname[:35])
        ref_dat = np.loadtxt('slug_shu_ref_10001.dat')
        label = r'$\ell/\Delta=%.1f$' % self.eldel
        ax.plot(data[:,0], data[:,1], ls='-', marker='x',label=label, color='r')
        ax.plot(ref_dat[:,0], ref_dat[:,1],'k')
        #ax.legend()
        ax.set_xlim(xlim[0],xlim[1])
        ax.set_ylim(ylim[0],ylim[1])
        
def make_plot24(R, N):
    xlim = [4.5,7.5]
    ylim = [2.5,5.]
    eldels = [6,12,24,32]
    sigdels = [1.5,2,3,4]#,5]
    eldels = [6,32]

    Nx = 4 #4 columns
    Ny = 2 #5 rows
    fig, axes = plt.subplots(Ny,Nx, sharex='all',sharey='all')
    for i in range(Nx):
        #loop over eldels
        title = r'$\sigma/\Delta = %.1f$' % sigdels[i]
        axes[0,i].set_title(title)
        for j in range(Ny):
            #loop over sigdels
            pars = shu_pars(eldels[j],sigdels[i],N,R)
            pars.add_plot_sig(axes[j,i],xlim,ylim)
    fig.subplots_adjust(hspace=0.,wspace=0.)
    title = 'GP-R%d; N=%d' % (R,N)
    plt.suptitle(title)
#    plt.setp(axes[0,0].get_yticklabels(),visible=False)
#    for ax in axes[:,0]:
#        plt.setp(ax.get_yticklabels(),visible=False)
    plt.setp([a.get_xticklabels() for a in axes[-1,1::2] ], visible=False)
    plt.setp([a.get_yticklabels() for a in axes[0::2, 0 ]], visible=False)
    for j in range(Ny):
        ax = axes[j,-1]
        label = r'$\ell/\Delta=%.1f$' % eldels[j]
        ax.text(1.,0.5,label,transform=ax.transAxes)
    fig.set_figheight(4.)
    fig.set_figwidth(8.)
    outname = 'shu_GPR%d_N%d' % (R,N)
    plt.savefig(outname,dpi=200)

def make_plot54(R, N):
    xlim = [4.5,7.5]
    ylim = [2.5,5.]
    eldels = [6,12,24,32]
    sigdels = [1.5,2,3,4]#,5]
    eldels = [6,32]

    Nx = 2 #4 columns
    Ny = 4 #5 rows
    fig, axes = plt.subplots(Ny,Nx, sharex='all',sharey='all')
    for i in range(Nx):
        #loop over eldels
        title = r'$\ell/\Delta = %.1f$' % eldels[i]
        axes[0,i].set_title(title)
        for j in range(Ny):
            #loop over sigdels
            pars = shu_pars(eldels[i],sigdels[j],N,R)
            pars.add_plot(axes[j,i],xlim,ylim)
    fig.subplots_adjust(hspace=0.,wspace=0.)
    title = 'GP-R%d; N=%d' % (R,N)
    plt.suptitle(title)
#    plt.setp(axes[0,0].get_yticklabels(),visible=False)
#    for ax in axes[:,0]:
#        plt.setp(ax.get_yticklabels(),visible=False)
    plt.setp([a.get_xticklabels() for a in axes[-1,1::2] ], visible=False)
    plt.setp([a.get_yticklabels() for a in axes[1::2, 0 ]], visible=False)
    fig.set_figheight(8.)
    fig.set_figwidth(4.)
    outname = 'shu_GPR%d_N%d' % (R,N)
    plt.savefig(outname,dpi=200)
#make_plot54(2,200)
#plt.show()

RR = [2,3]
NN = [200,400,800]
#RR = [2]
#NN = [200]
for R in RR:
    for N in NN:
        print R, N
        make_plot24(R,N)


