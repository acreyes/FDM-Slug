import numpy as np
import matplotlib.pyplot as plt

class data():
    def __init__(self, R):
        basen = 'slug_sod_gpr%d_' % R
        zeron = basen + 'zero_10001.dat'
        nonzeron = basen + 'nonzero_10001.dat'
        highn = basen+'drop_10001.dat'
        lown  = basen+'drop_high_10001.dat'
        bothn = basen+'drop_both_10001.dat'


        self.R = R
        if (R == 4) : self.R = 3
        self.zero = self.get_data(zeron)
        self.nonzero = self.get_data(nonzeron)
        self.high = self.get_data(highn)
        self.low = self.get_data(lown)
        self.both = self.get_data(bothn)

    def get_data(self, fname):
        fdata = np.loadtxt(fname)
        return(fdata[:,0:2])

    def add_plot(self, ax, data, **kwargs):
        ax.plot(self.zero[:,0], self.zero[:,1], 'k-',label='zero mean')
        ax.plot(data[:,0],data[:,1],'r+', **kwargs)


    def make_plot(self):
        fig, ax = plt.subplots(4, sharex=True, sharey=True)
        ax[0].set_title('GP-R%d' % self.R)
        self.add_plot(ax[0], self.nonzero, label='non-zero mean')
        self.add_plot(ax[1], self.high, label='drop largest')
        self.add_plot(ax[2], self.low, label='drop lowest')
        self.add_plot(ax[3], self.both, label='drop both')
        for axes in ax:
            axes.legend()
        fig.subplots_adjust(hspace=0)
        fig.set_figheight(10.)
        fig.set_figwidth(4.5)
        outname = 'sod_betas_GPR%d' % self.R
        plt.savefig(outname,dpi=200)
        return(fig, ax)

R2 = data(2)
R3 = data(4)

fig1, ax1 = R2.make_plot()
fig2, ax2 = R3.make_plot()
plt.show()
