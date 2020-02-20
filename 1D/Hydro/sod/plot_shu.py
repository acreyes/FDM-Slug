import numpy as np
import matplotlib.pyplot as plt

def add_plot(fname, ax, **kwargs):
    data = np.loadtxt(fname)
    ax.plot(data[:,0],data[:,2],**kwargs)

#fig, ax = plt.subplot(2,2)

eldel = [3,6,12,50]
eldel = [6,12]
colors = ['g','r','b']
eldel = [12]

for ldel in eldel:
    fig, axis = plt.subplots()
    fname = "slug_sod_lamb0_10001.dat"
    add_plot(fname,axis,label=r'$\lambda=0$',color='k',marker='^')
    fname = "slug_sod_lamb1_10001.dat"
    add_plot(fname,axis,label=r'$\lambda=1$',c='cyan',marker='+')
    fname = "slug_sod_lamb10_10001.dat" 
    add_plot(fname,axis,label=r'$\lambda=10$',c='r',marker='x')
    title = r'$\ell/\Delta=%d$' % ldel

    axis.set_title(title)
    plt.legend()
    #plt.xlim(4.5,7.5)
    #plt.ylim(2.5,5.)
    #savename = 'shu_zoom_R3_eldel%d.png' % ldel
    #plt.savefig(savename)


plt.show()


