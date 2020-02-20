import numpy as np
import matplotlib.pyplot as plt

def add_plot(fname, ax, **kwargs):
    data = np.loadtxt(fname)
    ax.plot(data[:,0],data[:,1],**kwargs)

#fig, ax = plt.subplot(2,2)

eldel = [3,6,12,50]
#eldel = [6,12]
colors = ['g','r','b']
eldel = [3]

for ldel in eldel:
    fig, axis = plt.subplots()

    fname = "slug_shu_ref_10001.dat"
    add_plot(fname,axis,label=r'REF',color='k')
    
    fname = "slug_shu_200_R2_sig4.0_el12.0_10001" 
    add_plot(fname,axis,color='r',marker='x',label='GP-R2')#,mec='k')
    fname = "slug_shu_weno_10001.dat"
    add_plot(fname,axis,c='cyan',marker='+',label='WENO-JS')#,mec='k')
    title = r'$\ell/\Delta=%d$' % ldel

    
    #axis.set_title(title)
    plt.legend()
#    plt.xlim(4.5,7.5)
#    plt.ylim(2.5,5.)
    savename = 'shu_R2-weno_eldel%d_full.png' % ldel
    plt.savefig(savename)


plt.show()


