import numpy as np
import matplotlib.pyplot as plt


cmaps = [('Perceptually Uniform Sequential', [
            'viridis', 'plasma', 'inferno', 'magma']),
         ('Sequential', [
            'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
            'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
            'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']),
         ('Sequential (2)', [
            'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
            'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
            'hot', 'afmhot', 'gist_heat', 'copper']),
         ('Diverging', [
            'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
            'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic']),
         ('Qualitative', [
            'Pastel1', 'Pastel2', 'Paired', 'Accent',
            'Dark2', 'Set1', 'Set2', 'Set3',
            'tab10', 'tab20', 'tab20b', 'tab20c']),
         ('Miscellaneous', [
            'flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
            'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg', 'hsv',
            'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral', 'gist_ncar'])]



nx=1600
ny=400
nFiles = 26

maxDens = 25.4411016001
minDens = 1.38882060577

"""


for i in range(0,nFiles,1):
    if i < 10:
        gpM= "slug_dmr_gp1d_12_0p7_16x4_1000"+str(i)+".dat"
    else:
        gpM= "slug_dmr_gp1d_12_0p7_16x4_100"+str(i)+".dat"

    print i
    data = np.loadtxt(gpM)
    z = data[:,2]
    
    if z.max() > maxDens:
        maxDens = z.max()
    if z.min() < minDens:
        minDens = z.min()

    print minDens,maxDens
    
    #contours = np.linspace(0., 25., 30)

'''
    plt.pcolormesh(X,Y,Z,cmap='PiYG')
    #plt.contour(X,Y,Z,contours)
    plt.axis('scaled')
    #plt.axis('equal')
    #plt.xlim(0.,4.)
    #plt.xlim(2.7,3.5)
    plt.ylim(0.,1.)
    plt.title('GP-DMR')
    #plt.show()
    plt.savefig(str(i)+'.png')
'''
"""

for i in range(0,nFiles,1):
    dt = 0.01
    if i < 10:
        #gpM= "slug_dmr_gp1d_12_0p7_16x4_1000"+str(i)+".dat"
        gpM="slug_weno_p7_16x4_1000"+str(i)+".dat"
    else:
        #gpM= "slug_dmr_gp1d_12_0p7_16x4_100"+str(i)+".dat"
        gpM="slug_weno_p7_16x4_100"+str(i)+".dat"
    print i
    data = np.loadtxt(gpM)

    if i == 0:
        x = data[:,0]
        y = data[:,1]
        X = np.reshape(x,(nx,ny))
        Y = np.reshape(y,(nx,ny))

        
    z = data[:,2]
    Z = np.reshape(z,(nx,ny))
    
    #contours = np.linspace(0., 25., 30)

    colorCode=cmaps[5][1][-1]
    colorCode=cmaps[2][1][4]
    ##PiYG = cmaps[3][1][0]
    plt.pcolormesh(X,Y,Z,cmap=colorCode,vmin=minDens,vmax=maxDens)
    #plt.colorbar()
    #plt.contour(X,Y,Z,contours)
    plt.axis('scaled')
    #plt.axis('equal')
    #plt.xlim(0.,4.)
    #plt.xlim(2.7,3.5)
    plt.ylim(0.,1.)
    t=i*dt
    #plotTitle='GP-DMR, time='+str(t)+' sec'
    plotTitle='WENO-DMR, time='+str(t)+' sec'
    plt.title(plotTitle)
    #plt.show()

    #plt.savefig(str(i)+'.png', dpi=1200)
    #plt.close(str(i)+'.png')

    plt.savefig(str(i)+'_weno.png', dpi=1200)
    plt.close(str(i)+'_weno.png')
