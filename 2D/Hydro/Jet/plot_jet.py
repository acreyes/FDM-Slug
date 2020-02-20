import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import sys
sys.path.insert(0, '../Python')
from slugger2D import *

N = 1
var = 0
# (0,1,2,3) (dens, u, v, P)
#for velocity arrows
skip = [10,10]
scale = 75
bckgrnd = [0., 0.]
#contour params
levels = np.linspace(1., 7., 15)


name = "jet_GP_224_1%04d" % N
data  = slug_data(name)
#data.data[:,:,0] = np.log(data.data[:,:,0])
fig, ax = plt.subplots()
cax = data.plot_var(var,ax)#,norm=clr.LogNorm(vmin=1.,vmax=7.))
#data.vel_arrows(ax,skip,scale,bckgrnd)
#CS = data.contours(var, ax, colors='k', levels=levels)
ax.set_aspect(aspect=1)
fig.colorbar(cax)
plt.title('FOG')
plt.show()
