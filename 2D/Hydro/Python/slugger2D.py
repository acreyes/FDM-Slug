#supress warning due to deprecated use of numpy inside of h5py
#needed as of h5py version 2.7
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import h5py as hdf
warnings.resetwarnings()
import numpy as np
import matplotlib.pyplot as plt

class slug_data():
    def __init__(self, name):
        self.fname = "%s.slug" % name
        self.data = []
        self.xCoord = []
        self.yCoord = []
        self.betas = []
        self.get_data()


    def get_data(self):
        f = hdf.File(self.fname, "r")
        self.data = np.array(f.get('prim_vars'))
        self.betas = np.array(f.get('betas'))
        self.xCoord = np.array(f.get('xCoord'))
        self.yCoord = np.array(f.get('yCoord'))
        f.close()
        return

    def plot_var(self,  var, ax, **kwargs):
        xx, yy = np.meshgrid(self.xCoord, self.yCoord)
        zz = self.data[:,:,var]
        cmap = ax.pcolormesh(xx,yy,zz,**kwargs)
        return(cmap)

    def contours(self, var, ax, **kwargs):
        xx, yy = np.meshgrid(self.xCoord, self.yCoord)
        zz = self.data[:,:,var]
        CS = ax.contour(xx, yy, zz, **kwargs)
        return(CS)

    def vel_arrows(self, ax, skip, scale, rel):
        i = skip[0]
        j = skip[1]
        xx, yy = np.meshgrid(self.xCoord[::i], self.yCoord[::j])
        u = self.data[::i,::j,1] - rel[0]
        v = self.data[::i,::j,2] - rel[1]
        ax.quiver(xx,yy,u,v,scale=scale,pivot='mid')

    def var_label(self,var):
        labels = [r'$\rho$',r'$u$',r'$v$',r'$P$']
        return(labels[var])
