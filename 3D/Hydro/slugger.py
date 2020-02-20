import numpy as np
import matplotlib.pyplot as plt
import mayavi.mlab as mvi
import h5py as hdf

class slug_data():
    def __init__(self, name):
        self.fname = "%s.slug" % name
        self.name = name
        self.data = []
        self.xCoord = []
        self.yCoord = []
        self.zCoord = []
        self.get_data()

    def get_data(self):
        f = hdf.File(self.fname,"r")
        self.data = np.array(f.get('prim_vars'))
        self.xCoord = np.array(f.get('xCoord'))
        self.yCoord = np.array(f.get('yCoord'))
        self.zCoord = np.array(f.get('zCoord'))
        f.close()
        return

    def plot_var(self, var):
        #obj = mvi.contour3d(self.data[:,:,:,var],contours=10,opacity=.4,colormap='blue-red')
        obj = mvi.pipeline.image_plane_widget(mvi.pipeline.scalar_field(self.data[:,:,:,var]),
                            plane_orientation='x_axes',
                            slice_index=10,
                        )
        obj = mvi.pipeline.image_plane_widget(mvi.pipeline.scalar_field(self.data[:,:,:,var]),
                            plane_orientation='y_axes',
                            slice_index=10,
                        )

name = "vortex_10000"
Sdat = slug_data(name)
#print (Sdat.data[:,:,10,0])

obj = Sdat.plot_var(0)
mvi.show()
