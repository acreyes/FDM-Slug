#supress warnings from vtk related to mayavi and h5py related to numpy
#as of mayavi 4.5
###### h5py   2.7
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import vtk
vtk.vtkObject.GlobalWarningDisplayOff()
import numpy as np
import matplotlib.pyplot as plt
import mayavi.mlab as mvi
import h5py as hdf


class slug_data():
    def __init__(self, name):
        self.fname = "%s.slug" % name
        self.name = name
        self.data = []
        self.N = []
        self.extent = []
        self.x = []
        self.y = []
        self.z = []
        self.xCoord = []
        self.yCoord = []
        self.zCoord = []
        self.bdry  = []
        self.get_data()

    def get_indx(self,xpt, ax):
        #return the index in axis ax of coordinate xpt
        axes = [self.x, self.y, self.z]
        Coord = axes[ax]
        indx = (np.abs(Coord-xpt)).argmin()
        return(indx)

    def make_mgrid(self, xx, yy, zz):
        Nx = len(xx)
        xmin = xx[0]
        xmax = xx[Nx-1]
        Ny = len(yy)
        ymin = yy[0]
        ymax = yy[Ny-1]
        Nz = len(zz)
        zmin = zz[0]
        zmax = zz[Nz-1]
        self.extent = [xmin,xmax,ymin,ymax,zmin,zmax]
        Nx = 1j*Nx
        Ny = 1j*Ny
        Nz = 1j*Nz
        #return(np.mgrid[zmin:zmax:Nz,ymin:ymax:Ny,xmin:xmax:Nx])
        return(np.mgrid[xmin:xmax:Nx,ymin:ymax:Ny,zmin:zmax:Nz])

    def get_data(self):
        f = hdf.File(self.fname,"r")
        
        x = np.array(f.get('xCoord'))
        y = np.array(f.get('yCoord'))
        z = np.array(f.get('zCoord'))
        self.N = np.array([len(x), len(y), len(z)])

        self.x, self.y, self.z = (x,y,z)
        self.xCoord, self.yCoord, self.zCoord = self.make_mgrid(x,y,z)
        
        
        data = np.array(f.get('prim_vars'))
        Nvars = len(data[0,0,0,:])
        for var in range(Nvars):
            #transpose to account for indexing difference between hdf5 & fortran
            data_var = data[:,:,:,var]
            self.data.append((data_var.transpose()))#[::-1,::-1,:])

        try:
            bdry = np.array(f.get('bdry_var'))
            self.bdry = bdry[:,:,:,0].transpose()
        except:
            pass
        f.close()
        return

    def iso_var(self, var, **kwargs):
        obj = mvi.contour3d(self.xCoord,self.yCoord,self.zCoord,self.data[var],**kwargs)
        return obj

    def iso_sec(self, var, **kwargs):
        N34 = self.N/2
        x = self.xCoord[N34[0]:,N34[1]:,N34[2]:]
        y = self.yCoord[N34[0]:,N34[1]:,N34[2]:]
        z = self.zCoord[N34[0]:,N34[1]:,N34[2]:]

        V = self.data[var][N34[0]:,N34[1]:,N34[2]:]

        obj = mvi.contour3d(x,y,z,V,**kwargs)
        return obj

    def bdry_surf(self, **kwargs):
        surf = mvi.contour3d(self.xCoord,self.yCoord,self.zCoord,self.bdry,contours=[1.],**kwargs)

    def cut_kwrd(self,ax):
        
        #default kwargs for plane cuts for use in self.cut_var
        axes = ['x_axes','y_axes','z_axes']
        kwargs = {}
        kwargs['plane_orientation'] = axes[ax]
        kwargs['slice_index'] = int(self.N[ax]/2)
#        kwargs['extent'] = self.extent
        return(kwargs)
    
    def cut_var(self, var, **kwargs):

        obj = mvi.pipeline.image_plane_widget(mvi.pipeline.scalar_field(self.xCoord,self.yCoord,self.zCoord,self.data[var]),**kwargs)
                  
        return obj

