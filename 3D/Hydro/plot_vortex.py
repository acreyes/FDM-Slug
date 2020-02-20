import numpy as np
import matplotlib.pyplot as plt
import mayavi.mlab as mvi
import sys
sys.path.insert(0, 'Python')
from slugger3D import *



N = 1
name = "vortex_1%04d" % N
var = 0
Sdat = slug_data(name)
'''
Sdat = slug_data(name)
#print (Sdat.data[:,:,10,0])
obj = Sdat.iso_var(var,contours=10,opacity=.4,colormap='blue-red')
bdry = Sdat.bdry_surf(color=(0,0,0))
'''

#plane cuts
kwargs = Sdat.cut_kwrd(1)
cuts = Sdat.cut_var(var,**kwargs)
kwargs = Sdat.cut_kwrd(0)
cuts = Sdat.cut_var(var,**kwargs)
mvi.show()
