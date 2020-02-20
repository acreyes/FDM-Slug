import numpy as np
import matplotlib.pyplot as plt
import mayavi.mlab as mvi
import sys
sys.path.insert(0, '../Python')
from slugger3D import *



N = 6
name = "windtunnel2D_1%04d" % N
var = 0

Sdat = slug_data(name)
#print (Sdat.data[:,:,10,0])
#obj = Sdat.cuts_var(var)
obj = Sdat.iso_var(var,contours=10,opacity=.2,colormap='blue-red')
mvi.show()
