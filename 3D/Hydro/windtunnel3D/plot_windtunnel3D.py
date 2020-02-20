import numpy as np
import matplotlib.pyplot as plt
import mayavi.mlab as mvi
import sys
sys.path.insert(0, '../Python')
from slugger3D import *



N = 7
name = "windtunnel3D_FOG_128_1%04d" % N
var = 0

Sdat = slug_data(name)

iso = True
if iso :
    obj = Sdat.iso_var(var,contours=10,opacity=.4,colormap='blue-red')
    bdry = Sdat.bdry_surf(color=(0,0,0))
else:
    #plane cuts
    kwargs = Sdat.cut_kwrd(1)
    cuts = Sdat.cut_var(var,**kwargs)
    kwargs = Sdat.cut_kwrd(0)
    cuts = Sdat.cut_var(var,**kwargs)
mvi.show()
