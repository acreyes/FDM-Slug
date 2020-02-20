import numpy as np
import matplotlib.pyplot as plt

#gpM = "slug_dmr_gpM-ld12-R2_512_10001.dat"
gpM= "slug_dmr_gp1d_12_0p7_16x4_10025.dat"
#weno = "slug_dmr_w5_512_10001.dat"

data = np.loadtxt(gpM)
x = data[:,0]
y = data[:,1]
z = data[:,2]

nx=1600
ny=400

X = np.reshape(x,(nx,ny))
Y = np.reshape(y,(nx,ny))
Z = np.reshape(z,(nx,ny))

contours = np.linspace(0., 25., 30)


plt.pcolormesh(X,Y,Z,cmap='PiYG')
#plt.contour(X,Y,Z,contours)
plt.axis('scaled')
#plt.axis('equal')
#plt.xlim(0.,4.)
#plt.xlim(2.7,3.5)
plt.ylim(0.,1.)
plt.title('GP-multiD')


plt.show()
