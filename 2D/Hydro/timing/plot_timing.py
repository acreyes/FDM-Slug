import numpy as np
import matplotlib.pyplot as plt

fnames = ['timings_weno.dat','timings_gpr1.dat','timings_gpr2.dat','timings_gpr3.dat']
markers = ['o','^','D','v']
labels  = ['WENO-JS', 'GP-R1', 'GP-R2', 'GP-R3']

fig, ax = plt.subplots()

for fname, marker, label in zip(fnames,markers,labels):
    data = np.loadtxt(fname)
    plt.plot(data[:,0],data[:,1],label=label,marker=marker)


plt.legend()
plt.ylabel(r'$L_1$ Error')
plt.xlabel('CPU Time (sec)')
plt.yscale('log')
plt.xscale('log')
plt.show()
