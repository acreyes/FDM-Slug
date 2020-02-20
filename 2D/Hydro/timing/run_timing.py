import numpy as np
import os
import h5py as hdf
import matplotlib.pyplot as plt

class slug_sim():
    def __init__(self,NN,order,R,iprocs,jprocs):
        self.NN    = NN
        self.order = order
        self.R     = R
        self.iprocs= iprocs
        self.jprocs= jprocs
        if order == 10 :
            self.flux = '.true.'
        else:
            self.flux = '.false.'
        self.times = np.zeros(len(NN))
        self.error = np.zeros(len(NN))

    def make_init(self, N):
        frd = open("timing.init", 'r')
        fwr = open("slug.init"  , 'w')

        for line in frd:
            fwr.write(line)
        line = "gr_nx %d\ngr_ny %d\n" % (N,N)
        fwr.write(line)
        line = "sim_order %d\n" % self.order
        fwr.write(line)
        line = "gp_radius %d\n" % self.R
        fwr.write(line)
        line = "bl_iProcs %d\nbl_jProcs %d\n" % (self.iprocs,self.jprocs)
        fwr.write(line)
        line = "sim_gpFlux %s\n" % self.flux
        fwr.write(line)

        fwr.close()
        frd.close()

    def get_data(self,N):
        fname = "vortex_10000.slug"
        f = hdf.File(fname,'r')
        ref_data = np.array(f.get('prim_vars'))
        f.close()
        
        fname = "vortex_10001.slug"
        f = hdf.File(fname,'r')
        time = f.get('eTime')
        delT = time[0]
        new_data = np.array(f.get('prim_vars'))
        f.close()

        ref_dens = ref_data[:,:,0]
        new_dens = new_data[:,:,0]

        dx = 10./N
        dxdy = dx*dx

        L1 = dxdy*np.sum(np.abs(ref_dens-new_dens))
        return delT, L1

    def run_sims(self,fname):
        for i in range(len(self.NN)):
            self.make_init(NN[i])
            cmnd = 'mpirun -np %d ../slugEuler1d' % (self.iprocs*self.jprocs)
            os.system(cmnd)
            self.times[i], self.error[i] = self.get_data(NN[i])
        np.savetxt(fname, np.transpose([self.times,self.error]))


#these set the number of procs to use in the simulations
#nprocs = iprocs*jprocs
iprocs = 2
jprocs = 2

#resolutions to be run for each set up
NN = [25,50,100,200]

#these set parameters for each set up
#N refers to the number of setups
N = 4
##file names datas will be saved to
fname = ['timings_weno.dat','timings_gpr1.dat','timings_gpr2.dat','timings_gpr3.dat']
order = [5,10,10,10]
R     = [2,1,2,3]

#plotting things: note that data is saved so plot can be recreated outside of this script
fig, ax = plt.subplots()
marker = ['o','^','D','v']
label  = ['WENO-JS', 'GP-R1', 'GP-R2', 'GP-R3']
for i in range(5):
    sim = slug_sim(NN,order[i],R[i],iprocs,jprocs)
    sim.run_sims(fname[i])
    ax.plot(sim.times,sim.error,label=label[i],marker=marker[i])

plt.legend()
plt.ylabel(r'$L_1$ Error')
plt.xlabel('CPU Time (sec)')
plt.yscale('log')
plt.xscale('log')
plt.show()
