import numpy as np
import os

class sim_pars():
    def __init__(self, N, order):
        self.N = N
        self.order = order
        self.cfl = 0.8
        self.name = 'vortex'
        self.red_cfl = False
        self.get_name()
        self.get_cfl()

    def get_name(self):
        if self.order == 10:
            schm = "GP-R1_RK4"
        else:
            schm = "WENO_RK4"
        self.name = "G_%s_%d_red" % (schm, self.N)
        return()

    def get_cfl(self):
        #assums rk4 + 5th order in space
        m = 5.
        n = 4.
        N0 = 25.
        if self.red_cfl :
            self.cfl = (0.8)*(N0/float(self.N))**(m/n)
        return

    def make_init(self):
        self.get_name()
        self.get_cfl()
        frd = open("conv.init", 'r')
        fwr = open("slug.init", 'w')
        for line in frd:
            fwr.write(line)
        line = "sim_name '%s'\n" % self.name
        fwr.write(line)
        line = "gr_nx %d\n" % (self.N)
        fwr.write(line)
        line = "sim_cfl %e\n" % self.cfl
        fwr.write(line)
        line = "sim_order %d \n" % self.order
        fwr.write(line)
        if self.order == 10:
            line = "sim_gpFlux .true.\n"
        else:
            line = "sim_gpFlux .false.\n"
        fwr.write(line)
        fwr.close()
        frd.close()

    def run_sim(self):
        self.make_init()
        print self.N
        os.system("../slugEuler1d -> dump.dat")


    
NN = [25,50,100,200,400]
order = 5
sim = sim_pars(NN[0], order)
sim.cfl = 0.8
sim.red_cfl = True
for N in NN:
    sim.N = N
    sim.run_sim()


