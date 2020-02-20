import numpy as np
import os

class shu_pars():
    def __init__(self, eldel, sigdel, nx, R):
        self.name = 'shu_default'
        self.eldel = eldel
        self.sigdel = sigdel
        self.nx = nx
        self.R = R
        self.get_name()

    def get_name(self):
        self.name = 'shu_%d_R%d_el%.1f_sig%.1f' % (self.nx, self.R, self.eldel, self.sigdel)

    def make_init(self):
        self.get_name()

        #load default pars
        fid = open('shu.init', 'r')
        fwr = open('slug.init', 'w')
        for line in fid:
            fwr.write(line)
        line = "\nsim_name '%s' \n" % self.name
        fwr.write(line)

        line = 'gp_eldel %.1f \n' % self.eldel
        fwr.write(line)
        line = 'gp_sigdel %.1f \n' % self.sigdel
        fwr.write(line)
        line = 'gp_radius %d \n' % self.R
        fwr.write(line)
        line = 'gr_nx %d \n' % self.nx
        fwr.write(line)

        fid.close()
        fwr.close()

    def run_shu(self):
        self.make_init()
        print self.name
        os.system('../../slugEuler1d -> dump.dat')


N = 200
eldel = 12
sigdel = 2
R = 2

NN = [200, 400, 800]
eldels = [6,12,24,32]
sigdels = [1.5,2,3,4,5]
RR = [2,3]

for R in RR:
    for N in NN:
        for eldel in eldels:
            for sigdel in sigdels:
                
                pars = shu_pars(eldel, sigdel, N, R)
                pars.run_shu()
            print ''
        print ''
    print ''

