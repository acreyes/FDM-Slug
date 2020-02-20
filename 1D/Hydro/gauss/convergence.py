import numpy as np
import matplotlib.pyplot as plt


class L1_vars(object):
    def __init__(self):
        self.dens = []
        self.velx = []
        self.vely = []
        self.pres = []

def get_fname(N, i):
    fname = "slug_G_GP-R2_RK4_%d_red_1%04d.dat" % (N, i)
    fname = "slug_G_WENO_RK3_%d_nored_1%04d.dat" % (N, i)
    return(fname)

def get_data(N):

    items = ['dens','velx','vely']
    dx = 1./float(N)
    dxdy = dx*dx
    L1 = {}
    #get initial data
    fname = get_fname(N, 0)
    init_dat   = np.loadtxt(fname)
    
    #get final data
    fname = get_fname(N, 1)
    finl_dat   = np.loadtxt(fname)

    nvars = [1,2,3]
    for nvar,item in zip(nvars, items):
        init = init_dat[:,nvar]
        finl = finl_dat[:,nvar]

        ERR = 0.
        tmax = 2.0
        u=1.
        v=1.
        gamm = 1.4
        beta = 5.
        maxL1 = 0.
        

        for i in range(len(init)):
            '''
            x = finl_dat[i,0] - (5.+tmax*u)
            y = finl_dat[i,1] - (5.+tmax*v)
            r2 = x**2+y**2
            T = 1. - (gamm-1.)*beta*beta*np.exp(1.-r2)/(8.*gamm*np.pi*np.pi)
            '''
            #ref = T**(1./(gamm-1.))
            ref = init[i]

            delta = abs(finl[i]-ref)
            ERR += (delta)
            if (delta > maxL1): maxL1 = (delta)

        L1[item] = ERR*dx
        if (item == 'dens'): print 'N= %d Linf = %e' % (N, maxL1)
        '''
    for item in items:
        init = init_dat[item][0,0,:,:]
        finl = finl_dat[item][0,0,:,:]

        #now calculate the L1-error
        ERR = 0.
        
        for i in range(N):
            for j in range(N):
                delta = abs(finl[i,j] - init[i,j])
                ERR += delta
        L1[item] = ERR*dxdy
        '''
    return(L1)

#NN = [ 32,64, 128, 256]#, 512]#, 1024]
#NN = [40,80,160]
NN = [25,50,100,200,400]
#HghOrd_F = L1_vars()
HghOrd_T = L1_vars()
items = ['dens','velx','vely']
for N in NN:
    HghOrd = '.True.'
    L1_T = get_data(N)
    HghOrd_T.dens.append(L1_T['dens'])
    HghOrd_T.velx.append(L1_T['velx'])
    HghOrd_T.vely.append(L1_T['vely'])
                     

#now lets write the data
#fname_T = 'GP-R2.dat'
#fname_T = 'WENO.dat'
#fname_T = 'WENO-Z.dat'
fname_T  = 'GPR2-RK4_red.dat'
fname_T = "WENO-RK3_nored.dat"
#fname_T = 'WENO_cntr.dat'
#fname_T = 'GPMD-R2.dat'
#fname_T = 'MD_cntr.dat'

fid_T = open(fname_T, 'w')
#fid_F = open(fname_F, 'w')

for i in range(len(NN)):
    line_T = '%d    %e   %e   %e \n' % (NN[i], HghOrd_T.dens[i], HghOrd_T.velx[i], HghOrd_T.vely[i])
    fid_T.write(line_T)
    #line_F = '%d    %e   %e   %e \n' % (NN[i], HghOrd_F.dens[i], HghOrd_F.velx[i], HghOrd_F.vely[i])
    #fid_F.write(line_F)

fid_T.close()
#fid_F.close()
