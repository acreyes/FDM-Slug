import numpy as np
import os
import matplotlib.pyplot as plt

def get_data(N):
    
    items = ['dens','velx','vely']
    dx = 1./float(N)

    #get initial data
    fname = "slug_gauss_10000.dat"
    init_dat   = np.loadtxt(fname)
    
    #get final data
    fname = "slug_gauss_10001.dat"
    finl_dat   = np.loadtxt(fname)

    nvar = 1

    init = init_dat[:,nvar]
    finl = finl_dat[:,nvar]

    ERR = 0.
    tmax = 2.0
    u=1.
    v=1.
    gamm = 1.4
    beta = 5.
    maxL1 = 0.
        
    L1 = 0.
    for i in range(len(init)):
        ref = init[i]

        delta = abs(finl[i]-ref)
        ERR += (delta)
        if (delta > maxL1): maxL1 = (delta)

        L1 = ERR*dx

        
    return(L1)

def make_init(ell,N):
    fname = "ell_test.init"
    frd = open(fname, 'r')
    fwr = open('slug.init', 'w')

    for line in frd:
        fwr.write(line)
    line = "gp_ell %e \n" % ell
    fwr.write(line)
    line = "gr_nx %d" % N
    fwr.write(line)
    fwr.close()
    frd.close()
    return

def error_ell(ell,N):
    make_init(ell,N)
    os.system("../slugEuler1d -> dump.dat")
    L1 = get_data(N)
    print ell, L1
    return(L1)

def run_ell(N):
    errors = []
    N_ell = 75
    Lmin = .01
    Lmax = 0.6
    LL = np.linspace(Lmin, Lmax, N_ell)
    for L in LL:
        L1 = (error_ell(L,N))
        errors.append(L1)
    fname = "ell_RK3_%d.dat" % N
    fout = open(fname, 'w')
    for i in range(len(LL)):
        line = "%e %e\n" % (LL[i], errors[i])
        fout.write(line)
    fout.close()

NN = [25,50,100,200]
N = 25
for N in NN:
    run_ell(N)
