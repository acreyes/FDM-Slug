import numpy as np
import os
import matplotlib.pyplot as plt
import scipy.interpolate as intp

def get_data(N):
    
    items = ['dens','velx','vely']
    dx = 9./float(N)

    #get initial data
    fname = "slug_shu_ref_10001.dat"
    ref_dat   = np.loadtxt(fname)
    
    #get final data
    fname = "slug_shu_10001.dat"
    finl_dat   = np.loadtxt(fname)

    nvar = 1

    f = intp.interp1d(ref_dat[:,0],ref_dat[:,nvar])
    x = finl_dat[:,0]
    xx = x[x<6.5]
    dens = finl_dat[:,1]
    ref  = f(xx)
    finl = dens[x<6.5]

    L1 = np.sum(np.abs(finl-ref))*dx

        
    return(L1)

def make_init(ell,N):
    fname = "ell_test.init"
    frd = open(fname, 'r')
    fwr = open('slug.init', 'w')

    for line in frd:
        fwr.write(line)
    line = "gp_eldel %e \n" % ell
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
    N_ell = 50
    Lmin = 1.
    Lmax = 10.
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
NN = [32,64,128,256]#,512]
for N in NN:
    run_ell(N)
