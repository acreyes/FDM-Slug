import numpy as np
import matplotlib.pyplot as plt
import os

def make_init(name, Nx, order, torder, rk, dt,sigma):
    sigdel = sigma
    '''
    sigdel = 0.
    if order != 5:
        sigdel = order
        order = 10
    radius = (order - 1)/2
    '''
    sigma = 0.
    fname = "init_conv.init"
    fread = open(fname, 'r')
    fwrite = open("slug.init", 'w')
    for line in fread:
        fwrite.write(line)
    #fname has now been copied to "slug.init"
    fwrite.write('\n')
    #fwrite.write("sim_name " + name + '\n')
    fwrite.write("sim_dt " + str(dt) + '\n')
    fwrite.write("gr_nx " + str(Nx) + '\n' )
    fwrite.write("sim_order " + str(order) + '\n')
    fwrite.write("sim_RK " + rk + '\n')
    fwrite.write("sim_Torder " + str(torder) + '\n')
    #fwrite.write("sim_sigdel " + str(sigdel) + '\n')
    fwrite.write("gr_radius " + str(sigdel))
    fwrite.write("\n")
    #fwrite.write('# IO frequency\nsim_ioTfreq 1.8 #real\nsim_ioNfreq -1    #positive integer; zero or negative if not used')
    #DONE

def L1(fname, N, dt):
    '''
    Nt = int(1./dt)
    t = Nt*dt
    if (abs(1.-t) > abs(1.-(t+dt))):
        t += dt
    elif (abs(1.-t) > abs(1.-(t-dt))):
        t -= dt
    #print t
    '''
    ref_dat = np.loadtxt('slug_gauss_sim_10000.dat')
    new_dat = np.loadtxt('slug_gauss_sim_10001.dat')
    

    Navg = 1.1024/N
    kbeg = 0
    kend = kbeg + Navg
    summ = 0.
    for i in range(len(new_dat[:,0])):
        x = new_dat[i,0]
        ref = sum(ref_dat[kbeg:kend,1])/Navg
        #delta = abs(anal(x,t) - new_dat[i,1])
        delta = abs(ref_dat[i,1] - new_dat[i,1])
        summ += delta
        kbeg += Navg
        kend += Navg
    return(summ/float(N))

def con(N, order, torder, rk, dt,sigma):
    #note that I don't need to bother with naming, but I'm just leaving this here because I'm lazy
    name ='sod_' + '_' + str(N)
    #make init file
    make_init(name, N, order, torder, rk, dt,sigma)
    #run simulation
    os.system("./slugeuler1d -> dump.dat")
    fname = 'slug_' + name + '_10001.dat'
    #calculate L1 error
    norm = L1(fname, N, dt)
    return(norm)

def logline(x, x0, y0, m):
    return(y0*(x/x0)**-m)

def conv_dt(dt0, N0, N, order, m):
    #calculates dt to use
    if order == 10:
        r = m*N/(N0)
    else:
        r = order
    dt = dt0*(float(N0)/float(N))**(float(r)/float(m))
    #this is dt reduction, make sure sim_fixdt is set to .true. in init_conv.init
    
    dt = dt0*(float(N0)/float(N))
    #print 'dt = ' + str(dt)
    return(dt)

def conv_data(order, torder, rk, sigma):
    print "spatial order: " + str(order) + "   torder: " + str(torder)
    dt0 = .0125#0.00625000 #gauss @ 32
    #dt0 = 0.02327189*16/32 #sod
    #dt0 = 0.09793092 #shu
    N0 = 16
    N = N0
    Nx = []
    data = []
    Ns = 6
    for i in range(Ns):
        dt = conv_dt(dt0, N0, N, order, torder)
        print str(i+1) + " of " + str(Ns)
        print dt
        norm = con(N, order, torder, rk, dt,sigma)
        data.append(norm)
        Nx.append(N)
        N = int(2*N)
    return(Nx, data)

#starting resolution and #of sims to run
#make sure that N0 and Ns match the same vars in conv_data function
##Their only use is for making the NN variable that is used for plotting loglines
N0 = 16
Ns = 6
Nf = N0*(2**Ns)
NN = np.linspace(N0, Nf, 500)
'''
#use this w/ dt reduction
rk_plm = [2]
rk_ppm = [2,3]
rk_wen = [2,3,4]
rk_gp  = [2,3,4]
'''
#use this one w/o dt reduction
#these are the rk-orders that will be run
rk_plm = [2,3,4]
rk_ppm = [2,3,4]
rk_wen = [2,3,4]
rk_gp  = [2,3,4]

rk_plm = [4]
rk_ppm = [2,4]
rk_wen = [4]
rk_gp  = [4]

rk     = [4]

plm = {}
ppm = {}
wen = {}
gp1 = {}
gp2 = {}
gp3 = {}
rk_plots = {}
rk_ax    = {}

for order in rk:
    #initialize all the figures and axis
    rk_plots[order] = plt.figure()
    rk_ax[order]    = rk_plots[order].add_subplot(111)
    rk_ax[order].set_xscale('log')
    rk_ax[order].set_yscale('log')
    #rk_ax[order].set_ylim([1.e-10,1e-1])
    tit = 'RK' + str(order)##
    plt.title(tit)

'''
for order in rk_plm:
    Nx, plm[order] = conv_data(2, order, ".true.", 2)
    rk_ax[order].plot(Nx, plm[order], 'ro', label='PLM')
    rk_ax[order].plot(NN, logline(NN, Nx[1], plm[order][1], 2), 'r--')

for order in rk_ppm:
    Nx, ppm[order] = conv_data(3, order, ".true.", 3)
    rk_ax[order].plot(Nx, ppm[order], 'b+', label='PPM')
    #rk_ax[order].plot(NN, logline(NN, Nx[1], ppm[order][1], 3), 'b--')
    #rk_ax[order].plot(NN, logline(NN, N0, ppm[order][0], 2), 'b--')

for order in rk_wen:
    Nx, wen[order] = conv_data(5, order, ".true.", 3)
    rk_ax[order].plot(Nx, wen[order], marker='o', linestyle='', color='k', label='WENO')
    #rk_ax[order].plot(NN, logline(NN, Nx[4], wen[order][4], 2), 'k--')
    #rk_ax[order].plot(NN, logline(NN, N0, wen[order][0], 4), 'g--')
'''
sigmas = [1, 2, 3]
markers = ['+','x','^','*']
for order in rk_gp:
    for sigma,mark in zip(sigmas,markers):
        Nx, gp1[order] = conv_data(10, order, ".true.", sigma)
        lab = r'R = ' + str(sigma)
        rk_ax[order].plot(Nx, gp1[order], marker=mark, linestyle='', label=lab)
    print Nx
    print gp1[order]
    #rk_ax[order].plot(Nx, gp2[order], marker='x', linestyle='-', color='b', label=r'$\frac{\sigma}{\Delta}=24$')
    #rk_ax[order].plot(Nx, gp3[order], marker='+', linestyle='-', color='k', label=r'$\frac{\sigma}{\Delta}=48$')
    #print gp[order]
    rk_ax[order].plot(NN, logline(NN, Nx[4], gp1[order][4], 1), 'k--')
    rk_ax[order].plot(NN, logline(NN, Nx[1], gp1[order][1], 2), 'g--')



for ax in rk_ax.itervalues():
    ax.legend()

plt.show()







