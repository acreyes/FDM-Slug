import numpy as np

def eoc(x1,x2,y1,y2):
    return(-(np.log(y2) -np.log(y1) )/(np.log(x2)-np.log(x1)))

def str_L1(err):
    num = '%.2e' % err
    splits = num.split('e')
    return('%s \\times 10^{%s}' % (splits[0], splits[1]))

def str_line(L1, OC):
    line = "$ %s $ & %s " % (str_L1(L1), OC)
    return (line)

fname = "GPR1-RK4_red.dat"
gpr1 = np.loadtxt(fname)

fname = "GPR2-RK4_red.dat"
gpr2 = np.loadtxt(fname)

fname = "GPR3-RK4_red.dat"
gpr3 = np.loadtxt(fname)


fname = "WENO-RK4_red.dat"
weno = np.loadtxt(fname)

N = len(gpr2[:,0])
var = 1
print 'GP'
for i in np.arange(0,N,1):
    if i == 0 :
        OC1 = '--'
        OC2 = '--'
        OC3 = '--'
        OCW = '--'
    else :
        OC1 = "%.2f" % eoc(gpr1[i-1,0],gpr1[i,0],gpr1[i-1,1],gpr1[i,1])
        OC2 = "%.2f" % eoc(gpr2[i-1,0],gpr2[i,0],gpr2[i-1,1],gpr2[i,1])
        OC3 = "%.2f" % eoc(gpr3[i-1,0],gpr3[i,0],gpr3[i-1,1],gpr3[i,1])
        OCW = "%.2f" % eoc(weno[i-1,0],weno[i,0],weno[i-1,1],weno[i,1])
    line1 =  str_line(gpr1[i,var], OC1)
    line2 =  str_line(gpr2[i,var], OC2)
    line3 =  str_line(gpr3[i,var], OC3)
    lineW =  str_line(weno[i,var], OCW)
    #line = "%d & $%s $ & %s & %d & $%s $ & %s \\\\" % ( int(gpr2[i,0]), str_L1(gpr2[i,1]), OC2,int(weno[i,0]), str_L1(weno[i,1]), OCW)
    line = "%d & %s & %s %s & %s \\\\" % (int(gpr2[i,0]), line1, line2, line3, lineW)
    print line

