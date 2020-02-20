#FC	= gfortran-mp-4.8 -O2 -g
FC	= gfortran

LDFLAGS = -framework accelerate

FFLAGS_OPT2 = -ggdb -g -O3 -fdefault-real-8 -fdefault-double-8\
	-ffree-line-length-none -Wuninitialized

# quadruple precision
FFLAGS_OPT4 = -ggdb -g -O3 -freal-4-real-16 -freal-8-real-16\
	-ffree-line-length-none -Wuninitialized

# double precision
FFLAGS_DEBUG2 = -ggdb  -g -fdefault-real-8 -fdefault-double-8\
#	-ffree-line-length-none -Wuninitialized

# quadruple precision
FFLAGS_DEBUG4 = -ggdb  -g -freal-4-real-16 -freal-8-real-16\
	-ffree-line-length-none -Wuninitialized

FFLAGS_DOUB = $(FFLAGS_OPT2)
FFLAGS_QUAD = $(FFLAGS_OPT4)

EXE_FILE = slugEuler1d
OBJS  = driver_euler1d.o \
	read_initFile.o\
	read_pars.o \
	sim_data.o  \
	sim_init.o \
	sim_initBlock.o \
	gp_data.o \
	sim_GPinit.o \
	gp_eigens.o \
	gp_Fluxinit.o \
	grid_data.o \
	grid_init.o \
	grid_finalize.o\
	io.o\
	eos.o\
	primconsflux.o \
	soln_ReconEvolveAvg.o \
	soln_RK4.o \
	soln_RK3.o \
	soln_RK2.o \
	soln_reconstruct.o \
	soln_getFlux.o \
	soln_gpFlux.o \
	hll.o \
	hllc.o \
	roe.o \
	sim_interfaces.o \
	soln_update.o \
	soln_FOG.o \
	soln_PLM.o \
	soln_PPM.o \
	soln_WENO.o \
	soln_GP.o \
	bc.o \
	cfl.o \
	eigensystem.o \
	averageState.o \
	char_tracing.o \
	GP.o \
	WENO.o \
	linalg.o \
	slopeLimiter.o

########################################################################################
#COMPILING AND LINKING USING GENERIC SUFFIX RULE FOR F90

$(EXE_FILE) : $(OBJS)
	@$(FC) $(FFLAGS_DOUB) $(OBJS) -o $(EXE_FILE) $(LDFLAGS)
	@echo "code is now linking..."

#some GP routines need quad precision
GP.o: %.o : %.F90
	$(FC) $(FFLAGS_QUAD) -c $<

sim_GPinit.o: %.o : %.F90
	$(FC) $(FFLAGS_QUAD) -c $<

gp_Fluxinit.o: %.o : %.F90
	$(FC) $(FFLAGS_QUAD) -c $<

linalg.o: %.o : %.F90
	$(FC) $(FFLAGS_QUAD) -c $<


#LET'S APPLY GENERIC SUFFIX RULE HERE FOR FORTRAN 90
.SUFFIXES : 
.SUFFIXES : .F90 .o

.F90.o:
	$(FC) $(FFLAGS_DOUB) -c $<

#######################################################################################
#SOME USEFUL COMMANDS
clean:
	@rm -f *.o *.mod *~ slugEuler1d

#######################################################################################
#LET'S DEFINE SOME MODULE DEPENDENCIES!
driver_euler1d.o: sim_data.o grid_data.o io.o bc.o eos.o

read_pars.o     : grid_data.o sim_data.o gp_data.o read_initFile.o

eos.o		: grid_data.o sim_data.o

grid_init.o	: grid_data.o read_initFile.o
grid_finalize.o : grid_data.o

hll.o		: grid_data.o primconsflux.o
roe.o		: grid_data.o primconsflux.o eigensystem.o

io.o		: grid_data.o sim_data.o


primconsflux.o  : grid_data.o eos.o

hll.o		: grid_data.o primconsflux.o
hllc.o          : grid_data.o primconsflux.o
roe.o		: grid_data.o primconsflux.o eigensystem.o

soln_getFlux.o  	: grid_data.o sim_data.o primconsflux.o sim_interfaces.o

sim_init.o	: sim_data.o read_initFile.o
sim_initBlock.o : sim_data.o grid_data.o primconsflux.o

gp_data.o               : GP.o
sim_GPinit.o 		: linalg.o GP.o gp_data.o
gp_Fluxinit.o           : linalg.o sim_data.o grid_data.o GP.o gp_data.o
soln_gpFlux.o           : gp_data.o sim_data.o
WENO.o                  : sim_data.o

soln_update.o		: grid_data.o primconsflux.o
soln_ReconEvolveAvg.o 	: grid_data.o sim_data.o primconsflux.o
soln_reconstruct.o 	: grid_data.o sim_data.o eigensystem.o primconsflux.o sim_interfaces.o char_tracing.o
soln_FOG.o		: grid_data.o
soln_PLM.o		: grid_data.o sim_data.o slopeLimiter.o eigensystem.o
soln_PPM.o              : grid_data.o sim_data.o slopeLimiter.o eigensystem.o char_tracing.o
soln_WENO.o             : grid_data.o sim_data.o eigensystem.o  char_tracing.o WENO.o
soln_GP.o		: grid_data.o sim_data.o eigensystem.o WENO.o gp_data.o

#######################################################################################
