module sim_data

#include "definition.h"

  implicit none

  !! numerics
  real, save :: sim_cfl, sim_tmax, sim_outputIntervalTime, sim_WENeps, sim_dt, sim_sigdel, sim_sigma
  integer, save :: sim_order, sim_nStep, sim_Torder, sim_mval
  character(len=MAX_STRING_LENGTH), save :: sim_name, sim_limiter,sim_riemann, sim_WENO, sim_quad
  logical, save :: sim_charLimiting, sim_RK, sim_fixDt, sim_nlim, sim_intFlux

  !! ICs
  real, save :: sim_gamma
  real, save :: sim_densL,sim_velxL,sim_presL
  real, save :: sim_densR,sim_velxR,sim_presR
  real, save :: sim_shockLoc
  real, save :: sim_smallPres
  character(len=MAX_STRING_LENGTH), save :: sim_icType

  !! BCs
  character(len=MAX_STRING_LENGTH), save :: sim_bcType

  !! IO
  integer, save :: sim_ioNfreq
  real,    save :: sim_ioTfreq

end module sim_data
