subroutine scal(vL,vR,Flux)

#include "definition.h"  

  use grid_data
  use primconsflux, only : prim2flux,prim2cons


  implicit none
  real, dimension(NUMB_VAR), intent(IN) :: vL,vR !prim vars
  real, dimension(NSYS_VAR), intent(OUT):: Flux
  
  Flux = 0.
  Flux = vL(DENS_VAR) !upwind scalar advection
  return
end subroutine scal
