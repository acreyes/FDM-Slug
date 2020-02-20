subroutine hll(vL,vR,Flux)

#include "definition.h"  

  use grid_data
  use primconsflux, only : prim2flux,prim2cons


  implicit none
  real, dimension(NUMB_VAR), intent(IN) :: vL,vR !prim vars
  real, dimension(NSYS_VAR), intent(OUT):: Flux 

  real, dimension(NSYS_VAR) :: FL,FR,uL,uR
  real :: sL,sR,aL,aR

  call prim2flux(vL,FL)
  call prim2flux(vR,FR)
  call prim2cons(vL,uL)
  call prim2cons(vR,uR)

  
  ! left and right sound speed a
  aL = sqrt(vL(GAMC_VAR)*vL(PRES_VAR)/vL(DENS_VAR))
  aR = sqrt(vR(GAMC_VAR)*vR(PRES_VAR)/vR(DENS_VAR))

  ! fastest left and right going velocities
  sL = min(vL(VELX_VAR) - aL,vR(VELX_VAR) - aR)
  sR = max(vL(VELX_VAR) + aL,vR(VELX_VAR) + aR)

  ! numerical flux
  if (sL >= 0.) then
     Flux(DENS_VAR:ENER_VAR) = FL(DENS_VAR:ENER_VAR)
  elseif ( (sL < 0.) .and. (sR >= 0.) ) then
     Flux(DENS_VAR:ENER_VAR) = (    sR*FL(DENS_VAR:ENER_VAR) &
                                   -sL*FR(DENS_VAR:ENER_VAR) &
                               +sR*sL*(uR(DENS_VAR:ENER_VAR) &
                                      -uL(DENS_VAR:ENER_VAR)))/(sR-sL)
  else
     Flux(DENS_VAR:ENER_VAR) = FR(DENS_VAR:ENER_VAR)
  endif
  return
end subroutine hll
