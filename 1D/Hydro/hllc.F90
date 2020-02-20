subroutine hllc(vL,vR,Flux)

#include "definition.h"  

  use grid_data
  use primconsflux, only : prim2flux,prim2cons


  implicit none
  real, dimension(NUMB_VAR), intent(IN) :: vL,vR !prim vars
  real, dimension(NSYS_VAR), intent(OUT):: Flux 

  real, dimension(NSYS_VAR) :: FL,FR,uL,uR,Uhll,UstarL,UstarR
  real :: sL,sR,aL,aR,uStar,pStar,pTot
  real :: numerL, denomL, numerR, denomR
  real :: dStarL, dStarR

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

  ! Get HLL states for later use
  if (sL > 0.) then
     Uhll(DENS_VAR:ENER_VAR) = uL(DENS_VAR:ENER_VAR)
  elseif ((sL <= 0.) .and. (sR >= 0.)) then
     Uhll(DENS_VAR:ENER_VAR) = &
          ( sR*uR(DENS_VAR:ENER_VAR) &
           -sL*uL(DENS_VAR:ENER_VAR) &
             - FR(DENS_VAR:ENER_VAR) &
             + FL(DENS_VAR:ENER_VAR)&
           )/(sR - sL)
  else
     Uhll(DENS_VAR:ENER_VAR) = uR(DENS_VAR:ENER_VAR)
  endif

  ! Get uStar
  uStar = vR(DENS_VAR)*vR(VELX_VAR)*(sR-vR(VELX_VAR)) &
         -vL(DENS_VAR)*vL(VELX_VAR)*(sL-vL(VELX_VAR)) &
         +vL(PRES_VAR)-vR(PRES_VAR)
  uStar = uStar/( vR(DENS_VAR)*(sR-vR(VELX_VAR)) &
                 -vL(DENS_VAR)*(sL-vL(VELX_VAR)))

  ! Convenient parameters
  numerL = sL-vL(VELX_VAR)
  denomL = sL-uStar
  numerR = sR-vR(VELX_VAR)
  denomR = sR-uStar

  ! Get pStar
  pStar = vL(DENS_VAR)*numerL*(uStar-vL(VELX_VAR)) + vL(PRES_VAR)


  ! density
  dStarL = uL(DENS_VAR)*numerL/denomL
  dStarR = uR(DENS_VAR)*numerR/denomR

  ! left and right star regions
  UstarL(DENS_VAR) = dStarL
  UstarL(MOMX_VAR) = dStarL*uStar
  UstarL(ENER_VAR) = uL(ENER_VAR)*numerL/denomL+&
       (pStar*uStar - vL(PRES_VAR)*vL(VELX_VAR))/denomL

  UstarR(DENS_VAR) = dStarR
  UstarR(MOMX_VAR) = dStarR*uStar
  UstarR(ENER_VAR) = uR(ENER_VAR)*numerR/denomR+&
       (pStar*uStar - vR(PRES_VAR)*vR(VELX_VAR))/denomR




  ! numerical flux
  if (sL >= 0.) then
     Flux(DENS_VAR:ENER_VAR) = FL(DENS_VAR:ENER_VAR)

  elseif ( (sL < 0.) .and. (uStar >= 0.) ) then
     Flux(DENS_VAR:ENER_VAR) = FL(DENS_VAR:ENER_VAR) &
          + sL*(UstarL(DENS_VAR:ENER_VAR) - uL(DENS_VAR:ENER_VAR))
  elseif ( (uStar < 0.) .and. (sR >= 0.) ) then
     Flux(DENS_VAR:ENER_VAR) = FR(DENS_VAR:ENER_VAR) &
          + sR*(UstarR(DENS_VAR:ENER_VAR) - uR(DENS_VAR:ENER_VAR))
  else
     Flux(DENS_VAR:ENER_VAR) = FR(DENS_VAR:ENER_VAR)
  endif

  return
end subroutine hllc
