subroutine soln_getFlux(V)

#include "definition.h"  

  use grid_data
  use sim_data
  use primconsflux
  use sim_interfaces

  implicit none
  real, dimension(NUMB_VAR,gr_imax), intent(IN) :: V

  
  procedure (rmn_slvr) :: hll, hllc, roe
  procedure (rmn_slvr), pointer :: RP
  integer :: i, var
  real    :: D2, D4, Fiph_fac, fac2, fac4
  real, dimension(2) :: coeff2
  real, dimension(3) :: coeff3
  real, dimension(4) :: coeff4
  real, dimension(5) :: coeff5
  real, dimension(NSYS_VAR,gr_imax) :: cntrFlux, intFlux

  fac2 = -1./24.
  fac4 = 3./640.  !this is in my notes and also in Chen et. al. and delZanna et al
!  fac4 = 7./5760. !this is from Jiang et al
  coeff2 = 4.*(/1., 1. /)
  coeff4 = (/8., -72., -72., 8. /)/3.
  Fiph_fac = 1. + fac2*4.*(-2.) + fac4*128./3.

  coeff3 = (/1., -2., 1. /)
  coeff5 = (/1., -4., 6., -4., 1. /)
  
  if (sim_riemann == 'hll') then
     RP => hll
  elseif (sim_riemann == 'hllc') then
     RP => hllc
  elseif (sim_riemann == 'roe') then
     RP => roe
  else
     RP => hll
  end if

  if (sim_intFlux) then
     !here we correct the fluxes using the interface fluxes
     do i = gr_ibeg-2, gr_iend+3
        call RP(gr_vR(  DENS_VAR:GAME_VAR,i-1),&
             gr_vL(  DENS_VAR:GAME_VAR,i  ),&
             intFlux(DENS_VAR:ENER_VAR,i)   )
     end do

     do i = gr_ibeg, gr_iend+1
        do var = DENS_VAR, ENER_VAR
           if (.false.) then
              call soln_gpFlux(intFlux(var,i-2:i+2),5,gr_flux(var,i))
           else
              D2 = dot_product(coeff3, intFlux(var,i-1:i+1))
              D4 = dot_product(coeff5, intFlux(var,i-2:i+2))
              gr_flux(var,i) = intFlux(var,i) +fac2*D2 + fac4*D4
           end if
        end do
        !     print *, gr_flux(:,i)
     end do

    
  else
     !correct fluxes using cell center fluxes
     
     !calculate interface fluxes
     do i = gr_ibeg, gr_iend+1
        call RP(gr_vR(  DENS_VAR:GAME_VAR,i-1),&
             gr_vL(  DENS_VAR:GAME_VAR,i  ),&
             intFlux(DENS_VAR:ENER_VAR,i)   )
     end do
     
     !calculate cell centered fluxes
     do i = gr_i0, gr_imax
        call prim2flux(V(:,i),cntrFlux(:,i))
     end do

     !now we need to correct interface fluxes to approximate numerical flux
     do i = gr_ibeg, gr_iend+1
        do var = DENS_VAR, ENER_VAR
           D2 = dot_product(coeff2, cntrFlux(var,i-1:i  ))
           D4 = dot_product(coeff4, cntrFlux(var,i-2:i+1))
           gr_flux(var,i) = Fiph_fac*intFlux(var,i) + fac2*D2 + fac4*D4
        end do
        !     print *, gr_flux(:,i)
     end do

  end if

  


  return
end subroutine soln_getFlux
