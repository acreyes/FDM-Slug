subroutine soln_update(dt)

#include "definition.h"
  
  use grid_data
  use primconsflux, only : cons2prim

  implicit none
  real, intent(IN) :: dt
  integer :: i
  real :: dtx

  dtx = dt/gr_dx

  !! update conservative vars
  do i = gr_ibeg, gr_iend
     gr_U(DENS_VAR:ENER_VAR,i) = gr_U(DENS_VAR:ENER_VAR,i) - &
          dtx*(gr_flux(DENS_VAR:ENER_VAR,i+1) - gr_flux(DENS_VAR:ENER_VAR,i))
  end do


  !! get updated primitive vars from the updated conservative vars
  do i = gr_ibeg, gr_iend
     ! Eos is automatically callled inside cons2prim
     call cons2prim(gr_U(DENS_VAR:ENER_VAR,i),gr_V(DENS_VAR:GAME_VAR,i))
     !gr_V(DENS_VAR:GAME_VAR,i) = gr_vR(:,i)
  end do
  

  return
end subroutine soln_update
