subroutine soln_RK2(dt, j, Vj, Flux)
  !performs the jth step of the RK3 algorithm
  !Total Flux = (k1 + k2 + 2k3)/6

#include "definition.h"

  use grid_data
  use primconsflux

  implicit none
  real, intent(IN) :: dt
  integer, intent(IN) :: j
  real, dimension(NUMB_VAR,gr_imax), intent(INOUT) :: Vj
  real, dimension(NSYS_VAR,gr_imax), intent(INOUT) :: Flux

  real :: dtx
  integer :: i
  real, dimension(NUMB_VAR,gr_imax) :: Uj

  call soln_reconstruct(dt, Vj)
  call soln_getFlux(Vj)

  dtx = dt/gr_dx

  if (j == 1) then
     do i = gr_ibeg, gr_iend
        !get the U1 state
        Uj(DENS_VAR:ENER_VAR,i) = gr_U(DENS_VAR:ENER_VAR,i) - &
             dtx *(gr_flux(DENS_VAR:ENER_VAR,i+1) - gr_flux(DENS_VAR:ENER_VAR,i))
        !print *, gr_flux(ENER_VAR,i+1) - gr_flux(ENER_VAR,i), i
        call cons2prim(Uj(DENS_VAR:ENER_VAR,i),Vj(DENS_VAR:GAME_VAR,i))
     end do
  end if

  Flux(DENS_VAR:ENER_VAR,gr_i0:gr_imax) = Flux(DENS_VAR:ENER_VAR,gr_i0:gr_imax) + &
       .5*gr_flux(DENS_VAR:ENER_VAR,gr_i0:gr_imax)
end subroutine soln_RK2
