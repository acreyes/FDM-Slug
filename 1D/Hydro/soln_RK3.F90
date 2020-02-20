subroutine soln_RK3(dt, j, Vj, Flux)
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

  real :: F, dtx
  integer :: i
  real, dimension(NUMB_VAR,gr_imax) :: Uj

  !F is the factor that multiplies the Kj flux
  if (j == 3) then
     F = 2./3.
  else
     F = 1./6.
  end if

  call soln_reconstruct(dt, Vj)
  call soln_getFlux(Vj)

  dtx = dt/gr_dx

  if (j == 1) then
     do i = gr_ibeg, gr_iend
        !getting the U1 state
        !U1 = U0 + k1
        Uj(DENS_VAR:ENER_VAR,i) = gr_U(DENS_VAR:ENER_VAR,i) - &
             dtx *(gr_flux(DENS_VAR:ENER_VAR,i+1) - gr_flux(DENS_VAR:ENER_VAR,i))
     end do
  elseif (j == 2) then
     do i = gr_ibeg, gr_iend
        !getting the U2 state
        !U2 = U0 + 1/4(k1 +k2)
        !k1 = 6.*Flux
        Uj(DENS_VAR:ENER_VAR,i) = gr_U(DENS_VAR:ENER_VAR,i) - &
             0.25*dtx*( 6.*(Flux(DENS_VAR:ENER_VAR,i+1) -Flux(DENS_VAR:ENER_VAR,i) ) +&
             (gr_flux(DENS_VAR:ENER_VAR,i+1) - gr_flux(DENS_VAR:ENER_VAR,i))  )
     end do
  end if

  !get updated prim vars if needed for j+1 step
  if (j .NE. 3) then
     do i = gr_ibeg, gr_iend
        ! Eos is automatically callled inside cons2prim
        call cons2prim(Uj(DENS_VAR:ENER_VAR,i),Vj(DENS_VAR:GAME_VAR,i))
     end do
  end if

  !finally add kj flux to total flux with weights
  Flux(DENS_VAR:ENER_VAR,gr_i0:gr_imax) = Flux(DENS_VAR:ENER_VAR,gr_i0:gr_imax) + &
       F*gr_flux(DENS_VAR:ENER_VAR,gr_i0:gr_imax)

end subroutine soln_RK3
