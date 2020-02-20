subroutine soln_ReconEvolveAvg(dt)

#include "definition.h"  

  use grid_data
  use sim_data
  use primconsflux
  use bc

  implicit none
  real, intent(IN) :: dt
  real, dimension(NUMB_VAR,gr_imax) :: Vj
  real, dimension(NSYS_VAR,gr_imax) :: Flux
  integer :: i,j

  ! conservative left and right states
  real, dimension(NSYS_VAR,gr_imax) :: uL, uR


!!$  !get cell-centered fluxes
!!$  do i = gr_i0, gr_imax
!!$     call prim2flux(gr_V(DENS_VAR:NUMB_VAR, i), Flux(DENS_VAR:ENER_VAR))
!!$  end do
  if (sim_Torder == 1) then
     !this is just forward euler
     !should also be called w/ sim_RK set to true
     call soln_reconstruct(dt, gr_V)
     call soln_getFlux(gr_V)
  elseif (sim_Torder == 2) then
     !second order in time
     if(sim_RK) then
        !do RK2
        Flux(DENS_VAR:ENER_VAR,gr_i0:gr_imax) = 0.
        !initialize Vj as V0
        Vj(DENS_VAR:NUMB_VAR,gr_i0:gr_imax) = gr_V(DENS_VAR:NUMB_VAR,gr_i0:gr_imax)
        !do RK2 steps
        do j = 1,2
           call soln_RK2(dt, j, Vj, Flux)
           call bc_apply(Vj)
        end do
        gr_flux(DENS_VAR:ENER_VAR,gr_i0:gr_imax) = Flux(DENS_VAR:ENER_VAR,gr_i0:gr_imax)
     else
        !this is just midpoint rule
        call soln_reconstruct(dt, gr_V)
        call soln_getFlux(gr_V)
     end if
     
  elseif (sim_Torder == 3) then
     !RK3 in time
     Flux(DENS_VAR:ENER_VAR,gr_i0:gr_imax) = 0.
     !initialize Vj as V0
     Vj(DENS_VAR:NUMB_VAR,gr_i0:gr_imax) = gr_V(DENS_VAR:NUMB_VAR,gr_i0:gr_imax)
     !do RK3 steps
     do j = 1,3
        call soln_RK3(dt, j, Vj, Flux)
        call bc_apply(Vj)
     end do
     gr_flux(DENS_VAR:ENER_VAR,gr_i0:gr_imax) = Flux(DENS_VAR:ENER_VAR,gr_i0:gr_imax)
     
  elseif (sim_Torder == 4) then
     !RK4 in time
     Flux(DENS_VAR:ENER_VAR,gr_i0:gr_imax) = 0.
     !initialize Vj as V0
     Vj(DENS_VAR:NUMB_VAR,gr_i0:gr_imax) = gr_V(DENS_VAR:NUMB_VAR,gr_i0:gr_imax)
     !do rk4 steps
     do j = 1, 4
        !each step adds the flux at the kj'th state with the appropriate weight
        !to the total flux
        call soln_RK4(dt, j, Vj, Flux)
        call bc_apply(Vj)
     end do
     gr_flux(DENS_VAR:ENER_VAR,gr_i0:gr_imax) = Flux(DENS_VAR:ENER_VAR,gr_i0:gr_imax)
     
  end if
     

  return
end subroutine soln_ReconEvolveAvg
