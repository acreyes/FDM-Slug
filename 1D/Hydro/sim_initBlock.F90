subroutine sim_initBlock()

#include "definition.h"
  
  use sim_data
  use grid_data, only : gr_V,gr_U,gr_i0,gr_imax,gr_xCoord,gr_dx,gr_ngc
  use primconsflux, only : prim2cons
  
  implicit none

  integer :: i
  real :: ekin, eint, x, gamm, del, dens
  
  ! generate x-coordinate
  do i = gr_i0,gr_imax
     gr_xCoord(i) = (real(i-gr_ngc)-0.5)*gr_dx
  end do

  do i = gr_i0,gr_imax
     if (sim_icType == 'shock') then
        if (gr_xCoord(i) < sim_shockLoc) then
           gr_V(DENS_VAR,i) = sim_densL
           gr_V(VELX_VAR,i) = sim_velxL
           gr_V(PRES_VAR,i) = sim_presL
        else
           gr_V(DENS_VAR,i) = sim_densR
           gr_V(VELX_VAR,i) = sim_velxR
           gr_V(PRES_VAR,i) = sim_presR
        end if

     elseif (sim_icType == 'blast') then
        !do blas IC
        if (gr_xCoord(i) <= 0.9 .AND. gr_xCoord(i) > .1) then
           !middle state
           gr_V(DENS_VAR,i) = 1.
           gr_V(VELX_VAR,i) = 0.
           gr_V(PRES_VAR,i) = .01
        elseif (gr_xCoord(i) <= .1) then
           gr_V(DENS_VAR,i) = 1.
           gr_V(VELX_VAR,i) = 0.
           gr_V(PRES_VAR,i) = 1000.
        else
           gr_V(DENS_VAR,i) = 1.
           gr_V(VELX_VAR,i) = 0.
           gr_V(PRES_VAR,i) = 100.
        endif
     elseif (sim_icType == 'isen_pulse') then
        x = gr_xCoord(i) - 0.5
        gamm = .5*(sim_gamma-1.)
        del = 0.1
        
        dens = 1. + 0.1*EXP(-x**2/(del)**2)
        gr_V(DENS_VAR,i) = dens
        gr_V(VELX_VAR,i) = SQRT(sim_gamma)/gamm * (dens**(gamm)-1.)
        gr_V(PRES_VAR,i) = 1.*dens**(sim_gamma)
           
     elseif (sim_icType == 'shu') then
        !do IC for shu-osher problem
        !transform the domain [0,9] onto the one give for the problem:[-4.5,4.5]
        x = gr_xCoord(i) - 4.5
        if (x < -4.) then
           !left state
           gr_V(DENS_VAR,i) = 3.857143
           gr_V(VELX_VAR,i) = 2.629369
           gr_V(PRES_VAR,i) = 10.33333
        else
           gr_V(DENS_VAR,i) = 1 + .2*SIN(5.*x)
           gr_V(VELX_VAR,i) = 0.
           gr_V(PRES_VAR,i) = 1.
        end if

     elseif (sim_icType == 'gauss') then
        x = gr_xCoord(i) - .5
        gr_V(DENS_VAR,i) = 1. + EXP(-100.*x**2)
        gr_V(VELX_VAR,i) = 1.
        gr_V(PRES_VAR,i) = 1./sim_gamma
     elseif (sim_icType == 'sin') then
        x = gr_xCoord(i)
        gr_V(DENS_VAR,i) = 1.5 + SIN(2.*PI*x)
        gr_V(VELX_VAR,i) = 1.
        gr_V(PRES_VAR,i) = 1./sim_gamma
        
     end if
     gr_V(GAMC_VAR,i) = sim_gamma
     gr_V(GAME_VAR,i) = sim_gamma
     gr_V(EINT_VAR,i) = gr_V(PRES_VAR,i)/(gr_V(GAME_VAR,i)-1.)/gr_V(DENS_VAR,i)
     
  end do

  ! also initialize conservative vars
  do i = gr_i0,gr_imax
     call prim2cons(gr_V(:,i), gr_U(DENS_VAR:ENER_VAR,i))
  end do

  
end subroutine sim_initBlock
