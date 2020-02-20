subroutine old_soln_WENO(dt, V)

#include "definition.h"

  use grid_data
  use sim_data
  use eigensystem
  use primconsflux

  implicit none
  real, intent(IN) :: dt
  real, dimension(NUMB_VAR,gr_imax), intent(IN) :: V

  integer :: i, s, var, split
  real :: RHS, delV, wsum
  logical :: conservative
  real, dimension(NUMB_WAVE) :: lambda
  real, dimension(NSYS_VAR,NUMB_WAVE) :: reig, leig

  real, dimension(3) ::  beta
  real, dimension(3, 2) :: omegas, omegab, poly, gamm
  real, dimension(3  , 3,   2) :: coeff
  real, dimension(NSYS_VAR, 5, 2) :: flux
  real, dimension(5          ) :: Stencil
  real, dimension(NSYS_VAR   ) :: vL, vR, Us, Uss, Fs, Fss
  real, dimension(NSYS_VAR, 2) :: int_vals
  real, dimension(NUMB_VAR)   :: Viph


  !initialize w/ zeros
  vL   = 0.; vR   = 0.
  
  beta  = 0.

  conservative = .true.
  
  gamm(1, 1) = .1
  gamm(2, 1) = .6
  gamm(3, 1) = .3

  gamm(1, 2) = .3
  gamm(2, 2) = .6
  gamm(3, 2) = .1

  coeff(1, 1:3, 1) = (/ 2., -7., 11. /)
  coeff(2, 1:3, 1) = (/-1.,  5., 2.  /)
  coeff(3, 1:3, 1) = (/ 2.,  5.,-1.  /)

  coeff(1, 1:3, 2) = (/-1.,  5., 2. /)
  coeff(2, 1:3, 2) = (/ 2.,  5.,-1. /)
  coeff(3, 1:3, 2) = (/11., -7., 2. /)

  coeff = coeff/6.

  do i = gr_ibeg-1, gr_iend+1
     call eigenvalues(V(DENS_VAR:GAME_VAR,i),lambda)
     call prim2cons(V(:,i), Us)
     call prim2cons(V(:,i+1), Uss)
     call cons2prim(0.5*(Us + Uss), Viph)

     call left_eigenvectors (Viph(DENS_VAR:GAME_VAR),conservative,leig)
     call right_eigenvectors(Viph(DENS_VAR:GAME_VAR),conservative,reig)

     !do flux-splitting
     do s = 1,5
        ![s] stencil
        call prim2cons(V(:,i+s-3), Us)
        call prim2flux(V(:,i+s-3), Fs)
        ![s'] stencil
        call prim2cons(V(:,i+s-2), Uss)
        call prim2flux(V(:,i+s-2), Fss)

        do var = 1, NUMB_WAVE
           flux(var, s, 1) = 0.5*dot_product(leig(:,var), Fs(:)  + gr_maxAlphas(var)*Us(:) )
           flux(var, s, 2) = 0.5*dot_product(leig(:,var), Fss(:)  - gr_maxAlphas(var)*Uss(:) )
        end do
     end do

     do var = 1, NUMB_WAVE
        do split = 1,2
           ! calculate smoothness indicators (betas)
           beta(1) = 13./12.*(flux(var, 1, split) - 2.*flux(var, 2, split) + flux(var, 3, split))**2 &
                + 0.25*(flux(var, 1, split) - 4.*flux(var, 2, split) + 3.*flux(var, 3, split))**2
           
           beta(2) = 13./12.*(flux(var, 2, split) - 2.*flux(var, 3, split) + flux(var, 4, split))**2 &
                + 0.25*(flux(var, 2, split) - flux(var, 4, split))**2

           beta(3) = 13./12.*(flux(var, 3, split) - 2.*flux(var, 4, split) + flux(var, 5, split))**2 &
                + 0.25*(3.*flux(var, 3, split) - 4.*flux(var, 4, split) + flux(var, 5, split))**2

           
           !now lets calculate weight factors
           do s = 1,3
              if (sim_WENO == '5') then
                 omegab(s, split) = gamm(s, split)/((sim_WENeps + beta(s))**sim_mval)
              elseif (sim_WENO == 'Z') then
                 omegab(s, split) = gamm(s, split)*(1 + (ABS(beta(1) - beta(3))/(sim_WENeps + beta(s)))**sim_mval)
              end if
           end do
           wsum = SUM(omegab(:, split))
           !normalize omegas
           omegas(:,split) = omegab(:,split)/wsum

           !calculate the interface values
           poly(1, split) = dot_product(flux(var, 1:3, split), coeff(1, 1:3, split))
           poly(2, split) = dot_product(flux(var, 2:4, split), coeff(2, 1:3, split))
           poly(3, split) = dot_product(flux(var, 3:5, split), coeff(3, 1:3, split))
        end do !split = 1,2
        vR(var) = dot_product(omegas(:, 1), poly(:, 1))
        vL(var) = dot_product(omegas(:, 2), poly(:, 2))
     end do !split flux vars
     !now all that is left is to recombine the split fluxes into the interface flux
     do var = DENS_VAR,PRES_VAR
        !print *, abs(vR(var)-flux(var,3,1))
        gr_flux(var,i+1) = dot_product( vL(:) + vR(:), reig(var,:))
     end do
  end do
  return
end subroutine old_soln_WENO
