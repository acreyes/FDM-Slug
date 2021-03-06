subroutine soln_PPM(dt, V)

#include "definition.h"  

  use grid_data
  use sim_data
  use eigensystem
  use primconsflux

  implicit none
  real, intent(IN) :: dt
  real, dimension(NUMB_VAR,gr_imax), intent(IN) :: V

  integer :: i, var, s
  real, dimension(NUMB_WAVE) :: lambda
  real, dimension(NSYS_VAR,NUMB_WAVE) :: reig, leig
  logical :: conservative
  real, dimension(NUMB_VAR) :: Viph
  real, dimension(NSYS_VAR) :: vL, vR, Fs, Fss, Us, Uss
  real, dimension(3, 2, NSYS_VAR) :: delV
  real, dimension(2, NSYS_VAR) :: C0, C1, C2
  real, dimension(5, 2, NSYS_VAR) :: stencil !stencil for split fluxes
  

  ! we need primitive eigenvectors
  conservative = .true.
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
           stencil(s, 1, var) = 0.5*dot_product(leig(:,var), Fs(:)  + gr_maxAlphas(var)*Us(:) )
           stencil(s, 2, var) = 0.5*dot_product(leig(:,var), Fss(:) - gr_maxAlphas(var)*Uss(:))
        end do
     end do
     !do interpolation on split fluxes to i+1/2 interface
     do s = 1,2
        !s now indexes the left or right split fluxes
        !get TVD slopes
        call del_V(stencil(:, s, :), delV(:, s, :))
        
        vL(DENS_VAR:PRES_VAR) = 0.5*(stencil(2, s, :) + stencil(3, s, :)) &
             + ( delV(1, s, :) - delV(2, s, :) )/6.
        vR(DENS_VAR:PRES_VAR) = 0.5*(stencil(3, s, :) + stencil(4, s, :)) &
             + ( delV(2, s, :) - delV(3, s, :) )/6.

!!$        !enforce monotonicity conditions
        do var = DENS_VAR, PRES_VAR
           if ( (vR(var) - stencil(3,s,var))*(stencil(3,s,var)-vL(var)) <= 0.) then
              vL(var) = stencil(3, s, var)
              vR(var) = stencil(3, s, var)
           end if
           if (6.*(vR(var)-vL(var))*(stencil(3,s,var)-.5*(vL(var)+vR(var))) > (vR(var)-vL(var))**2 ) then
              vL(var) = 3.*stencil(3, s, var) - 2.*vR(var)
           end if
           if ( 6.*(vR(var)-vL(var))*(stencil(3,s,var)-.5*(vL(var)+vR(var))) < -(vR(var)-vL(var))**2 ) then
              vR(var) = 3.*stencil(3, s, var) - 2.*vL(var)
           end if
        end do
        
        C2(s, DENS_VAR:PRES_VAR) = 6.*( .5*(vR(DENS_VAR:PRES_VAR) + vL(DENS_VAR:PRES_VAR)) &
             - stencil(3, s, :))
        C1(s, DENS_VAR:PRES_VAR) = vR(DENS_VAR:PRES_VAR) - vL(DENS_VAR:PRES_VAR)
        C0(s, DENS_VAR:PRES_VAR) = stencil(3, s, :) - C2(s, DENS_VAR:PRES_VAR)/12.
  end do
     

  vL(DENS_VAR:PRES_VAR) = C0(2, DENS_VAR:PRES_VAR) - .5*C1(2, DENS_VAR:PRES_VAR) &
       + .25*C2(2, DENS_VAR:PRES_VAR)
  vR(DENS_VAR:PRES_VAR) = C0(1, DENS_VAR:PRES_VAR) + .5*C1(1, DENS_VAR:PRES_VAR) &
       + .25*C2(1, DENS_VAR:PRES_VAR)

  !now all that is left is to recombine the split fluxes into the interface flux
  do var = DENS_VAR,PRES_VAR
     gr_flux(var,i+1) = dot_product( vL(:) + vR(:), reig(var,:))
  end do
     
  end do
  return
end subroutine soln_PPM

subroutine del_V(stencil, delV)

  use grid_data
  use sim_data
  use slopeLimiter

  implicit none

  real, dimension(5, NSYS_VAR), intent(IN) :: stencil
  real, dimension(3, NSYS_VAR), intent(OUT) :: delV

  real, dimension(NUMB_VAR)  :: delL,delR
  integer :: nVar, s

  do s = 2,4
     delL(DENS_VAR:PRES_VAR) = stencil(s  , :) - stencil(s-1, :)
     delR(DENS_VAR:PRES_VAR) = stencil(s+1, :) - stencil(s  , :)
     ! slope limiting
     do nVar = DENS_VAR,PRES_VAR
        if (sim_limiter == 'minmod') then
           call minmod(delL(nVar),delR(nVar),delV(s, nVar))
        elseif (sim_limiter == 'vanLeer') then
           call vanLeer(delL(nVar),delR(nVar),delV(s, nVar))
        elseif (sim_limiter == 'mc') then
           call mc(delL(nVar),delR(nVar),delV(s, nVar))
        endif
     end do
  end do

  
  return
end subroutine del_V

