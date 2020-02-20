subroutine soln_PLM(dt, V)

#include "definition.h"  

  use grid_data
  use sim_data
  use slopeLimiter
  use eigensystem
  use primconsflux

  implicit none
  real, intent(IN) :: dt
  real, dimension(NUMB_VAR,gr_imax), intent(IN) :: V
  integer :: i

  real, dimension(NUMB_WAVE) :: lambda
  real, dimension(NSYS_VAR,NUMB_WAVE) :: reig, leig
  logical :: conservative
  real, dimension(NSYS_VAR) :: Fs, Us, Uss, Fss, Ui, Uipo
  integer :: kWaveNum
  real :: delL, delR, delV
  real, dimension(NUMB_VAR)  :: Viph
  real, dimension(3, 2, NSYS_VAR) :: stencil
  real, dimension(2, NSYS_VAR) :: Fpm !the reconstructed split fluxes
  real, dimension(2) :: pm
  integer :: nVar, s, split
  

  pm = (/ 1., -1./)
  ! we need conservative eigenvectors
  conservative = .true.

  do i = gr_ibeg-1, gr_iend

     call eigenvalues(V(DENS_VAR:GAME_VAR,i),lambda)
     call prim2cons(V(:,i), Ui)
     call prim2cons(V(:,i+1), Uipo)
     call cons2prim(0.5*(Ui + Uipo), Viph)
     !Viph = 0.5*(V(:,i) + V(:,i+1))
     call left_eigenvectors (Viph(DENS_VAR:GAME_VAR),conservative,leig)
     call right_eigenvectors(Viph(DENS_VAR:GAME_VAR),conservative,reig)

     do s = 1,3
        !get pointwise flux & cons vars at stencil point
        ![s] stencil
        call prim2cons(V(:,i+s-2), Us)
        call prim2flux(V(:,i+s-2), Fs)
        ![s'] stencil
        call prim2cons(V(:,i+s-1), Uss)
        call prim2flux(V(:,i+s-1), Fss)
        !print *, i +s - 1, i + s - 2, i
        !do flux splitting
        do nvar = 1, NUMB_WAVE
           stencil(s, 1, nvar) = 0.5*dot_product(leig(:,nvar), Fs(:)  + gr_maxAlphas(nvar)*Us(:) )
           stencil(s, 2, nvar) = 0.5*dot_product(leig(:,nvar), Fss(:) - gr_maxAlphas(nvar)*Uss(:))
        end do
     end do

     
     do s = 1, 2 !loop over split fluxes
        ! here s is the left or right split flux
        do kWaveNum = 1, NUMB_WAVE
           delL = stencil(2, s, kWaveNum) - stencil(1, s, kWaveNum) !   i   - (i-1)
           delR = stencil(3, s, kWaveNum) - stencil(2, s, kWaveNum) ! (i+1) -   i
           !do slope limiting
           if (sim_limiter == 'minmod')then
              call minmod(delL, delR, delV)
           elseif (sim_limiter == 'vanLeer') then
              call vanLeer(delL, delR, delV)
           elseif (sim_limiter == 'mc') then
              call mc(delL, delR, delV)
           endif
           !finish interpolating inteface split-flux
           Fpm(s, kWaveNum) = stencil(2, s, kWaveNum) + pm(s)*0.5*delV
        end do
     end do

     !Put back together the interface flux from the right-eigenvector
     do nVar = DENS_VAR, PRES_VAR
        gr_flux(nVar, i+1) = dot_product( Fpm(1,:) + Fpm(2,:), reig(nVar,:))
        !print *, gr_flux(i+1,nVar)
     end do
     !print *,
  end do
  return
end subroutine soln_PLM
