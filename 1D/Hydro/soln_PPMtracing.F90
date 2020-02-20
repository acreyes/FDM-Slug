module char_tracing

#include "definition.h"

  use grid_data
  use sim_data
  use eigensystem

contains

subroutine soln_PPMtracing(dt, V, vL, vR, lambda, leig, reig)
  !subroutine to do characteristic tracing for a particular cell for PPM and put the edge states into vL and vR

  implicit none
  real, intent(IN) :: dt
  real, dimension(NSYS_VAR), intent(IN   ) :: lambda
  real, dimension(NSYS_VAR), intent(INOUT) :: vL, vR
  real, dimension(NSYS_VAR,NUMB_WAVE), intent(IN) :: leig, reig
  real, dimension(NUMB_VAR), intent(IN   ) :: V

  logical :: conservative
  integer :: kWaveNum
  real :: lambdaDtDx, delC1, delC2
  real, dimension(NSYS_VAR) :: sigL1, sigL2, sigR1, sigR2, vecL1, vecR1, vecL2, vecR2, C0, C1, C2

  !need primitive eigenvectors
!!$  conservative = .FALSE.
!!$  call eigenvalues(V(DENS_VAR:GAME_VAR), lambda)
!!$  call left_eigenvectors(V(DENS_VAR:GAME_VAR), conservative, leig)
!!$  call right_eigenvectors(V(DENS_VAR:GAME_VAR), conservative, reig)

  C0 = V(DENS_VAR:NSYS_VAR)
  C1 = vR - vL
  C2 = 2.*(vR+vL-C0)

  !set initial sum to 0
  sigL1(DENS_VAR:ENER_VAR) = 0.
  sigL2(DENS_VAR:ENER_VAR) = 0.
  sigR1(DENS_VAR:ENER_VAR) = 0.
  sigR2(DENS_VAR:ENER_VAR) = 0.

  do kWaveNum = 1, NUMB_WAVE
     lambdaDtDx = lambda(kWaveNum)*dt/gr_dx
     delC1 = gr_dx*     dot_product(leig(DENS_VAR:PRES_VAR, kWaveNum), C1(DENS_VAR:PRES_VAR))
     delC2 = (gr_dx**2)*dot_product(leig(DENS_VAR:PRES_VAR, kWaveNum), C2(DENS_VAR:PRES_VAR))

     if (sim_riemann == 'roe') then
        if (lambda(kWaveNum) > 0.) then
           vecR1(DENS_VAR:PRES_VAR) = (1. - lambdaDtDx)*delC1*reig(DENS_VAR:PRES_VAR, kWaveNum)
           sigR1(DENS_VAR:PRES_VAR) = sigR1(DENS_VAR:PRES_VAR) + 0.5*vecR1(DENS_VAR:PRES_VAR)

           vecR2(DENS_VAR:PRES_VAR) = (1. - 2.*lambdaDtDx + 4.*(lambdaDtDx**2)/3.)*&
                delC2*reig(DENS_VAR:PRES_VAR, kWaveNum)
           sigR2(DENS_VAR:PRES_VAR) = sigR2(DENS_VAR:PRES_VAR) + 0.25*vecR2(DENS_VAR:PRES_VAR)
        elseif (lambda(kWaveNum) < 0.) then
           vecL1(DENS_VAR:PRES_VAR) = (-1. - lambdaDtDx)*delC1*reig(DENS_VAR:PRES_VAR, kWaveNum)
           sigL1(DENS_VAR:PRES_VAR) = sigL1(DENS_VAR:PRES_VAR) + 0.5*vecL1(DENS_VAR:PRES_VAR)

           vecL2(DENS_VAR:PRES_VAR) = (1. + 2.*lambdaDtDx + 4.*(lambdaDtDx**2)/3.)*&
                delC2*reig(DENS_VAR:PRES_VAR, kWaveNum)
           sigL2(DENS_VAR:PRES_VAR) = sigL2(DENS_VAR:PRES_VAR) + .25*vecL2(DENS_VAR:PRES_VAR)
        endif
     elseif (sim_riemann == 'hll') then
        vecR1(DENS_VAR:PRES_VAR) = 0.5*(1. - lambdaDtDx)*delC1*reig(DENS_VAR:PRES_VAR, kWaveNum)
        sigR1(DENS_VAR:PRES_VAR) = sigR1(DENS_VAR:PRES_VAR) + vecR1(DENS_VAR:PRES_VAR)

        vecR2(DENS_VAR:PRES_VAR) = 0.25*(1. - 2.*lambdaDtDx + 4.*(lambdaDtDx**2)/3.)*&
             delC2*reig(DENS_VAR:PRES_VAR, kWaveNum)
        sigR2(DENS_VAR:PRES_VAR) = sigR2(DENS_VAR:PRES_VAR) + vecR2(DENS_VAR:PRES_VAR)

        vecL1(DENS_VAR:PRES_VAR) = 0.5*(1. + lambdaDtDx)*delC1*reig(DENS_VAR:PRES_VAR, kWaveNum)
        sigL1(DENS_VAR:PRES_VAR) = sigL1(DENS_VAR:PRES_VAR) - vecL1(DENS_VAR:PRES_VAR)

        vecL2(DENS_VAR:PRES_VAR) = 0.25*(1. + 2.*lambdaDtDx + 4.*(lambdaDtDx**2)/3.)*&
             delC2*reig(DENS_VAR:PRES_VAR, kWaveNum)
        sigL2(DENS_VAR:PRES_VAR) = sigL2(DENS_VAR:PRES_VAR) + vecL2(DENS_VAR:PRES_VAR)
     end if
  enddo
  
  !Now do PPM reconstruction for left and right states
  vL(DENS_VAR:PRES_VAR) = C0(DENS_VAR:PRES_VAR) + sigL1(DENS_VAR:PRES_VAR) + &
       sigL2(DENS_VAR:PRES_VAR)
  vR(DENS_VAR:PRES_VAR) = C0(DENS_VAR:PRES_VAR) + sigR1(DENS_VAR:PRES_VAR) + &
       sigR2(DENS_VAR:PRES_VAR)


end subroutine soln_PPMtracing

end module char_tracing
