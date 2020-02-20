subroutine soln_GP2(dt, V)
#include "definition.h"

  use grid_data, only: gr_radius, gr_GP_w2, gr_GPzk2
  use sim_data, only: sim_charLimiting, sim_mval, sim_RK
  use eigensystem
  use char_tracing
  use WENO

  implicit none
  real, intent(IN) :: dt
  real, dimension(NUMB_VAR,gr_imax), intent(IN) :: V

  integer :: i, var, s, N, k, R, M
  logical :: conservative

  !real :: vim2, vim1, vi, vip1, vip2, f0, delV, delL, delR, GL, GR, u_starL, u_starR, logP, sum_wbar, pow
  real :: sum_wbar
  real, dimension(NSYS_VAR) :: vecL, vecR, vL, vR, C0, C1, C2

  real, dimension(NUMB_WAVE) :: lambda
  real, dimension(NSYS_VAR,NUMB_WAVE) :: reig0, leig0
  
  real, dimension(2*gr_radius+1) :: G, un
  real, dimension(gr_radius+1) :: Gk, unk, vLk, vRk, beta_k
  real, dimension(2,gr_radius+1) :: weights, wbar
  
  !real, dimension(gr_radius+1,gr_radius+1,2) :: z_weights
  !real, dimension(gr_radius+1,2) :: lin_weights

!!$  real, dimension(5) :: G, un
!!$  real, dimension(3) :: Gk, unk, vLk, vRk, beta_k
!!$  real, dimension(2,3) :: weights, wbar
  !profiling
  real :: start, finish

  vlk = 0.
  vrk = 0.
  wbar = 0.
  weights = 0.
  
  conservative = .false.
  R = gr_radius
  M = R+1
  N = 2*R+1
  G = 0.
  !beta_k = 1.
  !wbar = 1.
  !allocate(weights(2,M))
  !allocate(wbar(2,M))
!!$  z_weights = gr_GPzk2
!!$  lin_weights = gr_gp_w2
  do i = gr_ibeg-1,gr_iend+1
     
     !copy cell-centered values to left and right states
     gr_vL(DENS_VAR:NUMB_VAR,i) = V(DENS_VAR:NUMB_VAR,i)
     gr_vR(DENS_VAR:NUMB_VAR,i) = V(DENS_VAR:NUMB_VAR,i)

     !get eigen information
     call eigenvalues(V(DENS_VAR:GAME_VAR,i), lambda)
     call left_eigenvectors(V(DENS_VAR:GAME_VAR,i), conservative, leig0)
     call right_eigenvectors(V(DENS_VAR:GAME_VAR,i), conservative, reig0)
     !this is the global stencil
     do var = DENS_VAR, PRES_VAR
        do s = 1, N
           if (sim_charLimiting) then
              G(s) = dot_product(V(DENS_VAR:PRES_VAR, i+s-M), leig0(DENS_VAR:PRES_VAR, var))
           else
              G(s) = V(var, i+(s-M))
           end if
        end do
        

        call gp_betas(G, R, beta_k)
        !call betas(G, R, beta_k)
!!$        do k = 1, R+1
!!$           beta_k(k) = dot_product(G(k:k+R), matmul(gr_GP_Kki2, G(k:K+R)))
!!$        end do

!!$        !original WENO5 weights
        do k = 1, M
           wbar(1, k) = gr_GP_w2(1,k)/(1.e-36+beta_k(k))**sim_mval
           wbar(2, k) = gr_GP_w2(2,k)/(1.e-36+beta_k(k))**sim_mval
        end do
!!$       
        sum_wbar = SUM(wbar(1,:))
        weights(1,:) = wbar(1,:)/sum_wbar
        sum_wbar = SUM(wbar(2,:))
        weights(2,:) = wbar(2,:)/sum_wbar
        
        do k = 1, M
           vLk(k) = dot_product(gr_GPzk2(1,1:M, k), G(k:k+R)) !- f0*dot_product(gr_GPzk2(1, :, k), unk)
           vRk(k) = dot_product(gr_GPzk2(2,1:M, k), G(k:k+R)) !- f0*dot_product(gr_GPzk2(2, :, k), unk)
        end do
!!$
        vecL(var) =  dot_product(weights(1,:), vLk(:))
        vecR(var) =  dot_product(weights(2,:), vRk(:))

!!$        
!!$        vecL(var) =  dot_product(gr_GP_w2(1:M,1), vLk(1:M))
!!$        vecR(var) =  dot_product(gr_GP_w2(1:M,2), vRk(1:M))
        
     end do !recon. vars
      
     do var = DENS_VAR, PRES_VAR
        if (sim_charLimiting) then
           !project char vars back onto prim vars
           vL(var) = dot_product(reig0(var,1:NUMB_WAVE),vecL(1:NUMB_WAVE))
           vR(var) = dot_product(reig0(var,1:NUMB_WAVE),vecR(1:NUMB_WAVE))
        else
           vL(var) = vecL(var)
           vR(var) = vecR(var)
        end if
     end do !var

     

     if (.not. sim_RK) then
        !do char tracing
!!!!!!!!!!!!!!! THis is PPM tracing !!!!!!!!!!!!!!!
        C2(DENS_VAR:PRES_VAR) = 6.*( .5*(vR(DENS_VAR:PRES_VAR) + vL(DENS_VAR:PRES_VAR)) &
             - V(DENS_VAR:PRES_VAR,i))
        C1(DENS_VAR:PRES_VAR) = vR(DENS_VAR:PRES_VAR) - vL(DENS_VAR:PRES_VAR)
        C0(DENS_VAR:PRES_VAR) = V(DENS_VAR:PRES_VAR,i) - C2(DENS_VAR:PRES_VAR)/12.
        call soln_PPMtracing(dt, V(DENS_VAR:GAME_VAR,i), C0(DENS_VAR:PRES_VAR), &
             C1(DENS_VAR:PRES_VAR), C2(DENS_VAR:PRES_VAR), vL, vR, lambda, leig0, reig0)
     end if


     gr_vL(DENS_VAR:PRES_VAR, i) = vL(DENS_VAR:PRES_VAR)
     gr_vR(DENS_VAR:PRES_VAR, i) = vR(DENS_VAR:PRES_VAR)
     
  end do !i

end subroutine soln_GP2
 !WENOZ weights
!!$        do k = 1, R+1
!!$           wbar(1, k) = gr_GP_w2(1, k) * (1. + ( ABS( beta_k(1) - beta_k(3) )/(1.e-36+beta_k(k)) )**sim_mval)
!!$           wbar(2, k) = gr_GP_w2(1,k) * (1. + (  ABS( beta_k(1) - beta_k(3) )/(1.e-36+beta_k(k)) )**sim_mval)
!!$        end do
!!$        !WENO-M weights
!!$        do k = 1, R+1
!!$           do s = 1, 2
!!$              wbar(s, k) = weights(s,k) * ( gr_GP_w2(s,k) + gr_GP_w2(s,k)**2 - &
!!$                   3.*gr_GP_w2(s,k)*weights(s,k) + weights(s,k)**2)/( gr_GP_w2(s,k)**2 + &
!!$                   weights(s,k)*(1.-2.*gr_GP_w2(s,k)) )
!!$           end do
!!$        end do
!!$
!!$        sum_wbar = SUM(wbar(1,:))
!!$        weights(1,:) = wbar(1,:)/sum_wbar
!!$        sum_wbar = SUM(wbar(2,:))
!!$        weights(2,:) = wbar(2,:)/sum_wbar
