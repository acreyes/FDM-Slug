subroutine soln_reconstruct(dt, V)

#include "definition.h"  

  use gp_data, only: gp_radius
  use grid_data
  use sim_data
  use eigensystem
  use primconsflux
  use sim_interfaces
  use char_tracing

  implicit none
  real, intent(IN) :: dt
  real, dimension(NUMB_VAR,gr_imax), intent(IN) :: V

  real, allocatable, dimension(:,:) :: stencil
  real, dimension(NUMB_WAVE)          :: lambda
  real, dimension(NSYS_VAR,NUMB_WAVE) :: reig, leig
  integer :: i, s, var, Nx, im, ip

  !reconstruction pointers
  procedure (recon1D)          :: soln_WENO, soln_GP, soln_FOG
  procedure (recon1D), pointer :: recon

  !figure out reconstruction to use
  select case(sim_order)
  case(1)
     !FOG
     Nx = 1
     im = 0
     ip = 0
     recon => soln_FOG
  case(5)
     !WENO
     Nx = 5
     im = -2
     ip =  2
     recon => soln_WENO
  case(10)
     Nx = 2*gp_radius +1
     im = -gp_radius
     ip =  gp_radius
     recon => soln_GP
  case default
     !default to WENO reconstruction
     Nx = 5
     im = -2
     ip =  2
     recon => soln_WENO
  end select

  allocate(stencil(Nx,NUMB_VAR))

  do i = gr_ibeg-3, gr_iend+3
     !make reconstruction stencil
     
     do var = 1, NUMB_VAR
        do s = im,ip,1
           stencil(s+ip+1, var) = V(var,i+s)
        end do
     end do
     !now reconstruct left and right states
     call recon(dt, Nx, stencil, gr_vL(:,i), gr_vR(:,i))

     !char tracing
     if (.not. sim_RK) then
        call eigenvalues(V(:,i),lambda)
        call left_eigenvectors (V(:,i), .false., leig)
        call right_eigenvectors(V(:,i), .false., reig)
        call soln_PPMtracing(dt, V(:,i), gr_vL(:,i), gr_vR(:,i), &
                             lambda, leig, reig)
     end if
  end do
deallocate(stencil)
!!$  !this is only needed for the regular FDM
!!$  do i = gr_ibeg, gr_iend
!!$     !need to get the max char wavespeeds in the whole domain
!!$     call eigenvalues(V(:,i), lambda)
!!$     do nvar = 1,NUMB_WAVE
!!$        gr_maxAlphas(nvar) = MAX(gr_maxAlphas(nvar), ABS(lambda(nvar)))
!!$     end do
!!$  end do
!!$
!!$  do i = gr_i0, gr_imax
!!$     call prim2flux(V(:,i),Flux(:,i))
!!$  end do
!!$
!!$  do i = gr_ibeg-1, gr_iend
!!$     Viph = 0.5*(V(:,i) + V(:,i+1))
!!$     call left_eigenvectors (Viph(DENS_VAR:GAME_VAR),.true.,leig)
!!$     call right_eigenvectors(Viph(DENS_VAR:GAME_VAR),.true.,reig)
!!$
!!$     !make stencil
!!$     
!!$     do s = 1, 5
!!$        !        print *, '+', i, i+s-3
!!$        call prim2cons(V(:,i+s-3), U)
!!$        do var = 1, NUMB_WAVE
!!$           stencil(s,var) = 0.5*dot_product(leig(:,var), Flux(:,i+s-3) + gr_maxAlphas(var)*U(:))
!!$           
!!$        end do
!!$     end do
!!$     !print *, '+'
!!$     call recon(dt, Nx, stencil, Fp)
!!$
!!$     do s = 1, 5
!!$        !        print *, '-', i, i+4-s
!!$        call prim2cons(V(:,i+4-s), U)
!!$        do var = 1, NUMB_WAVE
!!$           
!!$           stencil(s,var) = 0.5*(dot_product( leig(:,var), Flux(:,i+1+3-s) - gr_maxAlphas(var)*U(:) ))
!!$        end do
!!$     end do
!!$     !print *, '-'
!!$     call recon(dt, Nx, stencil, Fm)
!!$
!!$     !put back together the fluxes
!!$     do var = DENS_VAR, ENER_VAR
!!$        gr_flux(var,i+1) = dot_product(Fp + Fm, reig(var,:))
!!$     end do
!!$  end do
  
  return
end subroutine soln_reconstruct
