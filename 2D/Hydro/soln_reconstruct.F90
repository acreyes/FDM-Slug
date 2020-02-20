subroutine soln_reconstruct(dt, V)

#include "definition.h"  

  use gp_data, only: gp_radius, gp_Npts, gp_stencil
  use grid_data
  use sim_data
  use eigensystem
  use primconsflux
  use sim_interfaces
  use reconstruction


  implicit none
  real, intent(IN) :: dt
  real, dimension(NUMB_VAR,gr_imax(XDIM),gr_imax(YDIM)), intent(IN) :: V

  real, allocatable, dimension(:,:) :: recon_stencil
  real, dimension(NUMB_WAVE)          :: lambda
  real, dimension(NSYS_VAR,NUMB_WAVE) :: reig, leig
  integer :: i, j, si, sj, scntr, var, Nx, im, ip, dir, i0, j0, is, js, ig, jg
  integer, dimension(NDIM) :: ibeg, iend

  !reconstruction pointers
  procedure (recon1D)          :: soln_FOG , soln_WENO, soln_gpWENO, soln_gpAO
  procedure (recon1D), pointer :: recon
  procedure (reconMD)          :: soln_gpMD
  procedure (reconMD), pointer :: recon2D

  if (sim_reconMultiD) then
     select case(sim_order)
     case(10)
        Nx = gp_Npts
        allocate(recon_stencil(NDIM,gp_Npts))
        recon_stencil = gp_stencil
        recon2D => soln_gpMD
     case default
        Nx = gp_Npts
        allocate(recon_stencil(NDIM,gp_Npts))
        recon_stencil = gp_stencil
        recon2D => soln_gpMD
     end select
     call reconMD_ij(gr_ibeg(XDIM),gr_iend(XDIM),gr_ibeg(YDIM),gr_iend(YDIM), Nx, recon_stencil, V, recon2D)
     !x guard cells
     call reconMD_ij(gr_ibeg(XDIM)-3,gr_ibeg(XDIM)-1,gr_ibeg(YDIM),gr_iend(YDIM), Nx, recon_stencil, V, recon2D)
     call reconMD_ij(gr_iend(XDIM)+1,gr_iend(XDIM)+3,gr_ibeg(YDIM),gr_iend(YDIM), Nx, recon_stencil, V, recon2D)
     !y guard cells
     call reconMD_ij(gr_ibeg(XDIM),gr_iend(XDIM),gr_ibeg(YDIM)-3,gr_ibeg(YDIM)-1, Nx, recon_stencil, V, recon2D)
     call reconMD_ij(gr_ibeg(XDIM),gr_iend(XDIM),gr_iend(YDIM)+1,gr_iend(YDIM)+3, Nx, recon_stencil, V, recon2D)
  else

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
        if (sim_AO) then
           recon => soln_gpAO
        else
           recon => soln_gpWENO
        end if
     case default
        !default to WENO reconstruction
        Nx = 5
        im = -2
        ip =  2
        recon => soln_WENO
     end select

     do dir = XDIM, NDIM
        !reconstruction inside domain
        ibeg = gr_ibeg
        iend = gr_iend
        call recon1D_ij(ibeg, iend, im, ip, dir, V, recon)
        !reconstructions needed in GC regions
        !left most first
        ibeg(dir) = gr_ibeg(dir)-3
        iend(dir) = gr_ibeg(dir)-1
        call recon1D_ij(ibeg, iend, im, ip, dir, V, recon)
        !right most next
        ibeg(dir) = gr_iend(dir) + 1
        iend(dir) = gr_iend(dir) + 3
        call recon1D_ij(ibeg, iend, im, ip, dir, V, recon)
     end do
!!$     dir = XDIM
!!$     call recon1D_ij(gr_ibeg(XDIM)-3,gr_ibeg(XDIM)-1,gr_ibeg(YDIM),gr_iend(YDIM), im, ip, dir, V, recon)
!!$     call recon1D_ij(gr_iend(XDIM)+1,gr_iend(XDIM)+3,gr_ibeg(YDIM),gr_iend(YDIM), im, ip, dir, V, recon)
!!$     dir = YDIM
!!$     call recon1D_ij(gr_ibeg(XDIM),gr_iend(XDIM),gr_ibeg(YDIM)-3,gr_ibeg(YDIM)-1, im, ip, dir, V, recon)
!!$     call recon1D_ij(gr_ibeg(XDIM),gr_iend(XDIM),gr_iend(YDIM)+1,gr_iend(YDIM)+3, im, ip, dir, V, recon)

  end if

  
  return
end subroutine soln_reconstruct
