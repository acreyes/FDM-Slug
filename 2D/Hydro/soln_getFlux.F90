subroutine soln_getFlux(V)

#include "definition.h"  

  use grid_data
  use sim_data
  use primconsflux
  use sim_interfaces
  use gp_data, only: gp_radius

  implicit none
  real, dimension(NUMB_VAR,gr_imax(XDIM),gr_imax(YDIM)), intent(IN) :: V

  
  procedure (rmn_slvr) :: hll, hllc, roe
  procedure (rmn_slvr), pointer :: RP
  procedure (num_flux) :: soln_intFlux, soln_cntrFlux, soln_gpFlux
  procedure (num_flux), pointer :: NumFlux
  integer :: i, j, var, dir, l, m, BDRY_VEL
  integer :: ibeg, jbeg, iend, jend, im, jm, cntr, Npts
  real    :: D2, D4, Fiph_fac, fac2, fac4
  real, dimension(2) :: coeff2
  real, dimension(3) :: coeff3
  real, dimension(4) :: coeff4
  real, dimension(5) :: coeff5, flux_vec
  real, dimension(NSYS_VAR,gr_imax(XDIM),gr_imax(YDIM)) ::  intFlux
  real, dimension(NSYS_VAR,gr_imax(XDIM),gr_imax(YDIM)) ::  cntrFlux

  ! general left and right Riemann problem (rp) state vectors
  real, dimension(NUMB_VAR) :: rpL, rpR

  if (sim_intFlux) then
     if (sim_GPFlux .and. sim_order == 10) then
        NumFlux => soln_gpFlux
        Npts    =  5 !2*gp_radius + 1 
     else
        Npts = 5
        NumFlux => soln_intFlux
     end if
  else
     Npts = 5
     NumFlux => soln_cntrFlux
  end if
  
  if (sim_riemann == 'hll') then
     RP => hll
  elseif (sim_riemann == 'hllc') then
     RP => hllc
  elseif (sim_riemann == 'roe') then
     RP => roe
  else
     RP => hll
  end if

  if (sim_intFlux) then
     !here we correct the fluxes using the interface fluxes
     do dir = XDIM, NDIM
        select case(dir)
        case(XDIM)
           ibeg = gr_ibeg(XDIM)-2
           iend = gr_iend(XDIM)+3
           jbeg = gr_ibeg(YDIM)
           jend = gr_iend(YDIM)
           im = 1
           jm = 0
           BDRY_VEL = VELX_VAR
        case(YDIM)
           ibeg = gr_ibeg(XDIM)
           iend = gr_iend(XDIM)
           jbeg = gr_ibeg(YDIM)-2
           jend = gr_iend(YDIM)+3
           im = 0
           jm = 1
           BDRY_VEL = VELY_VAR
        end select
        do i = ibeg, iend
           do j = jbeg, jend

              !store left and right states
              rpL(DENS_VAR:GAME_VAR) = gr_vR(DENS_VAR:GAME_VAR,i-im, j-jm, dir)
              rpR(DENS_VAR:GAME_VAR) = gr_vL(DENS_VAR:GAME_VAR,i   , j   , dir)
              
#ifdef BDRY_VAR
              ! solid internal boundary
              ! Cell i and i-1:
              if (gr_V(BDRY_VAR,i,j) > 0.0 .and. gr_V(BDRY_VAR,i-im,j-jm) < 0.0) then
                 rpR = rpL
                 rpR(BDRY_VEL) = -rpL(BDRY_VEL)
!!$                 gr_vR(DENS_VAR:GAME_VAR,i,j,dir) = gr_vL(DENS_VAR:GAME_VAR,i,j,dir)
!!$                 gr_vR(VELX_VAR,i,j,dir) = -gr_vL(VELX_VAR,i,j,dir)
              end if

              if (gr_V(BDRY_VAR,i,j) < 0.0 .and. gr_V(BDRY_VAR,i-im,j-jm) > 0.0) then
                 rpL = rpR
                 rpL(BDRY_VEL) = -rpR(BDRY_VEL)
!!$                 gr_vL(DENS_VAR:GAME_VAR,i,j,dir) = gr_vR(DENS_VAR:GAME_VAR,i,j,dir)
!!$                 gr_vL(VELX_VAR,i,j,dir) = -gr_vR(VELX_VAR,i,j,dir)
              end if
#endif          
!!$              call RP(gr_vR  (DENS_VAR:GAME_VAR,i-im, j-jm, dir),&
!!$                         gr_vL  (DENS_VAR:GAME_VAR,i,       j,       dir),&
!!$                         intFlux(DENS_VAR:ENER_VAR,i,       j             ),dir)

              call RP(rpL(DENS_VAR:GAME_VAR), rpR(DENS_VAR:GAME_VAR), intFlux(DENS_VAR:ENER_VAR,i,j),dir)
           end do
        end do

        do i = gr_ibeg(XDIM), gr_iend(XDIM)+im
           stencilDo: do j = gr_ibeg(YDIM), gr_iend(YDIM)+jm
#ifdef BDRY_VAR
              do l = i-2*im, i+2*im
                 do m = j-2*jm, j + 2*jm
                    if (gr_V(BDRY_VAR,l,m) == 1.0) then
                    !if (.true.) then
                       gr_flux(:,i,j,dir) = intFlux(:,i,j)
                       cycle stencilDo
                    end if
                 end do
              end do
#endif              
              do var = DENS_VAR, ENER_VAR
                 cntr = 0
                 do l = i-2*im, i + 2*im
                    do m = j-2*jm, j + 2*jm
                       cntr = cntr + 1
                       flux_vec(cntr) = intFlux(var,l,m)
                    end do
                 end do
                 call NumFlux(flux_vec, Npts, gr_flux(var,i,j,dir))
              end do
           end do stencilDo
         end do
      end do

     !now we get y-fluxes
    
  else
     !correct fluxes using cell center fluxes
     !calculate interface fluxes
     do dir = XDIM, NDIM

        select case(dir)
        case(XDIM)
           ibeg = gr_ibeg(XDIM)
           iend = gr_iend(XDIM)+1
           jbeg = gr_ibeg(YDIM)
           jend = gr_iend(YDIM)
           im = 1
           jm = 0
        case(YDIM)
           ibeg = gr_ibeg(XDIM)
           iend = gr_iend(XDIM)
           jbeg = gr_ibeg(YDIM)
           jend = gr_iend(YDIM)+1
           im = 0
           jm = 1
        end select
        do i = ibeg, iend
           do j = jbeg, jend

              rpL(DENS_VAR:GAME_VAR) = gr_vR   (DENS_VAR:GAME_VAR,i-im,j-jm, dir)
              rpR(DENS_VAR:GAME_VAR) = gr_vL   (DENS_VAR:GAME_VAR,i,      j,       dir)
              
#ifdef BDRY_VAR
              ! solid internal boundary
              ! Cell j and j-1:
              if (gr_V(BDRY_VAR,i,j) > 0.0 .and. gr_V(BDRY_VAR,i,j-1) < 0.0) then
                 rpR = rpL
                 rpR(VELY_VAR) = -rpL(VELY_VAR)
!!$                 gr_vR(DENS_VAR:GAME_VAR,i,j,dir) = gr_vL(DENS_VAR:GAME_VAR,i,j,dir)
!!$                 gr_vR(VELY_VAR,i,j,dir) = -gr_vL(VELY_VAR,i,j,dir)
              end if

              if (gr_V(BDRY_VAR,i,j) < 0.0 .and. gr_V(BDRY_VAR,i,j-1) > 0.0) then
                 rpL = rpR
                 rpL(VELY_VAR) = -rpR(VELY_VAR)
!!$                 gr_vL(DENS_VAR:GAME_VAR,i,j,dir) = gr_vR(DENS_VAR:GAME_VAR,i,j,dir)
!!$                 gr_vL(VELY_VAR,i,j,dir) = -gr_vR(VELY_VAR,i,j,dir)
              end if
#endif    
!!$              call RP(gr_vR   (DENS_VAR:GAME_VAR,i-im,j-jm, dir),&
!!$                         gr_vL   (DENS_VAR:GAME_VAR,i,      j,       dir),&
!!$                         intFlux(DENS_VAR:ENER_VAR,i,       j             ),dir)
              
              call RP(rpL(DENS_VAR:GAME_VAR), rpR(DENS_VAR:GAME_VAR), intFlux(DENS_VAR:ENER_VAR,i,j),dir)
              
           end do
        end do

        do i = gr_i0(XDIM), gr_imax(XDIM)
           do j = gr_i0(YDIM), gr_imax(YDIM)
              call prim2flux(V(:,i,j),cntrFlux(:,i,j),dir)
           end do
        end do

        

        do i = gr_ibeg(XDIM), gr_iend(XDIM)+im
           do j = gr_ibeg(YDIM), gr_iend(YDIM)+jm
              do var = DENS_VAR, ENER_VAR
                 cntr = 0
                 do l = i-2*im, i + im
                    do m = j-2*jm, j + jm
                       cntr = cntr + 1
                       flux_vec(cntr) = cntrFlux(var,l,m)
                    end do
                 end do
                 flux_vec(5) = intFlux(var,i,j)
                 call NumFlux(flux_vec, Npts, gr_flux(var,i,j,dir))
              end do
           end do
        end do

        
        
     end do

  end if

  


  return
end subroutine soln_getFlux
