subroutine soln_RK2(dt, k, Vk, Flux)
  !performs the jth step of the RK3 algorithm
  !Total Flux = (k1 + k2 + 2k3)/6

#include "definition.h"

  use grid_data
  use primconsflux

  implicit none
  real, intent(IN) :: dt
  integer, intent(IN) :: k
  real, dimension(NUMB_VAR,gr_imax(XDIM),gr_imax(YDIM)), intent(INOUT) :: Vk
  real, dimension(NSYS_VAR,gr_imax(XDIM),gr_imax(YDIM),2), intent(INOUT) :: Flux !(var,i,j,ndim)

  real :: dtx,dty
  integer :: i,j
  real, dimension(NSYS_VAR,gr_imax(XDIM),gr_imax(YDIM)) :: Uk
  call soln_reconstruct(dt, Vk)
  call soln_getFlux(Vk)
  dtx = dt/gr_dx
  dty = dt/gr_dy
  Uk = 0.
  if (k == 1) then
     do i = gr_ibeg(XDIM), gr_iend(XDIM)
        do j = gr_ibeg(YDIM), gr_iend(YDIM)
           Uk(DENS_VAR:ENER_VAR,i,j) = gr_U(DENS_VAR:ENER_VAR,i,j)                                - &
                dtx*(gr_flux(DENS_VAR:ENER_VAR,i+1,j,XDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,XDIM)) - &
                dty*(gr_flux(DENS_VAR:ENER_VAR,i,j+1,YDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,YDIM))
           !print *, "x", gr_flux(:,i,j+1,2) - gr_flux(:,i,j,2), j
           call cons2prim(Uk(DENS_VAR:ENER_VAR,i,j),Vk(DENS_VAR:GAME_VAR,i,j))
        end do
     end do
     
  end if

  
  Flux(:,:,:,:) = Flux(:,:,:,:) + 0.5*gr_flux(:,:,:,:)
end subroutine soln_RK2
