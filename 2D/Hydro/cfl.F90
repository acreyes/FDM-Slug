subroutine cfl(dt)

#include "definition.h"
  
  use grid_data
  use sim_data, only: sim_cfl
  use eigensystem
  use block_data, only: bl_delT

  implicit none
  real, intent(OUT) :: dt
  !real, allocatable, save :: delT[:]
  integer :: i,j, var, step
  real :: maxSpeed, lambda, cs, u, v, h, temp
  real, dimension(NUMB_WAVE) :: lambdax, lambday
  !allocate(delT[*])
  maxSpeed = -1e30
  bl_delT = 1e30
!  gr_maxalphas = -1e30
  !! update conservative vars
  do i = gr_ibeg(XDIM), gr_iend(XDIM)
     do j = gr_ibeg(YDIM), gr_iend(YDIM)
#ifdef BDRY_VAR        
        if (gr_V(BDRY_VAR,i,j) == -1.0) then
           cs = sqrt(gr_V(GAMC_VAR,i,j)*gr_V(PRES_VAR,i,j)/gr_V(DENS_VAR,i,j))
           u = gr_V(VELX_VAR,i,j); v = gr_V(VELY_VAR,i,j)
           u = abs(u) + cs
           v = abs(v) + cs
           bl_delT = min(gr_dx/u, gr_dy/v, bl_delT)
        end if
#else
        cs = sqrt(gr_V(GAMC_VAR,i,j)*gr_V(PRES_VAR,i,j)/gr_V(DENS_VAR,i,j))
        u = gr_V(VELX_VAR,i,j); v = gr_V(VELY_VAR,i,j)
        u = abs(u) + cs
        v = abs(v) + cs
        bl_delT = min(gr_dx/u, gr_dy/v, bl_delT)
#endif           
     end do
  end do
  ! cfl
!!$  Sync All
!!$  call co_min(bl_delT)
  !dt = sim_cfl/bl_delT
  step = 2
  do while (step/2 .le. num_images())
     sync all
     if (this_image() + step/2 .le. num_images()) then
        temp = min(bl_delT, bl_delT[this_image()+step/2])
     else
        temp = bl_delT
     end if
     sync all
     bl_delT = temp
     step = step*2
  end do
  sync all
  if (this_image() .ne. 1) bl_delT = bl_delT[1]
  dt = sim_cfl*bl_delT
  !deallocate(delT)
  return

end subroutine cfl
