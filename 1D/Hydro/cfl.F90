subroutine cfl(dt)

#include "definition.h"
  
  use grid_data
  use sim_data, only: sim_cfl

  implicit none
  real, intent(OUT) :: dt
  integer :: i
  real :: maxSpeed, lambda, cs

  maxSpeed = -1.e30
  !! update conservative vars
  do i = gr_ibeg, gr_iend
     cs = sqrt(gr_V(GAMC_VAR,i)*gr_V(PRES_VAR,i)/gr_V(DENS_VAR,i))
     lambda=(abs(gr_V(VELX_VAR,i)) + cs)
     maxSpeed=max(maxSpeed,lambda)
  end do
  ! cfl
  dt = sim_cfl*gr_dx/maxSpeed

  return

end subroutine cfl
