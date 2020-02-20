subroutine soln_FOG(dt, Nx, V, vL, vR)

#include "definition.h"  

  implicit none
  real, intent(IN) :: dt
  integer, intent(IN) :: Nx
  real, dimension(Nx,NUMB_VAR), intent(IN) :: V
  real, dimension(NUMB_VAR), intent(INOUT) :: vL, vR

  vL(:) = V(1,:)
  vR(:) = V(1,:)

  return
end subroutine soln_FOG
