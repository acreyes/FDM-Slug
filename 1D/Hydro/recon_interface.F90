module recon_interface
#include "definition.h"
  abstract interface
     subroutine recon_subroutine(dt, Nx, stencil, F)
       
       integer, intent(IN) :: Nx
       real   , intent(IN) :: dt
       
       real, dimension(NSYS_VAR), intent(INOUT ) :: F
       real, dimension(Nx, NSYS_VAR), intent(IN) :: stencil
     end subroutine recon_subroutine
  end interface
  
end module recon_interface
