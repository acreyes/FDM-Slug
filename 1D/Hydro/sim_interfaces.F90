module sim_interfaces
#include "definition.h"
  abstract interface
     subroutine recon1D(dt, Nx, stencil, vL, vR)
       
       integer, intent(IN) :: Nx
       real   , intent(IN) :: dt
       
       real, dimension(NUMB_VAR), intent(INOUT ) :: vL, vR
       real, dimension(Nx, NUMB_VAR), intent(IN) :: stencil
     end subroutine recon1D
  end interface

  abstract interface
     subroutine rmn_slvr(vR, vL, flux)

       real, dimension(NUMB_VAR), intent(IN ) :: vL, vR
       real, dimension(NSYS_VAR), intent(OUT) :: flux

     end subroutine rmn_slvr
  end interface
  
end module sim_interfaces
