subroutine soln_gpFlux(F, Npts, Fiph)
  !get the numerical fluxes using GP
  !I think I'm only intending to use this with the interface fluxes
#include "definition.h"

  use gp_data, only: gpF_Z2, gpF_Z4, gp_radius
  use sim_data, only: sim_DongwookFlux
  
  implicit none

  integer, intent(IN) :: Npts
  real, dimension(Npts), intent(IN) :: F
  real, intent(OUT) :: Fiph

  logical :: switch
  integer :: R
  real :: fac2, fac4, D2, D4
  real, dimension(3) :: coeff3
  real, dimension(5) :: coeff5


  fac2 = -1./24.
  fac4 = 3./640.
  coeff3 = (/1., -2., 1. /)
  coeff5 = (/1., -4., 6., -4., 1. /)

  coeff3 = (/1.0041350655174979, -2.0082530792912037, 1.0041350655174979/)
  coeff5 = (/1.0208093285109237, -4.0827232307974528, 6.1238292148346121,&
            -4.0827232307974528,  1.0208093285109237 /)

!!$  if (sim_DongwookFlux) then
!!$     !doesn't work :(
!!$     !Fiph = dot_product(gpF_Z(1,:), F)
!!$  else
  R = 2
  !Adam's way, not done yet!!
  D2 = dot_product(gpF_Z2(:), F(2:4))
  D4 = dot_product(gpF_Z4(:), F)

!!$  D2 = dot_product(coeff3, F(2:4))
!!$  D4 = dot_product(coeff5, F)
  !print *, Fiph - (F(3) + fac2*D2 + fac4*D4) !difference with Finite Difference fluxes
  
  !end if

  Fiph = F(R+1) + fac2*D2 + fac4*D4

  

  return

end subroutine soln_gpFlux
