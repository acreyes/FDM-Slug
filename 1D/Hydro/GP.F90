module GP
#include "definition.h"

  abstract interface
     function gp_kernel1D(x, y, eldel) result(f)
       real, intent(IN) :: x, y, eldel
       real :: f
     end function gp_kernel1D
  end interface

contains

  function SE(x, y, eldel) result(f)
    real, intent(IN) :: x, y, eldel
    real :: f, r

    r = abs(x-y)/eldel
    
    f = EXP(-0.5*r**2)
    return
  end function SE

  function d2_SE(x, y, eldel) result(f)
    implicit none
    real, intent(IN) :: x, y, eldel
    real :: f, ell

    f = SE(x,y,eldel)/(eldel**2)*(((x-y)/eldel)**2-1.)
    return
    
  end function d2_SE

  function d4_SE(x, y, eldel) result(f)
    implicit none
    real, intent(IN) :: x, y, eldel
    real :: f, ell, xmy

    xmy = (x - y)/eldel

    f = SE(x,y,eldel)/(eldel**4)*( (xmy)**4 - 6.*(xmy)**2 + 3. )
    return
  end function d4_SE

end module GP
