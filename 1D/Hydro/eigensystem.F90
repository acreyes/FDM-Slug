module eigensystem

#include "definition.h"
  
  use grid_data

contains

  subroutine eigenvalues(V,lambda)
    implicit none

    real, dimension(NUMB_VAR), intent(IN)  :: V
    real, dimension(NUMB_WAVE),intent(OUT) :: lambda

    real :: a, u

    ! sound speed
    a = sqrt(V(GAMC_VAR)*V(PRES_VAR)/V(DENS_VAR))!;print*,a,V(GAMC_VAR),V(PRES_VAR),V(DENS_VAR)
    u = V(VELX_VAR)
    
    lambda(SHOCKLEFT) = u - a
    lambda(CTENTROPY) = u
    lambda(SHOCKRGHT) = u + a
    
    return
  end subroutine eigenvalues


  
  subroutine right_eigenvectors(V,conservative,reig)
    implicit none
    real, dimension(NUMB_VAR), intent(IN)  :: V
    logical :: conservative
    real, dimension(NSYS_VAR,NUMB_WAVE), intent(OUT) :: reig

    real :: a, u, d, g, ekin, hdai, hda, H, p
    
    ! sound speed, and others
    a = sqrt(V(GAMC_VAR)*V(PRES_VAR)/V(DENS_VAR))
    u = V(VELX_VAR)
    d = V(DENS_VAR)
    g = V(GAMC_VAR) - 1.
    ekin = 0.5*u**2
    hdai = 0.5*d/a
    hda  = 0.5*d*a
    p = V(PRES_VAR)
    H = ekin + a**2/g
    
    if (conservative) then
       !! Conservative eigenvector
!!$       reig(DENS_VAR,SHOCKLEFT) = 1.
!!$       reig(VELX_VAR,SHOCKLEFT) = u - a
!!$       reig(PRES_VAR,SHOCKLEFT) = ekin + a**2/g - a*u
!!$       reig(:,SHOCKLEFT) = -hdai*reig(:,SHOCKLEFT)
!!$
!!$       reig(DENS_VAR,CTENTROPY) = 1.
!!$       reig(VELX_VAR,CTENTROPY) = u
!!$       reig(PRES_VAR,CTENTROPY) = ekin
!!$       
!!$       reig(DENS_VAR,SHOCKRGHT) = 1.
!!$       reig(VELX_VAR,SHOCKRGHT) = u + a
!!$       reig(PRES_VAR,SHOCKRGHT) = ekin + a**2/g + a*u
!!$       reig(:,SHOCKRGHT) = hdai*reig(:,SHOCKRGHT)


       reig(DENS_VAR,SHOCKLEFT) = 1.
       reig(VELX_VAR,SHOCKLEFT) = u - a
       reig(PRES_VAR,SHOCKLEFT) = H - a*u

       reig(DENS_VAR,CTENTROPY) = 1.
       reig(VELX_VAR,CTENTROPY) = u
       reig(PRES_VAR,CTENTROPY) = ekin
       
       reig(DENS_VAR,SHOCKRGHT) = 1.
       reig(VELX_VAR,SHOCKRGHT) = u + a
       reig(PRES_VAR,SHOCKRGHT) = H + u*a
    else
       !! Primitive eigenvector
       reig(DENS_VAR,SHOCKLEFT) = -hdai
       reig(VELX_VAR,SHOCKLEFT) = 0.5
       reig(PRES_VAR,SHOCKLEFT) = -hda

       reig(DENS_VAR,CTENTROPY) = 1.
       reig(VELX_VAR,CTENTROPY) = 0.
       reig(PRES_VAR,CTENTROPY) = 0.

       reig(DENS_VAR,SHOCKRGHT) = hdai
       reig(VELX_VAR,SHOCKRGHT) = .5
       reig(PRES_VAR,SHOCKRGHT) = hda
              
    endif
    
    return
  end subroutine right_eigenvectors


  subroutine left_eigenvectors(V,conservative,leig)
    implicit none
    real, dimension(NUMB_VAR), intent(IN)  :: V
    logical :: conservative
    real, dimension(NSYS_VAR,NUMB_WAVE), intent(OUT) :: leig

    real :: a, u, d, g, gi, ekin, hdai, hda, ha2i
    
    ! sound speed, and others
    a = sqrt(V(GAMC_VAR)*V(PRES_VAR)/V(DENS_VAR))
    u = V(VELX_VAR)
    d = V(DENS_VAR)
    g = V(GAMC_VAR) - 1.
    gi = 1./g
    ekin = 0.5*u**2
    hdai = 0.5*d/a
    hda  = 0.5*d*a
    ha2i = 0.5/(a**2)
    
    if (conservative) then
       !! Conservative eigenvector
!!$       leig(DENS_VAR,SHOCKLEFT) = -ekin - a*u*gi
!!$       leig(VELX_VAR,SHOCKLEFT) = u+a*gi
!!$       leig(PRES_VAR,SHOCKLEFT) = -1.
!!$       leig(:,SHOCKLEFT) = g*leig(:,SHOCKLEFT)/(d*a)
!!$
!!$       leig(DENS_VAR,CTENTROPY) = d*(-ekin + gi*a**2)/a
!!$       leig(VELX_VAR,CTENTROPY) = d*u/a
!!$       leig(PRES_VAR,CTENTROPY) = -d/a
!!$       leig(:,CTENTROPY) = g*leig(:,CTENTROPY)/(d*a)
!!$       
!!$       leig(DENS_VAR,SHOCKRGHT) = ekin - a*u*gi
!!$       leig(VELX_VAR,SHOCKRGHT) = -u+a*gi
!!$       leig(PRES_VAR,SHOCKRGHT) = 1.
!!$       leig(:,SHOCKRGHT) = g*leig(:,SHOCKRGHT)/(d*a)

       !here is where the i
       leig(DENS_VAR,SHOCKLEFT) = ha2i*(g*ekin + u*a)
       leig(VELX_VAR,SHOCKLEFT) = -ha2i*(g*u+a)
       !leig(VEL2_VAR,SHOCKLEFT) = -ha2i*g*v2
       leig(PRES_VAR,SHOCKLEFT) = g*ha2i

       leig(DENS_VAR,CTENTROPY) = 1.-ha2i*g*ekin*2.
       leig(VELX_VAR,CTENTROPY) = g*u/(a**2)
       !leig(VEL2_VAR,CTENTROPY) = g*v2/(a**2)
       leig(PRES_VAR,CTENTROPY) = -g/(a**2)
       
       leig(DENS_VAR,SHOCKRGHT) = ha2i*(g*ekin-u*a)
       leig(VELX_VAR,SHOCKRGHT) = -ha2i*(g*u-a)
       !leig(VEL2_VAR,SHOCKRGHT) = -ha2i*g*v2
       leig(PRES_VAR,SHOCKRGHT) = ha2i*g
    else
       !! Primitive eigenvector
       leig(DENS_VAR,SHOCKLEFT) = 0.
       leig(VELX_VAR,SHOCKLEFT) = 1.
       leig(PRES_VAR,SHOCKLEFT) = -1./(d*a)

       leig(DENS_VAR,CTENTROPY) = 1.
       leig(VELX_VAR,CTENTROPY) = 0.
       leig(PRES_VAR,CTENTROPY) = -1./(a**2)

       leig(DENS_VAR,SHOCKRGHT) = 0.
       leig(VELX_VAR,SHOCKRGHT) = 1.
       leig(PRES_VAR,SHOCKRGHT) = 1./(d*a)
       

    endif
    
    return
  end subroutine left_eigenvectors


  
end module eigensystem
