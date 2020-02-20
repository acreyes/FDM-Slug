module primconsflux

#include "definition.h"
  
  use grid_data
  use sim_data, only : sim_gamma, sim_smallPres
  use eos, only : eos_cell
  
contains

  subroutine prim2cons(V,U)
    implicit none
    real, dimension(NUMB_VAR), intent(IN)  :: V
    real, dimension(NSYS_VAR), intent(OUT) :: U

    real :: ekin, eint

    U(DENS_VAR) = V(DENS_VAR)
    U(MOMX_VAR) = V(DENS_VAR)*V(VELX_VAR)

    ekin = 0.5*V(DENS_VAR)*V(VELX_VAR)**2
    eint = V(PRES_VAR)/(V(GAME_VAR)-1.)
    U(ENER_VAR) = ekin + eint
    
  end subroutine prim2cons


  subroutine cons2prim(U,V)
    implicit none
    real, dimension(NSYS_VAR), intent(IN)  :: U
    real, dimension(NUMB_VAR), intent(OUT) :: V
    real :: eint, ekin, pres
    
    V(DENS_VAR) = U(DENS_VAR)
    V(VELX_VAR) = U(MOMX_VAR)/U(DENS_VAR)
    ekin = 0.5*V(DENS_VAR)*V(VELX_VAR)**2
    eint = max(U(ENER_VAR) - ekin, sim_smallPres) !eint=rho*e
    eint = eint/U(DENS_VAR)
    ! get pressure by calling eos
    call eos_cell(U(DENS_VAR),eint,sim_gamma,pres)
    V(PRES_VAR) = pres
    V(EINT_VAR) = eint*U(DENS_VAR)
    V(GAMC_VAR) = sim_gamma
    V(GAME_VAR) = sim_gamma
    
  end subroutine cons2prim

  subroutine prim2flux(V,Flux)
    implicit none
    real, dimension(NUMB_VAR), intent(IN)  :: V
    real, dimension(NSYS_VAR), intent(OUT) :: Flux

    real :: ekin,eint,ener
    
    Flux(DENS_VAR) = V(DENS_VAR)*V(VELX_VAR)
    Flux(MOMX_VAR) = Flux(DENS_VAR)*V(VELX_VAR) + V(PRES_VAR)
    ekin = 0.5*V(VELX_VAR)*V(VELX_VAR)*V(DENS_VAR)
    eint = V(PRES_VAR)/(V(GAME_VAR)-1.)
    ener = ekin + eint
    Flux(ENER_VAR) = V(VELX_VAR)*(ener + V(PRES_VAR))
    
  end subroutine prim2flux

  subroutine cons2flux
    implicit none
    
  end subroutine cons2flux
  
end module primconsflux
