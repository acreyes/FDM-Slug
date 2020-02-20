subroutine soln_gpAO(dt, Nx, V, vL, vR, dir)

#include "definition.h"

  use gp_data
  use eigensystem
  use WENO
  use sim_data, only: sim_WENeps, sim_mval, sim_charLimiting

  implicit none
  real, intent(IN) :: dt
  integer, intent(IN) :: Nx, dir 
  real, dimension(Nx, NUMB_VAR), intent(IN   ) :: V
  real, dimension(NUMB_VAR),   intent(INOUT) :: vL, vR

  integer :: var, s, R
  real    :: w_norm, gamm_HI, gamm_LO, vML, vMR, v5L, v5R, nuM, nu5, ratio

  real, dimension(NSYS_VAR,NUMB_WAVE) :: reig, leig

  real, dimension(Nx) :: stencil
  real, dimension(5) :: smth_ind
  real, dimension(2) :: beta
  real, dimension(5,2) :: nonLin_w, ENO_intp
  real, dimension(5) :: lin_w, nonLin
  real, dimension(NSYS_VAR) :: tempL, tempR

  R = gp_radius
  vL = V(R+1,:)
  vR = V(R+1,:)

  gamm_HI = 0.85
  gamm_LO = 0.85

  lin_w(1) = 0.5*(1.-gamm_HI)*(1.-gamm_LO)
  lin_w(2) =     (1.-gamm_HI)*gamm_LO
  lin_w(3) = 0.5*(1.-gamm_HI)*(1.-gamm_LO)
  lin_w(4) = gamm_HI
  lin_w(5) = gamm_HI

  if (sim_charLimiting) then
     call eigenvectors(V(R+1,:), .false., reig, leig, dir)
  end if

  do var = 1, NSYS_VAR
     do s = 1, Nx
        if (sim_charLimiting) then
           stencil(s) = dot_product(leig(:,var), V(s,DENS_VAR:ENER_VAR))
        else
           stencil(s) = V(s,var)
        end if
     end do

     call gpA_betas(stencil, R, smth_ind)

     do s = 1, 5
        nonLin(s) = lin_w(s)/(sim_WENeps+smth_ind(s))**sim_mval
     end do
     w_norm = SUM(nonLin(1:3))
     nonLin_w(:,1) = nonLin(:)/(w_norm+nonLin(4)) !1 for 5pt
     nonLin_w(:,2) = nonLin(:)/(w_norm+nonLin(5)) !2 for 2R+1

     !calculate ENO_intp
     do s = 1, 3
        ENO_intp(s,1) = dot_product(gp_zk(1,1:3,s), stencil(R+1+s-3:R+1+s-1))
        ENO_intp(s,2) = dot_product(gp_zk(2,1:3,s), stencil(R+1+s-3:R+1+s-1))
     end do

     ENO_intp(4,1) = dot_product(gpA_z5(1,:),stencil(R+1-2:R+1+2))
     ENO_intp(4,2) = dot_product(gpA_z5(2,:),stencil(R+1-2:R+1+2))

     ENO_intp(5,1) = dot_product(gpA_z(1,:),stencil(:))
     ENO_intp(5,2) = dot_product(gpA_z(2,:),stencil(:))

     !here we calculate the hybridiztions w/ the 3 point stencils
     !5 to 3 hybridization
     ratio = nonLin_w(4,1)/lin_w(4)
     v5L = ratio * ENO_intp(4,1) + &
           dot_product(nonLin_w(1:3,1)-ratio*lin_w(1:3), ENO_intp(1:3,1))
     v5R = ratio * ENO_intp(4,2) + &
           dot_product(nonLin_w(1:3,1)-ratio*lin_w(1:3), ENO_intp(1:3,2))
     !2R+1 to 3 hybridization
     ratio = nonLin_w(5,2)/lin_w(5)
     vML = ratio * ENO_intp(5,1) + &
           dot_product(nonLin_w(1:3,2)-ratio*lin_w(1:3), ENO_intp(1:3,1))
     vMR = ratio * ENO_intp(5,2) + &
           dot_product(nonLin_w(1:3,2)-ratio*lin_w(1:3), ENO_intp(1:3,2))

     !now we do the recursive hybridization
     nuM = gamm_HI/(sim_WENeps+smth_ind(5))**sim_mval
     nu5 = (1.-gamm_HI)/(sim_WENeps+smth_ind(4))**sim_mval

     w_norm = nuM + nu5
     nuM = nuM/w_norm
     nu5 = nu5/w_norm

     !7,5,3 hybridization
     tempL(var) = nuM/gamm_HI * vML + (nu5-nuM*(1.-gamm_HI)/gamm_HI) * v5L
     tempR(var) = nuM/gamm_HI * vMR + (nu5-nuM*(1.-gamm_HI)/gamm_HI) * v5R
     
!!$     tempL(var) = v5L !just do 5,3 AO
!!$     tempR(var) = v5R !just do 5,3 AO

!!$     tempL(var) = vML !just do 7,3 AO
!!$     tempR(var) = vMR !just do 7,3 AO

  end do

  do var = DENS_VAR, PRES_VAR
     if (sim_charLimiting) then
        vL(var) = dot_product(tempL(:), reig(var,:))
        vR(var) = dot_product(tempR(:), reig(var,:))
     else
        vL(var) = tempL(var)
        vR(var) = tempR(var)
     end if
  end do

  

end subroutine soln_gpAO
