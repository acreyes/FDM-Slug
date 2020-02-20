subroutine soln_gpMD(dt, Npts, V, vL, vR)
  
#include "definition.h"

  use WENO
  use gp_data
  use eigensystem
  use sim_data
  
  implicit none

  real   , intent(IN) :: dt
  integer, intent(IN) :: Npts

  real, dimension(Npts    , NUMB_VAR), intent(IN)    :: V
  real, dimension(NUMB_VAR, NDIM    ), intent(INOUT) :: vL, vR

  real, dimension(NSYS_VAR,NUMB_WAVE) :: reig, leig
  real, dimension(NSYS_VAR) :: vecL, vecR
  real, dimension(Npts)     :: stencil
  real, dimension(5)        :: stencil1D
  real, dimension(4)        :: lin_w, smth_ind, non_lin_w, tempL, tempR
  
  integer :: var, dir, s, zi
  real    :: gammHI, gammLO, sumW

  !initilizations
  stencil1D = 0.
  stencil = 0.
  !linear weights for adaptive order hybridization w/ weno5 stencil gp prediction
  gammHI = 0.85
  gammLO = 0.85
  lin_w(1) = gammHI                      !multiD prediction
  lin_w(2) = 0.5*(1.-gammHI)*(1.-gammLO) !left biased 1Dstencil
  lin_w(3) = (1.-gammHI)*gammLO          !central ENO stencil
  lin_w(4) = 0.5*(1.-gammHI)*(1.-gammLO) !right bieased 1Dstencil
  !copy cell centered value to faces
  vL(:,XDIM) = V(gp_cntrPt,:)
  vR(:,XDIM) = V(gp_cntrPt,:)
  vL(:,YDIM) = V(gp_cntrPt,:)
  vR(:,YDIM) = V(gp_cntrPt,:)
  do dir = XDIM, NDIM
     select case(dir)
     case(XDIM)
        zi = 1
     case(YDIM)
        zi = 3
     end select
     if (sim_charLimiting) then
        call eigenvectors(V(gp_cntrPt,:),.false.,reig,leig,dir)
     end if
     
     do var = 1, NSYS_VAR
        if (sim_charLimiting) then
           do s = 1, Npts
              stencil(s) = dot_product(leig(:,var), V(s,DENS_VAR:ENER_VAR))
           end do
        else
           do s = 1, Npts 
              stencil(s) = V(s,var)
           end do
        end if
        do s = 1, 5
           stencil1D(s) = stencil(gp_1Dstencil(s,dir))
        end do
        
        !get smoothness indicators
        call gpM_beta(stencil, Npts, smth_ind(1))
        call betas(stencil1D, 2, smth_ind(2:4))

        !non-linear weights
        do s = 1, 4
           non_lin_w(s) = lin_w(s)/(sim_WENeps+smth_ind(s))**sim_mval
        end do
        sumW = SUM(non_lin_w)
        non_lin_w = non_lin_w/sumW
        !now we get predictions from each stencil
        tempL(1) = dot_product(gpM_z(zi  ,:), stencil(:))
        tempR(1) = dot_product(gpM_z(zi+1,:), stencil(:))
        do s = 1, 3
           tempL(s+1) = dot_product(gp_zk(1,1:3,s), stencil1D(s:s+2))
           tempR(s+1) = dot_product(gp_zk(2,1:3,s), stencil1D(s:s+2))
        end do
        !now we take convex combination
        vecL(var) = dot_product(non_lin_w,tempL)
        vecR(var) = dot_product(non_lin_w,tempR)
!!$        vecL(var) = tempL(3)
!!$        vecR(var) = tempR(3)
     end do
     do var = DENS_VAR, PRES_VAR
        if (sim_charLimiting) then
           vL(var,dir) = dot_product(vecL(:), reig(var,:))
           vR(var,dir) = dot_product(vecR(:), reig(var,:))
        else
           vL(var,dir) = vecL(var)
           vR(var,dir) = vecR(var)
        end if
     end do
  end do


  return


end subroutine soln_gpMD
