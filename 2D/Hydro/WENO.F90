module WENO

#include "definition.h"

contains

  subroutine gpA_betas(V, R, beta)
    use gp_data, only : gpA_Pvecs, gpA_Pvecs5, gp_det, gpA_det, gpA_det5

    integer,                   intent(IN   ) :: R
    real   , dimension(2*R+1), intent(IN   ) :: V
    real   , dimension(5    ), intent(INOUT) :: beta


    real, dimension(3) :: weno_beta
    
    integer :: i, N

    call gp_betas(V(R+1-2:R+1+2) , 2, beta(1:3))
    !    call betas(V(R+1-2:R+1+2) , 2, weno_beta(1:3))
    beta(1:3) = beta(1:3) + gp_det - 3.*LOG(2.*PI)

!!$    print *, beta(1:3)
!!$    print *, weno_beta
!!$    print *, 

    beta(4) = dot_product(V(R+1-2:R+1+2), gpA_Pvecs5(:,1))**2 + gpA_det5 - 5.*LOG(2.*PI)
    do i = 2, 5
       beta(4) = beta(4) + dot_product(V(R+1-2:R+1+2), gpA_Pvecs5(:,i))**2
    end do
    beta(4)=beta(4)
!!$    print *, V(2:6)
!!$    print *, beta(1:3), beta(4)
!!$    print *,
    N = 2*R+1
    beta(5) = dot_product(V, gpA_Pvecs(:,1))**2 + gpA_det - REAL(N)*LOG(2.*PI)
    do i = 2, N
       beta(5) = beta(5) + dot_product(V, gpA_Pvecs(:,i))**2
    end do
    beta(5)=beta(5)
    return
  end subroutine gpA_betas

  subroutine gpM_beta(V, Npts, beta)
    use gp_data, only: gpM_Pvecs
    implicit none

    integer              , intent(IN   ) :: Npts
    real, dimension(Npts), intent(IN   ) :: V
    real                 , intent(INOUT) :: beta

    integer :: i
    beta = 0.
    do i = 1, Npts
       beta = beta + dot_product(V(:), gpM_Pvecs(:,i))**2
    end do

    return
  end subroutine gpM_beta

  subroutine gp2_betas(V, R, beta)
    use gp_data , only:  gp_wr, gp_U2
    implicit none

    integer,                   intent(IN   ) :: R
    real   , dimension(2*R+1), intent(IN   ) :: V
    real   , dimension(  R+1), intent(INOUT) :: beta

    integer i,j, N, s
    real, dimension(R) :: data
    real :: Vi, P
    !real, dimension(R+1) :: gamms, kxp, f, gamms2
    !real, dimension(R+1) :: weno_beta
    !real :: gamm

    N = R+1
    Vi = V(R+1)
    !call gp_betas(V,R,beta)
    !loop over stencils
    do s = 1, N
       !make stencil leaving out Vi
       i = 1
       do j = 1, R+1-s
          data(i) = V(j+s-1)
          i = i + 1
       end do
       do j = R+1-s+2, R+1
          data(i) = V(j+s-1)
          i = i +1
       end do
       P = gp_U2(s)*(dot_product(gp_wr(:,s),data-vi))**2
       beta(s) = abs(P)!*(sum(vi-data)**2)!sqrt(2.*PI*gp_U2(s))*EXP(P)
    end do
    return
  end subroutine gp2_betas

  subroutine gp_betas(V, R, beta)
    use gp_data , only:  gp_Pvecs, gp_Zvecs
    implicit none

    integer,                   intent(IN   ) :: R
    real   , dimension(2*R+1), intent(IN   ) :: V
    real   , dimension(  R+1), intent(INOUT) :: beta

    integer i, N, s
    !real, dimension(R+1) :: gamms, kxp, f, gamms2
    real, dimension(R+1) :: stencil, fdotv
    !real :: gamm

    N = R+1
    !loop over stencils
    do s = 1 ,N
       !loop over eigen values
       stencil = V(s:s+R) !- V(R+1)
       beta(s) = (dot_product(stencil, gp_Pvecs(:,1))**2)
       do i = 2, N
          beta(s) = beta(s) + (dot_product(stencil, gp_Pvecs(:,i))**2)
          !fdotv(i) = dot_product(stencil, gp_Pvecs(:,1))**2
       end do
       !fdotv(maxloc(fdotv)) = 0.
       !beta(s) = sum(fdotv)
    end do

    return
  end subroutine gp_betas



  subroutine betas(V, R, beta)
    !subroutine to calculate the smoothness-indicators for a WENO scheme on a 2R+1 point stencil
    implicit none

    integer,                   intent(IN   ) :: R
    real   , dimension(2*R+1), intent(IN   ) :: V
    real   , dimension(  R+1), intent(INOUT) :: beta

    integer :: i

    i = R+1
    !print *, 'betas'

    select case(R)
    case(1)
       beta(1) = (V(i  )-V(i-1))**2
       beta(2) = (V(i+1)-V(i  ))**2
    case(2)
       beta(1) = 13./12.*(V(i-2) - 2.*V(i-1) + V(i  ) )**2 + 0.25*(   V(i-2) - 4.*V(i-1) + 3.*V(i  ) )**2
       beta(2) = 13./12.*(V(i-1) - 2.*V(i  ) + V(i+1) )**2 + 0.25*(   V(i-1)             -    V(i+1) )**2
       beta(3) = 13./12.*(V(i  ) - 2.*V(i+1) + V(i+2) )**2 + 0.25*(3.*V(i  ) - 4.*V(i+1) +    V(i+2) )**2
    case(3)
       beta(1) = V(i-3)*(  547.*V(i-3) -  3882.*V(i-2) + 4642.*V(i-1) - 1854.*V(i  )) + &
                 V(i-2)*( 7043.*V(i-2) - 17246.*V(i-1) + 7042.*V(i  )               ) + &
                 V(i-1)*(11003.*V(i-1) -  9402.*V(i  )                              ) + 2107.*V(i  )**2
       beta(2) = V(i-2)*(  267.*V(i-2) -  1642.*V(i-1) + 1602.*V(i  ) -  494.*V(i+1)) + &
                 V(i-1)*( 2843.*V(i-1) -  5966.*V(i  ) + 1922.*V(i+1)               ) + &
                 V(i  )*( 3443.*V(i  ) -  2522.*V(i+1)                              ) +  547.*V(i+1)**2
       beta(3) = V(i-1)*(  547.*V(i-1) -  2522.*V(i  ) + 1922.*V(i+1) -  494.*V(i+2)) + &
                 V(i  )*( 3443.*V(i  ) -  5966.*V(i+1) + 1602.*V(i+2)               ) + &
                 V(i+1)*( 2843.*V(i+1) -  1642.*V(i+2)                              ) +  267.*V(i+2)**2
       beta(4) = V(i  )*( 2107.*V(i  ) -  9402.*V(i+1) + 7042.*V(i+2) - 1854.*V(i+3)) + &
                 V(i+1)*(11003.*V(i+1) - 17246.*V(i+2) + 4642.*V(i+3)               ) + &
                 V(i+2)*( 7043.*V(i+2) -  3882.*V(i+3)                              ) +  547.*V(i+3)**2
    case DEFAULT
       beta = 1.
    end select

    return
  end subroutine betas

end module WENO
