subroutine gp_AOeigen()

#include "definition.h"

  use grid_data
  use gp_data

  implicit none

  real, allocatable, dimension(:,:) :: C, C5
  integer :: N, LDA, LWORK, INFO, i, j, k, m, LDA5, LWORK5

  real, allocatable, dimension(:) :: W, WORK, W5, WORK5

  real :: eldel

!!$  if (gp_el == 0.) then
!!$     gp_el = gr_dx*gp_eldel
!!$  elseif(gp_eldel == 0.) then
!!$     gp_eldel = gp_el/gr_dx
!!$  end if
!!$  print *, eldel
  eldel=gp_sigdel
  N = 2*gp_radius + 1
  LDA = N
  LWORK = 66*N
  allocate(C(N,N))
  allocate(W(N))
  allocate(WORK(LWORK))

  do i = 1, N
     do j = 1, N
        C(i,j) = SE(REAL(i), REAL(j), eldel)
     end do
  end do

  call DSYEV('V', 'L', N, C, LDA, W, WORK, LWORK, INFO)

  gpA_det = -LOG(PRODUCT(W))
  allocate(gpA_Pvecs(N,N))
  do i = 1, N
     gpA_Pvecs(:,i) = C(:,i)/sqrt(W(i))
  end do

  N = 5
  LDA5 = N
  LWORK5 = 66*N
  allocate(C5(N,N))
  allocate(W5(N))
  allocate(WORK5(LWORK5))

  do i = 1, N
     do j = 1, N
        C5(i,j) = SE(REAL(i), REAL(j), eldel)
     end do
  end do

  call DSYEV('V', 'L', N, C5, LDA5, W5, WORK5, LWORK5, INFO)

  gpA_det5 = -LOG(PRODUCT(W5))
  allocate(gpA_Pvecs5(N,N))
  do i = 1, N
     gpA_Pvecs5(:,i) = C5(:,i)/sqrt(W5(i))
  end do


  contains
  function SE(x, y, eldel) result(f)
    implicit none
    real, intent(IN) :: x, y, eldel
    real :: f, r
    r = abs(x-y)
    f = EXP( -0.5*(r/eldel)**2 )
    return
  end function SE

end subroutine gp_AOeigen
