subroutine gp_AOinit()

#include "definition.h"

  use grid_data, only: gr_dx
  use gp_data
  use linalg
  use GP

  implicit none


  real, dimension(2*gp_radius+1,2*gp_radius+1) :: C, L
  real, dimension(2,2*gp_radius+1) :: T, Z
  real, dimension(2*gp_radius+1)   :: stencil, u

  real, dimension(5,5) :: C5, L5
  real, dimension(2,5) :: T5, Z5
  real, dimension(5)   :: stencil5, u5

  integer :: i,j,N,R,m,LR,COL,ROW
  real :: Xdel



  

  if (gp_el == 0.) then
     gp_el = gr_dx*gp_eldel
  elseif(gp_eldel == 0.) then
     gp_eldel = gp_el/gr_dx
  end if

  Xdel = gp_el/gr_dx
  gp_Xdel = Xdel

  !initialize
  R = gp_radius
  N = 2*R+1
  C = 0.
  T = 0.
  u = 1.


  !This handles the global (2R+1) stencil
  do i = 1, N
     stencil(i) = REAL(i-gp_radius-1)
  end do

  !first thing is to calculate the covariance matrix according to eq. 15
  !since C is symmetric only bother with one side of the diaganol
  do i = 1,N
     do j = 1,N
        C(i,j) = SE(stencil(i),stencil(j),Xdel)
     end do
     T(1, i) = SE(-0.5, stencil(i), Xdel)
     T(2, i) = SE( 0.5, stencil(i), Xdel)
  end do

  !now we need to solve the linear eqns for v & Z (see eqs 30-32)
  call chol(C, N, L)
  call solve_CZT(C, Z, T, L, N)

  !now we do the 5 point stencil we wish to hybridize with the 2R+1 stencil


  R = 2
  N = 2*R+1
  C5 = 0.
  T5 = 0.
  u5 = 0.
  do i = 1, N
     stencil5(i) = REAL(i-R-1)
  end do

  do i = 1,N
     do j = 1,N
        C5(i,j) = SE(stencil5(i),stencil5(j),Xdel)
     end do
     T5(1, i) = SE(-0.5, stencil5(i), Xdel)
     T5(2, i) = SE( 0.5, stencil5(i), Xdel)
  end do

  !now we need to solve the linear eqns for v & Z (see eqs 30-32)
  call chol(C5, N, L5)
  call solve_CZT(C5, Z5, T5, L5, N)

  !truncate to double precision
  allocate(gpA_z (2,2*R+1)); gpA_z  = Z
  allocate(gpA_z5(2,5)    ); gpA_z5 = Z5


end subroutine gp_AOinit
  
