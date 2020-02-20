subroutine sim_GPinit2()
  !want to initialize all the GP parameters that do not depend on the data, and only on the grid geometries.
  !I am using sig/del = 2
  !This includes calculating the following:
  !****** the covariance matrix C 
  !****** the cross-correlation matrix T 
  !****** The corresponding vectors and matrices to be used in actual calculations (see eqs 30-32) v & Z
  !******  ****** C.v  = u  (gr_GPv)
  !******  ****** C.Z* = T* (gr_GPZ)

  !for the case of general WENO like candidate stencils for the global stencil radius R
  !the global stencil will be of size 2R+1
  !then each of the smaller stencils will be of size R+1 --> 3 for the case of WENO5
  !and there will be R+1 of these smaller stencils (l) so l = 1, R+1
  !each stencil S_l will extend from i-R+(l-1) to i+(l-1) --> R=1 for WENO5 and l
#include "definition.h"

  use grid_data
  use linalg
  use GP

  implicit none
  !integer (kind=4), parameter :: N = 2*gr_radius+1+1
  real, dimension(2*gr_radius+1,2*gr_radius+1) :: C, L, W, Vt, Sigmai
  real, dimension(2,2*gr_radius+1) :: T
  real, dimension(2*gr_radius+1,2) :: B
  real, dimension(2*gr_radius+1)   :: stencil, u, S

!!$  real, dimension(3, 3) :: Ck, Lk, Wk, Vtk, Sigmaik
!!$  real, dimension(2,3)  :: Tk
!!$  real, dimension(3, 2) :: Bk
!!$  real, dimension(3)    :: stencilk, uk, Sk
!!$  
!!$  real, dimension(6, 3) :: Zmat
!!$  real, dimension(6)    :: Zvec
!!$  real, dimension(3)    :: vec, ul
!!$  real, dimension(5)    :: un
!!$
!!$  real, dimension(3,3)  :: K_char, L_char
!!$  real, dimension(3)    :: un3, char_stencil

  real, dimension(gr_radius+1, gr_radius+1) :: Ck, Lk, Wk, Vtk, Sigmaik, Ki
  real, dimension(2,gr_radius+1)  :: Tk
  real, dimension(gr_radius+1, 2) :: Bk
  real, dimension(gr_radius+1)    :: stencilk, uk, Sk

  !I think I put this 6 here to be 2*R+1 + 1 equation for the mean function (!?!)
  real, dimension(2*gr_radius+2, gr_radius+1) :: Zmat
  real, dimension(2*gr_radius+2)    :: Zvec
  real, dimension(gr_radius+1)    :: vec, ul, Pk
  real, dimension(2*gr_radius+1)    :: un

  real, dimension(gr_radius+1,gr_radius+1)  :: K_char, L_char
  real, dimension(gr_radius+1)    :: unl, char_stencil, eig_stencil
  
  integer :: LR,i,j,N,M, ROW, COL, R
  real :: small, sigma, sigdel
  !for dgetrf/s

  
  !initialize
  R = gr_radius
  N = 2*R+1
  C = 0.
  T = 0.
  u = 1.

  do i = 1, N
     stencil(i) =  REAL(i - R - 1)       !0.5*REAL(2*(i-R)-1)
     !stencil(i) = (i-R-1)*gr_dx
  end do

  !first thing is to calculate the covariance matrix according to eq. 15
  !since C is symmetric only bother with one side of the diaganol
  do i = 1,N
     do j = 1,N
        C(i,j) = SE(stencil(i), stencil(j))     
     end do
     T(1,i) = SE(stencil(i), -0.5)
     T(2,i) = SE(stetncil(i), 0.5)
  end do
  !now we need to solve the linear eqns for v & Z (see eqs 30-32)
  
  
  call chol(C, N, L)
  call solve_Axb(C, gr_GPv, u, L, N)
  call solve_CZT(C, gr_GPZ, T, L, N)

  !this should only be needed for calculating the conditional probability i.e. I haven't found a need for it yet
  gr_U2_16(1) = K(0.5,0.5)   - dot_product(gr_GPZ(1,:), T(1,:))
  gr_U2_16(2) = K(-0.5,-0.5) - dot_product(gr_GPZ(2,:), T(2,:))
  print *, 'global stencil complete'
  !now lets solve for the gp vars over the local stencils (size = R+1)
  N = R+1
  Ck = 0.
  Tk = 0.
  uk = 1.

  !this is the stencil over the central cell
  !used for eigen decomp

  do m = 1,N
     !make the m-th stencil 
     !these are just the R+1 eno stencils
     !this should extend from -R+m-1 to m-1
     do i = 1, N
        stencil(i) = REAL( i-1 - R + m-1 )
     end do
!!$     print *, 'm/l = ', m
!!$     print *, stencil(1:N)
!!$     print *,
     !now lets calculate the covariance matrix

     do i = 1,N
        do j = 1,N
           Ck(i,j) = SE(stencil(i), stencil(j))
        end do
        Tk(1,i) = SE(stencil(i),-0.5)
        Tk(2,i) = SE(stencil(i), 0.5)
     end do

     !now we do cholesky
     call chol(Ck, N, Lk)
     call solve_Axb(Ck, gr_GPvk(:,m), uk, Lk, N)
     call solve_CZT(Ck, gr_GPZk(:,:,m), Tk, Lk, N)
     !z-vectors for eigen system
     do i = 1, N
        do j = 1, N
           Pk(j) = quad_cross(stencil(j), stencil(i))
        end do
        call solve_Axb(Ck, gr_gp_Zvecs(:,i,m), Pk(:), Lk, N)
     end do

  end do
  print *, 'ENO stencils complete'
  !now we have the Z vectors for all of our stencils
  !lets now compute the linear weights
  !we do this by solving Ax=b, using least squares for the overdetermined system
  !----where x are the linear weights
  !----A is going to be the matrix of Z_k vectors
  !----b is the Z vector for the larger stencil

  !lets first make Z matrix
  !Zmat = 0.
  ROW = 2*R+1
  COL = R+1
  ul = 1.
  un = 1.
  
  do LR = 1, 2
     Zmat = 0.
     print *, 
     do m = 1, COL
        Zmat(m:m+R,m) = gr_GPZk(LR,:, m)
        !this would be for non-zero prior mean
        Zmat(2*R+2,m) = dot_product(gr_GPZk(LR,:,m), ul)
     end do
     
     Zvec(1:ROW) = gr_GPZ(LR,:)
     Zvec(2*R+2)   = dot_product(gr_GPZ(LR,:),un)

     call LSTSQ(ROW+1, COL, Zmat, gr_GP_w(LR,:), Zvec)
!!$     gr_GP_w(1,:) = (/.3, .6, .1 /)
!!$     gr_GP_w(2,:) = (/.1, .6, .3 /)
     
  end do

  !last thing to do is calculate K-inverse for ths smaller stencils
  !we need this in order to compute the marginal-likelyhood later
  N = R+1
  sigma  = sim_sigma
  sigdel = sim_sigdel
  
!!$  sim_sigdel = 1.
!!$  sim_sigma  = sim_sigdel*gr_dx
  do i = 1,N
     do j = 1,N
        Ck(i,j) =  SE(stencil(i),stencil(j))!intg_kernel(stencil(i), stencil(j))
     end do
  end do
  do i = 1, N
     vec = 0.
     vec(i) = 1.
     call solve_Axb(Ck, Ki(:,i), vec, Lk, N)
     gr_gp_Kki = matmul(gr_gp_Zvecs(:,:,1), matmul(Ki(:,:), transpose(gr_gp_Zvecs(:,:,1))))

  end do

  print *, 'linear weights computed'

!!$  sim_sigma  = sigma
!!$  sim_sigdel = sigdel



  !experimental things for gp-char tracing
!!$  if (.not. sim_RK) then
!!$     !we are going to do GP characteristic tracing
!!$
!!$     K_char = 0.
!!$     char_stencil = (/0., -0.5, 0.5  /)
!!$     K_char(1,1) = quad_exact(0.,0.)
!!$     do i = 2,3
!!$        K_char(1,i) = quad_cross(0., char_stencil(i))
!!$        K_char(i,1) = K_char(1,i)
!!$     end do
!!$     do i = 2,3
!!$        do j = 2,3
!!$           K_char(i,j) = K(char_stencil(i), char_stencil(j))
!!$        end do
!!$     end do
!!$
!!$     !now we invert the matrix using the cholesky decomposition
!!$     call chol(K_char, 3, L_char)
!!$
!!$     do i = 1, 3
!!$        unl = 0.
!!$        unl(i) = 1.
!!$        call solve_Axb(K_char, gr_GP_char_Kki(:,i), unl, L_char, 3)
!!$     end do
!!$
!!$
!!$
!!$     gr_GP_char_Kki2 = gr_GP_char_Kki
!!$
!!$     
!!$  end if !.not. sim_RK
  
end subroutine sim_GPinit2

