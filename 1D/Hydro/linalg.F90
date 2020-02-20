module linalg

contains

  subroutine chol(A, N, L)
    !Calculate the cholesky decomposition (L) of matrix A that is NxN
    implicit none
    integer             , intent(IN)  :: N
    real, dimension(N,N), intent(IN)  :: A

    real, dimension(N,N), intent(OUT) :: L

    integer :: i, j, k

    !initialize L
    L = 0.
    L(1,1) = SQRT(A(1,1))

    !do decomposition
    do i = 2,N
       do j = 1,i-1
          !do entries below diaganol
          L(i,j) = A(i,j)
          do k = 1, j-1
             L(i,j) = L(i,j) - L(i,k)*L(j,k)
          end do
          L(i,j) = L(i,j)/L(j,j)
       end do
       !after finishing below diaganol, take care of diaganol elements
       L(i,i) = A(i,i)
       do k = 1, i-1
          L(i,i) = L(i,i) - L(i,k)**2
       end do
       L(i,i) = SQRT(L(i,i))
    end do

  end subroutine chol

  subroutine LSTSQ(M, N, A, x, b)
    !compute the least squares solution to the overdetermined system (ONLY!) Ax=b
    !using the QR factorization of A
    !done using Householder reflections
    implicit none
    integer, intent(IN) :: M, N
    real, dimension(M,N), intent(INOUT) :: A
    real, dimension(M)  , intent(INOUT) :: b

    real, dimension(N)  , intent(INOUT) :: x

    real, dimension(M) :: QTb
    real, dimension(M,1) :: v
    real, dimension(N,N) :: R

    integer :: i, j
    real    :: v2, gam

    R = 0.
    QTb = 0.
    V = 0.
    !we're going to compute the relevant items from the QR factorization
    !R is going to be stored in A
    !we are also going to calculate Q^T*b for the least squares solve
    do i = 1, N
       v = 0.
       v(i:M,1) = A(i:M,i)
       !vi = xi + sign(xi)*x^2
       v(i,1  ) = A(i,i) + SIGN( SQRT(dot_product( A(i:M,i),A(i:M,i) )), A(i,i) )

       v2 = SQRT(dot_product(v(i:M,1),v(i:M,1)))
       v(:,1) = v(:,1)/v2

       
       A(i:M,i:N) = A(i:M,i:N) - 2.*MATMUL( v(i:M,:), MATMUL( TRANSPOSE(v(i:M,:)), A(i:M,i:N) ) )

       
       gam = -2.*dot_product(v(i:M,1), b(i:M))
       b(i:M) = b(i:M) + gam*v(i:M,1)
       
    end do

    !now we just to do a back-substitution to solve Rx=Q^T*b
    do i = N, 1, -1
       x(i) = b(i)
       do j = i+1, N
          x(i) = x(i) - A(i,j)*x(j)
       end do
       x(i) = x(i)/A(i,i)
    end do


  end subroutine LSTSQ

  subroutine solve_Axb(A, x, b, L, N)
    !solve matrix equation Ax = b for vector x, A is NxN given the cholesky decomposition L of A
    implicit none
    integer, intent(IN) :: N
    real, dimension(N  ), intent(IN) :: b
    real, dimension(N,N), intent(IN) :: A, L

    real, dimension(N) , intent(OUT) :: x

    real, dimension(N)   :: y
    integer :: i, j

    !initialize
    x = 0.
    y = 0.
    
    !first we need to solve Ly=b using forward sub
    y(1) = b(1)/L(1,1)
    do i = 2,N
       y(i) = b(i)
       do j = 1, i-1
          y(i) = y(i) - L(i,j)*y(j)
       end do
       y(i) = y(i)/L(i,i)
    end do

    !now solve L*x = y
    x(N) = y(N)/L(N,N)
    do i = N-1, 1, -1
       x(i) = y(i)
       do j =i+1, N
          x(i) = x(i) - L(j,i)*x(j)
       end do
       x(i) = x(i)/L(i,i)
    end do

  end subroutine solve_Axb

  subroutine solve_CZT(C, Z, T, L, N)
    !subroutine to solve the matrix equation C.Z* = T* for the matrix Z
    !note that the matrix equation is for the transpose of Z and T
    !C shoul be NxN, while Z and T are 2xN
    implicit none
    integer,                 intent(IN ) :: N
    real   , dimension(N,N), intent(IN ) :: C, L
    real   , dimension(2,N), intent(IN ) :: T

    real, dimension(2,N), intent(OUT) :: Z

    !solve for the first column of Z*
    call solve_Axb(C, Z(1, 1:N), T(1, 1:N), L, N)
    !solve fore the second column of Z*
    call solve_Axb(C, Z(2, 1:N), T(2, 1:N), L, N)
  end subroutine solve_CZT

  function analytic_weights2(K, Ks) result(z)
    implicit none
    
    real, dimension(2,2), intent(IN) :: K
    real, dimension(2  ), intent(IN) :: Ks

    real, dimension(2) :: z

    z(1)=(Ks(1)*K(1,1) - Ks(2)*K(1,2))/(K(1,1)**2 - K(1,2)**2)

    z(2)=(Ks(2)*K(1,1) - Ks(1)*K(1,2))/(K(1,1)**2 - K(1,2)**2)

    return
  end function analytic_weights2


  function analytic_weights3(K, Ks) result(z)
    implicit none

    real, dimension(3,3), intent(IN) :: K
    real, dimension(3  ), intent(IN) :: Ks

    real, dimension(3) :: z

    z(1)=(Ks(3)*(-(K(1,1)*K(1,3)) + K(1,2)*K(2,3)) + Ks(2)*(-(K(1,1)*K(1,&
         2)) + K(1,3)*K(2,3)) + Ks(1)*(K(1,1)**2 - K(2,3)**2))/(K(1,1)**3 + 2*&
         K(1,2)*K(1,3)*K(2,3) - K(1,1)*(K(1,2)**2 + K(1,3)**2 + K(2,3)**2))

    z(2)=(Ks(2)*(K(1,1)**2 - K(1,3)**2) + Ks(3)*(K(1,2)*K(1,3) - K(1,1)*K&
         (2,3)) + Ks(1)*(-(K(1,1)*K(1,2)) + K(1,3)*K(2,3)))/(K(1,1)**3 + 2*K(1&
         ,2)*K(1,3)*K(2,3) - K(1,1)*(K(1,2)**2 + K(1,3)**2 + K(2,3)**2))

    z(3)=(Ks(3)*(K(1,1)**2 - K(1,2)**2) + Ks(2)*(K(1,2)*K(1,3) - K(1,1)*K&
         (2,3)) + Ks(1)*(-(K(1,1)*K(1,3)) + K(1,2)*K(2,3)))/(K(1,1)**3 + 2*K(1&
         ,2)*K(1,3)*K(2,3) - K(1,1)*(K(1,2)**2 + K(1,3)**2 + K(2,3)**2))

    
    return

  end function analytic_weights3

end module linalg
