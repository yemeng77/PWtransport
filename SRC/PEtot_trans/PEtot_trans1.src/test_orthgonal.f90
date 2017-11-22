PROGRAM test_orthogonal
  
  IMPLICIT NONE

  INTEGER, PARAMETER :: n_row=20, n_col=10
  DOUBLE COMPLEX, DIMENSION(n_row,n_col) :: A
  
  DOUBLE COMPLEX, DIMENSION(n_row,n_col) :: A_tmp, A_bak
  DOUBLE COMPLEX, DIMENSION(n_col,n_col) :: S, U, R
  DOUBLE PRECISION, DIMENSION(n_col) :: Lamda
  DOUBLE PRECISION, DIMENSION(n_col,n_col) :: Lamda_inv_half
  INTEGER :: i, j, k, ierr
  DOUBLE COMPLEX :: overlap

  DOUBLE PRECISION :: xr, xi, vol
  
  DOUBLE COMPLEX, DIMENSION(3*n_col) :: work
  INTEGER :: lwork, info
  DOUBLE PRECISION, DIMENSION(3*n_col) :: rwork

       
  CHARACTER(len=1) :: uplo, trans
  DOUBLE PRECISION :: alpha, beta

  ! Random matrix A
  DO i=1, n_row
     DO j=1, n_col
        CALL RANDOM_NUMBER(xr)
        CALL RANDOM_NUMBER(xi)
        A(i,j) = DCMPLX(xr,xi)
        WRITE(7,*) i,j,A(i,j)
     ENDDO
  ENDDO

  A_bak = A

  vol=1.d0

  ! Calculate overlap matrix. 
  ! S = A* . A
  DO i=1, n_col
     DO j=i, n_col
        overlap = (0.d0,0.d0)
        DO k=1, n_row
           overlap = overlap + CONJG(A(k,i))*A(k,j)                     
        ENDDO
        overlap = overlap*vol
        S(i,j) = overlap
        S(j,i) = CONJG(overlap)
     ENDDO
  ENDDO

  DO i=1, n_col
     DO j=1, n_col
        WRITE(71,*) i, j, S(i,j)
     ENDDO
  ENDDO
  

  ! Diagonalize overlap matrix.
  ! S = U* . Lamda . U
  lwork=3*n_col
  CALL zheev('V','U',n_col,S,n_col,Lamda,work,lwork,rwork,info)
  IF (info /= 0) THEN
     WRITE(6,*) "Somthing is wrong with diagonalization of S. info", info
     STOP
  ENDIF
  
  
  ! Check if S is positive definite.
  IF (.NOT. ALL(Lamda > 0.d0)) STOP "S is not positive definite."
  
  
  ! Orthogonalize A.
  ! R = U* . Lamda_inv_half
  ! A = A . R
  Lamda_inv_half = 0.d0
  DO k=1, n_col
     Lamda_inv_half(k,k) = 1.d0/DSQRT(Lamda(k))
  ENDDO

  R = MATMUL(S,Lamda_inv_half)
  A_tmp = (0.d0,0.d0)
  DO i=1, n_row
     DO j=1, n_col
        DO k=1, n_col
           A_tmp(i,j) = A_tmp(i,j) + A(i,k)*R(k,j)
        ENDDO
     ENDDO
  ENDDO
  
  A = A_tmp

  DO i=1, n_col
     DO j=1, n_col
        WRITE(81,*) i, j, DOT_PRODUCT(A(:,i),A(:,j))
     ENDDO
  ENDDO
  


  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Calculate overlap matrix.
  ! S = A* . A

  A = A_bak


  uplo = 'l'
  trans = 'c'
  alpha = 1.d0
  beta = 0.d0
  S = CMPLX(0.d0,0.d0)
  CALL zherk(uplo,trans,n_col,n_row,alpha,A,n_row,beta,S,n_col)
  
  U = S*vol
  DO i=1, n_col
     DO j=1, i-1
        U(j,i) = CONJG(U(i,j))
     ENDDO
  ENDDO

  DO i=1, n_col
     DO j=1, n_col
        WRITE(72,*) i, j, U(i,j)
     ENDDO
  ENDDO


  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Cholesky decomposition
  ! A = U**H * U
  uplo = 'U'
  CALL zpotrf(uplo,n_col,U,n_col,info)
  IF (info /= 0) THEN
     WRITE(6,*) "Somthing is wrong with cholesky decomposition of S. info", info
     STOP
  ENDIF
     
  DO i=1, n_col
     DO j=1, n_col
        WRITE(100,*) i, j, U(i,j)
     ENDDO
  ENDDO

  DO i=1, n_col
     U(i,1:i-1) = (0.d0,0.d0)
!     DO j=1, n_col
!        IF (i>j) U(i,j) = (0.d0,0.d0)
!     ENDDO
  ENDDO

  R = U

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Calculate the inverse of cholesky factor.
  ! U = U**(-1)
  CALL ztrtri(uplo,'N',n_col,U,n_col,info)
  IF (info /= 0) THEN
     WRITE(6,*) "Somthing is wrong with inversion of cholesky factor of U. info", info
     STOP
  ENDIF


  S = MATMUL(R,U)

  DO i=1, n_col
     DO j=1, n_col
        WRITE(101,*) i, j, S(i,j)
     ENDDO
  ENDDO



  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Orthogonalize A.
  ! A = A * U**(-1)
  A = MATMUL(A,U)

  DO i=1, n_col
     DO j=1, n_col
        WRITE(82,*) i, j, DOT_PRODUCT(A(:,i),A(:,j))
     ENDDO
  ENDDO





END PROGRAM test_orthogonal
