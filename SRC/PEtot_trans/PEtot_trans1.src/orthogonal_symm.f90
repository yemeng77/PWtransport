SUBROUTINE orthogonal(A,n_row,n_col)
  ! On entry, the matrix A (of size n_row X n_col) contains column vectors a1, a2, a3, ... an_col, that are not orthogonal each other.
  ! On exit, A contains column vectors that are orthognormal.

  USE fft_data
  USE load_data
  USE DATA
  
  IMPLICIT NONE

  INCLUDE 'mpif.h'
  INCLUDE 'param.escan_real_f90'
  
  INTEGER, INTENT(in) :: n_row, n_col
  DOUBLE COMPLEX, DIMENSION(n_row,n_col) :: A
  
  DOUBLE COMPLEX, DIMENSION(n_row,n_col) :: A_tmp
  DOUBLE COMPLEX, DIMENSION(n_col,n_col) :: S, U, R
  DOUBLE PRECISION, DIMENSION(n_col) :: Lamda
  DOUBLE PRECISION, DIMENSION(n_col,n_col) :: Lamda_inv_half
  INTEGER :: i, j, k, ierr
  DOUBLE COMPLEX :: overlap
  

  CHARACTER(len=1) :: uplo, trans
  DOUBLE PRECISION :: alpha, beta

  DOUBLE COMPLEX, DIMENSION(3*n_col) :: work
  INTEGER :: lwork, info
  DOUBLE PRECISION, DIMENSION(3*n_col) :: rwork

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Calculate overlap matrix.
  ! S = A* . A

!!$  ! Clumsy way.
!!$  DO i=1, n_col
!!$     DO j=i, n_col
!!$        overlap = DOT_PRODUCT(A(:,i),A(:,j))
!!$        overlap = overlap*vol
!!$        CALL global_sumc(overlap)
!!$        S(i,j) = overlap
!!$        S(j,i) = CONJG(overlap)
!!$     ENDDO
!!$  ENDDO

  ! BLAS3
  uplo = 'l'
  trans = 'c'
  alpha = 1.d0
  beta = 0.d0
  U = CMPLX(0.d0,0.d0)
  CALL zherk(uplo,trans,n_col,n_row,alpha,A,n_row,beta,U,n_col)
  
  CALL mpi_allreduce(U,S,n_col*n_col,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_K,ierr)

  S = S*vol
  DO i=1, n_col
     DO j=1, i-1
        S(j,i) = CONJG(S(i,j))
     ENDDO
  ENDDO



  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
  

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Orthogonalize A.
  ! R = U* . Lamda_inv_half
  ! A = A . R
  Lamda_inv_half = 0.d0
  DO k=1, n_col
     Lamda_inv_half(k,k) = 1.d0/DSQRT(Lamda(k))
  ENDDO

  R = MATMUL(S,Lamda_inv_half)
  
  A = MATMUL(A,R)
  


END SUBROUTINE orthogonal
