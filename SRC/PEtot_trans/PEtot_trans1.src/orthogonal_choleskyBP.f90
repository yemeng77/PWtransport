SUBROUTINE orthogonal_choleskyBP(A,n_row,n_col,ng_n,U)
  ! On entry, the matrix A (of size n_row X n_col) contains column vectors a1, a2, a3, ... an_col, that are not orthogonal each other.
  ! however, only 1:ng_n of each column 1:n_col is used as the vector
  ! On exit, A contains column vectors that are orthognormal.

  ! Decompose A = (Psi**H * Psi) = U**H * U.
  ! The orthogonal vectors Phi = Psi * U**(-1).
  ! This is equivalent to Gram-Schmidt orthogonalization progressing from the first column vector of Psi.

  USE fft_data
  USE load_data
  USE DATA
  
  IMPLICIT NONE

  INCLUDE 'mpif.h'
  INCLUDE 'param.escan_real_f90'
  
  INTEGER, INTENT(in) :: n_row, n_col
  COMPLEX*16, DIMENSION(n_row,n_col) :: A, A_tmp, A_tmp2, AA
  COMPLEX*16, ALLOCATABLE,DIMENSION(:,:) :: CC
  
!!$  COMPLEX*16, DIMENSION(mx,mx) :: S, S_t, U
  COMPLEX*16, DIMENSION(mx,mx) :: U
  complex*16, allocatable, dimension (:,:) :: S
  COMPLEX*16, DIMENSION(n_col,n_col) :: SS, SS_tmp
  INTEGER :: i, j, k, ierr, info, icolor_bb, icol, irow
  
  CHARACTER(len=1) :: uplo, trans
  COMPLEX*16, PARAMETER :: one=(1.d0,0.d0), minus_one=(-1.d0,0.d0), zero=(0.d0,0.d0)
  INTEGER :: status(MPI_STATUS_SIZE)
  INTEGER :: ireq(num_group_b)

  complex*16 :: sss
  integer :: ii, mk
  integer :: mem_cut,m1_band,mm1_band,m1_step,im

!!$  COMPLEX*16, DIMENSION(mg_nx,mx) :: ug_n1, ug_n2


!!$! This is old way without band-parellelization.
!!$  ug_n1 = (0.d0,0.d0)
!!$  ug_n1(:,band_dis(1):band_dis(2)) = A(:,:)
!!$  CALL mpi_allreduce(ug_n1,ug_n2,mg_nx*mx,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_B2,ierr)



!!$  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$  ! Calculate overlap matrix.
!!$  ! S = A* . A
!!$  uplo = 'u'
!!$  trans = 'c'
!!$
!!$  ! important
!!$  S = zero  
!!$  CALL zherk(uplo,trans,mx,ng_n,vol,ug_n2,mg_nx,zero,S,mx)
!!$  CALL mpi_allreduce(S,U,mx*mx,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_K,ierr)
!!$
!!$! the lower triangle of S and U is zero



!!    allocate(S(mx,mx))


  ! Send
  DO icolor_bb=0, num_group_b-1
     IF (icolor_bb /= icolor_b) THEN
        CALL mpi_isend(A,n_row*n_col,MPI_DOUBLE_COMPLEX, &
             icolor_bb,ikey_b,MPI_COMM_B2,ireq(icolor_bb+1),ierr)
     ENDIF
  ENDDO

  U = zero
  ! Receive
  DO icolor_bb=0, num_group_b-1
     IF (icolor_bb == icolor_b) THEN

        CALL zgemm('c','n',n_col,n_col,ng_n,one,A,n_row,A,n_row,zero,SS,n_col)
        DO i=1, nblock_band_mx
           irow = icolor_bb*nblock_band_mx + i
           DO j=1, nblock_band_mx
              icol = icolor_b*nblock_band_mx + j
!!              S(irow,icol) = SS(i,j)*vol
              U(irow,icol) = SS(i,j)*vol
           ENDDO
        ENDDO
        
     ELSE

        CALL mpi_recv(A_tmp,n_row*n_col,MPI_DOUBLE_COMPLEX, &
             icolor_bb,ikey_b,MPI_COMM_B2,status,ierr)
        CALL mpi_wait(ireq(icolor_bb+1),status,ierr)

        CALL zgemm('c','n',n_col,n_col,ng_n,one,A_tmp,n_row,A,n_row,zero,SS,n_col)
        DO i=1, nblock_band_mx
           irow = icolor_bb*nblock_band_mx + i
           DO j=1, nblock_band_mx
              icol = icolor_b*nblock_band_mx + j
!!              S(irow,icol) = SS(i,j)*vol
              U(irow,icol) = SS(i,j)*vol
           ENDDO
        ENDDO

     ENDIF
  ENDDO

     mem_cut=20.0E+6        ! 20 MB 
  if(mx*mx*16.lt.mem_cut) then
!!  CALL mpi_allreduce(S,U,mx*mx,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_K1,ierr)
  allocate(S(mx,mx))
  CALL mpi_allreduce(U,S,mx*mx,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_K1,ierr)
  U=S

  else      ! to avoid the possible large buffer needed in mpi_allreduce
   m1_band=mem_cut/(mx*16)
   if(m1_band.lt.1) m1_band=1
   allocate(S(mx,m1_band))
   m1_step=mx/m1_band
   do im=0,m1_step
    mm1_band=m1_band
    if(im.eq.m1_step) mm1_band=mx-m1_step*m1_band
    if(mm1_band.gt.0) then
!!  CALL mpi_allreduce(S(1,im*m1_band+1),U(1,im*m1_band+1),mm1_band*mx,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_K1,ierr)
    CALL mpi_allreduce(U(1,im*m1_band+1),S(1,1),mm1_band*mx,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_K1,ierr)
    U(1:mx,im*m1_band+1:im*m1_band+mm1_band)=S(1:mx,1:mm1_band)
    endif
   enddo
   endif

  deallocate(S)


  DO i=2, mx
     U(i,1:i-1) = zero
  ENDDO

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Cholesky decomposition
  ! A = U**H * U
  uplo = 'U'
  
  CALL zpotrf(uplo,mx,U,mx,info)

  IF (info /= 0) THEN
     WRITE(6,*) "Something is wrong with choleskyBP decomposition of S. info", info
     STOP
  ENDIF

! the lower triangle of U is still zero, never used


     

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Calculate the inverse of cholesky factor.
  ! U = U**(-1)


  CALL ztrtri(uplo,'N',mx,U,mx,info)

  IF (info /= 0) THEN
     WRITE(6,*) "Somthing is wrong with inversion of cholesky factor of U. info", info
     STOP
  ENDIF


  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Orthogonalize A.
  ! A = A * U**(-1)
  ! A = MATMUL(A,U)  
  ! unfortunately, half is wasted, because half of U is zero. 

!!$  CALL zgemm('N','N',ng_n,mx,mx,one,ug_n2,mg_nx,U,mx,zero,ug_n1,mg_nx)
!!$   A(:,:)=ug_n1(:,band_dis(1):band_dis(2))


  AA = zero

  ! Send
  DO icolor_bb=0, num_group_b-1
     IF (icolor_bb /= icolor_b) THEN
        CALL mpi_isend(A,n_row*n_col,MPI_DOUBLE_COMPLEX, &
             icolor_bb,ikey_b,MPI_COMM_B2,ireq(icolor_bb+1),ierr)
     ENDIF
  ENDDO


 
  DO i=1, nblock_band_mx
     irow = icolor_b*nblock_band_mx + i
     DO j=1, nblock_band_mx
        icol = icolor_b*nblock_band_mx + j
        SS(i,j) =  U(irow,icol)
     ENDDO
  ENDDO
  CALL zgemm('n','n',ng_n,n_col,n_col,one,A,n_row,SS,n_col,zero,A_tmp2,n_row)
  AA = AA + A_tmp2


  ! Receive
  DO icolor_bb=0, num_group_b-1
     IF (icolor_bb /= icolor_b) THEN

        DO i=1, nblock_band_mx
           irow = icolor_bb*nblock_band_mx + i
           DO j=1, nblock_band_mx
              icol = icolor_b*nblock_band_mx + j
              SS(i,j) =  U(irow,icol)
           ENDDO
        ENDDO

        CALL mpi_recv(A_tmp,n_row*n_col,MPI_DOUBLE_COMPLEX, &
             icolor_bb,ikey_b,MPI_COMM_B2,status,ierr)
        CALL mpi_wait(ireq(icolor_bb+1),status,ierr)
       
        CALL zgemm('n','n',ng_n,n_col,n_col,one,A_tmp,n_row,SS,n_col,zero,A_tmp2,n_row)
        Aa = AA + A_tmp2

     ENDIF
  ENDDO

  A = AA
  

 END SUBROUTINE orthogonal_choleskyBP
