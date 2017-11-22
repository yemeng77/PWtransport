SUBROUTINE orthogonal_projectionBP2(A,row_a,col_a,B,SB,row_b,col_b,n_mul)
  ! Orthogonalize column vectors a_i in matrice A to column vectors b_j in matrice B for all j .le. i.
  ! On entry, the matrice A contains complex*16 column vectors that are not orthogonal to B.
  ! On exit, the column vectors a_i are orthogonal to b_j (j .le. i).
  ! The vectors a_i are not normalized.
  !
  ! Things to test.
  ! 1. Ignoring of projections to higher wavefunctions is necessary? 
  !....... Makes little difference.
  !....... See all_band/GaAs64at_3.6_scf (Using the projection to all states.)
  ! 2. If so, we waste half of the computation by using level 3 BLAS.
  !    We might want to use level 2 BLAS.

  USE fft_data
  USE load_data
  USE DATA
  
  IMPLICIT NONE

  INCLUDE 'mpif.h'
  INCLUDE 'param.escan_real_f90'
  
  INTEGER, INTENT(in) :: row_a, col_a, row_b, col_b, n_mul
  COMPLEX*16, DIMENSION(row_a,col_a) :: A
  COMPLEX*16, DIMENSION(row_b,col_b), INTENT(in) :: B,SB
  COMPLEX*16, DIMENSION(row_b,col_b) :: B_tmp, B_tmp2
  
  COMPLEX*16, DIMENSION(mx,mx) :: S_t, SS_t
  COMPLEX*16, DIMENSION(col_b,col_a) :: S, SS
  COMPLEX*16, PARAMETER :: one=(1.d0,0.d0), minus_one=(-1.d0,0.d0), zero=(0.d0,0.d0)
  INTEGER :: status(MPI_STATUS_SIZE)
  INTEGER :: ireq(num_group_b)

  INTEGER :: ierr, i, j, irow, icol, icolor_bb

!!$  COMPLEX*16, DIMENSION(mg_nx,mx) :: ug_n1, ug_n2, p1, p2



!!$! Sanity check
!!$  A = B
!!$! Should lead to identity matrix.


  ! Send
  DO icolor_bb=0, num_group_b-1
     IF (icolor_bb /= icolor_b) THEN
        CALL mpi_isend(SB,row_b*col_b,MPI_DOUBLE_COMPLEX, &
             icolor_bb,ikey_b,MPI_COMM_B2,ireq(icolor_bb+1),ierr)
     ENDIF
  ENDDO

  S_t = zero
 
  ! Receive
  DO icolor_bb=0, num_group_b-1
     IF (icolor_bb == icolor_b) THEN

        CALL zgemm('c','n',col_b,col_a,n_mul,one,SB,row_b,A,row_a,zero,S,col_b)        
        DO i=1, nblock_band_mx
           irow = icolor_bb*nblock_band_mx + i
           DO j=1, nblock_band_mx
              icol = icolor_b*nblock_band_mx + j
              S_t(irow,icol) = S(i,j)*vol
           ENDDO
        ENDDO

     ELSE
        
        CALL mpi_recv(B_tmp,row_b*col_b,MPI_DOUBLE_COMPLEX, &
             icolor_bb,ikey_b,MPI_COMM_B2,status,ierr)
        CALL mpi_wait(ireq(icolor_bb+1),status,ierr)

        !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! Calculate overlap matrix.
        ! S = B**H * A
        ! S_ij is the projection of vector a_i to vector b_j.

        CALL zgemm('c','n',col_b,col_a,n_mul,one,B_tmp,row_b,A,row_a,zero,S,col_b)
        DO i=1, nblock_band_mx
           irow = icolor_bb*nblock_band_mx + i
           DO j=1, nblock_band_mx
              icol = icolor_b*nblock_band_mx + j
              S_t(irow,icol) = S(i,j)*vol
           ENDDO
        ENDDO
       
     ENDIF
  ENDDO

  CALL mpi_allreduce(S_t,SS_t,mx*mx,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_K1,ierr)

!cccccc if this is used, P(i) is only orthogonalized to ug(j) for j<=i, otherwise it is orth. to all ug
  DO icol=1, mx
     DO irow=icol+1, mx
        SS_t(irow,icol) = zero
     ENDDO
  ENDDO
!ccccccccccccccccccccccccccc



!!$
!!$
!!$  IF (inode_k==1) THEN
!!$     DO i=1, mx
!!$        WRITE(icolor_k+200,*) SS_t(i,:)
!!$     ENDDO
!!$  ENDIF
!!$  CALL system_flush(icolor_k+200)





!!$! Sanity check
!!$! This is old way without band-parellelization.
!!$  ug_n1 = (0.d0,0.d0)
!!$  ug_n1(:,band_dis(1):band_dis(2)) = B(:,:)
!!$  CALL mpi_allreduce(ug_n1,ug_n2,mg_nx*mx,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_B2,ierr) 
!!$
!!$  p1 = (0.d0,0.d0)
!!$  p1(:,band_dis(1):band_dis(2)) = A(:,:)
!!$  CALL mpi_allreduce(p1,p2,mg_nx*mx,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_B2,ierr)   
!!$
!!$  CALL zgemm('c','n',mx,mx,n_mul,one,ug_n2,mg_nx,p2,mg_nx,zero,SS_t,mx)
!!$  CALL mpi_allreduce(SS_t,S_t,mx*mx,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_K,ierr)
!!$
!!$  DO icol=1, mx
!!$     DO irow=1, icol
!!$        SS_t(irow,icol) = S_t(irow,icol)*vol
!!$     ENDDO
!!$
!!$     DO irow=icol+1, mx
!!$        SS_t(irow,icol) = zero
!!$     ENDDO
!!$  ENDDO
!!$ 
!!$  CALL zgemm('n','n',n_mul,mx,mx,minus_one,ug_n2,mg_nx,SS_t,mx,one,p2,mg_nx)
!!$  
!!$  A(:,:) = p2(:,band_dis(1):band_dis(2))
!!$! end of sanity check.





  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Orthogonalize A.
  ! A = A - B * S

  ! Send
  DO icolor_bb=0, num_group_b-1
     IF (icolor_bb /= icolor_b) THEN
        CALL mpi_isend(B,row_b*col_b,MPI_DOUBLE_COMPLEX, &
             icolor_bb,ikey_b,MPI_COMM_B2,ireq(icolor_bb+1),ierr)
     ENDIF
  ENDDO

  ! Receive
  DO icolor_bb=0, num_group_b-1

     DO i=1, nblock_band_mx
        irow = icolor_bb*nblock_band_mx + i
        DO j=1, nblock_band_mx
           icol = icolor_b*nblock_band_mx + j
           S(i,j) =  SS_t(irow,icol)
        ENDDO
     ENDDO

     ! B_tmp = - B * S
     IF (icolor_bb == icolor_b) THEN

        CALL zgemm('n','n',n_mul,col_a,col_b,minus_one,B,row_b,S,col_b,zero,B_tmp2,row_a)
        A(1:n_mul,:) = A(1:n_mul,:) + B_tmp2(1:n_mul,:)

     ELSE

        CALL mpi_recv(B_tmp,row_b*col_b,MPI_DOUBLE_COMPLEX, &
             icolor_bb,ikey_b,MPI_COMM_B2,status,ierr)
        CALL mpi_wait(ireq(icolor_bb+1),status,ierr)
       
        CALL zgemm('n','n',n_mul,col_a,col_b,minus_one,B_tmp,row_b,S,col_b,zero,B_tmp2,row_a)
        A(1:n_mul,:) = A(1:n_mul,:) + B_tmp2(1:n_mul,:)

     ENDIF
  ENDDO


  
END SUBROUTINE orthogonal_projectionBP2
