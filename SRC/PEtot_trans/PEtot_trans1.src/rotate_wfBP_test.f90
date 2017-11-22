SUBROUTINE rotate_wfBP_test(A,U,ng_n,sum_test)

  USE fft_data
  USE load_data
  USE DATA
  
  IMPLICIT NONE

  INCLUDE 'mpif.h'
  INCLUDE 'param.escan_real_f90'
  
  COMPLEX*16, DIMENSION(mg_nx,nblock_band_mx) :: A, A_tmp, A_tmp2, AA
  COMPLEX*16, DIMENSION(mx,mx), INTENT(in) :: U
!  INTEGER, INTENT(in) :: ng_n

  COMPLEX*16, DIMENSION(nblock_band_mx,nblock_band_mx) :: S, SS
  COMPLEX*16, PARAMETER :: one=(1.d0,0.d0), zero=(0.d0,0.d0)
  INTEGER :: status(MPI_STATUS_SIZE)
  INTEGER :: ireq(num_group_b)
  COMPLEX*16 :: sum_test

  INTEGER :: ierr, i, j, irow, icol, icolor_bb

  write(6,*) "test_rot1 ", inode_tot,A(1,1),A(1,2),A(2,1)
  
  AA = zero

  ! Send
  DO icolor_bb=0, num_group_b-1
     IF (icolor_bb /= icolor_b) THEN
        CALL mpi_isend(A,mg_nx*nblock_band_mx,MPI_DOUBLE_COMPLEX, &
             icolor_bb,ikey_b,MPI_COMM_B2,ireq(icolor_bb+1),ierr)
     ENDIF
  ENDDO

  write(6,*) "test_rot2 ", inode_tot,sum_test

  DO i=1, nblock_band_mx
     irow = icolor_b*nblock_band_mx + i
     DO j=1, nblock_band_mx
        icol = icolor_b*nblock_band_mx + j
        SS(i,j) =  U(irow,icol)
     ENDDO
  ENDDO

  write(6,*) "test_rot3 ", inode_tot,A(1,1),A(1,2),A(2,1),ng_n,nblock_band_mx,mg_nx
  call mpi_barrier(mpi_comm_world,ierr)

  CALL zgemm('n','n',ng_n,nblock_band_mx,nblock_band_mx,one,A,mg_nx,SS,nblock_band_mx,zero,A_tmp2,mg_nx)

  call mpi_barrier(mpi_comm_world,ierr)
  write(6,*) "test_rot4 ", inode_tot,ng_n,nblock_band_mx,mg_nx
  write(6,*) "test_rot4.1 ", inode_tot,sum_test

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

  write(6,*) "test_rot5 ", inode_tot,sum_test

        CALL mpi_recv(A_tmp,mg_nx*nblock_band_mx,MPI_DOUBLE_COMPLEX, &
             icolor_bb,ikey_b,MPI_COMM_B2,status,ierr)
        CALL mpi_wait(ireq(icolor_bb+1),status,ierr)

  write(6,*) "test_rot6 ", inode_tot,sum_test

        CALL zgemm('n','n',ng_n,nblock_band_mx,nblock_band_mx,one,A_tmp,mg_nx,SS,nblock_band_mx,zero,A_tmp2,mg_nx)
        AA = AA + A_tmp2

  write(6,*) "test_rot7 ", inode_tot,sum_test

     ENDIF
  ENDDO

  write(6,*) "test_rot8 ", inode_tot,sum_test,mg_nx,nblock_band_mx

  A = AA

  write(6,*) "test_rot9 ", inode_tot,sum_test
  
END SUBROUTINE rotate_wfBP_test
