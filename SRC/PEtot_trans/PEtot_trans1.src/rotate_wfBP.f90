SUBROUTINE rotate_wfBP(A,U,ng_nt,mg_nxt)

  USE fft_data
  USE load_data
  USE DATA
  
  IMPLICIT NONE

  INCLUDE 'mpif.h'
  INCLUDE 'param.escan_real_f90'
  INTEGER, INTENT(in) :: ng_nt,mg_nxt
  
  COMPLEX*16, DIMENSION(mg_nxt,nblock_band_mx) :: A, A_tmp, A_tmp2, AA
  COMPLEX*16, DIMENSION(mx,mx), INTENT(in) :: U

  COMPLEX*16, DIMENSION(nblock_band_mx,nblock_band_mx) :: S, SS
  COMPLEX*16, PARAMETER :: one=(1.d0,0.d0), zero=(0.d0,0.d0)
  INTEGER :: status(MPI_STATUS_SIZE)
  INTEGER :: ireq(num_group_b)

  INTEGER :: ierr, i, j, irow, icol, icolor_bb
  
  AA = zero

  ! Send
  DO icolor_bb=0, num_group_b-1
     IF (icolor_bb /= icolor_b) THEN
        CALL mpi_isend(A,mg_nxt*nblock_band_mx,MPI_DOUBLE_COMPLEX, &
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
  CALL zgemm('n','n',ng_nt,nblock_band_mx,nblock_band_mx,one,A,mg_nxt,SS,nblock_band_mx,zero,A_tmp2,mg_nxt)
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

        CALL mpi_recv(A_tmp,mg_nxt*nblock_band_mx,MPI_DOUBLE_COMPLEX, &
             icolor_bb,ikey_b,MPI_COMM_B2,status,ierr)
!cccc        CALL mpi_wait(ireq(icolor_bb+1),status,ierr)

        CALL zgemm('n','n',ng_nt,nblock_band_mx,nblock_band_mx,one,A_tmp,mg_nxt,SS,nblock_band_mx,zero,A_tmp2,mg_nxt)
        AA = AA + A_tmp2

     ENDIF
  ENDDO

  do icolor_bb=0, num_group_b-1
    if(icolor_bb.ne.icolor_b) then
        CALL mpi_wait(ireq(icolor_bb+1),status,ierr)
    endif
  enddo

  A = AA
  
END SUBROUTINE rotate_wfBP
