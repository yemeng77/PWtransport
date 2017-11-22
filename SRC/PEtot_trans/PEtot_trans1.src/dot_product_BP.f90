SUBROUTINE dot_product_BP(A,B,ng_n,U)

  USE fft_data
  USE load_data
  USE DATA
  
  IMPLICIT NONE

  INCLUDE 'mpif.h'
  INCLUDE 'param.escan_real_f90'
  
  COMPLEX*16, DIMENSION(mg_nx,nblock_band_mx) :: A, B, A_tmp, A_tmp2
!  INTEGER, INTENT(in) :: ng_n
  COMPLEX*16, DIMENSION(mx,mx) :: U
  complex*16, allocatable, dimension(:,:) :: S
  integer :: mem_cut,m1_band,mm1_band,m1_step,im


  COMPLEX*16, DIMENSION(nblock_band_mx,nblock_band_mx) :: SS
  COMPLEX*16, PARAMETER :: one=(1.d0,0.d0), zero=(0.d0,0.d0)
  INTEGER :: status(MPI_STATUS_SIZE)
  INTEGER :: ireq(num_group_b)

  INTEGER :: ierr, i, j, irow, icol, icolor_bb
  
!  S = zero
  U = zero
  CALL zgemm('c','n',nblock_band_mx,nblock_band_mx,ng_n,one,A,mg_nx,B,mg_nx,zero,SS,nblock_band_mx)

  DO i=1, nblock_band_mx
     irow = icolor_b*nblock_band_mx + i
     DO j=1, nblock_band_mx
        icol = icolor_b*nblock_band_mx + j
!        S(irow,icol) = SS(i,j)
        U(irow,icol) = SS(i,j)
     ENDDO
  ENDDO


 ! Send
  DO icolor_bb=0, num_group_b-1
     IF (icolor_bb /= icolor_b) THEN
        CALL mpi_isend(A,mg_nx*nblock_band_mx,MPI_DOUBLE_COMPLEX, &
             icolor_bb,ikey_b,MPI_COMM_B2,ireq(icolor_bb+1),ierr)
     ENDIF
  ENDDO

 
  ! Receive
  DO icolor_bb=0, num_group_b-1
     IF (icolor_bb /= icolor_b) THEN

        CALL mpi_recv(A_tmp,mg_nx*nblock_band_mx,MPI_DOUBLE_COMPLEX, &
             icolor_bb,ikey_b,MPI_COMM_B2,status,ierr)
        CALL mpi_wait(ireq(icolor_bb+1),status,ierr)

        CALL zgemm('c','n',nblock_band_mx,nblock_band_mx,ng_n,one,A_tmp,mg_nx,B,mg_nx,zero,SS,nblock_band_mx)
        DO i=1, nblock_band_mx
           irow = icolor_bb*nblock_band_mx + i
           DO j=1, nblock_band_mx
              icol = icolor_b*nblock_band_mx + j
!              S(irow,icol) = SS(i,j)
              U(irow,icol) = SS(i,j)
           ENDDO
        ENDDO

     ENDIF
  ENDDO

!  CALL mpi_allreduce(S,U,mx*mx,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_K1,ierr)

     mem_cut=20.0E+6        ! 20 MB
  if(mx*mx*16.lt.mem_cut) then
  allocate(S(mx,mx))
  CALL mpi_allreduce(U,S,mx*mx,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_K1,ierr)
  U=S

  else      ! to avoid the possible large buffer needed in mpi_allreduce, also to reduce S
   m1_band=mem_cut/(mx*16)
   if(m1_band.lt.1) m1_band=1
   allocate(S(mx,m1_band))
   m1_step=mx/m1_band
   do im=0,m1_step
    mm1_band=m1_band
    if(im.eq.m1_step) mm1_band=mx-m1_step*m1_band
    if(mm1_band.gt.0) then
    CALL mpi_allreduce(U(1,im*m1_band+1),S(1,1),mm1_band*mx,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_K1,ierr)
    U(1:mx,im*m1_band+1:im*m1_band+mm1_band)=S(1:mx,1:mm1_band)
    endif
   enddo
   endif

  deallocate(S)


  
END SUBROUTINE dot_product_BP
