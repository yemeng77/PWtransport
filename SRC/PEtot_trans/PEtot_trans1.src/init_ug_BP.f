      subroutine init_ug_BP(AL,iwg_in,workr_n,kpt,iranm)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************

****************************************
cccccc this is just for one kpt. 
****************************************
****  It stores the wavefunction in G space. 
******************************************
      use fft_data
      use load_data
      use data
      implicit double precision (a-h,o-z)

      include 'mpif.h'
      include 'param.escan_real'

      integer status(MPI_STATUS_SIZE)
******************************************
       real*8 AL(3,3),ALt(3,3)
***********************************************
       complex*16 workr_n(mr_n)


*************************************************
****  input ug from file  and check the consistency
*************************************************

       if(iwg_in.eq.1) then
       call read_ug()
       endif
************************************************
*** generate the initial wavefunction ug from random if iwg_in.eq.0
************************************************
      if(iwg_in.eq.0) then
         ng_n = ngtotnod(inode,kpt)
	do m=1,nblock_band_mx
           do ig=1, ng_n
              x1 = ran1(iranm)
              x2 = ran1(iranm)
              ug_n_bp(ig,m) = dcmplx(x1-0.5d0,x2-0.5d0)
           enddo
           do ig=ng_n+1, mg_nx
              ug_n_BP(ig,m) = dcmplx(0.d0,0.d0)
           enddo
	enddo
      endif
*************************************************
**** end generate the initial wavefunction from random
*************************************************
ccccccc it is not necessary to have the orthogonalization here
c        do m=1,mx
c        call orth_comp(ug_n(1,m),ug_n,m-1,1,kpt,
c     &   sug_n,swg,1,0,fnorm)     ! ipsp_all always be 1 here, since sug_n is not known yet
c	enddo
*************************************************
      return
      contains
*************************************************

       subroutine read_ug()

       complex*16,allocatable,dimension(:) :: ugtemp

        call mpi_barrier(MPI_COMM_K,ierr)
cccccc  file (11) should have been opened before the do loop for kpt=1,nkpt

ccccccccccccc check the header.
ccccc Note, here icolor_k.eq.0 within init_ug, but kpt runs for all kpts. 

        if(inode_k==nnodes_k.and.kpt.eq.1) then     ! only do this for kpt.eq.1

         read(11)n1t,n2t,n3t,mxt
         read(11)Ecutt
         read(11)ALt
         read(11)nnodes_o

         if(n1t.ne.n1.or.n2t.ne.n2.or.n3t.ne.n3) then
         write(6,*) "n1t,n2t,n3t changed, stop", n1t,n2t,n3t
         call mpi_abort(MPI_COMM_WORLD,ierr)
         endif

         if(Ecut.ne.Ecutt.or.mxt.ne.mx) then
         write(6,*) "Ecutt,mxt changed, stop", Ecutt,mxt
         call mpi_abort(MPI_COMM_WORLD,ierr)
         endif

         diff=0.d0
         do i=1,3
            do j=1,3
               diff=diff+dabs(AL(i,j)-ALt(i,j))
             enddo
         enddo

         if(diff.gt.0.001d0) then
            write(6,*) "the AL.ne.ALt, stop"
         call mpi_abort(MPI_COMM_WORLD,ierr)
         endif

         if(nnodes_o.ne.nnodes) then
           write(6,*) "nnodes_o changed, stop", nnodes_o,nnodes
         call mpi_abort(MPI_COMM_WORLD,ierr)
         endif

      endif     ! check the header, make sure it is the same
ccccccccccccccccccccccccccccccccccccccccc
         call mpi_barrier(MPI_COMM_K,ierr)

         if(inode_k.eq.nnodes_k) then
            allocate(ugtemp(mg_nx))

            do inode_tmp=1,nnodes_b
            do igroup_tmp=1,num_group_b

               do iwavefun=1,nblock_band_mx
                  read(11)ugtemp    ! assuming the mg_nx is the same as before
                   ug_n_BP(:,iwavefun)=ugtemp
               end do

          idest=(igroup_tmp-1)*nnodes_b+inode_tmp-1
          if(idest.ne.nnodes_k-1) then
          call mpi_send(ug_n_BP,mg_nx*nblock_band_mx,MPI_DOUBLE_COMPLEX,
     &      idest,102,MPI_COMM_K1,ierr)
          endif
            end do
            end do
c     Now read in the data for node nnodes
            deallocate(ugtemp)

         else   ! for other nodes

         call mpi_recv(ug_n_BP,mg_nx*nblock_band_mx,MPI_DOUBLE_COMPLEX,
     &    nnodes_k-1,102,MPI_COMM_K1,status,ierr)

         end if

        return
        end subroutine read_ug

        end

