      subroutine write_wg_BP(fwg_out,AL,islda,nkpt)

*************************************************************************
*** Written by Lin-Wang Wang, 2001
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************


ccccccccccccccccccccccccccccccccccccccccccc
ccccc This should only be called by one icolor group, it will write all the kpt into one wg file 

      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)

      include 'mpif.h'
                                                                                                          
      include 'param.escan_real'

       real*8 AL(3,3)
       character*20 fwg_out(2)

       complex*16,allocatable,dimension(:) :: ugtemp
       

       integer status(MPI_STATUS_SIZE)


       if(icolor_k.eq.0) then     ! only do this with the first k-point group



       do 280 iislda=1,islda
                                                                                                          
                                                                                                          
       if(inode_k.eq.1) then
         allocate(ugtemp(mg_nx))
         open(11,file=fwg_out(iislda),form='unformatted')
         rewind(11)
         write(11)n1,n2,n3,mx
         write(11)Ecut
         write(11)AL
         write(11)nnodes
       end if
                  
       do 180 kpt=1,nkpt
                  
ccccccccccccccccccccccccccccccccccccccccccccc
cccc need to change later, abolish the use of ug_n (for memory, etc)
c       call ugIO(ug_n,kpt,2,0,iislda)
       call ugIOBP(ug_n_BP,kpt,2,0,iislda,0,nkpt,islda)
ccccc the write and read will be done by icolor_k.eq.0 


       if(inode_k.eq.1) then

       do inode_tmp=1,nnodes_b
       do igroup_tmp=1,num_group_b
       idest=(igroup_tmp-1)*nnodes_b+inode_tmp-1
       if(idest.ne.0) then
       call mpi_recv(ug_n_BP,mg_nx*nblock_band_mx,MPI_DOUBLE_COMPLEX,
     &  idest,102,MPI_COMM_K1,status,ierr)
       endif

          do iwavefun=1,nblock_band_mx
          ugtemp(:)=ug_n_BP(:,iwavefun) 
          write(11) ugtemp
          enddo
        enddo
        enddo
     
        else   ! for other nodes

        call mpi_send(ug_n_BP,mg_nx*nblock_band_mx,MPI_DOUBLE_COMPLEX,
     &  0,102,MPI_COMM_K1,ierr)

        endif

180   continue    ! kpt=1,nkpt
              
      if (inode_k.eq.1) then
      close(11)
      deallocate(ugtemp)
      endif
280   continue   ! iislda=1,islda


      endif   ! icolor_k.eq.0

      return
      end
              


       

       
      
