      subroutine ugIOBP(ug_n_tmp,kpt,iflag,istep,iislda,
     &   iflag2,nkpt,islda)
******************************************
c     c     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
*     *  copyright (c) 2003, The Regents of the University of California,
*     *  through Lawrence Berkeley National Laboratory (subject to receipt of any
*     *  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************
cccccc     is all the icolor groups are writing or reading, 
cccccc     or only the group with kpt inside the [kpt_dis(1),kpt_dis(2)] are writing or reading ?

******************************************

****************************************
ccccc iflag=1, write, iflag=2, read
****************************************
******************************************
      use fft_data
      use load_data
      use data
      implicit double precision (a-h,o-z)

      include 'mpif.h'
      include 'param.escan_real'

      integer status(MPI_STATUS_SIZE)

      complex*16 ug_n_tmp(mg_nx,nblock_band_mx)
      character*8 fname
*************************************************
      if(nkpt.eq.1.and.islda.eq.1.and.istep.eq.0) return
ccccc don't do anything in this case
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       if(iflag2.lt.0) then         ! read ug according to icolor_k and kpt relationship
       if((iislda-1)*nkpt+kpt.lt.kpt_slda_dis(1).or.
     &    (iislda-1)*nkpt+kpt.gt.kpt_slda_dis(2)) return
       endif 

       if(iflag2.ge.0) then      ! read ug only for icolor_k.eq.iflag
       if(icolor_k.ne.iflag2) return
       endif 
      
*************************************************
      if(iflag.eq.2) then       ! read the wavefunction

      call mpi_barrier(MPI_COMM_K1,ierr)
 
c$$$         write(6,*) 'ugIOBP read', ikey_k, size(ug_n_tmp)
c$$$         call system_flush(6)
        
         if (ikey_k==nnodes_k-1) then
            
            fname="ugiofile"
            kpt1=mod(kpt,10)
            kpt2=mod((kpt-kpt1)/10,10)
            kpt3=mod((kpt-kpt1-kpt2*10)/100,10)
            kpt4=mod((kpt-kpt1-kpt2*10-kpt3*100)/1000,10)            
            open(10,file=
     &           fname//char(iislda+48)//char(istep+48)//
     &           char(kpt4+48)//char(kpt3+48)//char(kpt2+48)//
     &           char(kpt1+48),form="unformatted")
            rewind(10)
            
            do i=0,nnodes-1
               do iwavefun=1,mx
                  
                  icolor_bb = (iwavefun-1)/nblock_band_mx
                  ib = iwavefun - icolor_bb*nblock_band_mx

                  read(10) ug_n_tmp(:,ib)

                  ikey_bb = i
                  ikey_kk = icolor_bb*nnodes_b + ikey_bb

                  if (iwavefun==mx .or. ib==nblock_band_mx) then
                     if(ikey_kk==nnodes_k-1) then 
                                ! No broad casting. Do nothing
                     else
                        call mpi_send(ug_n_tmp,mg_nx*nblock_band_mx,
     &                       MPI_DOUBLE_COMPLEX,ikey_kk,102,
     &                       MPI_COMM_K1,ierr)
                     endif
                  endif
                  
               enddo
            enddo
            close(10)

         else                   ! for other nodes

            call mpi_recv(ug_n_tmp,mg_nx*nblock_band_mx,
     &           MPI_DOUBLE_COMPLEX,nnodes_k-1,102,MPI_COMM_K1,
     &           status,ierr)

         endif                  ! for the nodes
         call mpi_barrier(MPI_COMM_K1,ierr)

      endif                     ! for the iflag, read
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if(iflag.eq.1) then       ! write the wavefunction

         if (ikey_k==0) then
            fname="ugiofile"
            kpt1=mod(kpt,10)
            kpt2=mod((kpt-kpt1)/10,10)
            kpt3=mod((kpt-kpt1-kpt2*10)/100,10)
            kpt4=mod((kpt-kpt1-kpt2*10-kpt3*100)/1000,10)
            open(10,file=
     &           fname//char(iislda+48)//char(istep+48)//
     &           char(kpt4+48)//char(kpt3+48)//char(kpt2+48)//
     &           char(kpt1+48),form="unformatted")
            rewind(10)
         endif
      
         do i=0, nnodes-1
            do j=1, num_group_b
               icolor_bb = j-1
               ikey_bb = i
               ikey_kk = icolor_bb*nnodes_b + ikey_bb
                              
               if (ikey_kk==0) then
                  if (ikey_k==0) then
                     do m=1, nblock_band_mx
                        write(10) ug_n_tmp(:,m)
                     enddo 
                  endif
               else
                  if (ikey_k==ikey_kk) then
                     call mpi_send(ug_n_tmp,mg_nx*nblock_band_mx,
     &                    MPI_DOUBLE_COMPLEX,0,102,MPI_COMM_K1,ierr)
                  else if (ikey_k==0) then
                     call mpi_recv(ug_n_tmp,mg_nx*nblock_band_mx,
     &                    MPI_DOUBLE_COMPLEX,ikey_kk,102,
     &                    MPI_COMM_K1,status,ierr)
                     if (j==num_group_b) then
                        mend = mx - (j-1)*nblock_band_mx
                     else
                        mend = nblock_band_mx
                     endif
                     do m=1, mend
                        write(10) ug_n_tmp(:,m)
                     enddo 
                  endif
               endif
               call mpi_barrier(MPI_COMM_K1,ierr) 
            enddo    ! do j
         enddo      ! do i
            
         if(ikey_k==0) then
            rewind(10)
              do m=1, nblock_band_mx
              read(10) ug_n_tmp(:,m)   ! for the first node in each kpoint_group, read back, so the ug_n_tmp will not be changed after write
              enddo 
            close(10)
         endif


        call system_flush(6)

c$$$!         if (ikey_k==0 .and. iband==1) then
c$$$         if (ikey_k==0) then
c$$$            fname="ugiofile"
c$$$            kpt1=mod(kpt,10)
c$$$            kpt2=mod((kpt-kpt1)/10,10)
c$$$            kpt3=mod((kpt-kpt1-kpt2*10)/100,10)
c$$$            kpt4=mod((kpt-kpt1-kpt2*10-kpt3*100)/1000,10)
c$$$            open(10,file=
c$$$     &           fname//char(iislda+48)//char(istep+48)//
c$$$     &           char(kpt4+48)//char(kpt3+48)//char(kpt2+48)//
c$$$     &           char(kpt1+48),form="unformatted")
c$$$            rewind(10)
c$$$
c$$$            do m=1, nblock_band_mx
c$$$               write(10) ug_n_tmp(:,m)
c$$$               write(6,*) 'ugIOBP1 write', kpt, m, ikey_k, inode_tot,
c$$$     &              size(ug_n_tmp(:,m))
c$$$            enddo                        
c$$$         endif
c$$$      
c$$$
c$$$         do ikey_kk=1, nnodes_k-1
c$$$            icolor_bb = ikey_kk/nnodes_b
c$$$            ikey_bb = ikey_kk - icolor_bb*nnodes_b
c$$$            
c$$$            if(ikey_bb == nnodes_b-1) then
c$$$               mend = mx - (nnodes_b-1)*nblock_band_mx
c$$$            else
c$$$               mend = nblock_band_mx
c$$$            endif
c$$$            
c$$$            if (ikey_kk==ikey_k) then
c$$$               call mpi_send(ug_n_tmp,mg_nx*nblock_band_mx,
c$$$     &              MPI_DOUBLE_COMPLEX,0,102,MPI_COMM_K1,ierr)
c$$$            endif
c$$$            
c$$$            if(ikey_k==0) then
c$$$               call mpi_recv(ug_n_tmp,mg_nx*nblock_band_mx,
c$$$     &              MPI_DOUBLE_COMPLEX,ikey_kk,102,MPI_COMM_K1,
c$$$     &              status,ierr)
c$$$               do m=1, mend
c$$$                  write(10) ug_n_tmp(:,m)
c$$$                  write(6,*) 'ugIOBP2 write', kpt, m, ikey_k, ikey_kk,
c$$$     &                 size(ug_n_tmp(:,m))
c$$$               enddo
c$$$            endif
c$$$
c$$$         enddo
c$$$
c$$$         call mpi_barrier(MPI_COMM_K1,ierr)
c$$$         call system_flush(6)
c$$$         if(ikey_k==0) then
c$$$            close(10)
c$$$         endif

      endif                     ! for the iflag, write


      
      return
      end

