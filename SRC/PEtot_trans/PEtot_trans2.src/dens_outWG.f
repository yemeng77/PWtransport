      subroutine dens_outWG(AL,workr_n,m1_out,m2_out,
     &   fdens_out,kpt_dens)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
******************************************

****************************************
***** this subroutine output the charge density.
******************************************

      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)

      include 'mpif.h'
      include 'param.escan_real'
***********************************************
       complex*16 workr_n(mr_n)
       real*8,allocatable,dimension(:)  ::  temp_rho
       real*8 AL(3,3)
       character*20 fdens_out
       integer status(MPI_STATUS_SIZE)
       complex*16 cai

**********************************************
       common /com123b/m1,m2,m3,ngb,nghb,ntype,rrcut,msb

       if(inode.eq.1) then
       open(11,file=fdens_out,form="unformatted")
       rewind(11)
       write(11) n1,n2,n3
       write(11) AL
       write(11) nnodes
       endif

       cai=dcmplx(0.d0,1.d0)

      do 1000 m=m1_out,m2_out

      call d3fft_comp(ug_n(1,m),workr_n,-1,kpt_dens)
cccccccccccccccccccccccccccccccccccccccc
       do ii=1,nr_n
       iit=ii+(inode-1)*nr_n-1
       it=iit/(n3*n2)
       jt=(iit-it*n3*n2)/n3
       kt=iit-it*n3*n2-jt*n3
       xt=AL(1,1)*it/n1+AL(1,2)*jt/n2+AL(1,3)*kt/n3
       yt=AL(2,1)*it/n1+AL(2,2)*jt/n2+AL(2,3)*kt/n3
       zt=AL(3,1)*it/n1+AL(3,2)*jt/n2+AL(3,3)*kt/n3
       workr_n(ii)=workr_n(ii)*cdexp(cai*(
     & xt*akx(kpt_dens)+yt*aky(kpt_dens)+      ! special
     & zt*akz(kpt_dens)))
       enddo


      call mpi_barrier(MPI_COMM_WORLD,ierr)

       if(inode.eq.1) then
       write(11) (workr_n(i),i=1,nr_n)
       endif

       do i=1,nnodes-1
      call mpi_barrier(MPI_COMM_WORLD,ierr)
       if(inode==i+1) then
       call  mpi_send(workr_n,nr_n,MPI_DOUBLE_COMPLEX,0,
     &   100,MPI_COMM_WORLD,ierr)
       endif
       if(inode.eq.1) then
        call mpi_recv(workr_n,nr_n,MPI_DOUBLE_COMPLEX,i,
     &   100,MPI_COMM_WORLD,status,ierr)
       write(11) (workr_n(ii),ii=1,nr_n)
       endif
       enddo

1000   continue

       if(inode.eq.1)  close(11)
*********************************************

      return
      end

