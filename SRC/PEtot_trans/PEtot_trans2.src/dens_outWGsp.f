      subroutine dens_outWGsp(AL,workr_n,m1_out,m2_out,
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

       n1st=44
       n1fn=76
	if(inode.eq.1) then
	write(6,*) 'n1st,n1fn',n1st,n1fn
	endif

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
       sum=0.d0
       do ii=1,nr_n
       iit=ii+(inode-1)*nr_n-1
       it=iit/(n3*n2)
       jt=(iit-it*n3*n2)/n3
       kt=iit-it*n3*n2-jt*n3

       if(it.gt.n1st.and.it.lt.n1fn) then
       sum=sum+abs(workr_n(ii))**2
       endif
       enddo
       sum=sum*vol/(n1*n2*n3)
  
       call global_sumr(sum)

	write(36,*) vol/(n1*n2*n3)

       write(6,*) "m,sum=", m,sum
       write(16,*)  m,sum
      
1000   continue

       if(inode.eq.1)  close(11)
*********************************************

      return
      end

