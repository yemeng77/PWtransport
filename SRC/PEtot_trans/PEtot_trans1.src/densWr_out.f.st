      subroutine densWr_out(AL,workr_n,kpt_dens,
     & ispin_dens,iw_dens,fdens_out,idens_out,nkpt,islda)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

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
       character*25 file_tmp
       integer status(MPI_STATUS_SIZE)
       integer kpt_dens(2),ispin_dens(2),iw_dens(2)
       complex*16 cai
**********************************************
       common /com123b/m1,m2,m3,ngb,nghb,ntype,rrcut,msb

      if(icolor_k.ne.0) return

      do 100 ikpt=kpt_dens(1),kpt_dens(2)     ! this sub is called by icolor_k=0 only

      call gen_G_comp(ikpt,0)
      call fftprep_comp(n1,n2,n3)

      ikpt_100=ikpt/100
      ikpt_10=(ikpt-ikpt_100*100)/10
      ikpt_1=ikpt-ikpt_100*100-ikpt_10*10
      file_tmp=trim(fdens_out)//char(48+ikpt_100)//
     &  char(48+ikpt_10)//char(48+ikpt_1)


      if(inode_b.eq.1.and.icolor_b.eq.0) then
      open(11,file=file_tmp,form="unformatted")
      rewind(11)
      write(11) n1,n2,n3,nnodes,
     &  ispin_dens(1),ispin_dens(2),
     &  iw_dens(1),iw_dens(2)
      write(11) AL
      endif

      cai=dcmplx(0.d0,1.d0)

      do 80 iislda=ispin_dens(1),ispin_dens(2)


c      call ugIO(ug_n,ikpt,2,0,iislda)
      call ugIOBP(ug_n_BP,ikpt,2,0,iislda,0,nkpt,islda)


      do 80 m_all=iw_dens(1),iw_dens(2)


      call mpi_barrier(MPI_COMM_K1,ierr)

      if(m_all.ge.band_dis(1).and.
     &     m_all.le.band_dis(2)) then
      m=m_all-band_dis(1)+1


      call d3fft_comp(ug_n_BP(1,m),workr_n,-1,ikpt)

ccccccc for idens_out=: 12,22, add the cphase,    output phi_k(r)
ccccccc for idens_out=: 11,21, do not add cphase, output u_k(r)
      if(idens_out.eq.12.or.idens_out.eq.22) then   
      do ii=1,nr_n
       iit=ii+(inode_b-1)*nr_n-1
       it=iit/(n3*n2)
       jt=(iit-it*n3*n2)/n3
       kt=iit-it*n3*n2-jt*n3
       xt=AL(1,1)*it/n1+AL(1,2)*jt/n2+AL(1,3)*kt/n3
       yt=AL(2,1)*it/n1+AL(2,2)*jt/n2+AL(2,3)*kt/n3
       zt=AL(3,1)*it/n1+AL(3,2)*jt/n2+AL(3,3)*kt/n3
       workr_n(ii)=workr_n(ii)*cdexp(cai*(
     & xt*akx(ikpt)+yt*aky(ikpt)+zt*akz(ikpt)))
       enddo
       endif


        if(icolor_b.ne.0) then
        call mpi_send(workr_n,nr_n,MPI_DOUBLE_COMPLEX,0,
     &   101,MPI_COMM_B2,ierr)
ccccc put everything in the first icolor_b group (size nnodes_b)
        endif
        call mpi_barrier(MPI_COMM_B1,ierr)
        endif    !  m_all inside [band_dis(1),band_dis(2)], only one nnodes_b group is in there
ccccccccccccccccccccccccccccccccccccccccccccccc

        if(icolor_b.eq.0) then

         icolor_b_send=(m_all-1)/nblock_band_mx
         if(icolor_b_send.ne.0) then
        call mpi_recv(workr_n,nr_n,MPI_DOUBLE_COMPLEX,
     &   icolor_b_send,101,MPI_COMM_B2,status,ierr)
         endif


       if(inode_b.eq.1) then
       write(11) (dreal(workr_n(j)),j=1,nr_n),
     & (dimag(workr_n(j)),j=1,nr_n)
       endif

        
       do i=1,nnodes_b-1

      call mpi_barrier(MPI_COMM_B1,ierr)

       if(inode_b==i+1) then
       call  mpi_send(workr_n,nr_n,MPI_DOUBLE_COMPLEX,0,
     &   100,MPI_COMM_B1,ierr)
       endif

       if(inode_b.eq.1) then

        call mpi_recv(workr_n,nr_n,MPI_DOUBLE_COMPLEX,i,
     &   100,MPI_COMM_B1,status,ierr)

       write(11) (dreal(workr_n(j)),j=1,nr_n),
     & (dimag(workr_n(j)),j=1,nr_n)
       endif
       enddo

       endif     ! icolor_b.eq.0

80    continue
      if(inode_b.eq.1.and.icolor_b.eq.0) close(11)
100   continue

      return
      end


