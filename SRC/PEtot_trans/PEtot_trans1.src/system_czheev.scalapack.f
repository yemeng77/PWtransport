      subroutine system_czheev(w1,w2,n,z,m,EE,workx,
     &  lwork,workrx,info,MPI_COMM_K1)
*************************************************************************
*** Written by Lin-Wang Wang, 2001
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

ccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      include 'mpif.h'

      integer lworktmp
      parameter (lworktmp=500000)
      integer inode_tmp,nnodes_tmp,
     &     status(MPI_STATUS_SIZE),ierr
      integer MPI_COMM_K1

      character*1 w1,w2
      integer n,m,lwork,ncol
      integer inode_tmp2,inode_tmp3
      complex*16 workx(lwork),z(m,m),workrx(2000)
      complex*16 workxtmp(lworktmp)
      real*8 EE(n)
      complex*16,allocatable,dimension(:) :: workrx_m

      complex*16,allocatable,dimension(:,:) :: hh_m,ev_m
      integer ictxt,nprow,npcol,myrow,mycol,info

      integer dlen_,msec1,msec2
      parameter (dlen_=9,msec1=200,msec2=500)
      integer mb,nb,maxllda
      integer desca(dlen_),descz(dlen_)
      integer nsec(0:msec1),isec(2,msec2,0:msec1)
      integer irow1,irow2,icol,i,j,ib,ip,is,js,i1,i2
cc      data mb/32/,nb/32/
      integer unit

      integer numroc
      integer sys2blacs_handle
     

cccccccccccccccccccccccccccccccccccccc
c       mb=32
c       nb=32

ccccc for T3E, lapack routine (cheev is for complex*16 in T3E)
ccccc in Cray T3E, the single precision is defined as real*8 and complex*16
c      call cheev(w1,w2,n,z,m,EE,workx,
c     &  lwork,workrx,info)
cccccccccccccccccccccccccccccccccccccccccc

c      if(m.lt.200) then
      if(m.lt.100) then
      call zheev(w1,w2,n,z,m,EE,workx,
     &  lwork,workrx,info)

      else

      call mpi_comm_rank(MPI_COMM_K1,inode_tmp,ierr)
      call mpi_comm_size(MPI_COMM_K1,nnodes_tmp,ierr)
      nprow = int(sqrt(nnodes_tmp*1.d0)*1.00001d0)
      npcol = nprow
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      mb=n/nprow
      if(mb.gt.32) mb=32
      if(mb.lt.1) mb=1
      nb=mb
ccccccccccccccccccccccccccccccccccccccccccccccccc

      if(nprow.gt.msec1)then
        if(inode_tmp.eq.0)then
          write(6,*)'nprow.gt.msec1!',nprow,msec1
        endif
        call mpi_abort(MPI_COMM_WORLD,1,ierr)
      endif
      if(n/mb/nprow+1.gt.msec2)then
        if(inode_tmp.eq.0)then
          write(6,*)'n/mb/nprow+1.gt.msec2!',n,mb,nprow,n/mb/nprow+1,
     $      msec2
        endif
        call mpi_abort(MPI_COMM_WORLD,1,ierr)
      endif


ccccc get scalapack handler ictxt from MPI communicator
       ictxt=sys2blacs_handle(MPI_COMM_K1)
ccc input nprow,npcol, provide ictxt (a handle)
      call blacs_gridinit(ictxt,"R",nprow,npcol)
      call blacs_gridinfo(ictxt,nprow,npcol,myrow,mycol)

ccccc give the myrow, mycol index in the processor array [nprow,npcol] for this processor


cccc If a processor is not used, myrow, mycol will be negative
      if(myrow.ge.nprow.or.mycol.ge.npcol.
     &  or.myrow.lt.0.or.mycol.lt.0)then
cccc for this processor, there is nothing to do
        inode_tmp2=-1
        goto 10
      endif

cccc numroc calculates the number of columns or rows own by this processors 
cccc nb is the size of each block

      ncol=numroc(n,nb,mycol,0,npcol)
      maxllda=numroc(n,nb,myrow,0,nprow)


      call descinit(desca,n,n,mb,nb,0,0,ictxt,maxllda,info)
      call descinit(descz,n,n,mb,nb,0,0,ictxt,maxllda,info)
cccccc return the describer desca, descz
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc how can the subroutine knows how many processors (which MPI_COMM_XXX) are calling for the scalapack ?

      allocate(hh_m(maxllda,ncol))
      allocate(ev_m(maxllda,ncol))
      allocate(workrx_m(4*n))

ccccccccccccccccccccccccccccccccccccccccccccccccc
      if(myrow.ge.0.and.mycol.ge.0) then    ! only these processors are used
      nsec(:)=0
      do i=1,n,mb
        ib=int((i-1)/mb) 
        ip=mod(ib,nprow)
        nsec(ip)=nsec(ip)+1
        isec(1,nsec(ip),ip)=i
        isec(2,nsec(ip),ip)=min(n,i+mb-1)
      enddo

      icol=0
      do js=1,nsec(mycol)
        do j=isec(1,js,mycol),isec(2,js,mycol)
          icol=icol+1
          irow1=1
          do is=1,nsec(myrow)
            i1=isec(1,is,myrow)
            i2=isec(2,is,myrow)
            irow2=irow1+(i2-i1)
            hh_m(irow1:irow2,icol)=z(i1:i2,j)
            irow1=irow2+1
          enddo
        enddo
      enddo
       endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c  check if `lwork' is large enough
c if JOBZ='V', LWORK>=(NP0+NQ0+NB)*NB+3*N+N^2, where NP0=NUMROC(NN,NB,0,0,NPROW), NQ0=NUMROC(MAX(N,NB,2),NB,0,0,NPCOL)
      call pzheev(w1,w2,n,hh_m,1,1,desca,EE,ev_m,1,1,descz,workxtmp,
     &     -1,workrx_m,4*n,info)

      if(real(workxtmp(1)).gt.lworktmp)then
         if(inode_tmp.eq.0)then
            write(6,*)'increase lworktmp!',lworktmp,real(workxtmp(1))
         endif
         call mpi_abort(MPI_COMM_WORLD,2,ierr)
      endif

      call pzheev(w1,w2,n,hh_m,1,1,desca,EE,ev_m,1,1,descz,workxtmp,
     &     lworktmp,workrx_m,4*n,info)

      z(:,:)=0.d0

      icol=0
      do js=1,nsec(mycol)
         do j=isec(1,js,mycol),isec(2,js,mycol)
            icol=icol+1
            irow1=1
            do is=1,nsec(myrow)
               i1=isec(1,is,myrow)
               i2=isec(2,is,myrow)
               irow2=irow1+(i2-i1)
               z(i1:i2,j)=ev_m(irow1:irow2,icol)
               irow1=irow2+1
            enddo
         enddo
      enddo

      call zgsum2d(ictxt,'A'," ",n,n,z,n,-1,-1)

      deallocate(hh_m)
      deallocate(ev_m)
      deallocate(workrx_m)


      call blacs_gridexit(ictxt)
      inode_tmp2=inode_tmp
 10   continue


cccccccccccc There are some processors which didn't do pzheev, thus the results have to be passed to them

      call mpi_allreduce(inode_tmp2,inode_tmp3,1,
     & MPI_INTEGER,MPI_MAX,MPI_COMM_K1,ierr)

      call mpi_bcast(z,m*m,MPI_DOUBLE_COMPLEX,
     &   inode_tmp3, MPI_COMM_K1,ierr) 

      call mpi_bcast(EE,n,MPI_REAL8,
     &   inode_tmp3, MPI_COMM_K1,ierr) 

      endif    ! big endif

      return

      end subroutine
