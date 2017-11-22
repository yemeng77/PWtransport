       subroutine forcNLr(workr_n,fatom,
     &   occ,kpt,nkpt,iislda,islda,E_st)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************


       use fft_data
       use load_data
       use data

       implicit double precision (a-h,o-z)
       include 'param.escan_real'
       include "mpif.h"

       real*8 fatom(3,matom)
       real*8 fatom_tmp(3,natom),f_tmp(3)

       complex*16 workr_n(mr_n)

       integer nmap(matom)
       integer isNLa(9,matom)
       real*8 Dij0(32,32,mtype),Qij(32,32,mtype)
       real*8 y_tmp1(10000),y_tmp2(10000)
       real*8 sumy1(mref),sumy2(mref)


       real*8 occ_t(mtype)
       integer iiatom(mtype),icore(mtype),numref(matom),
     &   ityatom(matom),ipsp_type(mtype)
       real*8 occ(mst,nkpt,islda),E_st(mst,nkpt,islda)

       complex*16 s1(mref,matom),s2(mref,matom,3),stmp(mref,matom)
       complex*16 y

       complex*16, allocatable, dimension (:,:,:,:) :: ss,ss2
       complex*16, allocatable, dimension (:,:) :: ss_tmp,ss2_tmp
************************************************
*************************************************
       common /comisNLa/isNLa,Dij0,Qij,ipsp_all,ipsp_type
       common /comnmap/nmap
       common /comNL2/occ_t,iiatom,icore,numref,ityatom

cccccc only the appropriat icolor_k has called this routine for this kpt
cccccccccccccccccccccccccccccccccccccc
       do ia=1,natom
       fatom_tmp(:,ia)=0.d0
       enddo

       if(ipsp_all.eq.2) then
       allocate(ss(mref,mref,natom,3))
       allocate(ss2(mref,mref,natom,3))
       allocate(ss_tmp(mref,mref))
       allocate(ss2_tmp(mref,mref))
       ss=dcmplx(0.d0,0.d0)
       ss2=dcmplx(0.d0,0.d0)
       ss_tmp=dcmplx(0.d0,0.d0)
       ss2_tmp=dcmplx(0.d0,0.d0)
       endif

c       do 40 m=1,mx
       do 40 m=1,nblock_band_mx
       m_all=icolor_b*nblock_band_mx+m 

       call d3fft_comp(ug_n_BP(1,m),workr_n,-1,kpt)

       s1 = dcmplx(0.d0,0.d0)
       s2 = dcmplx(0.d0,0.d0)

       ico1=0
       ico2=0
       do ia=1,natom
       nref=numref(ia)

cccccccccccccccccccccccccccccc
c	do i=1,nmap(ia)
c         ico1=ico1+1
c         y=workr_n(indm(ico1))*cphase(ico1)  
c          do jj = 1,nref
c          ico2=ico2+1
c          s1(jj,ia)=s1(jj,ia)+wmask(ico2)*y
c          s2(jj,ia)=s2(jj,ia)+wmaskX(ico2,ixyz)*y
c          enddo
c        enddo
ccccccccccccccccccccccccccccccccccccccccc
        do i=1,nmap(ia)
        ico1=ico1+1
         y=workr_n(indm(ico1))*cphase(ico1)
         y_tmp1(i)=dreal(y)
         y_tmp2(i)=dimag(y)
        enddo
        
ccccccccccccccccccccccccccccc
cccccccccc wmask is real*8 not complex*16 !
ccccccccccc maybe combine this to a matrix*matrix multiplication ?
        
       call dgemv('N',nref,nmap(ia),1.d0,wmask(ico2+1),nref,
     & y_tmp1,1,0.d0,sumy1,1)
       call dgemv('N',nref,nmap(ia),1.d0,wmask(ico2+1),nref,
     & y_tmp2,1,0.d0,sumy2,1)

         if(nmap(ia).gt.0) then
        do jj=1,nref
        s1(jj,ia)=dcmplx(sumy1(jj),sumy2(jj))
        enddo
          endif

         do ixyz=1,3
       call dgemv('N',nref,nmap(ia),1.d0,wmaskX(ico2+1,ixyz),nref,
     & y_tmp1,1,0.d0,sumy1,1)
       call dgemv('N',nref,nmap(ia),1.d0,wmaskX(ico2+1,ixyz),nref,
     & y_tmp2,1,0.d0,sumy2,1)

         if(nmap(ia).gt.0) then
        do jj=1,nref
        s2(jj,ia,ixyz)=dcmplx(sumy1(jj),sumy2(jj))
        enddo
          endif
         enddo   ! ixyz

         ico2=ico2+nmap(ia)*nref
 
       enddo    ! ia=1,natom
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




       call mpi_barrier(MPI_COMM_b1,ierr)

       call mpi_allreduce(s1,stmp,natom*mref,
     $     MPI_DOUBLE_COMPLEX,MPI_SUM, MPI_COMM_b1,ierr)
       s1 = stmp

       do ixyz=1,3
       call mpi_allreduce(s2(1,1,ixyz),stmp,natom*mref,
     $     MPI_DOUBLE_COMPLEX,MPI_SUM, MPI_COMM_b1,ierr)
       s2(:,:,ixyz) = stmp(:,:)
       enddo

       s1=s1*vol/nr
       s2=s2*vol/nr

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       do 21 ixyz=1,3
       do 20 ia=1,natom
       iitype=ityatom(ia)
       nref=numref(ia)

       if(ipsp_type(iitype).eq.2) then
         do j2=1,nref
         do j1=1,nref
         ss(j1,j2,ia,ixyz)=ss(j1,j2,ia,ixyz)+dconjg(s1(j2,ia))*
     &                        s2(j1,ia,ixyz)*occ(m_all,kpt,iislda)
         ss2(j1,j2,ia,ixyz)=ss2(j1,j2,ia,ixyz)+dconjg(s1(j2,ia))*
     &    s2(j1,ia,ixyz)*occ(m_all,kpt,iislda)*E_st(m_all,kpt,iislda)
         enddo
         enddo

        else

         fx=0.d0
         do j1=1,nref
         fx=fx+Dij(j1,j1,ia,iislda)*dconjg(s1(j1,ia))*
     &                      s2(j1,ia,ixyz)*occ(m_all,kpt,iislda)
         enddo
         fatom_tmp(ixyz,ia)=fatom_tmp(ixyz,ia)+2*fx
        endif

20      continue
21      continue



40       continue

cccccccccccccccccccccccccccc



       do 50 ia=1,natom
       iitype=ityatom(ia)
       nref=numref(ia)

        if(ipsp_type(iitype).eq.1) then

         call mpi_allreduce(fatom_tmp(1,ia),f_tmp,3,MPI_REAL8,
     & MPI_SUM,MPI_COMM_b2,ierr)
       fatom(1:3,ia)=fatom(1:3,ia)+f_tmp(1:3)

        else

       do 51 ixyz=1,3

          call mpi_allreduce(ss(1,1,ia,ixyz),ss_tmp,mref*mref,
     & MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_b2,ierr)

          call mpi_allreduce(ss2(1,1,ia,ixyz),ss2_tmp,mref*mref,
     & MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_b2,ierr)

         fx=0.d0
         do j1=1,nref
         do j2=1,nref
         fx=fx+Dij(j1,j2,ia,iislda)*ss_tmp(j1,j2)
     &        -Qij(j1,j2,iitype)*ss2_tmp(j1,j2)
         enddo
         enddo
         fatom(ixyz,ia)=fatom(ixyz,ia)+2*fx
51     continue
       endif

50     continue
cccccccccccccccccccccccccccccccccccccc
       if(ipsp_all.eq.2) then
       deallocate(ss)
       deallocate(ss2)
       deallocate(ss_tmp)
       deallocate(ss2_tmp)
       endif

       return
       end
