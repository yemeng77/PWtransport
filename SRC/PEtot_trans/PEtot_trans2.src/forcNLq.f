       subroutine forcNLq(fatom,occ,kpt,nkpt)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
******************************************


       use fft_data
       use load_data
       use data

       implicit double precision (a-h,o-z)

       include 'param.escan_real'
       include "mpif.h"

       real*8 occ_t(mtype)
       integer iiatom(mtype),icore(mtype),numref(matom)

       real*8 fatom(3,matom)
       real*8 occ(mst,nkpt)
c Modificado por Txomin
c       complex*16 s(9),sx(9),sy(9),sz(9)
       complex*16 s(9),s_rec(9),sx(9),sx_rec(9)
       complex*16 sy(9),sy_rec(9),sz(9),sz_rec(9)
	complex*16 cc,cx,cy,cz,cai
*************************************************
       integer isNLa(9,matom)
       common /comisNLa/isNLa
       common /comNL2/occ_t,iiatom,icore,numref


       cai=dcmplx(0.d0,1.d0)
       s=dcmplx(0.d0,0.d0)
       sx=dcmplx(0.d0,0.d0)
       sy=dcmplx(0.d0,0.d0)
       sz=dcmplx(0.d0,0.d0)
c Modificado por Txomin
       s_rec=dcmplx(0.d0,0.d0)
       sx_rec=dcmplx(0.d0,0.d0)
       sy_rec=dcmplx(0.d0,0.d0)
       sz_rec=dcmplx(0.d0,0.d0)





       ng_n=ngtotnod(inode,kpt)

       do 40 ia=1,natom
       nref=numref(ia)

       do 40 m=1,mx

	s=dcmplx(0.d0,0.d0)
	sx=dcmplx(0.d0,0.d0)
	sy=dcmplx(0.d0,0.d0)
	sz=dcmplx(0.d0,0.d0)
c Modificado por Txomin
       s_rec=dcmplx(0.d0,0.d0)
       sx_rec=dcmplx(0.d0,0.d0)
       sy_rec=dcmplx(0.d0,0.d0)
       sz_rec=dcmplx(0.d0,0.d0)

	do i=1,ng_n
	cc=dconjg(ug_n(i,m))
	cx=cai*gkx_n(i,kpt)*cc
	cy=cai*gky_n(i,kpt)*cc
	cz=cai*gkz_n(i,kpt)*cc
          do j=1,nref
          s(j)=s(j)+wqmask(j,i,ia)*cc
          sx(j)=sx(j)+wqmask(j,i,ia)*cx
          sy(j)=sy(j)+wqmask(j,i,ia)*cy
          sz(j)=sz(j)+wqmask(j,i,ia)*cz
c Modificado por Txomin
          s_rec(j)=s_rec(j)+wqmask(j,i,ia)*cc
          sx_rec(j)=sx_rec(j)+wqmask(j,i,ia)*cx
          sy_rec(j)=sy_rec(j)+wqmask(j,i,ia)*cy
          sz_rec(j)=sz_rec(j)+wqmask(j,i,ia)*cz
	enddo
	enddo

*****************************
**** it used: mpi_allreduce(s,s,9,MPI_DOUBLE_COMPLEX,...)
**** That causes overflow, since the s(nref+1) might be used
**** in the previous atoms, and is not set to zero at this atom (in the previous version). 
**** Each time mpi_allreduce is called, it has been multiplied by node. 
**** And it will be done for mx times for each atom !

c Modificado por Txomin
c        call mpi_allreduce(s,s,nref,MPI_DOUBLE_COMPLEX,
        call mpi_allreduce(s,s_rec,nref,MPI_DOUBLE_COMPLEX,
     & MPI_SUM,MPI_COMM_WORLD,ierr)
     	s = s_rec
c        call mpi_allreduce(sx,sx,nref,MPI_DOUBLE_COMPLEX,
        call mpi_allreduce(sx,sx_rec,nref,MPI_DOUBLE_COMPLEX,
     & MPI_SUM,MPI_COMM_WORLD,ierr)
     	sx = sx_rec
c        call mpi_allreduce(sy,sy,nref,MPI_DOUBLE_COMPLEX,
        call mpi_allreduce(sy,sy_rec,nref,MPI_DOUBLE_COMPLEX,
     & MPI_SUM,MPI_COMM_WORLD,ierr)
     	sy = sy_rec
c        call mpi_allreduce(sz,sz,nref,MPI_DOUBLE_COMPLEX,
        call mpi_allreduce(sz,sz_rec,nref,MPI_DOUBLE_COMPLEX,
     & MPI_SUM,MPI_COMM_WORLD,ierr)
     	sz = sz_rec


	 fx=0.d0
	 fy=0.d0
	 fz=0.d0
	 do j=1,nref
	 fx=fx+dreal(dconjg(s(j))*sx(j))*isNLa(j,ia)
	 fy=fy+dreal(dconjg(s(j))*sy(j))*isNLa(j,ia)
	 fz=fz+dreal(dconjg(s(j))*sz(j))*isNLa(j,ia)
	 enddo

	 fatom(1,ia)=fatom(1,ia)+2*fx*vol**2*occ(m,kpt)
	 fatom(2,ia)=fatom(2,ia)+2*fy*vol**2*occ(m,kpt)
	 fatom(3,ia)=fatom(3,ia)+2*fz*vol**2*occ(m,kpt)
 40     continue
*****************************************
*** The factors are right. 
*** The 2 is from the derivative of (\psi*\psi)^2
*****************************************


       return
       end
