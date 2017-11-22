      subroutine get_atomW(AL,xatom,xatom_w,
     &   ntype,iatom,totNel,islda)

******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
******************************************


ccccc if ido_rho.eq.1, generate rho_n from atom, 
ccccc if ido_rho.eq.0, do not touch rho_n, rho_n is inputted.

      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)

      include 'param.escan_real'
      include 'mpif.h'
c Modificado por Txomin
c      real*8 xatom(3,matom),xatom_w(matom,2),w_tmp(matom)
      real*8 xatom(3,matom),xatom_w(matom,2),xatom_w_rec(matom,2)
      real*8 w_tmp(matom)

      real*8 AL(3,3),ALt(3,3)

      real*8 rho_atomIr(2000,mtype)

      real*8 occ_t(mtype)
c
      integer iiatom(mtype),iatom(matom),icore(mtype),numref(matom)
      integer itype_atom(matom)

      complex*16 cc

***************************************************
****  xatom(1),xatom(2),xatom(3) are the coord in unit of AL(3,3)
****  supercell edges
***************************************************

      common /comNL2/occ_t,iiatom,icore,numref
      common /comrho_atomIr/rho_atomIr

*******************************************************
*******************************************************
        do ia=1,natom
        do itype=1,ntype
        if(iiatom(itype).eq.iatom(ia)) itype_atom(ia)=itype
        enddo
        enddo
*******************************************************
        dx1cut=9/dsqrt(AL(1,1)**2+AL(2,1)**2+AL(3,1)**2)
        dx2cut=9/dsqrt(AL(1,2)**2+AL(2,2)**2+AL(3,2)**2)
        dx3cut=9/dsqrt(AL(1,3)**2+AL(2,3)**2+AL(3,3)**2)
        ddcut=4.5**2

c Modificado por Txomin
        xatom_w=0.d0
        xatom_w_rec=0.d0

        do 200 ii=1,nr_n

        jj=ii+(inode-1)*nr_n
        i=(jj-1)/(n2*n3)+1
        j=(jj-1-(i-1)*n2*n3)/n3+1
        k=jj-(i-1)*n2*n3-(j-1)*n3
        x1=(i-1.d0)/n1
        x2=(j-1.d0)/n2
        x3=(k-1.d0)/n3
        
        w_tmp=0.d0
        sum=0.d0

        do 100 ia=1,natom

        do 91 i1=-1,1
        dx1=xatom(1,ia)-x1+i1
        if(dabs(dx1).gt.dx1cut) goto 91
        do 92 i2=-1,1
        dx2=xatom(2,ia)-x2+i2
        if(dabs(dx2).gt.dx2cut) goto 92
        do 93 i3=-1,1
        dx3=xatom(3,ia)-x3+i3
        if(dabs(dx3).gt.dx3cut) goto 93
      
        dx=AL(1,1)*dx1+AL(1,2)*dx2+AL(1,3)*dx3
        dy=AL(2,1)*dx1+AL(2,2)*dx2+AL(2,3)*dx3
        dz=AL(3,1)*dx1+AL(3,2)*dx2+AL(3,3)*dx3

        dd=dx**2+dy**2+dz**2
        if(dd.gt.ddcut) goto 90
        d=dsqrt(dd)
        ir=d*200+1
        sum=sum+rho_atomIr(ir,itype_atom(ia))
        w_tmp(ia)=w_tmp(ia)+rho_atomIr(ir,itype_atom(ia))
90      continue
93      continue
92      continue
91      continue
100     continue

        do isp=1,islda
        sum1=rho_n(ii,isp)/sum
         do ia=1,natom
         if(w_tmp(ia).ne.0.d0) then
c Modificado por Txomin
c         xatom_w(ia,isp)=xatom_w(ia,isp)+sum1*w_tmp(ia)
         xatom_w_rec(ia,isp)=xatom_w(ia,isp)+sum1*w_tmp(ia)
         endif
         enddo
        enddo
         
200     continue
*****************************************
c Modificado por Txomin
c	call mpi_allreduce(xatom_w,xatom_w,matom*2,
	call mpi_allreduce(xatom_w,xatom_w_rec,matom*2,
     &  MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

        sum=0.d0
        do isp=1,islda
        do ia=1,natom
        sum=sum+xatom_w(ia,isp)
        enddo
        enddo
        sum=totNel/sum
        do isp=1,islda
        do ia=1,natom
        xatom_w(ia,isp)=xatom_w(ia,isp)*sum/occ_t(itype_atom(ia))
        enddo
        enddo

ccccccccccccccc
        if(inode.eq.1) then
        do isp=1,islda
        write(6,*) "***** iatom,xatom_w *******, isp=",isp
        do ia=1,natom
        write(6,3300) iatom(ia),xatom_w(ia,isp)
        enddo
        write(6,*) "**********************" 
        enddo
        endif
3300    format(i4,2x,f12.5)     


      return
cccccccccccccccccccccccccccccccccccccccc

      end
      
