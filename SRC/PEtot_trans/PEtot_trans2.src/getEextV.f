      subroutine getEextV(fatom,xatom,AL,
     &    ityatom,E_IVext)
***************************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
******************************************

ccccc This program needs some fine tunning.

       use fft_data
       use load_data
       use data

      implicit double precision (a-h,o-z)
      include 'param.escan_real'
      include "mpif.h"

      real*8 AL(3,3),ALI2(3,3),tmp(3)
      real*8 xatom(3,matom),fatom(3,matom)
      real*8 zatom(mtype)
      integer ityatom(matom)
      integer, allocatable, dimension(:) :: ncount_tmp

      common /comzatom/zatom

cccccccccccccccccccccccccccccccccccccccccccccc

      pi=4*datan(1.d0)
      ALI2=AL
      tmp=1.d0
      call gaussj(ALI2,3,3,tmp,1,1)
ccccccc Note, ALI2 is diff from ALI
ccccccc AL(i1,j)*ALI2(j,i2)=\delta_i1,i2

      
cccccc this is a slow, but a easy way to do the 
cccccc calculation. All the processor are doing the 
cccccc same thing. 
      E_IVext=0.d0
      do 1000 ia=1,natom

      x1=xatom(1,ia)
      y1=xatom(2,ia)
      z1=xatom(3,ia)
      if(x1.lt.0.d0) x1=x1+1.d0
      if(x1.ge.1.d0) x1=x1-1.d0
      if(y1.lt.0.d0) y1=y1+1.d0
      if(y1.ge.1.d0) y1=y1-1.d0
      if(z1.lt.0.d0) z1=z1+1.d0
      if(z1.ge.1.d0) z1=z1-1.d0
      i1=x1*n1      ! i1 between 0 and n1-1
      j1=y1*n2
      k1=z1*n3
      ii000=i1*n2*n3+j1*n3+k1+1
      ii100=mod(i1+1,n1)*n2*n3+j1*n3+k1+1
      ii010=i1*n2*n3+mod(j1+1,n2)*n3+k1+1
      ii001=i1*n2*n3+j1*n3+mod(k1+1,n3)+1
      ii110=mod(i1+1,n1)*n2*n3+mod(j1+1,n2)*n3+k1+1
      ii101=mod(i1+1,n1)*n2*n3+j1*n3+mod(k1+1,n3)+1
      ii011=i1*n2*n3+mod(j1+1,n2)*n3+mod(k1+1,n3)+1
      ii111=mod(i1+1,n1)*n2*n3+mod(j1+1,n2)*n3+
     &    mod(k1+1,n3)+1
      call bcast_vr_ii(vt000,ii000)
      call bcast_vr_ii(vt100,ii100)
      call bcast_vr_ii(vt010,ii010)
      call bcast_vr_ii(vt001,ii001)
      call bcast_vr_ii(vt110,ii110)
      call bcast_vr_ii(vt101,ii101)
      call bcast_vr_ii(vt011,ii011)
      call bcast_vr_ii(vt111,ii111)
      f1=i1+1-x1*n1
      f2=j1+1-y1*n2
      f3=k1+1-z1*n3
      v_av=vt000*f1*f2*f3+vt100*(1-f1)*f2*f3+
     &  vt010*f1*(1-f2)*f3+vt001*f1*f2*(1-f3)+
     &  vt110*(1-f1)*(1-f2)*f3+vt101*(1-f1)*f2*(1-f3)+
     &  vt011*f1*(1-f2)*(1-f3)+vt111*(1-f1)*(1-f2)*(1-f3)  

      dv1=(vt100-vt000)*f2*f3+(vt110-vt010)*(1-f2)*f3+
     & (vt101-vt001)*f2*(1-f3)+(vt111-vt011)*(1-f2)*(1-f3)

      dv2=(vt010-vt000)*f1*f3+(vt110-vt100)*(1-f1)*f3+
     & (vt011-vt001)*f1*(1-f3)+(vt111-vt101)*(1-f1)*(1-f3)

      dv3=(vt001-vt000)*f1*f2+(vt101-vt100)*(1-f1)*f2+
     & (vt011-vt010)*f1*(1-f2)+(vt111-vt110)*(1-f1)*(1-f2)

      dv1=dv1*n1
      dv2=dv2*n2
      dv3=dv3*n3

      dvx=dv1*ALI2(1,1)+dv2*ALI2(2,1)+dv3*ALI2(3,1)
      dvy=dv1*ALI2(1,2)+dv2*ALI2(2,2)+dv3*ALI2(3,2)
      dvz=dv1*ALI2(1,3)+dv2*ALI2(2,3)+dv3*ALI2(3,3)

      ch=zatom(ityatom(ia))
      E_IVext=E_IVext-ch*v_av   ! ch>0, and nuclien has negative charge in this convention 
      fatom(1,ia)=fatom(1,ia)-ch*dvx       ! fatom is the energy derivitive, points to Etot increase dir.
      fatom(2,ia)=fatom(2,ia)-ch*dvy
      fatom(3,ia)=fatom(3,ia)-ch*dvz
1000  continue

      return
      contains



      subroutine bcast_vr_ii(vt_tmp,ii_tmp)
      implicit double precision(a-h,o-z)
      real*8 vt_tmp
      integer ii_tmp

      isend_node=(ii_tmp-1)/nr_n
      if(isend_node.eq.inode-1) then
      vt_tmp=vext_n(ii_tmp-isend_node*nr_n)
      endif
      call mpi_bcast(vt_tmp,1,MPI_REAL8,isend_node,MPI_COMM_WORLD,ierr)
      return
      end subroutine bcast_vr_ii

      end



  
