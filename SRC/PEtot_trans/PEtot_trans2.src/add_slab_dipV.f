      subroutine add_slab_dipV(iflag,
     & xatom,ityatom,totNel,AL,islda,dv_jump,dv_jump_curr,
     & E_IVextau,E_rhoVextau,E_rhoVextau2,fatom,vion_p)
***************************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
******************************************

**************************************************
***   This program adds the automatically generated dipole field vextau into vion_p
***   It also calculates: E_rhoVextau=rho*vextau(current-rho)
***                       E_rhoVextau2=rho*vextau(update, output dv_jump,vion_p)
***                       E_IVextau= nuclei_ion*vextau(current_rho)
***                       add force term into fatom                       
***   They remove the previous additions in vion_p and fatom (from previous dv_jump)
***   iflag=0, prev. dv_jump=0, initial call
***   iflag=1, update dv_jump, perhaps using some mixing scheme

       use fft_data
       use load_data
       use data

      implicit double precision (a-h,o-z)
      include 'param.escan_real'
      include "mpif.h"

      real*8 vion_p(mr_n)
      real*8 AL(3,3),ALI2(3,3),tmp(3)
      real*8 xatom(3,matom),fatom(3,matom)
      real*8 zatom(mtype)
      integer ityatom(matom)
      integer, allocatable, dimension(:) :: ncount_tmp
      integer icount_jump

      common /comzatom/zatom

      common /comVext/ivext_in,ivext_dir,xvext_c,dv_mix,
     &   nite_mix,dv_jump_init

      save icount_jump
cccccccccccccccccccccccccccccccccccccccccccccc
      if(iflag.eq.0) then
      dv_jump=0.d0
      dv_jump_old=0.d0
      icount_jump=0
      endif

      pi=4*datan(1.d0)

      sum=0.d0
      sum1=0.d0
      do ii=1,nr_n
      jj=ii+(inode-1)*nr_n
      i=(jj-1)/(n2*n3)+1
      j=(jj-1-(i-1)*n2*n3)/n3+1
      k=jj-(i-1)*n2*n3-(j-1)*n3
      if(ivext_dir.eq.1)  xt=(i-1.d0)/n1
      if(ivext_dir.eq.2)  xt=(j-1.d0)/n2
      if(ivext_dir.eq.3)  xt=(k-1.d0)/n3
      if(xt.gt.xvext_c)   xt=xt-1.d0
      if(islda.eq.1) rho_t=rho_n(ii,1)
      if(islda.eq.2) rho_t=rho_n(ii,1)+rho_n(ii,2)
      sum=sum+rho_t*xt
      sum1=sum1+rho_t
      enddo

      call global_sumr(sum)
      call global_sumr(sum1)
     
      sum=sum/sum1*totNel

      do i=1,natom
      xt=xatom(ivext_dir,i)
      if(xt.gt.xvext_c)   xt=xt-1.d0
      sum=sum-zatom(ityatom(i))*xt
      enddo

      vol=al(3,1)*(al(1,2)*al(2,3)-al(1,3)*al(2,2))
     &     +al(3,2)*(al(1,3)*al(2,1)-al(1,1)*al(2,3))
     &     +al(3,3)*(al(1,1)*al(2,2)-al(1,2)*al(2,1))
       vol=dabs(vol)

       j=ivext_dir
       h=AL(1,j)*ALI(1,j)+AL(2,j)*ALI(2,j)+AL(3,j)*ALI(3,j)
       h=h/dsqrt(ALI(1,j)**2+ALI(2,j)**2+ALI(3,j)**2)
       Dipole=sum*h**2/vol

       dv_jump_old=dv_jump
       dv_jump=4*pi*Dipole

*******************************************
       icount_jump=icount_jump+1

       if(iflag.eq.0) then
       if(dv_jump_init.lt.1000) then
       dv_jump_update=dv_jump_init
       else
       dv_jump_update=dv_jump
       endif
       endif
       
       if(iflag.eq.1) then
       if(icount_jump.ge.nite_mix+2.and.
     &    mod(icount_jump-2,nite_mix).eq.0) then 
       dv_jump_update=dv_jump*dv_mix+dv_jump_old*(1-dv_mix)
       else
       dv_jump_update=dv_jump_old
       endif
       endif
       

*******************************************
       
       E_rhoVextau=0.d0
       E_rhoVextau2=0.d0
       do ii=1,nr_n
       jj=ii+(inode-1)*nr_n

       i=(jj-1)/(n2*n3)+1
       j=(jj-1-(i-1)*n2*n3)/n3+1
       k=jj-(i-1)*n2*n3-(j-1)*n3
       if(ivext_dir.eq.1)  xt=(i-1.d0)/n1
       if(ivext_dir.eq.2)  xt=(j-1.d0)/n2
       if(ivext_dir.eq.3)  xt=(k-1.d0)/n3
       if(xt.lt.xvext_c)  then
       vt=dv_jump*(xt-xvext_c+0.5d0)
       vt2=dv_jump_update*(xt-xvext_c+0.5d0)
       dvt=(dv_jump_update-dv_jump_old)*(xt-xvext_c+0.5d0)
       else
       vt=dv_jump*(xt-xvext_c-0.5d0)
       vt2=dv_jump_update*(xt-xvext_c-0.5d0)
       dvt=(dv_jump_update-dv_jump_old)*(xt-xvext_c-0.5d0)
       endif

       vion_p(ii)=vion_p(ii)+dvt
       if(islda.eq.1) rho_t=rho_n(ii,1)
       if(islda.eq.2) rho_t=rho_n(ii,1)+rho_n(ii,2)
       E_rhoVextau=E_rhoVextau+rho_t*vt
       E_rhoVextau2=E_rhoVextau2+rho_t*vt2
       enddo

       call global_sumr(E_rhoVextau)
       call global_sumr(E_rhoVextau2)
 
       E_rhoVextau=E_rhoVextau*vol/(n1*n2*n3)
       E_rhoVextau2=E_rhoVextau2*vol/(n1*n2*n3)



*********************************************************************

      ALI2=AL
      tmp=1.d0
      call gaussj(ALI2,3,3,tmp,1,1)

ccccccc Note, ALI2 is diff from ALI
ccccccc AL(i1,j)*ALI2(j,i2)=\delta_i1,i2

      
cccccc this is a slow, but a easy way to do the 
cccccc calculation. All the processor are doing the  same thing.

      E_IVextau=0.d0
      dv1=0.d0
      dv2=0.d0
      dv3=0.d0
      do 1000 ia=1,natom

       xt=xatom(ivext_dir,ia)
       if(xt.lt.xvext_c)  then
       vt=dv_jump*(xt-xvext_c+0.5d0)
       else
       vt=dv_jump*(xt-xvext_c-0.5d0)
       endif

      if(ivext_dir.eq.1) dv1=dv_jump-dv_jump_old
      if(ivext_dir.eq.2) dv2=dv_jump-dv_jump_old
      if(ivext_dir.eq.3) dv3=dv_jump-dv_jump_old

      dvx=dv1*ALI2(1,1)+dv2*ALI2(2,1)+dv3*ALI2(3,1)
      dvy=dv1*ALI2(1,2)+dv2*ALI2(2,2)+dv3*ALI2(3,2)
      dvz=dv1*ALI2(1,3)+dv2*ALI2(2,3)+dv3*ALI2(3,3)

      ch=zatom(ityatom(ia))
      E_IVextau=E_IVextau-ch*vt
      fatom(1,ia)=fatom(1,ia)-ch*dvx       ! fatom is the energy derivitive, points to Etot increase dir.
      fatom(2,ia)=fatom(2,ia)-ch*dvy
      fatom(3,ia)=fatom(3,ia)-ch*dvz
1000  continue


      if(inode.eq.1) then
      write(6,*) "The current actual dv_jump=",dv_jump
      write(6,*) "The updated        dv_jump=",dv_jump_update
      endif

      dv_jump_curr=dv_jump
      dv_jump=dv_jump_update
      

      return
      end

