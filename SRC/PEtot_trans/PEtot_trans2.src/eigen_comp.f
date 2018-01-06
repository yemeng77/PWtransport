      subroutine eigen_comp(ilocal,nline,mp,tol,
     &  vr,workr_n,kpt,Ewind,eigen,mxc)
****************************************
cc     Written by Meng Ye, December 28, 2017. 
cc     Copyright 2017 The Regents of the University of California
cc     The United States government retains a royalty free license in this work 
****************************************

****************************************
cc     Use the Chebyshev filter to calculate the interior eigen states near Eref
cc     See J. Comp. Phys. 325 (2016) 226–243
****************************************

      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)

      include 'mpif.h'

      include 'param.escan_real'
***********************************************
      integer status(MPI_STATUS_SIZE)

      real*8 Ewind(2),eigen(mst),Ebound(2)
c       complex*16 workr_n(mg_nx)
      complex*16 workr_n(*)   ! original workr_n is of mr_n which is larger, xwjiang
      integer mp,mxc,mstatus

      real*8 Emax,Emin,alpha,beta,c0,g0

      real*8, allocatable, dimension (:) :: cn,gn,E_st,err_st
      complex*16, allocatable, dimension (:,:) :: xg_n,yg_n,zg_n

      common /com123b/m1,m2,m3,ngb,nghb,ntype,rrcut,msb
      common /comEk/Ek

      ng_n=ngtotnod(inode,kpt)

      niter_lan=25
      call lanczos_comp(ilocal,niter_lan,vr,workr_n,kpt,Ebound)

      Emin=Ebound(1)
      Emax=Ebound(2)

****************************************
cc    map the x in [Emin,Emax] to [-1,1] using alpha*x+beta
cc    alpha=2/(Emax-Emin), beta=(Emin+Emax)/(Emin-Emax)
****************************************      
      alpha=2.0d0/(Emax-Emin)
      beta=(Emin+Emax)/(Emin-Emax)

      allocate(cn(mp))
      allocate(gn(mp))
      allocate(E_st(mx))
      allocate(err_st(mx))
      allocate(xg_n(ng_n,mx))
      allocate(yg_n(ng_n,mx))
      allocate(zg_n(ng_n,mx))

****************************************
cc    W(x) is a windows function, which is 1 in [Ewind(1),Ewind(2)] and 0 outside
cc    W(x)\approx\sum_0^mp gn*cn*Tn(alpha*x+beta)
cc    The filter: ug=W(H)ug=\sum_0^mp gn*cn*Tn(alpha*H+beta)ug
****************************************
      Emin=Ewind(1)
      Emax=Ewind(2)
      pi=4.0d0*datan(1.0d0)
      cn0=(dacos(alpha*Emax+beta)-dacos(alpha*Emin+beta))/pi
      gn0=1.0d0

      do n=1,mp
      cn(n)=2.0d0*(dsin(n*dacos(alpha*Emax+beta))
     &     -dsin(n*dacos(alpha*Emin+beta)))/(n*pi)
      gn(n)=(dsin(n*pi/(mp+1))/(n*pi/(mp+1)))**2
      enddo

      do 3000 nint=1,nline
****************************************
cc     xg=T0(alpha*H+beta)ug=ug
cc     yg=T1(alpha*H+beta)ug=(alpha*H+beta)ug=(alpha*H+beta)xg
****************************************
       xg_n=ug_n
       do m=1,mx
         call Hpsi_comp(xg_n(1,m),yg_n(1,m),ilocal,vr,workr_n,kpt)
         do i=1,ng_n
           yg_n(i,m)=alpha*yg_n(i,m)+beta*xg_n(i,m)
           ug_n(i,m)=gn0*cn0*xg_n(i,m)+gn(1)*cn(1)*yg_n(i,m)
         enddo
       enddo

       do n=2,mp
****************************************
cc     xg=T(n-2)(alpha*H+beta)ug
cc     yg=T(n-1)(alpha*H+beta)ug
cc     zg=Tn(alpha*H+beta)ug=2*(alpha*H+beta)T(n-1)ug-T(n-2)ug=2*(alpha*H+beta)yg-xg
cc     yg->xg zg->yg
****************************************
       do m=1,mx
         call Hpsi_comp(yg_n(1,m),zg_n(1,m),ilocal,vr,workr_n,kpt)
         do i=1,ng_n
           zg_n(i,m)=2.0d0*(alpha*zg_n(i,m)+beta*yg_n(i,m))-xg_n(i,m)
           ug_n(i,m)=ug_n(i,m)+gn(n)*cn(n)*zg_n(i,m)
         enddo
       enddo
       xg_n=yg_n
       yg_n=zg_n
       enddo

       do m=1,mx
         call orth_comp(ug_n(1,m),ug_n,m-1,1,kpt)
       enddo
       call diag_comp(ilocal,E_st,err_st,vr,workr_n,kpt)

****************************************
cc     Exit if err_st<=tol for all Ritz values in the target interval.
****************************************
       mxc=0
       mstatus=0

       do m=1,mx
       if (E_st(m).ge.Ewind(1).and.E_st(m).le.Ewind(2)) then
         mxc=mxc+1
         if (err_st(m).gt.tol) then
           mstatus=1
           goto 100
         endif
       endif
       enddo
100    continue

       if(mstatus.eq.(0).and.mxc.gt.(0)) goto 3001

3000  continue
3001  continue

      mxc=0
      zg_n=ug_n

      do m=1,mx
      if (E_st(m).ge.Ewind(1).and.E_st(m).le.Ewind(2)) then
      mxc=mxc+1
      eigen(mxc)=E_st(m)
      ug_n(:,mxc)=zg_n(:,m)
      endif
      enddo

      deallocate(cn)
      deallocate(gn)
      deallocate(E_st)
      deallocate(err_st)
      deallocate(xg_n)
      deallocate(yg_n)
      deallocate(zg_n)

***********************************
      if(inode.eq.1) then
       write(6,*) "*********************************"
       write(6,*) "eigen energies, in eV"
       write(6,101) (eigen(i)*27.211396d0, i=1,mxc)
       write(6,*) "*********************************"
      endif

101   format(5(f12.8,1x))

      return

      end

