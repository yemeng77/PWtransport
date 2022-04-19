      subroutine lanczos_esti(ilocal,niter,
     &  vr,workr_n,kpt,Eref)

****************************************
cc     Written by Meng Ye, December 28, 2017. 
cc     Copyright 2017 The Regents of the University of California
cc     The United States government retains a royalty free license in this work 
****************************************

****************************************
cc     Use niter-step Lanczos algorithm to estimate the bound of eigenvalues of H
****************************************

      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)

      include 'mpif.h'

      include 'param.escan_real'
***********************************************
      integer status(MPI_STATUS_SIZE)

c      complex*16 workr_n(mg_nx)
      complex*16 tmp,workr_n(*)   ! original workr_n is of mr_n which is larger, xwjiang
      real*8 Emax,Emin,Eref,ran1

      real*8 alpha(niter), beta(niter), diag(niter), subdiag(niter), work(2*niter-2)
      real*8 eigen_vec(niter, niter)
      real*8 prec(mg_nx)
      complex*16 pg(mg_nx),ug(mg_nx),ugh(mg_nx),ughh(mg_nx),qg(mg_nx)
      complex*16 pg_n(mg_nx, niter)

      common /com123b/m1,m2,m3,ngb,nghb,ntype,rrcut,msb
      common /comEk/Ek

      ng_n=ngtotnod(inode,kpt)

************************************************
**** generate random inital pg, then normalize it
**** qg = 0, beta = 1
************************************************
      do i=1,ng_n
         x=gkk_n(i,kpt)/Ek
         y=27.d0+x*(18.d0+x*(12.d0+x*8.d0))
         prec(i)=1.d0/dsqrt(1.d0+16.d0*x**4/y)
      enddo
 
      iranm=-2291-inode*3651
      x=ran1(iranm)
      do i=1,ng_n
      x=ran1(iranm)
      y=ran1(iranm)
      pg(i)=dcmplx(x-0.5d0,y-0.5d0)
      qg(i)=dcmplx(0.d0,0.d0)
      enddo

      s=0.d0
      do i=1,ng_n
      s=s+cdabs(pg(i))**2
      enddo
       
      call global_sumr(s)
      s=1.d0/dsqrt(s*vol)

      do i=1,ng_n
      pg(i)=s*pg(i)
      enddo

      beta_tmp=1.d0
      do 4000 m=1,niter

       if (m.ne.1) then
         do i=1,ng_n
           tmp=pg(i)
           pg(i)=qg(i)/beta_tmp
           qg(i)=-tmp*beta_tmp
          enddo
       endif
************************************************
**** qg = qg + H * pg
**** alpha = pg^* * qg
**** qg = qg - alpha * pg
**** beta = ||qg||
**** alpha and beta are the m-th element of diag and subdiag
************************************************
       do i=1,ng_n
         ug(i)=prec(i)*pg(i)
       enddo
       call Hpsi_comp(ug,ugh,ilocal,vr,workr_n,kpt)
       do i=1,ng_n
       ugh(i)=ugh(i)-Eref*ug(i)
       enddo
       call Hpsi_comp(ugh,ughh,ilocal,vr,workr_n,kpt)
       do i=1,ng_n
       ughh(i)=ughh(i)-Eref*ugh(i)
       qg(i)=qg(i)+prec(i)*ughh(i)
       enddo

       alpha_tmp=0.d0
       do i=1,ng_n
       alpha_tmp=alpha_tmp+dreal(dconjg(pg(i))*qg(i))
       enddo
       call global_sumr(alpha_tmp)
       alpha_tmp=alpha_tmp*vol

       do i=1,ng_n
       qg(i)=qg(i)-alpha_tmp*pg(i)
       enddo

       beta_tmp=0.d0
       do i=1,ng_n
       beta_tmp=beta_tmp+cdabs(qg(i))**2
       enddo
       call global_sumr(beta_tmp)
       beta_tmp=dsqrt(beta_tmp*vol)

       if (beta_tmp.eq.0.d0) stop

       alpha(m)=alpha_tmp
       beta(m)=beta_tmp
       do i=1,ng_n
       pg_n(i,m)=pg(i)
       enddo
       
4000  continue


***********************************************
      diag=alpha
      subdiag=beta
      call dstev('V',niter,diag,subdiag,eigen_vec,niter,work,info)
      if (info.ne.0) stop

      if(inode.eq.1) write(6,*) 'eigen,err'
      do m=1,niter
       E0=diag(m)
        do i=1,ng_n
         pg(i)=dcmplx(0.d0,0.d0)
         do ii=1,niter
         pg(i)=pg(i)+eigen_vec(ii,m)*pg_n(i,ii)
         enddo
        enddo

       res=0.d0
       do i=1,ng_n
         ug(i)=prec(i)*pg(i)
       enddo
       call Hpsi_comp(ug,ugh,ilocal,vr,workr_n,kpt)
       do i=1,ng_n
         ugh(i)=ugh(i)-Eref*ug(i)
       enddo
       call Hpsi_comp(ugh,ughh,ilocal,vr,workr_n,kpt)
       do i=1,ng_n
          ughh(i)=prec(i)*(ughh(i)-Eref*ugh(i))-E0*pg(i)
          res=res+cdabs(ughh(i))**2
       enddo
       call global_sumr(res)
       res=dsqrt(res*vol)
       if(inode.eq.1) write(6,*) E0,res
      enddo

      return

      end

