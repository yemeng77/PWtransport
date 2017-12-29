      subroutine lanczos_comp(ilocal,niter,
     &  vr,workr_n,kpt,Emax,Emin)

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
      real*8 Emax,Emin,ran1

      real*8, allocatable, dimension(:) :: diag,subdiag,work
      real*8, allocatable, dimension(:,:) :: eigen_vec
      complex*16, allocatable, dimension (:) :: pg,pgh,qg
      complex*16, allocatable, dimension (:,:) :: pg_n

      common /com123b/m1,m2,m3,ngb,nghb,ntype,rrcut,msb
      common /comEk/Ek

      ng_n=ngtotnod(inode,kpt)

      allocate(diag(niter))
      allocate(subdiag(niter))
      allocate(pg_n(ng_n,niter))
      allocate(pg(ng_n))
      allocate(pgh(ng_n))
      allocate(qg(ng_n))

************************************************
**** generate random inital pg, then normalize it
**** qg = 0, beta = 1
************************************************

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

      beta=1.d0


      do 4000 m=1,niter

       if (m.ne.1) then
         do i=1,ng_n
           tmp=pg(i)
           pg(i)=qg(i)/beta
           qg(i)=-tmp*beta
          enddo
       endif
************************************************
**** qg = qg + H * pg
**** alpha = pg^* * qg
**** qg = qg - alpha * pg
**** beta = ||qg||
**** alpha and beta are the m-th element of diag and subdiag
************************************************
       call Hpsi_comp(pg,pgh,ilocal,vr,workr_n,kpt)
       do i=1,ng_n
       qg(i)=qg(i)+pgh(i)
       enddo

       alpha=0.d0
       do i=1,ng_n
       alpha=alpha+dreal(dconjg(pg(i))*qg(i))
       enddo
       call global_sumr(alpha)
       alpha=alpha*vol

       do i=1,ng_n
       qg(i)=qg(i)-alpha*pg(i)
       enddo

       beta=0.d0
       do i=1,ng_n
       beta=beta+cdabs(qg(i))**2
       enddo
       call global_sumr(beta)
       beta=dsqrt(beta*vol)

       if (beta.eq.0.d0) stop

       diag(m)=alpha
       subdiag(m)=beta
       do i=1,ng_n
       pg_n(i,m)=pg(i)
       enddo

4000  continue


***********************************************
      deallocate(qg)
      allocate(eigen_vec(niter,niter))
      allocate(work(2*niter-2))

      call dstev('V',niter,diag,subdiag,eigen_vec,niter,work,info)

      deallocate(work)

      if (info.ne.0) stop

      Emax=-1.d20
      Emin=1.d20

      do m=1,niter
       E0=diag(m)
        do i=1,ng_n
         pg(i)=dcmplx(0.d0,0.d0)
         do ii=1,niter
         pg(i)=pg(i)+eigen_vec(ii,m)*pg_n(i,ii)
         enddo
        enddo

       call Hpsi_comp(pg,pgh,ilocal,vr,workr_n,kpt)
       res=0.d0
       do i=1,ng_n
       pgh(i)=pgh(i)-E0*pg(i)
       res=res+cdabs(pgh(i))**2
       enddo
       call global_sumr(res)
       res=dsqrt(res*vol)

       if (Emax.lt.E0+res) Emax=E0+res
       if (Emin.gt.E0-res) Emin=E0-res

      enddo

      if(inode.eq.1) then
      write(6,*) '*******************'
      write(6,*) 'Emax=',Emax*27.211396d0,'eV'
      write(6,*) 'Emin=',Emin*27.211396d0,'eV'
      write(6,*) '*******************'
      endif

      deallocate(diag)
      deallocate(subdiag)
      deallocate(pg)
      deallocate(pgh)
      deallocate(pg_n)
      deallocate(eigen_vec)

      return

      end

