      subroutine lanczos_esti(ilocal,niter,
     &  vr,workr_n,kpt,iislda,Ebound)

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

      real*8 vr(mr_n)
      complex*16 tmp,workr_n(mr_n),sumdum(mref,natom)
      real*8 Emax,Emin,Ebound(2),Etmp
      complex*16 pg(mg_nx),pgh(mg_nx),qg(mg_nx),swg(mg_nx)

      real*8, allocatable, dimension(:) :: diag,subdiag,work
      real*8, allocatable, dimension(:,:) :: eigen_vec
      complex*16, allocatable, dimension (:,:) :: pg_n,pgh_n

      real*8  Dij0(32,32,mtype),Qij(32,32,mtype)
      integer isNLa(9,matom),ipsp_type(mtype),ipsp_all

      common /comisNLa/isNLa,Dij0,Qij,ipsp_all,ipsp_type

      allocate(diag(niter))
      allocate(subdiag(niter))
      allocate(pg_n(niter,ng_n))
      allocate(pgh_n(niter,ng_n))

************************************************
**** pg is normalized random wavefunction
**** qg = 0, beta = 1
************************************************
      pg=ug_n_BP(1,1)
      qg=dcmplx(0.d0,0.d0)

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
       call Hpsi_comp_AllBandBP(pg,pgh,1,
     &      ilocal,vr,workr_n,kpt,1,
     &      swg,sumdum,iislda)
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

       if (beta.eq.0.d0) then
        if(inode_tot.eq.1) write(6,*) 'beta equals 0, stop, m=', m
       call mpi_abort(MPI_COMM_WORLD,ierr)
       endif

       diag(m)=alpha
       subdiag(m)=beta
       do i=1,ng_n
       pg_n(m,i)=pg(i)
       pgh_n(m,i)=pgh(i)
       enddo

4000  continue


***********************************************
      allocate(eigen_vec(niter,niter))
      allocate(work(2*niter-2))

      call dstev('V',niter,diag,subdiag,eigen_vec,niter,work,info)

      deallocate(work)

      if (info.ne.0) stop

      Emax=-1.d20
      Emin=1.d20

      do m=1,niter
        E0=diag(m)
        res=0.d0
        do i=1,ng_n
         tmp=dcmplx(0.d0,0.d0)
         do ii=1,niter
         tmp=tmp+eigen_vec(ii,m)*(pgh_n(ii,i)-E0*pg_n(ii,i))
         enddo
         res=res+cdabs(tmp)**2
        enddo
        call global_sumr(res)
        res=dsqrt(res*vol)

        Etmp=E0+res
        if (Emax.lt.Etmp) Emax=Etmp
        Etmp=E0-res
        if (Emin.gt.Etmp) Emin=Etmp

      enddo

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      call mpi_allreduce(Emin,Etmp,1,
     &  MPI_REAL8,MPI_MAX,MPI_COMM_B2,ierr)
      Ebound(1)=Etmp
      call mpi_allreduce(Emax,Etmp,1,
     &  MPI_REAL8,MPI_MIN,MPI_COMM_B2,ierr)
      Ebound(2)=Etmp

      deallocate(diag)
      deallocate(subdiag)
      deallocate(pg_n)
      deallocate(pgh_n)
      deallocate(eigen_vec)

      return

      end

