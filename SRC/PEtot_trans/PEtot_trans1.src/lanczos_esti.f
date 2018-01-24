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
      complex*16 workr_n(mr_n),sumdum(mref,natom),tmp_complex
      real*8 Emax,Emin,Ebound(2),Etmp
      complex*16 qg(mg_nx)

      real*8, allocatable, dimension(:) :: diag,subdiag,work
      real*8, allocatable, dimension(:,:) :: eigen_vec
      complex*16, allocatable, dimension (:,:) :: pg_n,pgh_n,sug_m

      real*8  Dij0(32,32,mtype),Qij(32,32,mtype)
      integer isNLa(9,matom),ipsp_type(mtype),ipsp_all

      common /comisNLa/isNLa,Dij0,Qij,ipsp_all,ipsp_type

      allocate(diag(niter))
      allocate(subdiag(niter))
      allocate(pg_n(mg_nx,niter))
      allocate(pgh_n(mg_nx,niter))
      if(ipsp_all.eq.1) then
      allocate(sug_m(1,niter))
      else
      allocate(sug_m(mg_nx,niter))
      endif

      ng_n=ngtotnod(inode,kpt)

************************************************
**** pg is normalized random wavefunction
**** qg = 0, beta = 1
************************************************
      pg_n(:,1)=ug_n_bp(:,1)
      qg=dcmplx(0.d0,0.d0)

      s=0.d0
      do i=1,ng_n
      s=s+cdabs(pg_n(i,1))**2
      enddo
      call global_sumr(s)
      s=1.d0/dsqrt(s*vol)

      do i=1,ng_n
      pg_n(i,1)=s*pg_n(i,1)
      enddo
      beta=1.d0

      do 4000 m=1,niter

       if (m.ne.1) then
         do i=1,ng_n
           pg_n(i,m)=qg(i)/beta
           qg(i)=-pg_n(i,m-1)*beta
          enddo
       endif
************************************************
**** qg = qg + H * pg
**** alpha = pg^* * qg
**** qg = qg - alpha * pg
**** beta = ||qg||
**** alpha and beta are the m-th element of diag and subdiag
************************************************
       call Hpsi_comp_AllBandBP(pg_n(1,m),pgh_n(1,m),1,
     &      ilocal,vr,workr_n,kpt,1,
     &      sug_m(1,m),sumdum,iislda)
       do i=1,ng_n
       qg(i)=qg(i)+pgh_n(i,m)
       enddo

       alpha=0.d0
       do i=1,ng_n
       alpha=alpha+dreal(dconjg(pg_n(i,m))*qg(i))
       enddo
       call global_sumr(alpha)
       alpha=alpha*vol

       do i=1,ng_n
       qg(i)=qg(i)-alpha*pg_n(i,m)
       enddo

       beta=0.d0
       do i=1,ng_n
       beta=beta+cdabs(qg(i))**2
       enddo
       call global_sumr(beta)
       beta=dsqrt(beta*vol)

       if (beta.eq.0.d0) then
        if(inode_tot.eq.1) write(6,*) 'beta equals 0, Lanczos algorithm 
     &     stop, m=', m
       call mpi_abort(MPI_COMM_WORLD,ierr)
       endif

       diag(m)=alpha
       subdiag(m)=beta

4000  continue


***********************************************
      allocate(eigen_vec(niter,niter))
      allocate(work(2*niter-2))
      call dstev('V',niter,diag,subdiag,eigen_vec,niter,work,info)
      if (beta.eq.0.d0) then
        if(inode_tot.eq.1) write(6,*) "Somthing is wrong with 
     &    diagonalization of tridiagonal matrix. info",info
        call mpi_abort(MPI_COMM_WORLD,ierr)
      endif
      deallocate(subdiag)
      deallocate(work)

      if (info.ne.0) stop

      Emax=-1.d20
      Emin=1.d20

      do m=1,niter
        s=0.d0
        if(ipsp_all.eq.2) then
          do i=1,ng_n
            tmp_complex=dcmplx(0.d0,0.d0)
            do ii=1,niter
            tmp_complex=tmp_complex+
     &              (pgh_n(i,ii)-diag(m)*sug_m(i,ii))*eigen_vec(ii,m)
            enddo
            s=s+cdabs(tmp_complex)**2
          enddo
        else 
          do i=1,ng_n
            tmp_complex=dcmplx(0.d0,0.d0)
            do ii=1,niter
            tmp_complex=tmp_complex+
     &              (pgh_n(i,ii)-diag(m)*pg_n(i,ii))*eigen_vec(ii,m)
            enddo
            s=s+cdabs(tmp_complex)**2
          enddo
        endif
        call global_sumr(s)
        s=dsqrt(s*vol)

        Etmp=diag(m)+s
        if (Emax.lt.Etmp) Emax=Etmp
        Etmp=diag(m)-s
        if (Emin.gt.Etmp) Emin=Etmp

      enddo

      deallocate(pg_n)
      deallocate(pgh_n)
      deallocate(sug_m)
      deallocate(diag)
      deallocate(eigen_vec)

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      call mpi_allreduce(Emin,Etmp,1,
     &     MPI_REAL8,MPI_MAX,MPI_COMM_B2,ierr)
      Ebound(1)=Etmp
      call mpi_allreduce(Emax,Etmp,1,
     &     MPI_REAL8,MPI_MIN,MPI_COMM_B2,ierr)
      Ebound(2)=Etmp


      return

      end

