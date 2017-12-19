      subroutine eigen_comp(ilocal,nline,
     &  vr,workr_n,kpt,Eref,AL,Ewind
     &  eigen,mxc)
****************************************
cc     Written by Meng Ye, December 20, 2017. 
cc     Copyright 2017 The Regents of the University of California
cc     The United States government retains a royalty free license in this work 
****************************************

****************************************
cc     First use Lanczos algorithm to make T=P^dagger * (H-Eref)**2 * P, where P is a ng_n*mx matrix with orthonormal columns, and T is a mx*mx real symmetric tridiagonal matrix.
cc     Then compute the eigen value and vector of T, to produce the effective eigen states in (Eref-Ewind,Eref+Ewind).
cc     These eigen states will store in ug_n, and their eigen value (H-Eref)**2 will store in eigen. The number of states in mxc.
****************************************

      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)

      include 'mpif.h'

      include 'param.escan_real'
***********************************************
      integer status(MPI_STATUS_SIZE)

      real*8 AL(3,3)
c       complex*16 workr_n(mg_nx)
      complex*16 t,workr_n(*)   ! original workr_n is of mr_n which is larger, xwjiang
      integer mxc,info

      real*8, allocatable, dimension(:) :: diag,subdiag,work
      real*8, allocatable, dimension(:,:) :: eigen_vec
      complex*16, allocatable, dimension (:) :: pg,pgh,pghh,qg
      complex*16, allocatable, dimension (:,:) :: pg_n

      common /com123b/m1,m2,m3,ngb,nghb,ntype,rrcut,msb
      common /comEk/Ek

      real*8 ran1

      ng_n=ngtotnod(inode,kpt)

      allocate(diag(mx))
      allocate(subdiag(mx))
      allocate(pg_n(ng_n,mx))
      allocate(pg(ng_n))
      allocate(pgh(ng_n))
      allocate(pghh(ng_n))
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


      do 4000 m=1,mx

       if (m.ne.1) then
         do i=1,ng_n
           t=pg(i)
           pg(i)=qg(i)/beta
           qg(i)=-t*beta
          enddo
        endif
************************************************
**** qg = qg + (H-Eref)**2 * pg
**** alpha = pg^* * qg
**** qg = qg - alpha * pg
**** beta = ||qg||
**** alpha and beta are the m-th element of diag and subdiag
************************************************
       call Hpsi_comp(pg,pgh,ilocal,vr,workr_n,kpt)
       do i=1,ng_n
       pgh(i)=pgh(i)-Eref*pg(i)
       enddo

       call Hpsi_comp(pgh,pghh,ilocal,vr,workr_n,kpt)
       do i=1,ng_n
       pghh(i)=pghh(i)-Eref*pgh(i)
       qg(i)=qg(i)+pghh(i)
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
       deallocate(pghh)
       deallocate(qg)
       allocate(eigen_vec(mx,mx))
       allocate(work(2*mx-2))

       call dstev('V',mx,diag,subdiag,eigen_vec,mx,work,info)

       deallocate(work)

       if (info.ne.0) stop

       mxc=0

       do m=1,mx
          Ediff=dsqrt(diag(m))
c         if (Ediff.le.Ewind) then
           do i=1,ng_n
             pg(i)=dcmplx(0.d0,0.d0)
             do ii=1,mx
               pg(i)=pg(i)+eigen_vec(ii,m)*pg_n(i,ii)
             enddo
           enddo

           call Hpsi_comp(pg,pgh,ilocal,vr,workr_n,kpt)
           do i=1,ng_n
           pgh(i)=pgh(i)-Eref*pg(i)
           enddo

           E0=0.d0
           s=0.d0
           do i=1,ng_n
             E0=E0+dreal(dconjg(pg(i))*pgh(i))
             s=s+cdabs(pg(i))**2
           enddo

           call global_sumr(E0)
           call global_sumr(s)
           E0=E0*vol
           s=s*vol

           E0=dabs(E0/s)

           if(inode.eq.1) then
           write(6,*) '*******************'
           write(6,*) 's=',s
           write(6,*) 'E0,eigen,diff=',E0,Ediff,dabs(E0-Ediff)
           write(6,*) 'coff mx,mx-1=',eigen_vec(mx,m),eigen_vec(mx-1,m)
           endif

c           mxc=mxc+1
c         endif
       enddo

       deallocate(diag)
       deallocate(subdiag)
       deallocate(pg)
       deallocate(pgh)
       deallocate(pg_n)
       deallocate(eigen_vec)

      return

      end

