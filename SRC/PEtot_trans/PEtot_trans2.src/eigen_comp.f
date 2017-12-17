      subroutine eigen_comp(ilocal,nline,
     &  vr,workr_n,kpt,Eref,AL,Ewind
     &  eigen,mxc)
****************************************
cc     Written by Meng Ye, December 20, 2017. 
cc     Copyright 2017 The Regents of the University of California
cc     The United States government retains a royalty free license in this work 
****************************************

****************************************
cc     First use Lanczos algorithm to make T=V^dagger * (H-Eref)**2 * V, where V is a ng_n*mx matrix with orthonormal columns, and T is a mx*mx real symmetric tridiagonal matrix.
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
      complex*16 workr_n(*)   ! original workr_n is of mr_n which is larger, xwjiang
      integer mxc,info

      real*8, allocatable, dimension(:) :: diag,subdiag,work
      real*8, allocatable, dimension(:,:) :: eigen_vec
      complex*16, allocatable, dimension (:) :: vg,vgh,wg
      complex*16, allocatable, dimension (:,:) :: vg_n

      common /com123b/m1,m2,m3,ngb,nghb,ntype,rrcut,msb
      common /comEk/Ek

      real*8 ran1

      ng_n=ngtotnod(inode,kpt)

      allocate(diag(mx))
      allocate(subdiag(mx))
      allocate(vg_n(ng_n,mx))
      allocate(vg(ng_n))
      allocate(vgh(ng_n))
      allocate(wg(ng_n))

************************************************
**** generate random inital vg, then normalize it 
************************************************
      iislda=1 

      iranm=-4513*(iislda-1)-1371*kpt-5616*inode
      do i=1,ng_n
      x=ran1(iranm)
      y=ran1(iranm)
      vg(i)=dcmplx(x-0.5d0,y-0.5d0)
      enddo

      s=0.d0
      do i=1,ng_n
      s=s+cdabs(vg(i))**2
      enddo
       
      call global_sumr(s)
      s=1.d0/dsqrt(s*vol)

      do i=1,ng_n
      vg(i)=s*vg(i)
      enddo
      

      do 4000 m=1,mx

       if (m.gt.1) then
************************************************
**** beta = ||wg||
**** vg = wg / beta
**** beta is the m-1 th element of subdiag
**** vg is the m th row of vg_n
************************************************
       beta=0.d0

       do i=1,ng_n
       beta=beta+cdabs(wg(i))**2
       enddo

       call global_sumr(beta)

       if (beta.eq.0.d0) stop

       beta=dsqrt(beta*vol)
       subdiag(m-1)=beta
       beta=1.d0/beta

       do i=1,ng_n
       vg(i)=beta*wg(i)
       enddo
       endif

       do i=1,ng_n
       vg_n(i,m)=vg
       enddo

************************************************
**** wg = (H-Eref)**2 * vg
************************************************
       call Hpsi_comp(vg,vgh,ilocal,vr,workr_n,kpt)
       do i=1,ng_n
       vgh(i)=vgh(i)-Eref*vg(i)
       enddo

       call Hpsi_comp(vgh,wg,ilocal,vr,workr_n,kpt)
       do i=1,ng_n
       wg(i)=wg(i)-Eref*vgh(i)
       enddo


************************************************
**** alpha = wg^dagger * vg, and it is the m th element of diag
**** wg = wg - alpha * vg, when m=1
**** wg = wg - alpha * vg - beta * vg(prev), when m>1
************************************************
       alpha=0.d0
       
       do i=1,ng_n
       alpha=alpha+dreal(dconjg(wg(i))*vg(i))
       enddo
       call global_sumr(alpha)
       alpha=alpha*vol
       diag(m)=alpha

       if (m.eq.1) then
          do i=1,ng_n
            wg(i)=wg(i)-alpha*vg(i)
          enddo
       else
         do i=1,ng_n
           wg(i)=wg(i)-alpha*vg(i)-beta*vg_n(i,m-1)
         enddo
       endif

4000  continue


***********************************************

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
             vg(i)=dcmplx(0.d0,0.d0)
             do ii=1,mx
               vg(i)=vg(i)+eigen_vec(m,ii)*vg_n(i,ii)
             enddo
           enddo

           call Hpsi_comp(vg,vgh,ilocal,vr,workr_n,kpt)

           E0=0.d0
           s=0.d0
           do i=1,ng_n
             E0=E0+dreal(dconjg(vgh(i))*vg(i))
             s=s+cdabs(vg(i))**2
           enddo

           call global_sumr(E0)
           call global_sumr(s)
           E0=E0*vol
           s=s*vol

           E0=dabs(E0/s-Eref)

           if(inode.eq.1) then
           write(6,*) '*******************'
           write(6,*) 's=',s*vol
           write(6,*) 'E0,eigen,diff=',E0,Ediff,dabs(E0-Ediff)
           write(6,*) 'coff mx,mx-1=',eigen_vec(m,mx),eigen_vec(m,mx-1)
           endif

c           mxc=mxc+1
c         endif
       enddo

       deallocate(diag)
       deallocate(subdiag)
       deallocate(vg)
       deallocate(vgh)
       deallocate(wg)
       deallocate(eigen_vec)

      return

      end

