      subroutine MINRESQLP_linear(ilocal,nline,tol,
     &  wgp_n0,vr,workr_n,kpt,Eref,AL,eigen,
     &  err_st,mxc,mstateT)
****************************************
cc     Use the quasi-minimal residual (QMR) method to slove linear equation (H-E)x=w_l
cc     doi: 10.1016/0168-9274(95)00089-5, Algorithm 5.1
cc     P(H-E)P is Hermitian, so J=I and q_n=p_n here (P^2 is the preconditoner).
cc     Written by Lin-Wang Wang, March 30, 2001. 
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
****************************************

      use fft_data
      use load_data
      use data

      use zminresqlpModule, only: ZMINRESQLP, ZSYMORTHO

      implicit double precision (a-h,o-z)

      include 'mpif.h'

      include 'param.escan_real'
***********************************************
      integer status(MPI_STATUS_SIZE)

       complex*16 wgc_n(mg_nx),wgp_n0(mg_nx,mstateT) ! wgc_n-->x, wgp_n0-->w_l
       complex*16, allocatable, dimension (:) :: wgp_n, ugh! wgp_n-->P*w_l

       real*8 vr(mr_n) ! poteantial
       real*8 prec(mg_nx), invM(mg_nx) ! preconditioner
       real*8 AL(3,3) ! lattice parameter
c       complex*16 workr_n(mg_nx)
       complex*16 workr_n(*)   ! original workr_n is of mr_n which is larger, xwjiang
       complex*16 Zcoeff(mx),zfac,cai
**********************************************
**** if iopt=0 is used, pghh_old can be deleted
**********************************************
       integer lin_st(mst),m_max(mstateT),mxc
       real*8 E_st(mst),err_st(mst),eigen(mst)
       real*8 Eref,coeff(mst)
       common /com123b/m1,m2,m3,ngb,nghb,ntype,rrcut,msb
       common /comEk/Ek

       complex*16 zdotc

       integer istop, itn
       real*8 rnorm, Arnorm, xnorm, Anorm, Acond
       real*8 maxxnorm, trancond, Acondlim

       maxxnorm = 1.d7
       trancond = 1
       Acondlim = 1.d15

       ng_n=ngtotnod(inode,kpt)

       cai=dcmplx(0.d0,1.d0)

       allocate(wgp_n(mg_nx))
       allocate(ugh(mg_nx))

cccccccccccccccccccccccc
       ! calculate the preconditioner
       x_0=dble(ng_n)
       y_0=0.d0
       do i=1,ng_n
       prec(i)=gkk_n(i,kpt)/Ek
       y_0=y_0+prec(i)
       enddo
       call global_sumr(x_0)
       call global_sumr(y_0)
       y_0=y_0/x_0
       do i=1,ng_n
       prec(i)=1.d0/dsqrt(prec(i)+y_0)
       enddo
       !prec=1.d0
        
       do i=1,ng_n
        x_0=gkk_n(i,kpt)/Ek
        y_0=27.d0+x_0*(18.d0+x_0*(12.d0+x_0*8.d0))
        prec(i)=1.d0/dsqrt(1.d0+16.d0*x_0**4/y_0)
       enddo
       do i=1,ng_n
         invM(i)=prec(i)**2
       enddo

       do 4000 iii=1,mstateT ! for all index 'l'
cccccccccccccccccccccccccccccccccccccccccccc
       do i=1,ng_n
         wgp_n(i)=wgp_n0(i,iii)
       enddo

       call ZMINRESQLP( ng_n, Aprod, wgp_n, nline, vol, Eref,
     &               Msolve, 99, .false., 1.d-1*tol, maxxnorm,
     &               trancond, Acondlim, wgc_n, istop, itn,
     &               rnorm, Arnorm, xnorm, Anorm, Acond )

       
       call Hpsi_comp(wgc_n,ugh,ilocal,vr,workr_n,kpt)

       ugh = ugh - Eref * wgc_n - wgp_n
       err = real(zdotc(ng_n, ugh, 1, ugh, 1))
       call global_sumr(err)
       err = sqrt(err*vol)

      if(inode.eq.1) then
      write(6,888) nint2,Eref*27.211396d0,E0,err
888   format(i4,"  Eref=",f12.5,"  Emin=",E10.4,1x,
     &  "  err=",E10.2)
      write(6,*) 'nint2,istop',nint2,istop
      endif
cccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccjjjjjjjjjj
cccc output wgc_n in real space, including the phase factor 
ccccccccccccccccccccccccccccccccccc
!      do i=1,ng_n
!       wgc_n(i)=prec(i)*wgc_n(i)
!       enddo

!      call d3fft_comp(wgc_n,workr_n,-1,kpt)

ccccccccccccccccccccccccccccccccccc
cc debug, write out size, xwjiang
c       if(inode.eq.1) then
c         write(6,*) "nr_n = ", nr_n
c         write(6,*) "size(workr_n) = ", size(workr_n)
c       endif

       do ii=1,nr_n
       iit=ii+(inode-1)*nr_n-1
       it=iit/(n3*n2)
       jt=(iit-it*n3*n2)/n3
       kt=iit-it*n3*n2-jt*n3
       xt=AL(1,1)*it/n1+AL(1,2)*jt/n2+AL(1,3)*kt/n3
       yt=AL(2,1)*it/n1+AL(2,2)*jt/n2+AL(2,3)*kt/n3
       zt=AL(3,1)*it/n1+AL(3,2)*jt/n2+AL(3,3)*kt/n3

       workr_n(ii)=workr_n(ii)*cdexp(cai*(
     & xt*akx(kpt)+yt*aky(kpt)+      ! special
     & zt*akz(kpt)))
       enddo

       call mpi_barrier(MPI_COMM_WORLD,ierr)

       if(inode.eq.1) then
       write(17) (workr_n(i),i=1,nr_n)
       endif

       do i=1,nnodes-1
       call mpi_barrier(MPI_COMM_WORLD,ierr)
       if(inode==i+1) then
       call  mpi_send(workr_n,nr_n,MPI_DOUBLE_COMPLEX,0,
     &   100,MPI_COMM_WORLD,ierr)
       endif
       if(inode.eq.1) then
        call mpi_recv(workr_n,nr_n,MPI_DOUBLE_COMPLEX,i,
     &   100,MPI_COMM_WORLD,status,ierr)
       write(17) (workr_n(ii),ii=1,nr_n)
       endif
       enddo

4000  continue

      deallocate(wgp_n)
***********************************************

       return
       
       contains
       subroutine Aprod(n, x, y)
        integer(4), intent(in)    :: n
        complex(8), intent(in)    :: x(n)
        complex(8), intent(out)   :: y(n)
        call Hpsi_comp(x,y,ilocal,vr,workr_n,kpt)
       end subroutine Aprod

       subroutine Msolve(n, x, y)
        integer(4), intent(in)    :: n
        complex(8), intent(in)    :: x(n)
        complex(8), intent(out)   :: y(n)
        do i=1,n
          y(i) = invM(i) * x(i)
        enddo
       end subroutine Msolve

      end
