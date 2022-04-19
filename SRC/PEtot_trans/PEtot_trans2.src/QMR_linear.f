      subroutine QMR_linear(ilocal,nline,tol,
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

      implicit double precision (a-h,o-z)

      include 'mpif.h'

      include 'param.escan_real'
***********************************************
      integer status(MPI_STATUS_SIZE)

       complex*16 ugh(mg_nx) ! ugh-->P(w_l-(H-E)x)
       complex*16 dg(mg_nx),dgh(mg_nx) ! dg-->d, dgh-->P(H-E)P*d
       complex*16 pg(mg_nx),pg1(mg_nx),pgh(mg_nx) ! pg-->p, pg1-->P*p, pgh-->P(H-E)P*p
       complex*16 vg(mg_nx) ! vg-->v

       complex*16 wgc_n(mg_nx),wgp_n0(mg_nx,mstateT) ! wgc_n-->x, wgp_n0-->w_l
       complex*16, allocatable, dimension (:) :: wgp_n ! wgp_n-->P*w_l

       real*8 vr(mr_n) ! poteantial
       real*8 prec(mg_nx) ! preconditioner
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

       ng_n=ngtotnod(inode,kpt)

       cai=dcmplx(0.d0,1.d0)

       allocate(wgp_n(mg_nx))

cccccccccccccccccccccccc
       ! calculate the preconditioner
       x=dble(ng_n)
       y=0.d0
       do i=1,ng_n
       prec(i)=gkk_n(i,kpt)/Ek
       y=y+prec(i)
       enddo
       call global_sumr(x)
       call global_sumr(y)
       y=y/x
       do i=1,ng_n
       prec(i)=1.d0/dsqrt(prec(i)+y)
       enddo
       !prec=1.d0
        
       do i=1,ng_n
        x=gkk_n(i,kpt)/Ek
        y=27.d0+x*(18.d0+x*(12.d0+x*8.d0))
        prec(i)=1.d0/dsqrt(1.d0+16.d0*x**4/y)
       enddo

       do 4000 iii=1,mstateT ! for all index 'l'
cccccccccccccccccccccccccccccccccccccccccccc
       do i=1,ng_n
       wgp_n(i)=prec(i)*wgp_n0(i,iii) ! wgp_n-->P*w_l
       enddo
       wgc_n=dcmplx(0.d0,0.d0) ! x_0=0
       call Hpsi_comp(wgc_n,ugh,ilocal,vr,workr_n,kpt)

       tau=0.d0
       do i=1,ng_n
       ugh(i)=wgp_n(i)-prec(i)*(ugh(i)-Eref*wgc_n(i)) ! ugh-->P*w_l-P(H-E)P*x_0
       vg(i)=ugh(i) ! v_1=P[w_l-(H-E)Px_0]
       tau=tau+cdabs(vg(i))**2
       enddo
       call global_sumr(tau)
       tau=tau*vol ! tau_1 = ||v_1||^2
       pg=vg ! p_1=v_1
       rho=tau ! rho_1=tau_1
       theta=0.d0 ! theta_0=0

       dg=dcmplx(0.d0,0.d0) ! d_0=0
       call Hpsi_comp(dg,dgh,ilocal,vr,workr_n,kpt)
       do i=1,ng_n
       dgh(i)=prec(i)*(dgh(i)-Eref*dg(i)) ! P(H-E)P*d_0
       enddo

       do 3000 nint2=1,nline
       do i=1,ng_n
       pg1(i)=prec(i)*pg(i) ! pg1-->P*p_n
       enddo
       call Hpsi_comp(pg1,pgh,ilocal,vr,workr_n,kpt)
       sigma=0.d0
       do i=1,ng_n
       pgh(i)=prec(i)*(pgh(i)-Eref*pg1(i)) ! t_n=P(H-E)P*p_n
       sigma=sigma+dreal(dconjg(pg(i))*pgh(i))
       enddo
       call global_sumr(sigma)
       sigma=sigma*vol ! sigma_n=<p_n|P(H-E)P|p_n>
       err3=dsqrt(dabs(sigma))
       if(err3.lt.tol) goto 3001 ! if sigma_n=0, stop
       alpha=rho/sigma ! alpha_n=rho_n/sigma_n
       do i=1,ng_n
       vg(i)=vg(i)-alpha*pgh(i) ! v_{n+1}=v_n-alpha_n*t_n
       enddo
       rho_old=rho ! rho_old-->rho_n
       rho=0.d0
       do i=1,ng_n
       rho=rho+cdabs(vg(i))**2
       enddo
       call global_sumr(rho)
       rho=rho*vol ! rho_{n+1}=||v_{n+1}||^2
       theta_old=theta ! theta_old-->theat_{n-1}
       theta=rho/tau ! theta_n=rho_{n+1}/tau_n
       psi=1.d0/(1.d0+theta) ! psi_n=1/(1+theta_n)
       tau=tau*theta*psi ! tau_{n+1}=tau_n*theta_n*psi_n
       temp1=psi*theta_old ! temp1=psi_n*theat_{n-1}
       temp2=psi*alpha ! teamp2=psi_n*alpha_n
       err=0.d0
       err2=0.d0
       do i=1,ng_n
       dg(i)=temp1*dg(i)+temp2*pg(i) ! d_n=(psi_n*theat_{n-1})d_{n-1}+(psi_n*alpha_n)p_n
       dgh(i)=temp1*dgh(i)+temp2*pgh(i) ! dgh-->P(H-E)P*d_n
       wgc_n(i)=wgc_n(i)+dg(i) ! x_n=x_{n-1}+d_n
       ugh(i)=ugh(i)-dgh(i) ! ugh-->P*w_l-P(H-E)P*x_n
       err=err+cdabs(ugh(i)/prec(i))**2
       err2=err2+cdabs(ugh(i))**2
       enddo
       call global_sumr(err)
       err=dsqrt(err*vol) ! err=||w_l-(H-E)P*x_n||
       call global_sumr(err2)
       err2=dsqrt(err2*vol) ! err2=||P*w_l-P(H-E)P*x_n||
       if(err.lt.tol) goto 3001 ! if converged, stop
       beta=rho/rho_old ! beta_n=rho_{n+1}/rho_n
          
       do i=1,ng_n
       pg(i)=vg(i)+beta*pg(i) ! p_{n+1}=v_{n+1}+beta_n*p_n
       enddo

      if(inode.eq.1) then
      write(6,777) nint2,err,err2
      endif
777   format(i8,2(E20.12,2x))

**********************************************
**** do 3000, is for the nline line minimization
**********************************************
3000  continue
3001  continue


      if(inode.eq.1) then
      write(6,888) nint2,Eref*27.211396d0,E0,err
888   format(i4,"  Eref=",f12.5,"  Emin=",E10.4,1x,
     &  "  err=",E10.2)
      write(6,*) 'nint2,E1',nint2,E1
      endif
cccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccjjjjjjjjjj
cccc output wgc_n in real space, including the phase factor 
ccccccccccccccccccccccccccccccccccc
       do i=1,ng_n
       wgc_n(i)=prec(i)*wgc_n(i)
       enddo

       call d3fft_comp(wgc_n,workr_n,-1,kpt)

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
      end
