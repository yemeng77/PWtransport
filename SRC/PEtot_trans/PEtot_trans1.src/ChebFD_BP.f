      subroutine ChebFD_BP(ilocal,nline,tol,
     &     E_st,err_st,vr,workr_n,kpt,iislda,
     &     Ewind,Np,Nlan,ikernel)
****************************************
cc     Written by Meng Ye, December 28, 2017. 
cc     Copyright 2017 The Regents of the University of California
cc     The United States government retains a royalty free license in this work 
****************************************

****************************************
cc     Use the Chebyshev filter diagonalization to calculate the interior eigen states in Ewind
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

      real*8 vr(mr_n)
      complex*16 workr_n(mr_n)

      real*8 Ebound(2),Ewind(2)
      real*8 E_st(mst),err_st(mst),err2_st(mst)
      real*8 E_m(nblock_band_mx)
      real*8 err_m(nblock_band_mx)
      real*8 tmp_real(nblock_band_mx),dumm(mst)
      integer Np,Nlan,ikernel,info

      real*8 alpha,beta,c0,g0

      complex*16 xg_n_bp(mg_nx,nblock_band_mx)
      complex*16 yg_n_bp(mg_nx,nblock_band_mx)
      complex*16 zg_n_bp(mg_nx,nblock_band_mx)
      complex*16 sumdum_m(mref,natom,nblock_band_mx)

      real*8, allocatable, dimension (:) :: cn,gn
      complex*16, allocatable, dimension (:,:) :: sug_m,hh

      parameter (lwork=60000)
      complex*16 workx(lwork)
      real*8 workrx(3*mst)

      real*8  Dij0(32,32,mtype),Qij(32,32,mtype)
      integer isNLa(9,matom),ipsp_type(mtype),ipsp_all

      common /comisNLa/isNLa,Dij0,Qij,ipsp_all,ipsp_type

      ng_n=ngtotnod(inode,kpt)

      if(Nlan.lt.20) Nlan=20
      call lanczos_esti(ilocal,Nlan,vr,workr_n,kpt,iislda,Ebound)

      if(inode_tot.eq.1) then
      write(6,*) "*********************************"
      write(6,*) "estimated energy bound, in eV"
      write(6,101) (Ebound(i)*27.211396d0, i=1,2)
      write(6,*) "*********************************"
      endif

      if(Ewind(1).eq.Ewind(2)) then
       if(inode_tot.eq.1) write(6,*)
     &   "The width of energy window is 0, cannot use ChebFD, stop"
       call mpi_abort(MPI_COMM_WORLD,ierr)
      endif
      if(Ewind(1).lt.Ebound(1)) Ewind(1)=Ebound(1)
      if(Ewind(2).gt.Ebound(2)) Ewind(2)=Ebound(2)

      if(inode_tot.eq.1) then
      write(6,*) "*********************************"
      write(6,*) "energy windows to use the ChebFD, in eV"
      write(6,101) (Ewind(i)*27.211396d0, i=1,2)
      write(6,*) "The degree of ChebFD:",Np
      write(6,*) "*********************************"
      endif

****************************************
cc    map the x in [Ebound(1),Ebound(2)] to [-1,1] using alpha*x+beta
cc    alpha=2/(Ebound(2)-Ebound(1)), beta=(Ebound(1)+Ebound(2))/(Ebound(1)-Ebound(2))
****************************************      
      alpha=2.0d0/(Ebound(2)-Ebound(1))
      beta=(Ebound(1)+Ebound(2))/(Ebound(1)-Ebound(2))

      if(mst.ne.mx) then
      if(inode_tot.eq.1) write(6,*) "mst.ne.mx,stop"
      call mpi_abort(MPI_COMM_WORLD,ierr)
      endif

      sumdum_m=0.d0
      if(ipsp_all.eq.1) then
      allocate(sug_m(1,nblock_band_mx))
      else
      allocate(sug_m(mg_nx,nblock_band_mx))
      endif
      sug_m=0.d0
      allocate(hh(mst,mst))

****************************************
cc    W(x) is a windows function, which is 1 in [Ewind(1),Ewind(2)] and 0 outside
cc    W(x)\approx\sum_0^mp gn*cn*Tn(alpha*x+beta)
cc    The filter: W(H)ug=\sum_0^mp gn*cn*Tn(alpha*H+beta)ug
****************************************
      allocate(cn(Np))
      allocate(gn(Np))
      pi=4.0d0*datan(1.0d0)
      c0=(dacos(alpha*Ewind(2)+beta)-dacos(alpha*Ewind(1)+beta))/pi
      g0=1.0d0

      do n=1,Np
      cn(n)=2.0d0*(dsin(n*dacos(alpha*Ewind(2)+beta))
     &     -dsin(n*dacos(alpha*Ewind(1)+beta)))/(n*pi)
      if(ikernel.gt.0) then
        gn(n)=(dsin(n*pi/(Np+1))/(n*pi/(Np+1)))**ikernel !Lanczos kernel
      else if(ikernel.eq.0) then
        gn(n)=((Np-n)*dcos(n*pi/Np)+dsin(n*pi/Np)/dtan(pi/Np))/Np !Jackson kernel
      else if(ikernel.eq.0) then
        gn(n)=dble(Np-n+1)/dble(Np+1) !Fej\'er kernel
      else
        if(inode_tot.eq.1) write(6,*) 
     &   "warning: ikernel is not recognized, no kernel is used."
        gn(n)=1.d0
      endif
      enddo

      do 3000 nint=1,nline
****************************************
cc     xg=T0(alpha*H+beta)ug=ug
cc     yg=T1(alpha*H+beta)ug=(alpha*H+beta)ug=(alpha*H+beta)xg
****************************************
       xg_n_bp=ug_n_bp
       call Hpsi_comp_AllBandBP(xg_n_bp,yg_n_bp,
     &      nblock_band_mx,ilocal,vr,workr_n,kpt,1,
     &      sug_m,sumdum_m,iislda)

       do m=1,nblock_band_mx
         do i=1,ng_n
           yg_n_bp(i,m)=alpha*yg_n_bp(i,m)+beta*xg_n_bp(i,m)
           ug_n_bp(i,m)=g0*c0*xg_n_bp(i,m)+gn(1)*cn(1)*yg_n_bp(i,m)
         enddo
       enddo

       do n=2,Np
****************************************
cc     xg=T(n-2)(alpha*H+beta)ug
cc     yg=T(n-1)(alpha*H+beta)ug
cc     zg=Tn(alpha*H+beta)ug=2*(alpha*H+beta)T(n-1)ug-T(n-2)ug=2*(alpha*H+beta)yg-xg
cc     yg->xg zg->yg
****************************************
       call Hpsi_comp_AllBandBP(yg_n_bp,zg_n_bp,
     &      nblock_band_mx,ilocal,vr,workr_n,kpt,1,
     &      sug_m,sumdum_m,iislda)

       do m=1,nblock_band_mx
         do i=1,ng_n
           zg_n_bp(i,m)=2.0d0*(alpha*zg_n_bp(i,m)+
     &      beta*yg_n_bp(i,m))-xg_n_bp(i,m)
           ug_n_bp(i,m)=ug_n_bp(i,m)+gn(n)*cn(n)*zg_n_bp(i,m)
         enddo
       enddo
       xg_n_bp=yg_n_bp
       yg_n_bp=zg_n_bp
       enddo

       call mpi_barrier(MPI_COMM_WORLD,ierr)

****************************************
cc     orthonormalize ug using Cholesky decompostion
cc     ug^T*ug=S=U^T*U, then ug(U^-1) is orthonormalized vectors.
***************************************
C        if(ipsp_all.eq.1) then
C        call orthogonal_choleskyBP(ug_n_bp,mg_nx,nblock_band_mx,
C      &      ng_n,hh)    
C        else
C        call orthogonal_choleskyBP2(ug_n_bp,sug_m,mg_nx,
C      &      nblock_band_mx,ng_n,hh,1)      ! iflag=1, sug_m is rotated    
 
C        endif
       call dot_product_BP(ug_n_bp,ug_n_bp,ng_n,hh)
       do i=1,mx
        hh(i,1:i-1)=dcmplx(0.0d0,0.0d0)
        hh(i,i:mx)=hh(i,i:mx)*vol
       enddo

       call zpotrf('U',mx,hh,mx,info)
       if(info.ne.0) then
        if(inode_tot.eq.1) write(6,*) "Something is wrong with 
     &    Cholesky decomposition of S. info", info
       call mpi_abort(MPI_COMM_WORLD,ierr)
       endif

       call ztrtri('U','N',mx,hh,mx,info)
       if(info.ne.0) then
        if(inode_tot.eq.1) write(6,*) "Somthing is wrong with inversion 
     &    of Cholesky factor of S. info", info
       call mpi_abort(MPI_COMM_WORLD,ierr)
       endif

       call rotate_wfBP(ug_n_bp,hh,ng_n,mg_nx)

****************************************
cc     subspace diagonalization
***************************************
       call Hpsi_comp_AllBandBP(ug_n_bp,xg_n_bp,
     &      nblock_band_mx,ilocal,vr,workr_n,kpt,1,
     &      sug_m,sumdum_m,iislda)
       call dot_product_BP(ug_n_bp,xg_n_bp,ng_n,hh)

       call system_czheev('V','U',mx,hh,mst,E_st,workx,
     &      lwork,workrx,info,MPI_COMM_K1)
       if(info.ne.0) then
        if(inode_tot.eq.1) write(6,*) "Somthing is wrong with 
     &    diagonalization of H. info", info
       call mpi_abort(MPI_COMM_WORLD,ierr)
       endif

       E_m(1:nblock_band_mx)=E_st(band_dis(1):band_dis(2))*vol

       call rotate_wfBP(ug_n_bp,hh,ng_n,mg_nx)
       call rotate_wfBP(xg_n_bp,hh,ng_n,mg_nx)
c       sumdum_m=dconjg(sumdum_m)
c       call rotate_wfBP(sumdum_m,hh,mref*natom,mref*natom)
c       sumdum_m=dconjg(sumdum_m)
       if(ipsp_all.eq.2) call rotate_wfBP(sug_m,hh,ng_n,mg_nx)

****************************************
cc     check convergency
***************************************
       err_m=0.d0
       do m=1,nblock_band_mx
       if(ipsp_all.eq.2) then
         do i=1,ng_n
         err_m(m)=err_m(m)
     &           +cdabs(xg_n_bp(i,m)-E_m(m)*sug_m(i,m))**2
         enddo            
       else
         do i=1,ng_n
         err_m(m)=err_m(m)
     &           +cdabs(xg_n_bp(i,m)-E_m(m)*ug_n_bp(i,m))**2
         enddo            
       endif
       enddo
       
       call mpi_allreduce(err_m,tmp_real,
     &      nblock_band_mx,MPI_REAL8,MPI_SUM,MPI_COMM_K,ierr)
       do m=1,nblock_band_mx
       err_m(m)=dsqrt(dabs(tmp_real(m)*vol))
       enddo

       E_st=0.d0
       E_st(band_dis(1):band_dis(2))=E_m(1:nblock_band_mx)
       call mpi_allreduce(E_st,dumm,mst,MPI_REAL8,
     &      MPI_SUM,MPI_COMM_B2,ierr)
       E_st=dumm

       err_st=0.d0
       err_st(band_dis(1):band_dis(2))=err_m(1:nblock_band_mx)
       call mpi_allreduce(err_st,dumm,mst,MPI_REAL8,
     &      MPI_SUM,MPI_COMM_B2,ierr)
       err_st=dumm
       err2_st=dumm

       if(inode_tot.eq.1) then
       write(6,*) "*********************************"
       write(6,*) "****** kpt=",kpt
       write(6,*) "err of each states, A.U"
       write(6,102) (err_st(i), i=1,mx)
       write(6,*) "eigen energies, in eV"
       write(6,103) (E_st(i)*27.211396d0, i=1,mx)
       write(6,*) "*********************************"
       endif

****************************************
cc     Exit if err_st<=tol for all Ritz values in the target interval.
****************************************
       do m=1,mst
       if (E_st(m).lt.Ewind(1).or.E_st(m).gt.Ewind(2)) err2_st(m)=0.d0
       enddo
       if(all(err2_st.lt.tol)) goto 3001

3000  continue
3001  continue

      deallocate(cn)
      deallocate(gn)
      deallocate(sug_m)
      deallocate(hh)

***********************************
      do m=1,mst
       if(err_st(m).gt.tol) then
        E_st(m)=(E_st(m)-E_st(m))/(E_st(m)-E_st(m))
        err_st(m)=0.d0
       endif
      enddo

      if(inode_tot.eq.1) then
       write(6,*) "*********************************"
       write(6,*) "****** kpt=",kpt
       write(6,*) "err of each states, A.U"
       write(6,102) (err_st(i), i=1,mx)
       write(6,*) "eigen energies, in eV"
       write(6,103) (E_st(i)*27.211396d0, i=1,mx)
       write(6,*) "*********************************"
      endif

101   format("[",2(f13.8,1x),"]")
102   format(5(E10.4,3x))
103   format(5(f13.8,1x))

      return

      end

