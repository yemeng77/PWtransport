c      subroutine DIIS_comp(ilocal,nline,tol,
c     &  E_st,err_st,ave_line,vr,workr_n,mCGbad,kpt,iislda)
      subroutine DIIS_comp(ilocal,mline,tol,
     &  E_st,err_st,ave_line,vr,workr_n,mCGbad,kpt,iislda)
cccc The original name of DIIS_comp is CG_new
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************
******************************************
ccccc NOTE:
ccccc  we should not keep track and rotate sumdum_m when ipsp_all.eq.1, to save time

cccc some memory over write problem

      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)

      include 'mpif.h'
      include 'param.escan_real'

***********************************************
c      parameter (mline=20)
      parameter (lwork=6000)

      complex*16 workx(lwork)
      real*8 workrx(3*(mst+500))

      complex*16 sumdum(mref,natom)
      complex*16, allocatable,dimension(:,:,:) :: sumdum_all,
     &  sumdum_m
      complex*16, allocatable,dimension(:,:,:,:) :: sumdum_all_BP


      complex*16 ugh(mg_nx),pg(mg_nx)
      complex*16 ugh_m(mg_nx,nblock_band_mx)
      complex*16,allocatable,dimension(:,:) :: sug_m

c      complex*16 Rg(mg_nx,mline+1),Rgh(mg_nx,mline+1)
c      complex*16 SRg(mg_nx,mline+1)
      complex*16,allocatable,dimension(:,:,:) :: Rg_BP,Rgh_BP,
     &   SRg_BP

       complex*16 beta_psi_tmp(nref_tot,mx/nnodes+1)
      
**** first, use the stupid but easy to program way

       real*8 vr(mr_n),workr_n(mr_n)

       real*8 prec(mg_nx)

       real*8 Dij0(32,32,mtype),Qij(32,32,mtype)     ! these two variable still defined as (32,32)
       integer isNLa(9,matom),ipsp_type(mtype)

c       complex*16 hhx(mst,mst),chh(mst,mst)
c       complex*16 U_cholesky(mx,mx)

        complex*16, allocatable, dimension(:,:) :: hhx

c       complex*16 hh(mline,mline),Z(mline,mline)
c       complex*16 hs(mline,mline),ss(mline,mline)
c       real*8 EE(mline),EE1(mline),w(mline)
c       real*8 aux(3*mline)
       complex*16,allocatable,dimension(:,:,:) :: hh_BP,z_BP,
     &  hs_BP,ss_BP
       real*8,allocatable,dimension(:,:) :: EE_BP,EE1_BP,w_BP,
     &  aux_BP
       real*8,allocatable,dimension(:) :: fnorm_BP,Eref_BP,err_BP

       real*8 occ_t(mtype)
       integer iiatom(mtype),icore(mtype),numref(matom),ityatom(mtype)
**********************************************
        integer lin_st(mst)
	real*8 E_st(mst),err_st(mst),dumm(mst),E_st_tmp(mst)
        real*8 eigen_m(nblock_band_mx),err_m(nblock_band_mx)
	complex*16 c,cc,c1,c2,c3
        complex*16 c_one,c_zero


       common /comisNLa/isNLa,Dij0,Qij,ipsp_all,ipsp_type
       common /com123b/m1,m2,m3,ngb,nghb,ntype,rrcut,msb
       common /comEk/Ek
       common /comNL2/occ_t,iiatom,icore,numref,ityatom

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccc  

       if(mst.ne.mx) then
       write(6,*) "mst.ne.mx, stop"
       stop
       endif

       allocate(hhx(mst,mst))

       amem=500     ! 500 MB for temperary memory 
       mband=amem*1000000/(mg_nx*(mline+1)*8*3)
       if(mband.gt.nblock_band_mx) mband=nblock_band_mx
       nband_iter_n=nblock_band_mx*0.9999999d0/mband+1
       if(mband.eq.nblock_band_mx) nband_iter_n=1

        allocate(Rg_BP(mg_nx,mband,mline+1))
        allocate(Rgh_BP(mg_nx,mband,mline+1))
        allocate(SRg_BP(mg_nx,mband,mline+1))
        allocate(hh_BP(mline,mline,mband))
        allocate(z_BP(mline,mline,mband))
        allocate(hs_BP(mline,mline,mband))
        allocate(ss_BP(mline,mline,mband))
        allocate(EE_BP(mline,mband))
        allocate(EE1_BP(mline,mband))
        allocate(w_BP(mline,mband))
        allocate(aux_BP(mline,mband))
        allocate(fnorm_BP(mband))
        allocate(Eref_BP(mband))
        allocate(err_BP(mband))
cccccccccc it might not be necessary to keep the track of sumdum_all_BP, it will be recalculated
       allocate(sumdum_m(mref,natom,nblock_band_mx))       
       allocate(sumdum_all(mref,natom,mline))
       allocate(sumdum_all_BP(mref,natom,mline,mband))
       sumdum_m=0.d0
       sumdum_all=0.d0
       sumdum_all_BP=0.d0

       if(ipsp_all.eq.2) then
        allocate(sug_m(mg_nx,nblock_band_mx))
        sug_m=0.d0
        endif


       ng_n=ngtotnod(inode,kpt)

c       if(nline.gt.mline) then
c       write(6,*) "nline > mline, stop", nline
c       stop
c       endif

       do i=1,ng_n
       x=gkk_n(i,kpt)/Ek
       y=27.d0+x*(18.d0+x*(12.d0+x*8.d0))
       prec(i)=y/(y+16.d0*x**4)
       enddo


ccccccccccccccccccc  subspace diagonalization at the beginning of the CG_AllBand
 
       c_one=dcmplx(1.d0,0.d0)
       c_zero=dcmplx(0.d0,0.d0)

       ugh_m=dcmplx(0.d0,0.d0)

 
            call Hpsi_comp_AllBandBP(ug_n_bp,ugh_m,
     &           nblock_band_mx,ilocal,vr,workr_n,kpt,1,
     &           sug_m,sumdum_m,iislda)


 
      call dot_product_BP(ugh_m,ug_n_bp,ng_n,hhx)

      call system_czheev('V','U',mx,hhx,mst,E_st_tmp,workx,
     &     lwork,workrx,info,MPI_COMM_K1)

      E_st=0.d0
      E_st(band_dis(1):band_dis(2)) =
     &     E_st_tmp(band_dis(1):band_dis(2))*vol

       call rotate_wfBP(ug_n_bp,hhx,ng_n,mg_nx)

       call rotate_wfBP(ugh_m,hhx,ng_n,mg_nx)

       if(ipsp_all.eq.2) then
       call rotate_wfBP(sug_m,hhx,ng_n,mg_nx)
       endif


c       do 1000 m=1,mx
c       do 1000 m=1, nblock_band_mx
       do 1000 iter=1,nband_iter_n

       mmband=mband
       if(iter.eq.nband_iter_n) 
     &   mmband=nblock_band_mx-(nband_iter_n-1)*mband
        if(mmband.lt.1) then
         write(6,*) "mmband.lt.1,stop",mmband
         stop
         endif

***************************************************
       mstart=(iter-1)*mband
***************************************************
       do m=1,mmband
       do i=1,ng_n
       Rg_BP(i,m,1)=ug_n_BP(i,mstart+m)
       enddo
       enddo
***************************************************
c       do 200 nt=1,nline
       do 200 nt=1,mline


ccccccc  the input Rg(1,nt) might not be normalized, since SRg is not known !

      if(nt.eq.1) then
       do m=1,mmband
       do i=1,ng_n
       Rgh_BP(i,m,1)=ugh_m(i,mstart+m)
       enddo
       enddo

         if(ipsp_all.eq.2) then
       do m=1,mmband
       do i=1,ng_n
       SRg_BP(i,m,1)=sug_m(i,mstart+m)
       enddo
       enddo
         endif

      else

c       call Hpsi_comp(Rg(1,nt),Rgh(1,nt),
c     &    ilocal,vr,workr_n,kpt,1,SRg(1,nt),sumdum,iislda) 


        call Hpsi_comp_AllBandBP(Rg_BP(1,1,nt),Rgh_BP(1,1,nt),
     &  mmband,ilocal,vr,workr_n,kpt,1,SRg_BP(1,1,nt),sumdum_m,iislda)

       endif


ccccccc  now, normalize Rg(1,nt)

c       call orth_comp(Rg(1,nt),Rg,0,1,kpt,
c     &  SRg,SRg(1,nt),ipsp_all,1,1,fnorm)     ! normalize Rg and SRg

       call orth_comp_DIIS(Rg_BP(1,1,nt),Rg_BP,0,1,kpt,
     &  SRg_BP,SRg_BP(1,1,nt),mmband,ipsp_all,1,1,fnorm_BP)     ! normalize Rg and SRg

        do m=1,mmband
        do i=1,ng_n
        Rgh_BP(i,m,nt)=fnorm_BP(m)*Rgh_BP(i,m,nt)
        enddo
        sumdum_m(:,:,m)=fnorm_BP(m)*sumdum_m(:,:,m)
        sumdum_all_BP(:,:,nt,m)=sumdum_m(:,:,m)
        enddo
cccccccccccccccccccccccccccccccccccccc
       if(nt.eq.1) then
       do m=1,mmband
       s=0.d0
       do i=1,ng_n
       s=s+dreal(Rgh_BP(i,m,nt)*dconjg(Rg_BP(i,m,nt)))
       enddo
       
       call global_sumr(s)
       s=s*vol
       Eref_BP(m)=s
       enddo  ! m

       endif

ccccccccccccccccccccccccccccccccccccccccccccccccc
       do 50 nt1=1,nt

       do 51 m=1,mmband

       c1=dcmplx(0.d0,0.d0)
       do i=1,ng_n
       c1=c1+Rgh_BP(i,m,nt)*dconjg(Rgh_BP(i,m,nt1))
       enddo

       c2=dcmplx(0.d0,0.d0)

       if(ipsp_all.eq.1) then
       do i=1,ng_n
       c2=c2+Rgh_BP(i,m,nt)*dconjg(Rg_BP(i,m,nt1))
       enddo
       else
       do i=1,ng_n
       c2=c2+Rgh_BP(i,m,nt)*dconjg(SRg_BP(i,m,nt1))+
     &       SRg_BP(i,m,nt)*dconjg(Rgh_BP(i,m,nt1)) 
       enddo
       c2=c2/2
       endif

       call global_sumc(c1)
       call global_sumc(c2)
       c1=c1*vol
       c2=c2*vol
       hh_BP(nt,nt1,m)=dconjg(c1)
       hh_BP(nt1,nt,m)=c1
       hs_BP(nt,nt1,m)=dconjg(c2)
       hs_BP(nt1,nt,m)=c2

       if(ipsp_all.eq.2) then
       c3=dcmplx(0.d0,0.d0)
       do i=1,ng_n
       c3=c3+SRg_BP(i,m,nt)*dconjg(SRg_BP(i,m,nt1))
       enddo
       call global_sumc(c3)
       c3=c3*vol
       ss_BP(nt,nt1,m)=dconjg(c3)
       ss_BP(nt1,nt,m)=c3
       endif

 51    continue
 50    continue
       
cccccccccccccccccccccccccccccccccccccccccccc
       err_max=0.d0
       do 60 m=1,mmband

       do i1=1,nt
       do i2=1,nt
       z_BP(i1,i2,m)=hh_BP(i1,i2,m)-2*Eref_BP(m)*hs_BP(i1,i2,m)
       if(ipsp_all.eq.2) then
       z_BP(i1,i2,m)=z_BP(i1,i2,m)+Eref_BP(m)**2*ss_BP(i1,i2,m)
       endif
       enddo
       if(ipsp_all.eq.1) then
       z_BP(i1,i1,m)=z_BP(i1,i1,m)+Eref_BP(m)**2
       endif
       enddo


c       call system_czheev('v','U',nt,z_BP(1,1,m),mline,
c     &   EE_BP(1,m),workx,
c     &  lwork,workrx,info,MPI_COMM_K1)


       call zheev('v','U',nt,z_BP(1,1,m),mline,
     &   EE_BP(1,m),workx,
     &  lwork,workrx,info)


       do i=1,ng_n
       c=dcmplx(0.d0,0.d0)
       c1=dcmplx(0.d0,0.d0)
       do nt1=1,nt
       c=c+z_BP(nt1,1,m)*Rg_BP(i,m,nt1)
       c1=c1+z_BP(nt1,1,m)*Rgh_BP(i,m,nt1)
       enddo
       Rg_BP(i,m,nt+1)=c
       Rgh_BP(i,m,nt+1)=c1
       enddo



       if(ipsp_all.eq.2) then
       do i=1,ng_n
       c=dcmplx(0.d0,0.d0)
       do nt1=1,nt
       c=c+z_BP(nt1,1,m)*SRg_BP(i,m,nt1)
       enddo
       SRg_BP(i,m,nt+1)=c
       enddo
       endif

************************************************
*** update E_ref according to Rayleigh quotient
************************************************
       s=0.d0
       do i=1,ng_n
       s=s+dreal(Rgh_BP(i,m,nt+1)*dconjg(Rg_BP(i,m,nt+1)))
       enddo

       call global_sumr(s)
       s=s*vol
       Eref_BP(m)=s
**************************************************

       err_BP(m)=dsqrt(dabs(EE_BP(1,m)))
       if(err_BP(m).gt.err_max) err_max=err_BP(m)

60     continue

       if(nt.eq.mline.or.err_max.lt.tol) then

       do m=1,mmband
       lin_st(mstart+m)=nt
       do i=1,ng_n
       ug_n_BP(i,mstart+m)=Rg_BP(i,m,nt+1)
       enddo

       if(ipsp_all.eq.2) then
       do i=1,ng_n 
       sug_m(i,mstart+m)=SRg_BP(i,m,nt+1)
       enddo
       endif

       sumdum_m(:,:,m)=dcmplx(0.d0,0.d0)
       do nt1=1,nt
       sumdum_m(:,:,m)=sumdum_m(:,:,m)+dconjg(z_BP(nt1,1,m))*           ! sumdum_m is conjg(ug_n)
     &      sumdum_all_BP(:,:,nt1,m)
       enddo

       enddo    ! do m

       goto 201

       else

cccccccccccccccccccccccccccccccccccccccccc

c       call orth_comp(Rgh(1,nt+1),Rg,nt+1,3,kpt,
d     &   SRg,SRg(1,nt+1),ipsp_all,0,-1,fnorm)    ! SRg(1,nt+1) not used


       call orth_comp_DIIS(Rgh_BP(1,1,nt+1),Rg_BP,
     & nt+1,3,kpt,SRg_BP,SRg_BP(1,1,nt+1),mmband,ipsp_all,0,-1,
     & fnorm_BP)   


       do m=1,mmband
       do i=1,ng_n
       Rg_BP(i,m,nt+1)=Rgh_BP(i,m,nt+1)*prec(i)
       enddo
       enddo

c       call orth_comp(Rg(1,nt+1),Rg,nt,2,kpt,
c     &   SRg,SRg(1,nt+1),ipsp_all,0,1,fnorm)     ! SRg(1,nt+1) not used, not normalized 

        call orth_comp_DIIS(Rg_BP(1,1,nt+1),Rg_BP,
     &  nt,2,kpt,SRg_BP,SRg_BP(1,1,nt+1),mmband,ipsp_all,0,1,
     &  fnorm_BP)     ! SRg(1,nt+1) not used, not normalized


       endif

200    continue

201    continue

       do m=1,mmband
       err_m(m+mstart)=err_BP(m)
       eigen_m(m+mstart)=Eref_BP(m)
       enddo


1000   continue

cccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccc

       if(ipsp_all.eq.1) then
c       call orthogonal_choleskyBP(ug_n_bp,mg_nx,nblock_band_mx,
c     &     ng_n,U_cholesky)
       call orthogonal_choleskyBP(ug_n_bp,mg_nx,nblock_band_mx,
     &     ng_n,hhx)
       else
c       call orthogonal_choleskyBP2(ug_n_bp,sug_m,mg_nx,
c     &  nblock_band_mx,ng_n,U_cholesky,0)                   ! iflag=0, sug_m is not rotated for orthogonalization
       call orthogonal_choleskyBP2(ug_n_bp,sug_m,mg_nx,
     &  nblock_band_mx,ng_n,hhx,0)                   ! iflag=0, sug_m is not rotated for orthogonalization
       endif 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Subspace diagonalizaton.
 
      c_one=dcmplx(1.d0,0.d0)
      c_zero=dcmplx(0.d0,0.d0)

          call Hpsi_comp_AllBandBP(ug_n_bp,ugh_m,
     &           nblock_band_mx,ilocal,vr,workr_n,kpt,0,        ! isij=0, sug_m is not calculated
     &           sug_m,sumdum_m,iislda)

      call dot_product_BP(ugh_m,ug_n_bp,ng_n,hhx)

      call system_czheev('V','U',mx,hhx,mst,E_st_tmp,workx,
     &     lwork,workrx,info,MPI_COMM_K1)

      E_st=0.d0 

       E_st(1:mst) = E_st_tmp(1:mst)*vol

       call rotate_wfBP(ug_n_bp,hhx,ng_n,mg_nx)

c        chh=dconjg(hhx)
c        call rotate_wfBP(sumdum_m,chh,mref*natom,mref*natom)   

        sumdum_m=dconjg(sumdum_m)
        call rotate_wfBP(sumdum_m,hhx,mref*natom,mref*natom)   
        sumdum_m=dconjg(sumdum_m)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       beta_psi_tmp=dcmplx(0.d0,0.d0)
        mx_n=mx/nnodes+1
       do m1=1,nblock_band_mx
       m=m1+icolor_b*nblock_band_mx
       if((m-1)/mx_n+1.eq.inode) then
       im=m-(inode-1)*mx_n
       iref_t=0
       do ia=1,natom
       do iref=1,numref(ia)
       iref_t=iref_t+1
       beta_psi_tmp(iref_t,im)=sumdum_m(iref,ia,m1)
       enddo
       enddo
       endif
       enddo

       nsize_beta_psi=nref_tot*(mx/nnodes+1)
        call mpi_allreduce(beta_psi_tmp,beta_psi,nsize_beta_psi,
     &     MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_B2,ierr)
ccccccccccccccccccccccccccccccccccccccccccccccc

!     End of subspace diagonalization.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       deallocate(sumdum_all)
       deallocate(sumdum_all_BP)


       ave_line=0
       do i=1,nblock_band_mx
       ave_line=ave_line+lin_st(i)
       enddo
*****************************************************
****************************************************
       err_st=0.d0
       err_st(band_dis(1):band_dis(2))
     &    =err_m(1:nblock_band_mx) 

       call mpi_allreduce(err_st,dumm,mst,MPI_REAL8,
     &        MPI_SUM,MPI_COMM_B2,ierr)

        err_st = dumm


       if(inode_tot.eq.1) then
       write(6,*) "*********************************"
       write(6,*) "****** kpt=",kpt
       write(6,*) "err of each states, A.U"
       write(6,102) (err_st(i), i=1,mx)
       write(6,*) "eigen energies, in eV"
       write(6,103) (E_st(i)*27.211396d0, i=1,mx)
       write(6,*) "*********************************"
       write(6,*) "sum eigen(dbl occ, Ryd)=", Esum*2*2
       write(6,*) "*********************************"
       endif

101   format(5(i6,7x))
102   format(5(E10.4,3x))
103   format(5(f13.8,1x))
***********************************************
*** for metal only
***********************************************
        deallocate(Rg_BP)
        deallocate(Rgh_BP)
        deallocate(SRg_BP)
        deallocate(hh_BP)
        deallocate(z_BP)
        deallocate(hs_BP)
        deallocate(ss_BP)
        deallocate(EE_BP)
        deallocate(EE1_BP)
        deallocate(w_BP)
        deallocate(aux_BP)
        deallocate(fnorm_BP)
        deallocate(Eref_BP)
        deallocate(err_BP)
        deallocate(sumdum_m)
        deallocate(hhx)

       if(ipsp_all.eq.2) then
        deallocate(sug_m)
        endif

       call mpi_allreduce(ave_line,ave_line_tmp,1,MPI_REAL8,
     & MPI_SUM,MPI_COMM_B2,ierr)
       ave_line=ave_line_tmp/mx


      return
      end

