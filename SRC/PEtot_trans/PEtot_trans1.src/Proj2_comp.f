      subroutine Proj2_comp(ilocal,nline,tol,
     &     E_st,err_st,ave_line,vr,workr_n,mbad,ave_err,kpt,iislda)
****************************************
c     c     Written by Lin-Wang Wang, March 30, 2001. 
c     c     Band parallel version is written by Byounghak Lee, 2007.
*************************************************************************
*     *  copyright (c) 2003, The Regents of the University of California,
*     *  through Lawrence Berkeley National Laboratory (subject to receipt of any
*     *  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************
****************************************
****  The main conjugate gradient routine
****  mbad=0, normal iCG=1 run, ave_err not used
****  mbad>0, called from CG_new, fix the states: mx-mbad,mx to ave_err
******************************************
******************************************
      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)

      include 'mpif.h'
      include 'param.escan_real'
***********************************************
      complex*16 pg_m(mg_nx,nblock_band_mx)
      complex*16 pg_old_m(mg_nx,nblock_band_mx)
      complex*16 pgh_m(mg_nx,nblock_band_mx)
      complex*16 ugh_m(mg_nx,nblock_band_mx)
      complex*16 overlap_mat(nblock_band_mx,nblock_band_mx)
      complex*16 U_cholesky(mx,mx), W_cholesky(mx,mx)           
      complex*16, allocatable, dimension(:,:) :: C_tmp
      real*8 prec_m(mg_nx)
      real*8 s_m(nblock_band_mx)
      real*8 err_m(nblock_band_mx), err2_m(nblock_band_mx)
      real*8 E_ug_m(nblock_band_mx), eigen_m(nblock_band_mx)
      real*8 rr0_m(nblock_band_mx), rr00_m(nblock_band_mx)
      real*8 rr1_m(nblock_band_mx)
      real*8 beta_m(nblock_band_mx)
      real*8 E_upg_m(nblock_band_mx), E_pg_m(nblock_band_mx)
      real*8 theta_m(nblock_band_mx), cos_th_m(nblock_band_mx) 
      real*8 sin_th_m(nblock_band_mx)
      real*8 pred_E_m(nblock_band_mx),tmp_real(nblock_band_mx)
      real*8 dumm(mst)

      real*8 Dij0(32,32,mtype),Qij(32,32,mtype)
      integer isNLa(9,matom),ipsp_type(mtype)

      real*8 vr(mr_n)
      real*8 prec(mg_nx)
      complex*16 workr_n(mr_n)
      complex*16 cc, c_one, c_zero
**********************************************
****  if iopt=0 is used, pghh_old can be deleted
**********************************************
      integer lin_st(mst), idumm(mst)
      real*8 E_st(mst),err_st(mst),err2_st(mst), E_st_tmp(mst)
      real*8 Ef
      complex*16 sumdum(mref,natom),sumdum2(mref,natom)
      complex*16 sumdum_m(mref,natom,nblock_band_mx)
      complex*16 sumdum2_m(mref,natom,nblock_band_mx)
      real*8 occ_t(mtype)
      integer iiatom(mtype),icore(mtype),numref(matom),ityatom(mtype)

      common /com123b/m1,m2,m3,ngb,nghb,ntype,rrcut,msb
      common /comEk/Ek
      common /comisNLa/isNLa,Dij0,Qij,ipsp_all,ipsp_type
      common /comNL2/occ_t,iiatom,icore,numref,ityatom

      complex*16 sss

      complex*16 h(mst,mst),hh(mst,mst),chh(mst,mst)
      parameter (lwork=60000)
      complex*16 workx(lwork)
      real*8 workrx(3*mst)
      integer info
      complex*16, allocatable, dimension(:,:) :: ug_n_tmp, ugh_m_tmp
      complex*16, allocatable, dimension(:,:) :: pg_m_tmp, pg_m_tmp2
      complex*16,allocatable,dimension(:,:) :: sug_m,spg_m
      complex*16 beta_psi_tmp(nref_tot,mx/nnodes+1)

      if(mst.ne.mx) then
      write(6,*) "mst.ne.mx,stop"
      stop
      endif

      sumdum_m = 0.d0
      sumdum2_m = 0.d0

      lin_st = 0

      ng_n=ngtotnod(inode,kpt)

      if(ipsp_all.eq.1) then
      allocate(spg_m(1,nblock_band_mx))
      allocate(sug_m(1,nblock_band_mx))
      else
      allocate(spg_m(mg_nx,nblock_band_mx))
      allocate(sug_m(mg_nx,nblock_band_mx))
      sug_m=0.d0
      spg_m=0.d0
      endif

      err2_m=1.d0
      pg_m = (0.0d0,0.0d0)
      pg_old_m = (0.0d0,0.0d0)


ccccccccccccccccccc  subspace diagonalization at the beginning of the CG_AllBand
cccccc assume the input wave functions are orthogonal

       c_one=dcmplx(1.d0,0.d0)
       c_zero=dcmplx(0.d0,0.d0)

            call Hpsi_comp_AllBandBP(ug_n_bp,ugh_m,
     &           nblock_band_mx,ilocal,vr,workr_n,kpt,1,
     &           sug_m,sumdum_m,iislda)


            s_m = 0.d0
            do m=1, nblock_band_mx
               if(ipsp_all.eq.1) then
               do i=1,ng_n
                  s_m(m)=s_m(m)+cdabs(ug_n_bp(i,m))**2
               enddo
               else
               do i=1,ng_n
                  s_m(m)=s_m(m)+real(ug_n_bp(i,m)*dconjg(sug_m(i,m)))
               enddo
               endif
            enddo

            call mpi_allreduce(s_m,tmp_real,
     &           nblock_band_mx,MPI_REAL8,MPI_SUM,MPI_COMM_K,ierr)
            s_m=tmp_real

            do m=1, nblock_band_mx
               s_m(m)=dsqrt(1.d0/(s_m(m)*vol))
               sumdum_m(:,1:natom,m)=s_m(m)*sumdum_m(:,1:natom,m)
            enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      beta_psi_tmp = (0.d0,0.d0)
      mx_n=mx/nnodes+1

      do m1=1,nblock_band_mx
      m=m1+icolor_b*nblock_band_mx
         if((m-1)/mx_n+1.eq.inode) then
            im=m-(inode-1)*mx_n
            iref_t=0
            do ia=1,natom 
               do iref=1,numref(ia)
                  iref_t=iref_t+1
                  beta_psi_tmp(iref_t,im)
     &                 =sumdum_m(iref,ia,m1)
               enddo
            enddo
         endif
      enddo
      nsize_beta_psi = nref_tot*(mx/nnodes+1)

      call mpi_allreduce(beta_psi_tmp,beta_psi,nsize_beta_psi,
     &     MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_B2,ierr)
cccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(ipsp_all.eq.2) then
        deallocate(sug_m)
        deallocate(spg_m)
        endif

***********************************************

      return
      end

