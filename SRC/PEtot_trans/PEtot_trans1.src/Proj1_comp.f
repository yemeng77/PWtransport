      subroutine Proj1_comp(ilocal,nline,tol,
     &  E_st,err_st,ave_line,vr,workr_n,mbad,ave_err,kpt,iislda)
****************************************
cc     Written by Lin-Wang Wang, March 30, 2001. 
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************
****************************************
****   The main conjugate gradient routine
****   mbad=0, normal iCG=1 run, ave_err not used
****   mbad>0, called from CG_new, fix the states: mx-mbad,mx to ave_err
******************************************
******************************************

      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)

      include 'param.escan_real'
***********************************************
       complex*16 pg(mg_nx),ugh(mg_nx),pgh(mg_nx)
       complex*16 pg_old(mg_nx),ughh_old(mg_nx)
       complex*16,allocatable,dimension(:)  ::  spg

       real*8 Dij0(32,32,mtype),Qij(32,32,mtype)
       integer isNLa(9,matom),ipsp_type(mtype)

       real*8 vr(mr_n)
       real*8 prec(mg_nx)
       complex*16 workr_n(mr_n)
       complex*16 cc
**********************************************
**** if iopt=0 is used, pghh_old can be deleted
**********************************************
        integer lin_st(mst)
	real*8 E_st(mst),err_st(mst),err2_st(mst)
	real*8 Ef
        complex*16 sumdum(mref,matom),sumdum2(mref,matom)
        real*8 occ_t(mtype)
       integer iiatom(mtype),icore(mtype),numref(matom),ityatom(mtype)

       common /com123b/m1,m2,m3,ngb,nghb,ntype,rrcut,msb
       common /comEk/Ek
       common /comisNLa/isNLa,Dij0,Qij,ipsp_all,ipsp_type
       common /comNL2/occ_t,iiatom,icore,numref,ityatom

       ng_n=ngtotnod(inode,kpt)

************************************************
**** ugh = (H-Eref) * ug
************************************************
cccc This subroutine just go through all the wave functions, 
cccc calculate the project sumdum, store in beta_psi, to be output
        do 2000 m=1,mx

        call Hpsi_comp(ug_n_bp(1,m),ugh,ilocal,vr,workr_n,kpt,
     &   1,sug_n(1,m),sumdum,iislda)

ccccccccccccccc normalization
      s=0.d0
      if(ipsp_all.eq.1) then
      do i=1,ng_n
      s=s+cdabs(ug_n_bp(i,m))**2
      enddo
      else
      do i=1,ng_n
      s=s+real(ug_n_bp(i,m)*dconjg(sug_n(i,m)))
      enddo
      endif

      call global_sumr(s)
      s=dsqrt(1.d0/(s*vol))
      fnorm=s
ccccccccccccccccccccccccccccccccccccccc
        sumdum=fnorm*sumdum

cccccc  store sumdum in beta_psi
       mx_n=mx/nnodes+1
       if((m-1)/mx_n+1.eq.inode) then
       im=m-(inode-1)*mx_n
       iref_t=0
       do ia=1,natom 
       do iref=1,numref(ia)
       iref_t=iref_t+1
       beta_psi(iref_t,im)=sumdum(iref,ia)
       enddo
       enddo
       endif
**********************************************
2000  continue


***********************************************

      return
      end

