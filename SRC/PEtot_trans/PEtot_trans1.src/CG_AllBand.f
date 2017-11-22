      subroutine CG_AllBand(ilocal,nline,tol,
     &     E_st,err_st,ave_line,vr,workr_n,mbad,ave_err,kpt,iislda,
     &     imth)
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

c      complex*16 U_cholesky(mx,mx)
c      complex*16 h(mst,mst),hh(mst,mst),chh(mst,mst)
      complex*16,allocatable,dimension(:,:) :: hh

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

      allocate(hh(mst,mst))

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


      do i=1,ng_n
         x=gkk_n(i,kpt)/Ek
         y=27.d0+x*(18.d0+x*(12.d0+x*8.d0))
         prec_m(i)=y/(y+16.d0*x**4)   
      enddo
      rr0_m = 1.d+40
ccccccccccccccccccc  subspace diagonalization at the beginning of the CG_AllBand
cccccc assume the input wave functions are orthogonal

       c_one=dcmplx(1.d0,0.d0)
       c_zero=dcmplx(0.d0,0.d0)



            call Hpsi_comp_AllBandBP(ug_n_bp,ugh_m,
     &           nblock_band_mx,ilocal,vr,workr_n,kpt,1,
     &           sug_m,sumdum_m,iislda)


      call dot_product_BP(ugh_m,ug_n_bp,ng_n,hh)


      call system_czheev('V','U',mx,hh,mst,E_st_tmp,workx,
     &     lwork,workrx,info,MPI_COMM_K1)


cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      E_st(band_dis(1):band_dis(2)) =
     &     E_st_tmp(band_dis(1):band_dis(2))*vol

       call rotate_wfBP(ug_n_bp,hh,ng_n,mg_nx)
        call rotate_wfBP(ugh_m,hh,ng_n,mg_nx)
cccc  sumdum_m is proportional to dconjg(ug_n_bp)!  
c        chh=dconjg(hh)
c        call rotate_wfBP(sumdum_m,chh,mref*natom,mref*natom)
         sumdum_m=dconjg(sumdum_m)
        call rotate_wfBP(sumdum_m,hh,mref*natom,mref*natom)
         sumdum_m=dconjg(sumdum_m)

           if(ipsp_all.eq.2) then 
           call rotate_wfBP(sug_m,hh,ng_n,mg_nx)
           endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      
      nltot = 0
      do 3000 nint2=1,nline
         nltot=nltot+1

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccc normalize  ug_n_bp

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
               do i=1,ng_n
                  ug_n_bp(i,m)=s_m(m)*ug_n_bp(i,m)
                  ugh_m(i,m)=s_m(m)*ugh_m(i,m)
               enddo
               if(ipsp_all.eq.2) then
               do i=1,ng_n
                  sug_m(i,m)=s_m(m)*sug_m(i,m)
               enddo
               endif
            enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         s_m=0.d0
         do m=1, nblock_band_mx
            do i=1,ng_n
               s_m(m)=s_m(m)+dreal(ugh_m(i,m)*dconjg(ug_n_bp(i,m)))
            enddo
         enddo
         
         call mpi_allreduce(s_m,tmp_real,
     &        nblock_band_mx,MPI_REAL8,MPI_SUM,MPI_COMM_K,ierr)
         s_m=tmp_real

         
         do m=1, nblock_band_mx
            s_m(m)=s_m(m)*vol
            E_ug_m(m)=s_m(m)
            eigen_m(m)=s_m(m)          
            err_m(m)=0.d0
            if(ipsp_all.eq.1) then
            do i=1,ng_n
               err_m(m)=err_m(m)
     &              +cdabs(ugh_m(i,m)-s_m(m)*ug_n_bp(i,m))**2
               pg_m(i,m)=ugh_m(i,m)-s_m(m)*ug_n_bp(i,m)
            enddo            
            else
            do i=1,ng_n
               err_m(m)=err_m(m)
     &              +cdabs(ugh_m(i,m)-s_m(m)*sug_m(i,m))**2
               pg_m(i,m)=ugh_m(i,m)-s_m(m)*sug_m(i,m)
            enddo            
            endif

         enddo                  ! do m

         call mpi_allreduce(err_m,tmp_real,
     &        nblock_band_mx,MPI_REAL8,MPI_SUM,MPI_COMM_K,ierr)
         err_m=tmp_real



         do m=1, nblock_band_mx
            err_m(m)=dsqrt(dabs(err_m(m)*vol))
            err2_m(m)=0.d0
            do i=1,ng_n
               err2_m(m)=err2_m(m)+cdabs(pg_m(i,m))**2
            enddo
         enddo

         call mpi_allreduce(err2_m,tmp_real,
     &        nblock_band_mx,MPI_REAL8,MPI_SUM,MPI_COMM_K,ierr)
         err2_m=tmp_real


         do m=1, nblock_band_mx
            err2_m(m)=dsqrt(dabs(err2_m(m)*vol))
************************************************
****  begin conjugate gradient
****  should I use pg**2 formula, or pg*s*pg formula here !!??
****  This is still not clear. 
************************************************

            rr1_m(m)=0.d0
            rr00_m(m)=0.d0
            do i=1,ng_n
               rr00_m(m)=rr00_m(m)+cdabs(pg_m(i,m))**2*prec_m(i)
               rr1_m(m)=rr1_m(m)+cdabs(pg_m(i,m))**2*prec_m(i)
            enddo
         enddo                  ! do m

         call mpi_allreduce(rr00_m,tmp_real,
     &        nblock_band_mx,MPI_REAL8,MPI_SUM,MPI_COMM_K,ierr)
         rr00_m=tmp_real

         call mpi_allreduce(rr1_m,tmp_real,
     &        nblock_band_mx,MPI_REAL8,MPI_SUM,MPI_COMM_K,ierr)
         rr1_m=tmp_real


**********************************************
****  here, pg == ughh
**********************************************

         do m=1, nblock_band_mx
            beta_m(m)=rr1_m(m)/rr0_m(m)
            rr0_m(m)=rr00_m(m)
         enddo     
**********************************************
***** pg(i) is the line minimization direction
**********************************************
***** preconditioning should be done before this
**********************************************


         do m=1, nblock_band_mx
            do i=1,ng_n
               pg_m(i,m)=-pg_m(i,m)*prec_m(i)
     &              +beta_m(m)*pg_old_m(i,m)
            enddo
         enddo



         if(ipsp_all.eq.1) then
         call orthogonal_projectionBP(pg_m,mg_nx,nblock_band_mx,
     &        ug_n_bp,mg_nx,nblock_band_mx,ng_n)
         else
         call orthogonal_projectionBP2(pg_m,mg_nx,nblock_band_mx,
     &        ug_n_bp,sug_m,mg_nx,nblock_band_mx,ng_n)
         endif


         pg_old_m = pg_m



**********************************************
***** pgh = H * pg
**********************************************

         call Hpsi_comp_AllBandBP(pg_m,pgh_m,
     &        nblock_band_mx,ilocal,vr,workr_n,kpt,1,
     &        spg_m,sumdum2_m,iislda) 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         s_m = 0.d0
         do m=1, nblock_band_mx
            if(ipsp_all.eq.1) then
            do i=1,ng_n
               s_m(m)=s_m(m)+cdabs(pg_m(i,m))**2
            enddo
            else
            do i=1,ng_n
               s_m(m)=s_m(m)+real(pg_m(i,m)*dconjg(spg_m(i,m)))
            enddo
            endif
         enddo

         call mpi_allreduce(s_m,tmp_real,nblock_band_mx,
     $        MPI_REAL8,MPI_SUM,MPI_COMM_K,ierr)
         s_m=tmp_real

         do m=1, nblock_band_mx
            s_m(m)=dsqrt(1.d0/(s_m(m)*vol))
            sumdum2_m(:,1:natom,m)=s_m(m)*sumdum2_m(:,1:natom,m)
            do i=1,ng_n
               pg_m(i,m)=s_m(m)*pg_m(i,m)
               pgh_m(i,m)=s_m(m)*pgh_m(i,m)
            enddo
            if(ipsp_all.eq.2) then
            do i=1,ng_n
               spg_m(i,m)=s_m(m)*spg_m(i,m)
            enddo
            endif
         enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
**********************************************
         do m=1, nblock_band_mx
            E_upg_m(m)=0.d0
            E_pg_m(m)=0.d0

            do i=1,ng_n
               E_upg_m(m)=E_upg_m(m)
     &              + dreal(pg_m(i,m)*dconjg(ugh_m(i,m)))
               E_pg_m(m)=E_pg_m(m)
     &              + dreal(pg_m(i,m)*dconjg(pgh_m(i,m)))
            enddo
         enddo                  ! do m

         call mpi_allreduce(E_upg_m,tmp_real,
     &        nblock_band_mx,MPI_REAL8,MPI_SUM,MPI_COMM_K,ierr)
         E_upg_m=tmp_real

         call mpi_allreduce(E_pg_m,tmp_real,
     &        nblock_band_mx,MPI_REAL8,MPI_SUM,MPI_COMM_K,ierr)
         E_pg_m=tmp_real

         do m=1, nblock_band_mx
            E_upg_m(m)=E_upg_m(m)*vol
            E_pg_m(m)=E_pg_m(m)*vol
**********************************************
****  calculate the theta from E_ug,E_pg,E_upg:
****  ug_new = ug*cos(th)+pg*sin(th)
****  E=E_ug*cos(th)^2+2*E_upg*cos(th)*sin(th)+E_pg*sin(th)^2
**********************************************

            theta_m(m)=0.5d0
     &           *dabs(datan(2.d0*E_upg_m(m)/(E_ug_m(m)-E_pg_m(m))))
            cos_th_m(m)=dcos(theta_m(m))
            sin_th_m(m)=dsin(theta_m(m))

            pred_E_m(m)=E_ug_m(m)*cos_th_m(m)**2
     &           +2*E_upg_m(m)*cos_th_m(m)*sin_th_m(m)
     &           +E_pg_m(m)*sin_th_m(m)**2

**********************************************
****  update ug using theta
**********************************************

            do i=1,ng_n
               ug_n_bp(i,m)=ug_n_bp(i,m)*cos_th_m(m)
     &              +pg_m(i,m)*sin_th_m(m)
               ugh_m(i,m)=ugh_m(i,m)*cos_th_m(m)      ! for the last step, this has not been updated yet
     &           + pgh_m(i,m)*sin_th_m(m)
            enddo

               sumdum_m(:,1:natom,m)=
     &              sumdum_m(:,1:natom,m)*cos_th_m(m)
     &              + sumdum2_m(:,1:natom,m)*sin_th_m(m)

            if(ipsp_all.eq.2) then
            do i=1,ng_n
               sug_m(i,m)=sug_m(i,m)*cos_th_m(m)
     &              +spg_m(i,m)*sin_th_m(m)
            enddo
            endif

         enddo                ! do m

************************************************
******check convergency
************************************************
         lin_st = 0
         E_st = 0.d0
         err_st = 0.d0
         err2_st = 0.d0
         lin_st(band_dis(1):band_dis(2))=nint2 
         E_st(band_dis(1):band_dis(2))
     &        =eigen_m(1:nblock_band_mx)
         err_st(band_dis(1):band_dis(2))
     &        =err_m(1:nblock_band_mx)
         err2_st(band_dis(1):band_dis(2))
     &        =err2_m(1:nblock_band_mx)
         
         call mpi_allreduce(lin_st,idumm,mst,MPI_INTEGER,
     &        MPI_SUM,MPI_COMM_B2,ierr)
         lin_st = idumm 


         call mpi_allreduce(E_st,dumm,mst,MPI_REAL8,
     &        MPI_SUM,MPI_COMM_B2,ierr)
         E_st = dumm


         call mpi_allreduce(err_st,dumm,mst,MPI_REAL8,
     &        MPI_SUM,MPI_COMM_B2,ierr)
         err_st = dumm


         call mpi_allreduce(err2_st,dumm,mst,MPI_REAL8,
     &        MPI_SUM,MPI_COMM_B2,ierr)
         err2_st = dumm
        

**********************************************
****  do 3000, is for the nline line minimization
**********************************************
c 3000 continue
c 3001 continue
      c_one=dcmplx(1.d0,0.d0)
      c_zero=dcmplx(0.d0,0.d0)


    

      if(imth.eq.3.or.all(err2_st.lt.tol).or.
     & nint2.eq.nline) then   

      if(ipsp_all.eq.1) then
      call orthogonal_choleskyBP(ug_n_bp,mg_nx,nblock_band_mx,
     &     ng_n,hh)    
       else
      call orthogonal_choleskyBP2(ug_n_bp,sug_m,mg_nx,
     &     nblock_band_mx,ng_n,hh,1)      ! iflag=1, sug_m is rotated    
 
       endif

      call rotate_wfBP(ugh_m,hh,ng_n,mg_nx)

      sumdum_m=dconjg(sumdum_m)
      call rotate_wfBP(sumdum_m,hh,mref*natom,mref*natom)
      sumdum_m=dconjg(sumdum_m)

      endif

      if(all(err2_st.lt.tol)) goto 3001

3000  continue
3001  continue



      call dot_product_BP(ugh_m,ug_n_bp,ng_n,hh)

      call system_czheev('V','U',mx,hh,mst,E_st_tmp,workx,
     &     lwork,workrx,info,MPI_COMM_K1)

      E_st(1:mst) =  E_st_tmp(1:mst)*vol

       call rotate_wfBP(ug_n_bp,hh,ng_n,mg_nx)

       sumdum_m=dconjg(sumdum_m)
       call rotate_wfBP(sumdum_m,hh,mref*natom,mref*natom)
       sumdum_m=dconjg(sumdum_m)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
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



!     End of subspace diagonalization.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(ipsp_all.eq.2) then
        deallocate(sug_m)
        deallocate(spg_m)
        endif



      Esum=0.d0
      ave_line=0.d0
      ave_err2=0.d0
      do i=1, mx
         Esum=Esum+E_st(i)
         ave_line=ave_line+lin_st(i)
         ave_err2=ave_err2+err2_st(i)
      enddo
      ave_line=ave_line/mx
      ave_err2=ave_err2/mx

c      if(nline.ne.0.and.inode_k.eq.1) then
c         write(icolor_k+100,*) "*********************************"
c         write(icolor_k+100,*) "**** kpt=", kpt
c         write(icolor_k+100,*) "the number of line minimization used"
c         write(icolor_k+100,101) (lin_st(i),i=1,mx)
c         write(icolor_k+100,104) ave_err2
c 104     format("the err of each states, A.U.    Aver err2=", E10.4)
c         write(icolor_k+100,102) (err_st(i), i=1,mx)
c         write(icolor_k+100,*) "the eigen energies, in eV"
c         write(icolor_k+100,103) (E_st(i)*27.211396d0, i=1,mx)
c         write(icolor_k+100,*) "*********************************"
c         write(icolor_k+100,*) "sum eigen(dbl occ, Ryd)=", Esum*2*2
c         write(icolor_k+100,*) "*********************************"
c      endif
c      call system_flush(icolor_k+100)

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


 101  format(5(i6,7x))
 102  format(5(E10.4,3x))
 103  format(5(f13.8,1x))

      deallocate(hh)


***********************************************

      return
      end

