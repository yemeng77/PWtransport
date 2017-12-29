      subroutine Etranscalc(xatom,fatom,workr_n,E_tot,
     & iforce_cal,ido_rho,ido_vr,tolug,tolE,niter,nline,
     &  iCGmth,iscfmth,FermidE,mCGbad,E_st,err_st,AL,
     &  nkpt,ntype,convergE,islda,igga,
     &  n1w,dV,dE_evan,dk_same)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
******************************************

******************************************
cccccc this program reads the initial wavefunction from  ugIO(ug_n,kpt,2,0,iislda)
cccccc and write the converged wavefunction to           ugIO(ug_n,kpt,1,0,iislda)
******************************************

      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)

      include 'mpif.h'

      include 'param.escan_real'
******************************************
       real*8 AL(3,3),ALt(3,3)
********************************************
       real*8 E_st(mst,nkpt,islda),err_st(mst,nkpt,islda) 
       real*8 occ(mst,nkpt,islda)
       real*8 eigen(mst),Ewind(2)


       real*8,allocatable,dimension(:)  :: xyzmaptmp
       real*8,allocatable,dimension(:,:)  :: wmasktmp
       complex*16,allocatable,dimension(:)  :: v_dv
c       real*8,allocatable,dimension(:)  :: v_dv
       integer,allocatable,dimension(:) :: indmtmp

       real*8 xatom(3,matom)
       real*8 fatom(3,matom)

       real*8 occ_t(mtype)
       integer iiatom(mtype),icore(mtype),numref(matom)
       integer is_ref(mtype),ip_ref(mtype),id_ref(mtype)

       complex*16,allocatable,dimension(:,:) :: wgp_n

       integer nmap(matom)

       real*8 zatom(mtype)

       integer iCGmth(100),iscfmth(100)
       real*8 FermidE(100)
cccccccccccccccccccccccc

       integer iatom(matom),ityatom(matom)
       real*8  totNel
       integer smatr(3,3,48),nrot,ilocal
       real*8 Ealpha(mtype)
*************************************************
       character*20 vwr_atom(mtype)
       character*20 fwg_in,fwg_out,fdens_in,fdens_out
       character*20 file_dV
*************************************************

       common /comNL2/occ_t,iiatom,icore,numref
       common /comispd_ref/is_ref,ip_ref,id_ref
       common /comEk/Ek
       common /comzatom/zatom
       common /comnmap/nmap

       common /comMainEtot/iatom,ityatom,totNel,mxlow,
     &     smatr,nrot,ilocal,Ealpha,Ealphat

       common /comVext/ivext_in,ivext_dir,xvext_c,dv_mix,
     &   nite_mix,dv_jump_init
c

c saved arrays for mch_pulay
c
       real*8 AA(nmax,nmax)


       complex*16 workr_n(mr_n)

**************************************************
c initialize mpi and number each node
c
       pi=4.0d0*datan(1.0d0)

       E_NSC0 = 0.0d0  !amc
       TS0= 0.0d0
       n33 = 0 
       E_COR0 = 0.0d0 
       E_extV0= 0.0d0
       E_HXC0 = 0.0d0
       E_TOT0 = 0.0d0
       DVE0 = 0.0d0
       E_IVext=0.d0
       E_rhoVext=0.d0
       E_IVextau=0.d0
       E_rhoVextau=0.d0
       E_rhoVextau2=0.d0
       E_Vext0=0.d0

       Ek=0.5d0
       occ=1.d0

       call getewald(fatom,xatom,AL,ityatom,ewald)

       if(ivext_in.eq.1.or.ivext_in.eq.3) then
       call getEextV(fatom,xatom,AL,ityatom,E_IVext,vext_n)
       endif


       if(iforce_cal.eq.0) fatom=0.d0

       call getVrho(AL,vion_n,vionT_n,xatom,ntype,iatom,
     &  rhocr_n,totNel,ido_rho,workr_n,islda)

        if(ivext_in.eq.1.or.ivext_in.eq.3) then
        vion_n=vion_n+vext_n
        endif

*********************************************************
****  at this step, rho_n has already been obtained in any case of ido_rho
*********************************************************
      if(ivext_in.eq.2.or.ivext_in.eq.3) then
      iflag=0
      call add_slab_dipV(iflag,
     & xatom,ityatom,totNel,AL,islda,dv_jump,dv_jump_curr,
     & E_IVextau,E_rhoVextau,E_rhoVextau2,fatom,vion_n)

      if(inode.eq.1) then
      write(6,419) dv_jump
      write(22,419) dv_jump 
      write(22,*) "----------------------------------------"
      endif
419   format("initial dv_jump=", f14.6)

****  calculate dv_jump, and add it onto vion_n
      endif
*********************************************************
***  ilocal.eq.3, q space nolocal Kleimen_Bylander
*********************************************************
       if(inode.eq.1) write(6,*) "into getwmask, or getwq"

       if(ilocal.eq.3) then
       do kpt=1,nkpt
       call getwq(AL,ntype,iatom,ityatom,xatom,kpt)
       call wqIO(nkpt,kpt,1)     ! 1: write, 2: read
       enddo
       endif

*********************************************************
***  ilocal.eq.2, r space nolocal Kleimen_Bylander
*********************************************************
       if(ilocal.eq.2) then

ccccccc this formula sometime is not right, if many atoms 
ccccccc are located inside one processor

       mrb2=2*(4*pi/3*rcut**3)/(vol/(n1*n2*n3))
ccccc special version, for not equally distributed systems
       mrb2=3*mrb2

       mrb2_matom_node=mrb2*natom/nnodes

       allocate(wmask(9*mrb2_matom_node))
       allocate(xyzmap(3*mrb2_matom_node))
       allocate(cphase(mrb2_matom_node))
       allocate(indm(mrb2_matom_node))
       allocate(wmasktmp(9,mrb2))
       allocate(xyzmaptmp(3*mrb2))
       allocate(indmtmp(mrb2))

       iatsum = 0
       iatsum2=0 
       do 20 ia=1,natom
       iitype=ityatom(ia)
       call getwmask(xatom(1,ia),nmap(ia),indmtmp,
     &  ityatom(ia),wmasktmp,xyzmaptmp,AL,workr_n,mrb2,
     &  is_ref(iitype),ip_ref(iitype),id_ref(iitype),nref)
        numref(ia)=nref
       if(iatsum+nmap(ia).gt.mrb2_matom_node) then
       write(6,*) "iatsum.gt.mbr2_matom_node, stop",
     & ia,iatsum+nmap(ia),mbr2_matom_node 
       write(6,*) "atom might not be equally distributed"
       call mpi_abort(MPI_COMM_WORLD,ierr)
       endif

       i00=3*iatsum
         do i=1,nmap(ia)
         do j=1,nref
         iatsum2=iatsum2+1
         wmask(iatsum2)=wmasktmp(j,i)
         enddo
         enddo
         do i=1,nmap(ia)
         indm(i+iatsum)=indmtmp(i)
         enddo
         do i=1,3*nmap(ia)
         xyzmap(i+i00)=xyzmaptmp(i)
         enddo
       iatsum = iatsum + nmap(ia)
       
20     continue

       s = dfloat(iatsum)
       call global_sumr(s)
       avenmap = s/dfloat(natom)

       if(inode.eq.1) then
       write(6,*) "ave nmap=",avenmap
       endif

       endif !ilocal.eq.2


       call mpi_barrier(MPI_COMM_WORLD,ierr)

**************************************************
**** initial potential vr_in
**************************************************
       if(ido_vr.eq.1) then
       if(islda.eq.1.and.igga.eq.0) then
       call getpot2(rho_n(1,1),vion_n,rhocr_n,vr_in_n(1,1),
     &        v0,E_Hxc,workr_n,E_coul,E_ion)
       endif
       if(islda.eq.2.and.igga.eq.0) then
       call getpot3(rho_n,vion_n,rhocr_n,vr_in_n,v0,
     &           E_Hxc,workr_n,E_coul,E_ion)   
       endif

       if(islda.eq.1.and.igga.eq.1) then
       call getpot4(rho_n,vion_n,rhocr_n,vr_in_n,v0,
     &           E_Hxc,workr_n,E_coul,E_ion)   
       endif

       if(islda.eq.2.and.igga.eq.1) then
       call getpot5(rho_n,vion_n,rhocr_n,vr_in_n,v0,
     &           E_Hxc,workr_n,E_coul,E_ion)   
       endif
       
       endif
**************************************************
**** end initialize 
**************************************************

       if(inode.eq.1) write(6,*) "start iterations"

       nscf=0
       do 2000 nint=1,niter

******************************************************
**** tmperarily use vr_out as a work array
         do iislda=1,islda
         do i=1,nr_n
         vr_out_n(i,iislda)=rho_n(i,iislda)
         enddo
         enddo
******************************************************
**** wavefunction update calculations
**** vr_in & mx, are used in CG_real, CG_new, diag_comp.
****  CG_real, CG_new output the rho(r) of mx states.
******************************************************
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


       ave_line=0.d0
c       do 201 iislda=1,islda
c       do 200 kpt=1,nkpt
       iislda=1
       kpt=1

       if(ilocal.eq.2) then
       call getcphase()
       endif

       call gen_G_comp(kpt,0) 
       call fftprep_comp(n1,n2,n3)

       call ugIO(ug_n,kpt,2,0,iislda)

       if(ilocal.eq.3) then
       call wqIO(nkpt,kpt,2)
       endif



       dE=E2_trans-E1_trans
       if(N_trans.gt.1) then
       dE=dE/(N_trans-1)
       else
       dE=0.d0
       endif


       open(16,file="eigen_wg0")
       rewind(16)
       read(16,*) (eigen(i),i=1,mx)
       close(16)
       do i=1,mx
       eigen(i)=eigen(i)/27.211396d0
       enddo


       if(inode.eq.1) then
       open(17,file="wr_real.E",form="unformatted")
       rewind(17)
       write(17) N_trans
       write(17) n1,n2,n3,nnodes
       endif

       call mpi_barrier(MPI_COMM_WORLD,ierr)


       do 6000 iii=1,N_trans
       Eref=(E1_trans+dE*(iii-1))/27.211396d0


cccc determine mstate
cccccc write out the wave function first in file, then 
cccccc read it. It is a bit stupid, but it is okay
cccccccc this will be the minimum change from the current file. 
        if(inode.eq.1) then
        Ew=Eref*27.211396d0    ! change back to eV
        call find_wrl(Ew,n1w,n1,n2,n3,nnodes,dV,dE_evan,
     &  dk_same,AL,mstateT)
        write(17) iii,Eref*27.211396d0,mstateT
        endif

        call mpi_barrier(MPI_COMM_WORLD,ierr)
        call mpi_bcast(mstateT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        mstate=mstateT



        allocate(wgp_n(mg_nx,mstate))
	do i1=1,mstate
	do i2=1,mg_nx
	wgp_n(i2,i1)=dcmplx(0.d0,0.d0)
	enddo
	enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccc  This following is to read in the perturbation 
ccccc  wave function wl, and smooth it using the Fourier filter
ccccc  there will be mstate wl

ccccccccccccccccccccccccccccccccccccccccccccc


       do istate=1,mstate

       ii1=istate/10
       ii2=istate-ii1*10

       file_dV="wr_test."//char(48+ii1)//char(48+ii2)

       call rhoIO_comp(AL,workr_n,mr_n,2,file_dV)

       call d3fft_comp(wgp_n(1,istate),workr_n,1,kpt)

       ng_n=ngtotnod(inode,kpt)

ccccccccc  Ecut is passed in from param.escan_real
       beta=0.5

       do i=1,ng_n
       if(gkk_n(i,kpt).lt.Ecut*beta) then
       wgp_n(i,istate)=wgp_n(i,istate)
       else
       wgp_n(i,istate)=0.d0
       endif
       enddo 

       call d3fft_comp(wgp_n(1,istate),workr_n,-1,kpt)


ccccccccccccccccccc
       i0=n1-12+1
       dq=dsqrt(2*Ecut)-dsqrt(2*Ecut*beta)
       alpha=dq**2/16*
     &      (AL(1,1)**2+AL(2,1)**2+AL(3,1)**2)/n1**2

       alpha=alpha*2.d0


       do ii=1,nr_n
       jj=ii+(inode-1)*nr_n
       i=(jj-1)/(n2*n3)+1
       j=(jj-1-(i-1)*n2*n3)/n3+1
       k=jj-(i-1)*n2*n3-(j-1)*n3

       if(i.lt.n1/2+1) i=i+n1
       workr_n(ii)=workr_n(ii)*dexp(-(i-i0)**2*alpha)
       enddo
       call d3fft_comp(wgp_n(1,istate),workr_n,1,kpt)
ccccccccccccccccccccccccccccccccccccccccccccccccccc
       enddo    ! istate=1,mstate
ccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccc Inside CG_linear, the wr_real.E has already been written in real space.
       mp=2000
       dE=1.0/27.211396d0
       Ewind(1)=Eref-dE
       Ewind(2)=Eref+dE

       call eigen_comp(ilocal,nline,mp,tolE,
     &  vr_in_n(1,iislda),workr_n,kpt,Ewind,
     &  mxc,eigen)

c       call CG_linear(ilocal,nline,tolug,
c     &   wgp_n,vr_in_n(1,iislda),workr_n,
c     &   kpt,Eref,AL,eigen,mxlow,mstate)

       call mpi_barrier(MPI_COMM_WORLD,ierr)

       call flush(17)

        deallocate(wgp_n)

6000   continue

       if(inode.eq.1) then
       close(17)
       endif

       stop
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       call ugIO(ug_n,kpt,1,0,iislda)

       ave_line=ave_line+ave_linek
200    continue
201    continue

           ave_line=ave_line/nkpt/islda

           errmax=-100.d0
           errave=0.d0
           do iislda=1,islda
           do kpt=1,nkpt
           do m=1,mx
           if(err_st(m,kpt,iislda).gt.errmax) 
     &      errmax=err_st(m,kpt,iislda)
           errave=errave+err_st(m,kpt,iislda)
           enddo
           enddo
           enddo
           errave=errave/(mx*nkpt*islda)

********************************************************
         if(iscfmth(nint).le.0)  goto 99    ! non-self-consistent
           nscf=nscf+1

           call occup(FermidE(nint),E_st,totNel,nkpt,
     &       occ,E_NSC,TS,workr_n,Ef,islda)
          if(inode.eq.1) then
          write(6,*) "Ef=", Ef
          endif

    
         if(nrot.gt.1) then
         do iislda=1,islda
         call symmop(rho_n(1,iislda),workr_n)
         enddo
         endif

******************************************************
**************************************************
***** scf calc., updating vr_in 
**************************************************
      if(ivext_in.eq.2.or.ivext_in.eq.3) then
      iflag=1
      call add_slab_dipV(iflag,
     & xatom,ityatom,totNel,AL,islda,dv_jump,dv_jump_curr,
     & E_IVextau,E_rhoVextau,E_rhoVextau2,fatom,vion_n)

**** substract the old dv_jump from vion_n and fatom
**** add in the newly calculated dv_jump in vion_n and fatom
      endif
**************************************************
***  calculate the total energy E (or say E_trial)
***  and dvE=\int (v_out-v_in)^2 d^3r
**************************************************

       E_cor=0.d0
       drho=0.d0
       E_rhoVext=0.d0

       do iislda=1,islda
       do i=1,nr_n
       E_cor=E_cor+rho_n(i,iislda)*
     &              (vion_n(i)-vr_in_n(i,iislda))
       drho=drho+dabs(rho_n(i,iislda)-vr_out_n(i,iislda))
       E_rhoVext=E_rhoVext+rho_n(i,iislda)*Vext_n(i)
       enddo
       enddo

       call global_sumr(E_cor)
       call global_sumr(drho)
       call global_sumr(E_rhoVext)

       E_cor=E_cor*vol/nr
       drho=drho*vol/nr
       E_rhoVext=E_rhoVext*vol/nr
       E_cor=E_cor-E_rhoVext-E_rhoVextau2    ! E_cor is really vion(Z only)-vtot

       if(islda.eq.1.and.igga.eq.0) then
       call getpot2(rho_n(1,1),vion_n,rhocr_n,vr_out_n(1,1),
     &       v0,E_Hxc,workr_n,E_coul,E_ion)
       endif
       if(islda.eq.2.and.igga.eq.0) then
       call getpot3(rho_n,vion_n,rhocr_n,vr_out_n,v0,
     &           E_Hxc,workr_n,E_coul,E_ion)   
       endif
       if(islda.eq.1.and.igga.eq.1) then
       call getpot4(rho_n,vion_n,rhocr_n,vr_out_n,v0,
     &           E_Hxc,workr_n,E_coul,E_ion)   
       endif
       if(islda.eq.2.and.igga.eq.1) then
       call getpot5(rho_n,vion_n,rhocr_n,vr_out_n,v0,
     &           E_Hxc,workr_n,E_coul,E_ion)   
       endif

******  E_extV is the energy due to all the external potential corrections

       E_extV=E_IVext+0.5d0*E_IVextau+E_rhoVext+
     &                               0.5d0*E_rhoVextau

       E_tot=E_NSC+E_cor+E_Hxc+ewald+Ealphat+TS
       E_tot=E_tot+E_extV

       dvE=0.d0
       do iislda=1,islda
       do i=1,nr_n
       dvE=dvE+(vr_in_n(i,iislda)-vr_out_n(i,iislda))**2
       enddo
       enddo

       call global_sumr(dvE)

       dvE=dvE/nr
       dv_ave=dsqrt(dvE)

      if(islda.eq.2) then
      scharge1=0.d0
      scharge2=0.d0
      do ikpt=1,nkpt
      do m=1,mx
      scharge1=scharge1+occ(m,ikpt,1)
      scharge2=scharge2+occ(m,ikpt,2)
      enddo
      enddo
      dcharge=0.d0
      do i=1,nr_n
      dcharge=dcharge+dabs(rho_n(i,1)-rho_n(i,2))
      enddo
      call global_sumr(dcharge)
      dcharge=dcharge*vol/nr
      endif

       if(inode.eq.1) then

       write(6,*) "***************************************************"
       write(6,*) "iter=", nint
       write(6,*) "---------------------------------------------------"
       if(islda.eq.2) then
       write(6,411) scharge1,scharge2,dcharge
       write(6,*) "---------------------------------------------------"
       endif
       write(6,401) dvE, dvE-dvE0
       write(6,402) dv_ave,drho
       write(6,403) erraved,errmaxd
       write(6,404) errave,errmax
       write(6,*) "---------------------------------------------------"
       write(6,398) ewald
       write(6,399) Ealphat
       write(6,412) E_extV,E_extV-E_extV0 
       write(6,405) E_NSC, E_NSC-E_NSC0
       write(6,406) E_cor, E_cor-E_cor0
       write(6,407) E_Hxc, E_Hxc-E_Hxc0
       write(6,408) TS, TS-TS0
       write(6,409) E_tot, E_tot-E_tot0
       write(22,*) "-------------------------------------------"
       write(22,410) nint,ave_line,errave,erraved
       if(islda.eq.2) then
       write(22,411) scharge1,scharge2,dcharge
       endif
       if(ivext_in.eq.2.or.ivext_in.eq.3) then
       write(22,417) dv_jump,dv_jump_curr
       endif
       write(22,409) E_tot, E_tot-E_tot0
       call system_flush(22)
411   format("charges of spin up,down and loc_diff ", 3(f18.10,2x))
410   format("iter=",i4,"  ave_lin=", f5.1, "  ugerr CG=", E8.1,
     &         "  ugerr Diag=", E8.1 )
417   format("dv_jump(update,rho-calc)=",2(f12.5,1x))
 
      endif

       convergE=dabs(E_tot-E_tot0)
       if(convergE.lt.tolE) goto 2001     ! finished

       if(nint.ne.niter) then
       E_NSC0=E_NSC
       E_cor0=E_cor
       E_extV0=E_extV
       E_Hxc0=E_Hxc
       TS0=TS
       E_tot0=E_tot
       dvE0=dvE
       endif
**************************************************
*** mch_pulay: input the current vr_in->vr_out,
*** using previous such relation, and make a linear combination
*** output: a new vr_in->vr_out, with smaller difference
**************************************************
*** mch_kerk, Thomas3: input the vr_in -> vr_out, 
*** Using non linear-combination prediction (preconditioning)
*** output: a new vr_in->vr_out, with smaller difference.
***  the vr_in is used in the next iteration
**************************************************

******************************************************
       call mch_pulay(vr_in_n,vr_out_n,nscf,AA,nreset,islda)

       do iislda=1,islda
       if(iscfmth(nint).eq.1) then
       call mch_kerk(vr_in_n(1,iislda),vr_out_n(1,iislda),
     &     workr_n)
       else
       call Thomas3(vr_in_n(1,iislda),vr_out_n(1,iislda),
     &      Ef,totNel,workr_n,islda)
       endif
       enddo

***************************************************
*** end selfconsistent updating vr_in
***************************************************

99     continue      ! jumping point for non-self-consistency



       do 198 iislda=1,islda
       do 199 kpt=1,nkpt

       if(ilocal.eq.2) then
       call getcphase()
       endif

       call gen_G_comp(kpt,0) 
       call fftprep_comp(n1,n2,n3)

       call ugIO(ug_n,kpt,2,0,iislda)
 
       if(ilocal.eq.3) then
       call wqIO(nkpt,kpt,2)
       endif

       call diag_comp(ilocal,E_st(1,kpt,iislda),
     &  err_st(1,kpt,iislda),vr_in_n(1,iislda),workr_n,kpt)

       call ugIO(ug_n,kpt,1,0,iislda)

199    continue
198    continue

          errmaxd=-100.d0
          erraved=0.d0
          do iislda=1,islda
          do kpt=1,nkpt
          do m=1,mx
          if(err_st(m,kpt,iislda).gt.errmaxd)
     &                 errmaxd=err_st(m,kpt,iislda)
          erraved=erraved+err_st(m,kpt,iislda)
          enddo
          enddo
          enddo
          erraved=erraved/(mx*nkpt*islda)


2000   continue  
******************************************************
2001   continue    ! jump out point for convergE.lt.tolE

      if(inode.eq.1.and.nscf.gt.0) then
     
       write(22,*) "---------------------------------------------------"
       write(22,*) "E_Fermi(eV)=", Ef*27.211396d0
       write(22,*) "---------------------------------------------------"
       if(islda.eq.2) then
       write(22,411) scharge1, scharge2,dcharge
       write(22,*) "---------------------------------------------------"
       endif
       write(22,401) dvE, dvE-dvE0
401   format(" dvE, dvE(n)-dvE(n-1) = ", 2(E10.4,1x))
       write(22,402) dv_ave,drho
402   format(" dv_ave, drho_tot     = ", 2(E10.4,1x))
       write(22,403) erraved,errmaxd
403   format(" ug err,diag,[ave,max]= ", 2(E10.4,1x))
       write(22,404) errave,errmax
404   format(" ug err, CG ,[ave,max]= ", 2(E10.4,1x))
       write(22,*) "---------------------------------------------------"
       write(22,398) ewald
398   format(" Ewald        = ", E20.14)
       write(22,399) Ealphat
399   format(" Alpha        = ", E20.14)
       write(22,412) E_extV,E_extV-E_extV0
412   format(" E_extV       = ", E20.14,4x,E10.4)
       write(22,405) E_NSC, E_NSC-E_NSC0
405   format(" E_NSC        = ", E20.14,4x,E10.4)  
       write(22,406) E_cor, E_cor-E_cor0
406   format(" E[vion-vtot] = ", E20.14,4x,E10.4)  
       write(22,407) E_Hxc, E_Hxc-E_Hxc0
407   format(" E_Hxc        = ", E20.14,4x,E10.4)  
       write(22,408) TS, TS-TS0
408   format(" TS           = ", E20.14,4x,E10.4)  
       write(22,409) E_tot, E_tot-E_tot0
409   format(" E_tot        = ", E20.14,4x,E10.4)  
       write(22,*) "---------------------------------------------------"
       write(22,396) E_coul,E_Hxc-E_coul,E_ion
       write(22,413) E_rhoVext,E_IVext
       write(22,418) E_rhoVextau,E_IVextau
       write(22,*)"E_extV=E_rhoVext+E_IVext+0.5*(E_rhoVextau+E_IVextau)"
       write(22,*) "---------------------------------------------------"
       if(ivext_in.eq.2.or.ivext_in.eq.3) then
       write(22,414) ivext_dir,xvext_c,dv_jump
       write(22,*) "Dipole/unit_area = dv_jump/4pi"
       write(22,*) "---------------------------------------------------"
       endif
396   format(" E_Hart,E_xc,E_ion     =", 3(E16.10,2x))
413   format(" E_rhoVext,E_IVext     =", 2(E16.10,2x))   
418   format(" E_rhoVextau,E_IVextau =", 2(E16.10,2x))   
414   format(" ivext_dir,xvext_c,dv_jump=",i3,1x,f10.6,2x,E15.5)


c       do kpt=1,nkpt
c       write(22,*) "kpt= ", kpt
c       write(22,*) "err of each states, A.U"
c       write(22,102) (err_st(i,kpt), i=1,mx)
c       write(22,*) "eigen energies, in eV"
c       write(22,103) (E_st(i,kpt)*27.211396d0, i=1,mx)
c       write(22,*) "*********************************"
c       enddo

        endif
101   format(5(i6,7x))
102   format(5(E10.4,3x))
103   format(5(f12.8,1x))



******************************************************
*** calculate forces on each atom
****************************************************
**** using vr_out_n(i,1) and vionT_n(i) as working array, 
**** their contents will be destroyed. Bad practice
******************************************************
      if(iforce_cal.eq.1) then

      if(islda.eq.1.and.igga.eq.0) then
      do i=1,nr_n
      vr_out_n(i,1)=rho_n(i,1)
      vionT_n(i)=UxcCA(dabs(rho_n(i,1)+rhocr_n(i)),uxc2)
      enddo
      endif

      if(islda.eq.2.and.igga.eq.0) then
      do i=1,nr_n
      vr_out_n(i,1)=rho_n(i,1)+rho_n(i,2)
      call UxcCA2(dabs(rho_n(i,1)+rhocr_n(i)*0.5d0),
     &      dabs(rho_n(i,2)+rhocr_n(i)*0.5d0),
     &    vxc1,vxc2,uxc1,uxc2)
      vionT_n(i)=(vxc1+vxc2)/2
      enddo
      endif

      if(islda.eq.1.and.igga.eq.1) then
      s=0.d0
      do i=1,nr_n
      vr_out_n(i,1)=rho_n(i,1)
      s=s+dabs(rhocr_n(i))
      enddo
      call global_sumr(s)
      s=s*vol/nr
        if(s.gt.1.D-5) then
        call getpot4_force(rho_n,vionT_n,rhocr_n,workr_n)
        else
        vionT_n=0.d0
        endif
      endif
      

      if(islda.eq.2.and.igga.eq.1) then
      s=0.d0
      do i=1,nr_n
      vr_out_n(i,1)=rho_n(i,1)+rho_n(i,2)
      s=s+dabs(rhocr_n(i))
      enddo
      call global_sumr(s)
      s=s*vol/nr
        if(s.gt.1.D-5) then
        call getpot5_force(rho_n,vionT_n,rhocr_n,workr_n)
        else
        vionT_n=0.d0
        endif
      endif

        call forcLC(AL,vr_out_n(1,1),vionT_n,
     &   xatom,ntype,iatom,fatom,workr_n)


      if(ilocal.eq.3) then
      do iislda=1,islda
      do kpt=1,nkpt
       
      call ugIO(ug_n,kpt,2,0,iislda)
      call wqIO(nkpt,kpt,2)

      call forcNLq(fatom,occ(1,1,iislda),kpt,nkpt)
ccccccc the resulting force is only correct after symmforc, since
ccccccc only the unsymmetrized k-points are used

      enddo
      enddo
      endif

cccccccccccccccccccccccccccccccccccc
      if(ilocal.eq.2) then

      allocate(wmaskX(9*mrb2_matom_node))

      do 300 ixyz=1,3

       iatsum2=0
       do ia=1,natom
       iitype=ityatom(ia)
         call getwmaskX(xatom(1,ia),nmap(ia),indmtmp,
     &   ityatom(ia),wmasktmp,AL,workr_n,ixyz,mrb2,
     &    is_ref(iitype),ip_ref(iitype),id_ref(iitype),nref)

         if(nref.ne.numref(ia)) then
         write(6,*) "nref.ne.numref,getwmaskX,stop",nref,numref(ia)
         call mpi_abort(MPI_COMM_WORLD,ierr)
         endif

         do i=1,nmap(ia)
         do j=1,nref
         iatsum2=iatsum2+1
         wmaskX(iatsum2)=wmasktmp(j,i)
         enddo
         enddo
       enddo  ! natom 

       do 301 iislda=1,islda
       do 302 kpt=1,nkpt

       call getcphase()
       call gen_G_comp(kpt,0) 
       call fftprep_comp(n1,n2,n3)

       call ugIO(ug_n,kpt,2,0,iislda)

       call forcNLr(workr_n,fatom,
     &  occ(1,1,iislda),kpt,nkpt,ixyz)

302    continue
301    continue
300    continue
       deallocate(wmaskX)

      endif
cccccccccccccccccccccccccccccccccccc

       if(nrot.gt.1) then
       call symmopf(smatr,nrot,AL,fatom,xatom,iatom)
       endif


      endif   ! iforce_cal.eq.1
**************************************************

  

      if(ilocal.eq.2) then
      deallocate(wmask)
      deallocate(xyzmap)
      deallocate(cphase)
      deallocate(indm)
      deallocate(wmasktmp)
      deallocate(xyzmaptmp)
      deallocate(indmtmp)
      endif

      
      return

      contains
      
*****************************************************
*****************************************************
       subroutine getcphase()
       implicit double precision (a-h,o-z)

       complex*16 cai
    
       cai=dcmplx(0.d0,1.d0)

       ico1=0
       do ia=1,natom
        do i=1,nmap(ia)
        ico1=ico1+1
        cphase(ico1)=cdexp(cai*(xyzmap(ico1*3-2)*akx(kpt)+
     &   xyzmap(ico1*3-1)*aky(kpt)+xyzmap(ico1*3)*akz(kpt)))
        enddo
       enddo

       return
       end subroutine getcphase

*****************************************************

        end
      
