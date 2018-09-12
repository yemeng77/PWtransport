      program analyWaveM
      implicit double precision (a-h,o-z)

      include "mpif.h"
      integer status(MPI_STATUS_SIZE)

      parameter (nm=1200)
      parameter (neven=4000)
      parameter (nstm=200)
      parameter (nwellm=400)

      real*8 AL(3,3)
      real*8 E_line(nm,nm),E_linew(nm,nm)
      real*8 dE_dk1(nstm),dE_dk2(nstm)

      integer  numw(nm),ist_linew(nm,nm),ikpt_linew(nm,nm)
      integer ist_evan(nevan),ikpt_evan(nevan),iGX_evan(nevan),
     &        iband_evan(nevan)
      real*8 E_evan(nevan),E_evanC(nevan)

      integer ind(nm,nm),iflag(nm,nm),jj_st(nm)
      real*8 dE(nm,nm),overlap_max(nm,nm),E0(nm,nm)

      integer i1_st1_ind(nstm),i1_st2_ind(nstm),
     &  imx(10,nstm),num_mx(nstm),nline_w1(nstm),nline_w2(nstm)
      integer i1_st1_ind0(nstm),i1_st1_ind1(nstm)
      integer i1_st2_ind0(nstm),i1_st2_ind1(nstm)
      integer ifl_evan1(nstm),ifl_evan2(nstm)
      complex*16 cphase_ucw1(nstm),cphase_ucw2(nwellm)
      real*8 ak_w1(nstm),ak_w2(nstm),E_w1(nstm),E_w2(nstm)
     
      integer idble1(nstm),idble2(nstm),iposit1(nstm),iposit2(nstm)

      complex*16 ccy2_st1(nwellm,nwellm),ccy2_st2(nwellm,nwellm)
      real*8 weight1(nwellm),weight2(nwellm),weight_st(nwellm)
      integer ind_Well(nwellm)
      complex*16 cc_R1(nwellm,nwellm),cc_R2(nwellm,nwellm)

      complex*16 cc,cc1,cc2,cc3,cc_st

      character*7 fileh

      complex*16, allocatable, dimension (:,:,:) :: uc_test,uc_test2,
     &   uc_test3
      complex*16, allocatable, dimension (:,:,:,:) :: uc
      complex*16, allocatable, dimension (:,:,:,:) :: ucw1,ucw2
      complex*16, allocatable, dimension (:,:,:) :: cphase
      complex*16, allocatable, dimension (:,:) :: cc_matrix,cc_matrix0
      complex*16, allocatable, dimension (:,:) :: cc_matrix0_tmp
      complex*16, allocatable, dimension (:) :: cc_y,cwork
      complex*16, allocatable, dimension (:) :: uc_tmp,uc_tmp2
      real*8, allocatable, dimension (:) :: ss,rwork
      real*8, allocatable, dimension (:,:,:) :: phase
      real*8, allocatable, dimension (:) :: sum_w1,sum_w2
      real*8, allocatable, dimension (:,:) :: sum_w4
      complex*16,allocatable,dimension(:) :: cc_w1,cc_w2
      integer,allocatable,dimension (:) :: ipiv2

      complex*16, allocatable, dimension (:,:) :: cc_s1,cc_s2
      complex*16, allocatable, dimension (:,:) :: U_s1,U_s2
      complex*16, allocatable, dimension (:) :: cc_y_tmp
      complex*16,allocatable,dimension(:) :: x_tmp,y_tmp,cwork_tmp
      real*8,allocatable,dimension(:) :: s_scale,rwork_tmp

      integer num_wr_dis(2)
      complex*16 tmp_complex(nwellm,nwellm)
      real*8 tmp_real(nwellm)

      call mpi_init(ierr)
      call mpi_comm_size(MPI_COMM_WORLD,nnodes,ierr)
      call mpi_comm_rank(MPI_COMM_WORLD,inode,ierr)      
      inode=inode+1
cccccccccccccccccccccccccccccccccccccccccccccc
cccc preprocessing for the system wavefunction information
ccccccc  n1,n2,n3 are the grid points for the system
ccccccc  n1w,n2,n3 are the grid points for one electrode unit cell
ccccccc the current run in n1,n1w direction
      if(inode.eq.1) write(6,*) "cccccc Start preprocessing  ccccc"

      open(10,file="analyWave.input")
      rewind(10)
      read(10,*) dV
      read(10,*) n1,n2,n3,n1w
      read(10,*) nnposit1,nnposit2
      read(10,*) ist_init,ist_final
      read(10,*) AL(1,1)
      read(10,*) num_iter_evan,dE_evanL,dE_evanR
      close(10)
      if(inode.eq.1) write(6,*) "Bias voltage in eV = ", dV
      a11=AL(1,1) 
      if(num_iter_evan.eq.1.and.inode.eq.1.) then
      write(6,*) 'ONLY RUNNING WAVES INCLUDED!!'
      endif
ccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccc
cccc preprocessing for the electrode wavefunction information
cccc These are the information for the electrode states, their energy, index in the wavefunction file, etc.
      open(14,file="E_line_W2K2")
      read(14,*) nstw,nkptw,nintep
      do j=1,nstw
      read(14,*) j1,numw(j)
      read(14,*) (E_linew(i,j),i=1,numw(j))
      read(14,*) (ist_linew(i,j),i=1,numw(j))
      read(14,*) (ikpt_linew(i,j),i=1,numw(j))
      enddo
      read(14,*)
      read(14,*) num_evan
      do ii=1,num_evan
       read(14,*) iit,ikpt_evan(ii),ist_evan(ii),iGX_evan(ii),
     & iband_evan(ii),E_evan(ii),E_evanC(ii)
       enddo
      close(14)

cccccccccccccccccccccccccccccccccc

      pi=4*datan(1.d0)
      allocate(cphase(n1w,n2,n3))
      allocate(phase(n1w,n2,n3))
      allocate(ucw1(n1w,n2,n3,nstm))
      allocate(ucw2(n1w,n2,n3,nstm))

      do k=1,n3
      do j=1,n2
      do i=1,n1w
      phase(i,j,k)=pi*(i-1.d0)/(nkptw-1.d0)/n1w
      cphase(i,j,k)=cdexp(-dcmplx(0.d0,phase(i,j,k)))
      enddo 
      enddo
      enddo
      if(inode.eq.1) write(6,*) "cccccc End of preprocessing ccccc"
cccccccccccccccccccccccccccccccccccccccc
cccc end of preprocessing
cccccccccccccccccccccccccccccccccccccccc


      allocate(uc(n1,n2,n3,10))
cccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccc
      if(inode.eq.1) write(6,*) "cccccc Start reading wr_real.E ccccc"

      open(21,file="wr_real.E",form="unformatted")
      rewind(21)
      read(21) N_trans
      read(21) n1t,n2t,n3t,mnodes

      if(n1.ne.n1t.or.n2.ne.n2t.or.n3.ne.n3t) then
      if(inode.eq.1) write(6,*) "n1,n2,n3 changed in wr_real.E,stop",
     &  n1,n2,n3,n1t,n2t,n3t
      call mpi_abort(MPI_COMM_WORLD,1,ierr)
      endif

      if(ist_init.gt.N_trans) then
      if(inode.eq.1) write(6,*) "ist_init.gt.N_trans,stop",
     &  ist_init,N_trans
      call mpi_abort(MPI_COMM_WORLD,1,ierr)
      endif

      mr=n1*n2*n3
      mr_n=mr/mnodes
      allocate(uc_tmp(mr_n))

      if(inode.eq.1) then
      open(17,file="T.report")
      rewind(17)
      write(17,607)
607   format("  iE  band   ak         E(eV)             T        ", 
     & "        R              T-R           fit_L,fit_R(evan,/norm)",
     & "  fit_L,fit_R(No_ev,/norm)  er_eq    er_L    er_R    er_tot",
     & "  N_runW  N_sys_st_used"  )

      open(42,file="scatte_st_3D.out",form="unformatted")
      rewind(42)
      allocate(uc_tmp2(mr_n))
      numE=ist_final-ist_init+1
      write(42) n1,n2,n3,mnodes
      write(42) numE

      open(66,file="coefficients.out")
      rewind(66)
      write(66,*) 'Tranwellmission and coefficients of system'//
     & 'states for each scattering state'
      
      endif
ccccccccc  roll over to the starting position in the file 21
      do ist=1,ist_init-1
      read(21) ist_test,E_wave,mstate
      do ist1=1,mstate
      do iread=1,mnodes
      read(21)
      enddo
      enddo
      enddo
ccccccccccccccccccccccccccccccccccccccccccccc

      ist=ist_init-1
6000  continue
      ist=ist+1
      if(ist.gt.N_trans.or.ist.gt.ist_final) then
      if(inode.eq.1) then
      close(17)
      close(21)
      close(42)
      close(66)
      deallocate(uc_tmp2)
      endif
      deallocate(uc_tmp)
      deallocate(cphase)
      deallocate(phase)
      deallocate(ucw1)
      deallocate(ucw2)
      deallocate(ucL)
      deallocate(ucR)
      call mpi_abort(MPI_COMM_WORLD,1,ierr)
      endif

cccccccc  This subroutine reads in the system (to be also called "Well") wavefunction. 
cccccccc  These wavefuncitons in real space are the output from the special PEtot calculation

      read(21) ist_test,E_wave,mstate
      if(ist_test.ne.ist) then
      if(inode.eq.1) write(6,*) "ist_test.ne.ist,stop",ist_test,ist
      call mpi_abort(MPI_COMM_WORLD,1,ierr)
      endif
      num_stWell=mstate
      E=E_wave

      if(inode.eq.1) then
      write(6,*) "*********************"
      write(6,*) "*********************"
      write(6,*) "ist,num_stWell,E_wave=",ist,num_stWell,E_wave
      write(6,*) "*********************"
      write(6,*) "*********************"
      endif

      n_tmp1=num_stWell/nnodes
      n_tmp2=num_stWell-nnodes*n_tmp1
      if(inode.le.n_tmp2) then
        num_wr_dis(1)=(inode-1)*(n_tmp1+1)+1
        num_wr_dis(2)=inode*(n_tmp1+1)
      else
        num_wr_dis(1)=(inode-1)*n_tmp1+n_tmp2+1
        num_wr_dis(2)=inode*n_tmp1+n_tmp2
      endif
      nst1=num_wr_dis(2)-num_wr_dis(1)+1

      if(nst1.gt.10) then                       ! by xwjiang
      write(6,*) "nst.gt.10, stop",nst1
      call mpi_abort(MPI_COMM_WORLD,1,ierr)
      endif

      call wave_system(uc,n1,n2,n3,mnodes,mstate,num_wr_dis)


cccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccc  subroutine wave_decomp decomposes each well wavefunctions uc 
cccccccc  by the electrode states (both running waves and evanescence states). 
cccccccc  On return,  ccy2_st
cccccccc  ccy2_st(i_electrode_st,j_Well_st) is the complex coeff. of the
cccccccc  ith electrode state on the jth Well state.
cccccccc  Note, i_electrode_st  ucw and ucw^* are in i*2-1, and i*2
cccccccc  the j_Well_st  well states uc and uc^* are in j, and j+num_stWell
cccccccc
cccccccc also on return, the weight1(j_Well_st) indicates how well is the decomposition (sum of the 
cccccccc amplitudes) for each Well states. 1 means good, less than 1 means bad.  
cccccccc
cccccccc num_stw1 is the total number of electrode states used in the decomposition (fitting). 
cccccccc For uw and uw^*, it counts only as one. However, its information is stored in idble1. Note
cccccccc a evanescence state (not at Gamma point) also has idble1=2. At Gamma point, idble1=1. 
cccccccc
cccccccc num_run1 is the total number of running wave electrode states in thedecomposition. It
cccccccc counts only uw (not uw^*).
cccccccccccccccccccccccccccccccccc
cccccccc
cccccccc  First: the left hand side electrode, the results are  ..ww1 
cccccccc  The relative electrode energy is E+dV/2
      if(inode.eq.1) then
      write(6,*) "***********************************************"
      write(6,*) "**** electrode decomposition for left hand side"
      endif
      Ew0=E+dV/2
      ccy2_st1=dcmplx(0.d0,0.d0)
      cc_R1=dcmplx(0.d0,0.d0)
      weight1=0.d0
      call wave_decomp(ccy2_st1(1,1),Ew0,num_stw1,
     &  num_run1,idble1,iposit1,nnposit1,ucw1,
     &  n1w,n1,n2,n3,cphase_ucw1,dE_dk1,ak_w1,nline_w1,
     &  numw,E_linew,ist_linew,ikpt_linew,nstw,nkptw,num_mx,nintep,imx,
     &  cphase,phase,num_stWell,num_wr_dis,uc,E_evan,E_evanC,ist_evan,
     &  ikpt_evan,iGX_evan,num_evan,weight1,a11,num_iter_evan,dE_evanL,
     &  cc_R1,inode)
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call mpi_allreduce(ccy2_st1,tmp_complex,nwellm*nwellm,
     $     MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
      ccy2_st1=tmp_complex
      call mpi_allreduce(cc_R1,tmp_complex,nwellm*nwellm,
     $     MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
      cc_R1=tmp_complex
      call mpi_allreduce(weight1,tmp_real,nwellm,
     $     MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      weight1=tmp_real

cccccccc  Second: the right hand side electrode, the results are  ..ww2 
cccccccc  The relative electrode energy is E-dV/2
cccccccc  The current flow from the left hand side to right hand side. 
      if(inode.eq.1) then
      write(6,*) "************************************************"
      write(6,*) "**** electrode decomposition for right hand side"
      endif
      Ew0=E-dV/2
      ccy2_st2=dcmplx(0.d0,0.d0)
      cc_R2=dcmplx(0.d0,0.d0)
      weight2=0.d0
      call wave_decomp(ccy2_st2(1,1),Ew0,num_stw2,
     &  num_run2,idble2,iposit2,nnposit2,ucw2,
     &  n1w,n1,n2,n3,cphase_ucw2,dE_dk2,ak_w2,nline_w2,
     &  numw,E_linew,ist_linew,ikpt_linew,nstw,nkptw,num_mx,nintep,imx,
     &  cphase,phase,num_stWell,num_wr_dis,uc,E_evan,E_evanC,ist_evan,
     &  ikpt_evan,iGX_evan,num_evan,weight2,a11,num_iter_evan,dE_evanR,
     &  cc_R2,inode)
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call mpi_allreduce(ccy2_st2,tmp_complex,nwellm*nwellm,
     $     MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
      ccy2_st2=tmp_complex
      call mpi_allreduce(cc_R2,tmp_complex,nwellm*nwellm,
     $     MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
      cc_R2=tmp_complex
      call mpi_allreduce(weight2,tmp_real,nwellm,
     $     MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      weight2=tmp_real
      if(inode.eq.1) then
      write(6,*) "************************************************"
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccc
      
cccccccccccccccccccccccccccccccccc
      call select_electrode()
cccccccccccccccccccccccccccccccccccccccccccccc
cccc select_electrode takes ccy2_st1 as input, select only those running waves and 
cccc evanescence waves which has some composition values for some well states. 
cccc As a result, the eletrode states are re-indexed. 
cccc Now, there are total num_w1 electrode states (counting both uw, and uw^*). 
cccc There are num_wr1 running wave states (counting only uw). 
cccc
cccc  In the index ii=1,num_w1. If idble=2, then ii=i*2-1 is uw, ii=i*2 is uw^*,
cccc  So for ii.le.2*num_wr1, the ii is a running wave. 
cccc
cccc  The index mapping from the selected electrode state ii=1,num_w1, to 
cccc  the original state index are: ii=(1,num_w1)
ccccc iii=i1_st1_ind(ii), for original full index (both uw,uw^*), iii=1,[num_stw1,idble(iii)]. 
cccc                      to be used in ccy2_st1(iii,j_well_st)
cccc  ip=i1_st1_ind0(ii),for index ip=1,num_stw1 (only uw, no uw^*). 
cccc                      to be used in dE_dk1(ip), cphase_ncw1(ip). 
cccc  ib=i1_st1_ind1(ii),for idble(iii), either 1 (uw=uw^*) or 2 (both uw,uw^*)
cccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccc
      num_stWell0=num_stWell     ! the original number of Well states. 

      call select_Wellstate()
ccccccccccccccccccccccccccccccccccccccccccccccc
cccc  select_Wellstate, select the num_stWell system states according to weight1(ii)+weight2(ii)
cccc  It can reduce the num_stWell (this number can be changed). It also sorts the order of
cccc  the states according to the weight in a descending order. 
cccc  The mapping between the new index and the old index is:
ccccc
cccc  jj1=ind_Well(j1), j1=1,2*num_stWell(new), jj1=1,2*num_stWell0(old)
cccc
cccc  Note, in j1, uc,uc^* are i1*2-1,i1*2. But in jj1 (original), uc,uc^* are in jj1,jj1+num_stWell0
cccc
ccccc  It also output the weight_st(ind_Well(j1)). This is the weight to be used for the 
ccccc  Well state.  Basically, if there are enough num_stWell to satisfy all the linear equations (icase=2), 
ccccc  we request:  sum_j1 |cc_x(j1)/weight_st(ind_Well(j1))|^2  to be minimum
ccccc
ccccccccccccccccccccccccccccccccccccccccccccccc
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call mpi_bcast(ind_Well,nwellm,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(inode.eq.1) then
      write(6,*) "n_stWell;n_wr1,n_we1(L);n_wr2,n_we2(R)",
     &     num_stWell,";", num_wr1,num_w1-num_wr1*2, ";",
     &     num_wr2,num_w2-num_wr2*2
      endif

      if(2*num_stWell.lt.num_wr2+num_wr1) then
      if(inode.eq.1) write(6,*) "have no solution"
      goto 401
      endif
cccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccc
ccccc  num_st-1 will be the number of running wave and evanescence states to be eliminated 
ccccc  from the linear combination of the Well states, the other one equation is for the 
ccccc  incoming wave.   
ccccc   num_st is the number of linear equation need to be satisfied 
cccccccccccccccccccccccccccccccccccc

      call set_icase()

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if(inode.eq.1) write(6,*) "has a solution, icase=",icase

      n_tmp=2*num_stWell
      
      allocate(cc_y(num_st+2*n_tmp))
      allocate(cc_y_tmp(num_st+2*n_tmp))
      allocate(cc_matrix(num_st,2*num_stWell))
      allocate(cc_matrix0(num_st+2*n_tmp,n_tmp))
      allocate(cc_matrix0_tmp(num_st+2*n_tmp,n_tmp))

      lwork=20*(num_st+2*n_tmp)
      allocate(cwork(lwork))
      allocate(rwork(lwork))
      
      allocate(ss(num_st+2*n_tmp))
      allocate(ipiv2(num_st))
      allocate(sum_w1(n1))
      allocate(sum_w2(n1))
      allocate(sum_w4(n1,40))
      allocate(cc_w1(n1))
      allocate(cc_w2(n1))

ccccc Those arrays are deallocated!! They are allocated here again. 

      allocate(cc_s1(n_tmp,n_tmp))
      allocate(cc_s2(n_tmp,n_tmp))
      allocate(U_s1(n_tmp,n_tmp))
      allocate(U_s2(n_tmp,n_tmp))
      allocate(s_scale(n_tmp))
      allocate(y_tmp(n_tmp))
      allocate(x_tmp(n_tmp))
      allocate(cwork_tmp(2*n_tmp))
      allocate(rwork_tmp(2*n_tmp))

      if(inode.eq.1) then
      write(42) ist,num_wr2
      endif


      do 400 ii1=1,num_wr2    ! num_wr2 is the number of running waves from the right hand side
      if(inode.eq.1) then
      write(6,*) "**********************************"
      write(6,*) "***** solving for incoming running wave:",ii1 
      write(6,*) "**********************************"
cccccccc solve for the scattering state for this (ii1) incoming Bloch state from the right hand side
cccccccccccccccccc
ccccc jmx is the num_stWell state which has maximum ccy2_st2 for this ii1 scattering state
ccccc for this state, the coeff. cc=1.d0. 
ccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccc
**********************************************************************
**********************************************************************
****   setting up the linear equation cc_matrix(i,j)*cc_x(j)=cc_y(i) matrix 
****   elements cc_matrix. i=1,num_st, j=1,2*num_stWell
**********************************************************************
**********************************************************************
ccccc the index for ccy2_st1 is ind_Well(j1), index for cc_matrix is j1. 
******  index j1 is straight forward. 
******  index i2 runs through num_wr1,evenascence state 1, num_wr2, evanescence state 2. 
      do j1=1,2*num_stWell    

      i2=0
      do i1=1,num_wr1        ! left hand side, running state
      i2=i2+1
      if(dE_dk1(i1_st1_ind0(i1*2)).gt.0.d0) then     ! i1*2-1 forward scattering: dE_dk1.gt.0, 
      cc_matrix(i2,j1)=ccy2_st1(i1_st1_ind(i1*2),ind_Well(j1))       ! backward scattering equals zero
      else
      cc_matrix(i2,j1)=ccy2_st1(i1_st1_ind(i1*2-1),ind_Well(j1))       ! backward scattering equals zero
      endif
      enddo

      do i1=num_wr1*2+1,num_w1   ! left hand side, evanescence state
      if(ifl_evan1(i1).eq.1) then
      i2=i2+1 
      cc_matrix(i2,j1)=ccy2_st1(i1_st1_ind(i1),ind_Well(j1))
      endif
      enddo

      do i1=1,num_wr2        ! right hand side, running state
      i2=i2+1
      if(dE_dk2(i1_st2_ind0(i1*2)).gt.0.d0) then
      cc_matrix(i2,j1)=ccy2_st2(i1_st2_ind(i1*2-1),ind_Well(j1))  ! incoming scattering state equals zero or1 for i1.eq.ii1 
      else
      cc_matrix(i2,j1)=ccy2_st2(i1_st2_ind(i1*2),ind_Well(j1))  ! incoming scattering state equals zero or 1for i1.eq.ii1
      endif
      enddo

      do i1=num_wr2*2+1,num_w2   ! right hand side, evanescence state
      if(ifl_evan2(i1).eq.1) then
      i2=i2+1
      cc_matrix(i2,j1)=ccy2_st2(i1_st2_ind(i1),ind_Well(j1))
      endif
      enddo

       if(i2.ne.num_st) then
       write(6,*) "i2.ne.num_st, stop"
       stop
       endif

      enddo      !  do j1
**********************************************************************
**********************************************************************
cccc set up the cc_s(2*num_stWell,2*num_stWell) matrix
      do j1=1,2*num_stWell
      do j2=1,2*num_stWell
      cc_s1(j1,j2)=cc_R1(ind_Well(j1),ind_Well(j2))
      cc_s2(j1,j2)=cc_R2(ind_Well(j1),ind_Well(j2))
      enddo
      enddo
cccccccccccccccccccccccc
cccc Do a Cholesky factorization to cc_s1=U_s1^H*U_s1,cc_s2=U_s2^H*U_s2
cccc so: cc^* cc_s1 cc=|U_s1*cc_s1|^2
cccccccccccccccccccccccccc
      n_tmp=2*num_stWell
      x_tmp=dcmplx(1.d0,0.d0)
      s_scale=1.d0

      U_s1=cc_s1
      do i=2,n_tmp
      U_s1(i,1:i-1)=dcmplx(0.d0,0.d0)
      enddo
      call zpotrf('U',n_tmp,U_s1,n_tmp,info)
cccccccccccccccccccccccccccccccccccccccc

      x_tmp=dcmplx(1.d0,0.d0)

      U_s2=cc_s2
      do i=2,n_tmp
      U_s2(i,1:i-1)=dcmplx(0.d0,0.d0)
      enddo
      call zpotrf('U',n_tmp,U_s2,n_tmp,info)
**********************************************************************
**********************************************************************
***   Now, set up cc_y(i), i=1,num_st
**********************************************************************
**********************************************************************
      cc_y=dcmplx(0.d0,0.d0)

      i2=0
      do i1=1,num_wr1      ! left hand side
      i2=i2+1
      enddo

      do i1=num_wr1*2+1,num_w1
      if(ifl_evan1(i1).eq.1) then
      i2=i2+1
      endif
      enddo

      do i1=1,num_wr2     ! right hand side
      i2=i2+1
ccc this is the resulting coeff.  of the linear combination for the ii1 electrode state
ccc Note that,  if(dE_dk2(i1_st2_ind0(ii1*2)).gt.0.d0), then
ccc  this incoming state is i1_st2_ind(i1*2-1), otherwise, it is i1_st2_ind(i1*2). 
ccc  but in both case, the corresponding i2 index in cc_y(i2) is i2
      if(i1.eq.ii1) then       ! ii1 is the index of the incoming wave
      cc_y(i2)=dcmplx(1.d0,0.d0)    
      endif
      enddo

      do i1=num_wr2*2+1,num_w2   ! right hand side, evanescence state
      if(ifl_evan2(i1).eq.1) then
      i2=i2+1
      endif
      enddo

      if(i2.ne.num_st) then
      write(6,*) "i2.ne.num_st, stop"
      stop
      endif
ccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccc
cccc   Now, solve the equation:  cc_matrix*cc_x=cc_y
ccccc  \sum_(j=1,2*num_stWell) cc_matrix(i,j)*cc_x(j) = cc_y(i),   here i=1,num_st
cccccc  when 2*num_stWell > num_st, solve this equation by minimize ||cc_x(j)||
cccccccccccccccc   put a weight in cc_matrix(i,j)*weight_st(j)*cc_x(j)/weight_st(j)
cccccccc
       do j1=1,n_tmp       ! system states (number of coeff)
       do i1=1,num_st      ! electrode states
       cc_matrix0(i1,j1)=cc_matrix(i1,j1)
       enddo
       enddo

cccccccccccccccccc   put a factor here. We need the above num_st equation to be satisfied exactly, 
cccccccccccccccccc   Then, we will worry about the following states
       fact_left=1.D-7      ! tested, it doesn't make much differences, whether it is 1.D-5, or 1.D-8
       fact_right=1.D-7
cccccccc might want to change them to 1.D-5

       do j1=1,n_tmp
       do i1=1,n_tmp
       cc_matrix0(num_st+i1,j1)=U_s1(i1,j1)*fact_left
       cc_matrix0(num_st+n_tmp+i1,j1)=U_s2(i1,j1)*
     &    fact_right
       enddo
       enddo

       cc_matrix0_tmp=cc_matrix0
       cc_y_tmp=cc_y      ! cc_y is zero for j1>num_st
cccccccccccccccccccccccccccccc
      call zgelss(num_st+2*n_tmp,n_tmp,1,cc_matrix0,
     &  num_st+2*n_tmp,cc_y,
     &  num_st+2*n_tmp,ss,-1.d0,irank,cwork,lwork,rwork,info)

       fit_err1=0.d0
       fit_err2=0.d0
       fit_err3=0.d0
       do i1=1,num_st+2*n_tmp
       cc=dcmplx(0.d0,0.d0)
       do j1=1,n_tmp
       cc=cc+cc_matrix0_tmp(i1,j1)*cc_y(j1)
       enddo
       if(i1.le.num_st) then
       fit_err1=fit_err1+abs(cc_y_tmp(i1)-cc)**2
       endif
       if(i1.gt.num_st.and.i1.le.num_st+n_tmp) then
       fit_err2=fit_err2+abs(cc)**2/fact_left**2
       endif
       if(i1.gt.num_st+n_tmp) then
       fit_err3=fit_err3+abs(cc)**2/fact_right**2
       endif
       enddo
cccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccc  test, test
       fit_err22=0.d0
       fit_err33=0.d0
       do j1=1,n_tmp
       do j2=1,n_tmp
       fit_err22=fit_err22+dconjg(cc_y(j1))*cc_s1(j1,j2)*cc_y(j2)
       fit_err33=fit_err33+dconjg(cc_y(j1))*cc_s2(j1,j2)*cc_y(j2)
       enddo
       enddo
       endif

      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call mpi_bcast(cc_y,num_st+2*n_tmp,MPI_DOUBLE_COMPLEX,0,
     &  MPI_COMM_WORLD,ierr)
c        write(6,799) fit_err1,fit_err2,fit_err3 
c799    format("scatt state fit err(eq,left,right)",3(E11.3,1x))
c        write(6,798) fit_err1,fit_err22,fit_err33 
c798    format("scatt state fit err22(eq,left,right)",3(E11.3,1x))
ccccc fit_err is the overall error of the fit. It should be used to decide how good is this
ccccc scattering state
cccccccccccccc   construct the scattering wavefunction state uc_test in the whole space
cccccccccccccc   from the Well states
      uc_test=dcmplx(0.d0,0.d0)

      do j1=1,2*num_stWell
      iist=ind_Well(j1)
      if(iist.le.num_stWell0) then
      if(iist.ge.num_wr_dis(1).and.iist.le.num_wr_dis(2)) then
      iiist=iist-num_wr_dis(1)+1
      do k=1,n3
      do j=1,n2
      do i=1,n1
      uc_test(i,j,k)=uc_test(i,j,k)+uc(i,j,k,iiist)*cc_y(j1)
      enddo
      enddo
      enddo
      endif
      else
      iist=iist-num_stWell0
      if(iist.ge.num_wr_dis(1).and.iist.le.num_wr_dis(2)) then
      iiist=iist-num_wr_dis(1)+1
      do k=1,n3
      do j=1,n2
      do i=1,n1
      uc_test(i,j,k)=uc_test(i,j,k)+dconjg(uc(i,j,k,iiist)*cc_y(j1)     ! uc_test is the scattering state
      enddo
      enddo
      enddo
      endif
      endif     
      enddo    ! do j1

      do iread=1,mnodes
       do ii=1,mr_n
       jj=ii+(iread-1)*mr_n
       i=(jj-1)/(n2*n3)+1
       j=(jj-1-(i-1)*n2*n3)/n3+1
       k=jj-(i-1)*n2*n3-(j-1)*n3
       uc_tmp(ii)=uc_test(i,j,k)
       enddo
       call mpi_barrier(MPI_COMM_WORLD,ierr)
       call mpi_reduce(uc_tmp,uc_tmp2,mr_n,MPI_DOUBLE_COMPLEX,
     &  MPI_SUM,0,MPI_COMM_WORLD,ierr)
       if(inode.eq.1) then
         do ii=1,mr_n
         jj=ii+(iread-1)*mr_n
         i=(jj-1)/(n2*n3)+1
         j=(jj-1-(i-1)*n2*n3)/n3+1
         k=jj-(i-1)*n2*n3-(j-1)*n3
         uc_test(i,j,k)=uc_tmp2(ii)
         enddo
       endif
      enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call output_result()
400   continue
cccccccccccccccccccccccccccccccccccc
      deallocate(cc_matrix)
      deallocate(cc_matrix0)
      deallocate(cc_y)
      deallocate(cc_y_tmp)
      deallocate(cc_matrix0_tmp)
      deallocate(sum_w4)
      deallocate(cc_s1)
      deallocate(cc_s2)
      deallocate(u_s1)
      deallocate(u_s2)
      deallocate(s_scale)
      deallocate(y_tmp)
      deallocate(x_tmp)
      deallocate(cwork_tmp)
      deallocate(rwork_tmp)
      deallocate(cwork)
      deallocate(rwork)
      deallocate(ss)
      deallocate(ipiv2)
      deallocate(sum_w1)
      deallocate(sum_w2)
      deallocate(cc_w1)
      deallocate(cc_w2)
      deallocate(uc_test)
401   continue 
      goto 6000
      call mpi_finalize(ierr)

      contains


******************************************************************
******************************************************************
******************************************************************

      subroutine select_electrode()
      implicit double precision (a-h,o-z)

ccccc This subroutine selects the electrode (wire) state (running ane evanescence) to be used
ccccc in the fitting. (For example, many evanescence state will not be used). 
ccccc the selected wire states are: num_w1(include running and evanescence, both directions)
ccccc                               num_wr1 (running wave only, only one direction)
ccccc The mapping index: original iii(1,num_stw1), new i1(1,num_w1):  iii=i1_st1_ind(i1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      num_w1=0
      num_wr1=0

      iii=0
      do i1=1,num_stw1
      sum=0.d0
      do ip=1,idble1(i1)
      iii=iii+1
      do j1=1,num_stWell
      sum=sum+abs(ccy2_st1(iii,j1))**2
      enddo
      enddo
      sum=dsqrt(sum)

c      if(sum.gt.1.D-3) then      ! take the electrode states in the system_wave
      if(sum.gt.-1) then      ! take all the electrode state
        num_w1=num_w1+idble1(i1)
        i1_st1_ind(num_w1)=iii      ! original index with idble, both psi_wire, psi_wire^*
        i1_st1_ind0(num_w1)=i1      ! original index without idble, count psi_wire,psi_wire^* as one
        i1_st1_ind1(num_w1)=1       ! 1, or 2 for psi_wire, or psi_wire^*
        if(idble1(i1).eq.2) then
        i1_st1_ind(num_w1-1)=iii-1
        i1_st1_ind0(num_w1-1)=i1
        i1_st1_ind1(num_w1-1)=1
        i1_st1_ind1(num_w1)=2
        endif
        if(i1.le.num_run1) num_wr1=num_wr1+1
      endif
      enddo
cccccccccccccccccccccccccccccccccccccccccccccc
      num_w2=0
      num_wr2=0


      iii=0
      do i1=1,num_stw2

      sum=0.d0
      do ip=1,idble2(i1)
      iii=iii+1
      do j1=1,num_stWell
      sum=sum+abs(ccy2_st2(iii,j1))**2
      enddo
      enddo
      sum=dsqrt(sum)
cc      if(sum.gt.1.D-3) then
      if(sum.gt.-1) then     ! take all the electrode states


        num_w2=num_w2+idble2(i1)
        i1_st2_ind(num_w2)=iii
        i1_st2_ind0(num_w2)=i1
        i1_st2_ind1(num_w2)=1
        if(idble2(i1).eq.2) then
        i1_st2_ind(num_w2-1)=iii-1
        i1_st2_ind0(num_w2-1)=i1
        i1_st2_ind1(num_w2-1)=1
        i1_st2_ind1(num_w2)=2
        endif
        if(i1.le.num_run2) num_wr2=num_wr2+1

      endif
      enddo
ccccccccccccccccccccccccccccccccccccccccc
      return

      end subroutine select_electrode

******************************************************************
******************************************************************
******************************************************************

      subroutine select_Wellstate()
ccccccccccccccccccccccccccccccccccccccccccccc
ccccc This program also gives a sorted order for the num_stWell system states, 
ccccc according to their overall fitting magnitudes (weights, with evanescence states)
ccccc It can also reduce the number of num_stWell, for the usable system states !
ccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do i1=1,2*num_stWell
      weight_st(i1)=weight1(i1)+weight2(i1)
      enddo

      num_tmp=0
      num_tmp2=0

      do i1=1,2*num_stWell
      weight_max=-1.d0
      do i2=1,2*num_stWell
      if(weight_st(i2).gt.weight_max) then
      weight_max=weight_st(i2)
      i2m=i2
      endif
      enddo

      num_tmp=num_tmp+1
      ind_Well(num_tmp)=i2m
      if(weight_st(i2m).gt.0.2) num_tmp2=num_tmp2+1    !  a criterior to select num_stWell
      weight_st(i2m)=-2.d0
      enddo
      
      if(num_tmp.ne.2*num_stWell) then
      if(inode.eq.1) write(6,*) "num_tmp.ne.num_stWell, something wrong"
      stop
      endif


      num_stWell0=num_stWell


      if(num_tmp2.lt.2*num_stWell) then
c      write(6,*) "2*num_stWell reduced ", 2*num_stWell, num_tmp2 
      num_stWell=num_tmp2/2
      endif
cccc the index is stalled in ind_Well()
c	write(6,*) ' after reduce 2*num_stWell',2*num_stWell
ccccccccccccccccccccccccccc

      do i1=1,2*num_stWell0
      weight_st(i1)=1.d0/(1.d0+abs(weight1(i1)-1.d0)/0.001+
     &    abs(weight2(i1)-1.d0)/0.001)     ! 0.001 can be re-adjusted
      enddo
ccccc weight_st denotes the importance of this well state, it will determine how much
ccccc this state should participate in the linear combination process. This will allow
ccccc all the 2*num_stWell states to participate in the linear combination. 

      return
      end subroutine select_Wellstate
      
******************************************************************
******************************************************************
******************************************************************

      subroutine set_icase()
      implicit double precision (a-h,o-z)
  
      icase=0
******************************************************************
*****   case one
******************************************************************
cccc ignore the evanescence states, just eliminate the running waves
      num_st=num_wr2+num_wr1       
      ifl_evan1=0
      ifl_evan2=0
      icase=1
      return     ! only request the running wave eq. to be satisfied. 
ccccc Lin-Wang Wang, Oct.4, 2010
cccc This is modified, so for all the cases, it is icase=1, 
cccc which means the evanescence states coeff. are not forced to be zero in 
cccc the scattering states. Forcing them to be zero is two stringent, which 
cccc will cause problem (because in real case, they might not be zero).
cccc But note, the evanescence states (or its approximations) are used in 
cccc the well state decompositions
******************************************************************
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      return
      end subroutine set_icase

******************************************************************
******************************************************************
******************************************************************
      subroutine output_result()
      implicit double precision (a-h,o-z)
      real*8 T_store(100),abcc_store(100)
cccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccc
cccc  construct the fitted states uc_test2 from the electrode states
cccc  store the lateral averaged wavefunction square in sum_w1,sum_w2, to be output for plotting
ccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccc

      allocate(uc_test2(n1w,n2,n3))
      allocate(uc_test3(n1w,n2,n3))

      aL_forward=0.d0
      aL_backward=0.d0

cccc  actually, we don't need to do all lltot to get the current, this is just for plotting

      sum_ave=0.d0
      do 111 ll=0,0      
      uc_test2=dcmplx(0.d0,0.d0)
      uc_test3=dcmplx(0.d0,0.d0)

      do 101 i1=1,num_w1
      ip1=i1_st1_ind(i1)

      cc1=dcmplx(0.d0,0.d0)
      do j1=1,2*num_stWell
      cc1=cc1+ccy2_st1(ip1,ind_Well(j1))*cc_y(j1)
      enddo

      if(ll.eq.0) then     ! this is same for all ll
      if(i1.le.num_wr1*2) then
        if(mod(i1,2).eq.1) T1=abs(cc1)**2*dE_dk1(i1_st1_ind0(i1))
        if(mod(i1,2).eq.0) T1=-abs(cc1)**2*dE_dk1(i1_st1_ind0(i1))
        T_store(i1)=T1
        abcc_store(i1)=abs(cc1)
      endif
      if(i1.gt.num_wr1*2) then
      T1=0.d0
        T_store(i1)=0.d0
        abcc_store(i1)=abs(cc1)
      endif

      if(T1.gt.0.d0) then
      aL_forward=aL_forward+T1
      else
      aL_backward=aL_backward+T1
      endif

      endif
ccccccccccccccccccccccccccccccccccccccccccccccc

      ip=i1_st1_ind0(i1)
      cc=cphase_ucw1(ip)**0

      if(i1.le.num_wr1*2) then
      fact=1.d0
      else
      fact=0.d0   ! for evanescent states
      endif


      if(i1_st1_ind1(i1).eq.1) then
      do k=1,n3
      do j=1,n2
      do i=1,n1w
      uc_test2(i,j,k)=uc_test2(i,j,k)+
     &   cc1*ucw1(i,j,k,ip)*cc
      uc_test3(i,j,k)=uc_test3(i,j,k)+
     &   cc1*ucw1(i,j,k,ip)*cc*fact
      enddo
      enddo
      enddo 
      else
      do k=1,n3
      do j=1,n2
      do i=1,n1w
      uc_test2(i,j,k)=uc_test2(i,j,k)
     &  +cc1*dconjg(ucw1(i,j,k,ip)*cc)
      uc_test3(i,j,k)=uc_test3(i,j,k)
     &  +cc1*dconjg(ucw1(i,j,k,ip)*cc)*fact
      enddo
      enddo
      enddo 
      endif

101   continue
cccccccccccccccccccccccccccccccccccccc

cccc uc_test is the scattering state constructed directly from the system states
cccc uc_test2 is the unit length scattering state constructred from electrode states (running+evan) using coeff cc1
cccc uc_test2 is the unit length scattering state constructred from electrode states (running only) using coeff cc1
      sum1=0.d0
      sum2=0.d0
      sum11=0.d0
      diff=0.d0
      diff2=0.d0
      sum00=0.d0
      do i=1,n1w
      sum=0.d0
      do k=1,n3
      do j=1,n2
      sum1=sum1+abs(uc_test2(i,j,k))**2
      sum11=sum11+abs(uc_test3(i,j,k))**2
      sum2=sum2+abs(uc_test(i+nnposit1,j,k))**2
      diff=diff+abs(uc_test(i+nnposit1,j,k)-uc_test2(i,j,k))**2
      diff2=diff2+abs(uc_test(i+nnposit1,j,k)-uc_test3(i,j,k))**2
      sum=sum+abs(uc_test2(i,j,k))**2
      sum00=sum00+abs(ucw1(i,j,k,1))**2
      enddo
      enddo
      sum_w1(i+nnposit1)=sum
      cc_w1(i+nnposit1)=uc_test2(i,2,4)
      enddo

       fit_left=diff/sum2
       fit_left2=diff2/sum2
cccccccccccccccccccccccccccccccc

       sum2_nnposit=sum2

      if(ll.eq.0) then
       sum2_end=sum2
      endif
      sum_ave=sum_ave+sum2
111   continue
      sum_ave=sum_ave/14.d0
       amplitude_left=sum2_nnposit/sum_ave
       ratio_left=sum2_nnposit/sum2_end

       write(66,*) "num_w1=",num_w1
       write(66,887) (T_store(i1),i1=1,num_w1)
       write(66,886) (abcc_store(i1),i1=1,num_w1)
887    format(" Left eltrd channels, T=:", 20(E8.2,1x))
886    format(" Left eltrd channels, c=:", 20(E8.2,1x))

cccccccccccccccccccccccccccccccccccccccccccccccccc
      aR_forward=0.d0
      aR_backward=0.d0

ccccc  again, we don't need to do all lltot to get the current
      do 112  ll=0,0   
      uc_test2=dcmplx(0.d0,0.d0)
      uc_test3=dcmplx(0.d0,0.d0)

      do 102 i1=1,num_w2
      ip1=i1_st2_ind(i1)

      cc1=dcmplx(0.d0,0.d0)
      do j1=1,2*num_stWell
      cc1=cc1+ccy2_st2(ip1,ind_Well(j1))*cc_y(j1)
      enddo

      if(ll.eq.0) then      ! this is same for all lltot
      if(i1.le.num_wr2*2) then
       if(mod(i1,2).eq.1) T1=abs(cc1)**2*dE_dk2(i1_st2_ind0(i1))
       if(mod(i1,2).eq.0) T1=-abs(cc1)**2*dE_dk2(i1_st2_ind0(i1))
       T_store(i1)=T1
       abcc_store(i1)=abs(cc1)
      endif
      if(i1.gt.num_wr2*2) then
      T1=0.d0
       T_store(i1)=0.d0
       abcc_store(i1)=abs(cc1)
      endif

      if(T1.gt.0.d0) then
      aR_forward=aR_forward+T1
      else
      aR_backward=aR_backward+T1
      endif

      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ip=i1_st2_ind0(i1) 

      cc=cphase_ucw2(ip)**0

      if(i1.le.num_wr2*2) then
      fact=1.d0
      else
      fact=0.d0
      endif


      if(i1_st2_ind1(i1).eq.1) then
      do k=1,n3
      do j=1,n2
      do i=1,n1w
      uc_test2(i,j,k)=uc_test2(i,j,k)+
     &   cc1*ucw2(i,j,k,ip)*cc
      uc_test3(i,j,k)=uc_test3(i,j,k)+
     &   cc1*ucw2(i,j,k,ip)*cc*fact
      enddo
      enddo
      enddo 
      else
      do k=1,n3
      do j=1,n2
      do i=1,n1w
      uc_test2(i,j,k)=uc_test2(i,j,k)
     &  +cc1*dconjg(ucw2(i,j,k,ip)*cc)
      uc_test3(i,j,k)=uc_test3(i,j,k)
     &  +cc1*dconjg(ucw2(i,j,k,ip)*cc)*fact
      enddo
      enddo
      enddo 
      endif

102   continue

      sum1=0.d0
      sum11=0.d0
      sum2=0.d0
      diff=0.d0
      diff2=0.d0
      sum00=0.d0
      do i=1,n1w
      sum=0.d0
      do k=1,n3
      do j=1,n2
      sum1=sum1+abs(uc_test2(i,j,k))**2
      sum11=sum11+abs(uc_test3(i,j,k))**2
      sum2=sum2+abs(uc_test(i+nnposit2,j,k))**2
      diff=diff+abs(uc_test(i+nnposit2,j,k)-uc_test2(i,j,k))**2
      diff2=diff2+abs(uc_test(i+nnposit2,j,k)-uc_test3(i,j,k))**2
      sum=sum+abs(uc_test2(i,j,k))**2
      sum00=sum00+abs(ucw2(i,j,k,1))**2
      enddo
      enddo
      sum_w2(i+nnposit2)=sum
      cc_w2(i+nnposit2)=uc_test2(i,2,4)
      enddo

       fit_right=diff/sum2
       fit_right2=diff2/sum2


       sum2_nnposit=sum2
      if(ll.eq.lltot-1) sum2_end=sum2
112   continue
       amplitude_right=sum2_nnposit/sum_ave
       ratio_right=sum2_nnposit/sum2_end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       do iread=1,mnodes

       do ii=1,mr_n

       jj=ii+(iread-1)*mr_n

       i=(jj-1)/(n2*n3)+1
       j=(jj-1-(i-1)*n2*n3)/n3+1
       k=jj-(i-1)*n2*n3-(j-1)*n3

       uc_tmp2(ii)=uc_test(i,j,k)

       enddo
       write(42) uc_tmp2
       enddo

cccccccccccccccccccccccccccccccccccccccccccc



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       write(66,*) "num_w2=",num_w2
       write(66,885) (T_store(i1),i1=1,num_w2)
       write(66,884) (abcc_store(i1),i1=1,num_w2)
885    format("Right eltrd channels, T=:", 20(E8.2,1x))
884    format("Right eltrd channels, c=:", 20(E8.2,1x))

cccccccccccccccccccccccccccccccc
      deallocate(uc_test2)
      deallocate(uc_test3)

      aL_forward=aL_forward/aR_forward
      aL_bacward=aL_backward/aR_forward
      aR_backward=aR_backward/aR_forward


      write(6,406) ist,nline_w2(i1_st2_ind0(ii1*2)),
     &  (ak_w2(i1_st2_ind0(ii1*2))-1)/(nkptw-1),
     &  Ew0,aL_forward,aR_backward,
     &   aL_forward-aR_backward,
     &  fit_left,fit_right,fit_left2,fit_right2,
     &  fit_err1,fit_err2,fit_err3,
     &  fit_err1+fit_err2*fact_left**2+fit_err3*fact_right**2,
     &  num_st,2*num_stWell

ccccc fit_err1: the running wave boundary condition error, should be very small
ccccc fit_err2: the absolut left residual error  (exclude the evanescent states) (different from fit_left,not divided by norm_left)
ccccc fit_err1: the absolut right residual error (exclude the evanescent states) (different from fit_right,not divided by norm_right)
ccccc fit_err1+fit_err2*fact_left**2+fit_err3*fact_right**2, the fitted weighted total error
ccccc num_st: the number of equations (boundary condition requirements) used. In this program, it is the number of running waves
ccccc 2*num_stWell: the number of system states used. 



      write(17,406) ist,nline_w2(i1_st2_ind0(ii1*2)),
     &  (ak_w2(i1_st2_ind0(ii1*2))-1)/(nkptw-1),
     &  Ew0,aL_forward,aR_backward,
     &   aL_forward-aR_backward,
     &  fit_left,fit_right,fit_left2,fit_right2,
     &  fit_err1,fit_err2,fit_err3,
     &  fit_err1+fit_err2*fact_left**2+fit_err3*fact_right**2,
     &  num_st,2*num_stWell


405   format(2(i4,1x),f10.6,2x,f12.6,2x,3(E17.10,1x),3x,
     &  2(E13.6,1x),3x,4(E8.2,1x))
406   format(2(i4,1x),f10.6,1x,f12.6,2x,3(E17.10,1x),2x,
     &  2(E9.2,1x),2x,2(E9.2,1x),3x,4(E8.2,1x),1x,2(i3,1x))


       return
       end subroutine output_result
******************************************************************
******************************************************************
******************************************************************


        end






