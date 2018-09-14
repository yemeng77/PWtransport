      subroutine wave_decomp(ccy2_st,Ew,num_st,
     & num_run,idble,iposit,nnposit,ucw,
     & n1w,n1,n2,n3,cphase_ucw,dE_dk,ak_w,nline_w,
     & numw,E_linew,ist_linew,ikpt_linew,nstw,nkptw,
     & num_mx,nintep,imx,cphase,phase,num_stWell,num_wr_dis,uc,
     & E_evan,E_evanC,ist_evan,ikpt_evan,iGX_evan,num_evan,weight,a11,
     & num_iter_evan,dE_evan,cc_R)

ccccccccc  this subroutine generates the electrod state
ccccccccc  coefficients "ccy2_st" for "num_stWell" given input system wavefunction "uc". 
cccccccc  ccy2_st(i_electrode_st,j_Well_st) is the complex coeff. of the 
cccccccc  ith electrode state on the jth Well state. 
cccccccc  Note, i_electrode_st  ucw and ucw^* are in i*2-1, and i*2
cccccccc  the j_Well_st  well states uc and uc^* are in j, and j+num_stWell
cccccccc

      implicit double precision (a-h,o-z)
      include 'mpif.h'

      integer status(MPI_STATUS_SIZE)

      parameter (nm=1200)
      parameter (nevan=4000)
      parameter (nstm=200)
      parameter (nwellm=400)

      complex*16 uc(n1,n2,n3,10)
      real*8 E_linew(nm,nm)
      integer numw(nm),ist_linew(nm,nm),ikpt_linew(nm,nm)
      real*8 E_evan(nevan),E_evanC(nevan)
      integer ist_evan(nevan),ikpt_evan(nevan),iGX_evan(nevan),i
     &     used_evan(nevan)

      real*8 phase(n1w,n2,n3)
      complex*16 cphase(n1w,n2,n3)
      complex*16 ucw(n1w,n2,n3,nstm)
      real*8 dE_dk(nstm),ak_w(nstm),aI_tmp(nstm)
      integer nline_w(nstm)
      integer num_mx(nstm),imx(10,nstm)
      complex*16 cphase_ucw(nstm)

      real*8 kk1,kk2,kk3,k1k2,k2k3,k2k1,k3k1,k3k2,k1k3

      complex*16, allocatable, dimension (:,:,:) :: uc_test,
     &   uc_test2
      complex*16, allocatable, dimension (:,:,:,:) :: uc_R
      complex*16, allocatable, dimension (:) :: uc_tmp,uc_tmp2

      real*8 weight(nwellm)
      complex*16 ccy2_st(nwellm,nwellm),cc_R(nwellm,nwellm)
      complex*16 ccm(nwellm,nwellm),ccA(nwellm,nwellm)
      complex*16 ccy(nwellm,10),ccy2(nwellm),ccy3(nwellm)
      integer  ipiv(nstm),iposit(nstm),idble(nstm),ievan(nstm)

      complex*16 cc,cc1,cc2,cc3,cc_st

      character*7 fileh

      integer num_wr_dis(2)
      complex*16 tmp_complex(nwellm,nwellm)
      real*8 tmp_real(nwellm)
      integer inode,nnodes

      common /mpi_data/inode,nnodes

      fileh="wr.new."

      if(inode.eq.1) then
      open(77,file=trim(fileh)//"001",form="unformatted")
      rewind(77)
      read(77) n1t,n2t,n3t,mnodes
      close(77)
      if(n1w.ne.n1t.or.n2.ne.n2t.or.n3.ne.n3t) then
      write(6,*) "n1w,n2,n3.ne.n1t,n2t,n3t,stop",
     &   n1w,n2,n3,n1t,n2t,n3t
      call mpi_abort(MPI_COMM_WORLD,1,ierr)
      endif
      endif

      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call mpi_bcast(mnodes,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      pi=4*datan(1.d0)

      call find_state(Ew,num_tot0,num_st0,num_run,
     & idble,iposit,ievan,ucw,n1w,n2,n3,cphase_ucw,dE_dk,ak_w,nline_w,
     & numw,E_linew,ist_linew,ikpt_linew,nstw,nkptw,nintep,
     & num_mx,imx,cphase,phase,fileh,mnodes,
     & E_evan,E_evanC,ist_evan,ikpt_evan,iGX_evan,num_evan,dE_evan)

      ccy2_st1=dcmplx(0.d0,0.d0)
      cc_R1=dcmplx(0.d0,0.d0)
      weight1=0.d0

      ccm=dcmplx(0.d0,0.d0)
      ccA=dcmplx(0.d0,0.d0)
      ccy=dcmplx(0.d0,0.d0)
      ccy2=dcmplx(0.d0,0.d0)
      ccy3=dcmplx(0.d0,0.d0)
      ipiv=0
      
      allocate(uc_test(n1w,n2,n3))
      allocate(uc_test2(n1w,n2,n3))
      allocate(uc_R(n1w,n2,n3,num_stWell))
      uc_R=dcmplx(0.d0,0.d0)

      num_iter=0
      num_st_old=1
ccccccccccccccccccccccccccccccccccc
      do 2000 num_st=num_run,num_st0
      num_iter=num_iter+1
      do ii1=num_st_old,num_st
      do ii2=1,ii1
      cc1=dcmplx(0.d0,0.d0)
      cc2=dcmplx(0.d0,0.d0)
      do k=1,n3
      do j=1,n2
      do i=1,n1w
      cc1=cc1+ucw(i,j,k,ii1)*dconjg(ucw(i,j,k,ii2))
      cc2=cc2+dconjg(ucw(i,j,k,ii1))*dconjg(ucw(i,j,k,ii2))
      enddo
      enddo
      enddo

      ccm(iposit(ii2),iposit(ii1))=cc1             ! increase the overlap matrix
      ccm(iposit(ii1),iposit(ii2))=dconjg(cc1)
      if(idble(ii1).eq.2) then
      ccm(iposit(ii2),iposit(ii1)+1)=cc2
      ccm(iposit(ii1)+1,iposit(ii2))=dconjg(cc2)
      endif
      if(idble(ii2).eq.2) then
      ccm(iposit(ii1),iposit(ii2)+1)=cc2
      ccm(iposit(ii2)+1,iposit(ii1))=dconjg(cc2)
      endif
      if(idble(ii2).eq.2.and.idble(ii1).eq.2) then
      ccm(iposit(ii2)+1,iposit(ii1)+1)=dconjg(cc1)
      ccm(iposit(ii1)+1,iposit(ii2)+1)=cc1
      endif 

      enddo
      enddo
ccccccccccccccccccccccccccccccccccccccccccccccc
      num_st_old=num_st
      num_bad=0
      err_max=-1.d0

      do 600 num1=num_wr_dis(1),num_wr_dis(2)
      ist=num1-num_wr_dis(1)+1
      do ii1=num_st_old+1,num_st
      cc1=dcmplx(0.d0,0.d0)
      cc2=dcmplx(0.d0,0.d0)
      do k=1,n3
      do j=1,n2
      do i=1,n1w
      cc1=cc1+uc(i+nnposit,j,k,ist)*dconjg(ucw(i,j,k,ii1))
      cc2=cc2+uc(i+nnposit,j,k,ist)*ucw(i,j,k,ii1)
      enddo
      enddo
      enddo
      ccy(iposit(ii1),num1)=cc1
      if(idble(ii1).eq.2) then        ! not at the Gamma point, k.ne.0
      ccy(iposit(ii1)+1,num1)=cc2
      endif
      enddo
ccccccccccccccccccccccccccccccccccccccccc
ccc solve for: ccm(ii1,ii2)*ca(ii2)=ccy(ii1)
ccc then, uc= ucw(ii1)*ca(ii1)
cccccccccccccccccccccccccccccccccccccc
cccc lapack
      ccA=ccm
      num_tot=iposit(num_st)+idble(num_st)-1
      do i=1,num_tot            ! num_tot, number of electrode states, counting both phi_i and phi_i^*
      ccy2(i)=ccy(i,num1)
      enddo

      call zgesv(num_tot,1,ccA,nwellm,ipiv,ccy2,nwellm,info)

      do i=1,num_tot
      cc=dcmplx(0.d0,0.d0)
      do j=1,num_tot
      cc=cc+ccm(i,j)*ccy2(j)
      enddo
      ccy3(i)=cc
      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
       do i=1,num_tot
       ccy2_st(i,num1)=ccy2(i)
       enddo
       
       do i=1,num_st
       if(idble(i).eq.2) then
       ccy2_st(iposit(i),num1+num_stWell)=dconjg(ccy2(iposit(i)+1))         ! num1, from 1 to num_stWell, system states
       ccy2_st(iposit(i)+1,num1+num_stWell)=dconjg(ccy2(iposit(i)))
       endif
       if(idble(i).eq.1) then
       ccy2_st(iposit(i),num1+num_stWell)=dconjg(ccy2(iposit(i)))           ! this could be wrong for X point evanescent states, because it is not real
       endif

       enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      uc_test=dcmplx(0.d0,0.d0)
      uc_test2=dcmplx(0.d0,0.d0)

      aIcurrent=0.d0
      do ii=1,num_st
      if(ievan(ii).eq.0) then
      fact2=1.d0
      else
      fact2=0.d0
      endif
      if(idble(ii).eq.2) then
      aIcurrent=aIcurrent+(cdabs(ccy2(iposit(ii)))**2-
     &   cdabs(ccy2(iposit(ii)+1))**2)*dE_dk(ii)
      do k=1,n3
      do j=1,n2
      do i=1,n1w
      uc_test(i,j,k)=uc_test(i,j,k)+ccy2(iposit(ii))*ucw(i,j,k,ii)+
     &  ccy2(iposit(ii)+1)*dconjg(ucw(i,j,k,ii))
      uc_test2(i,j,k)=uc_test2(i,j,k)+fact2*(ccy2(iposit(ii))*
     &  ucw(i,j,k,ii)+ccy2(iposit(ii)+1)*dconjg(ucw(i,j,k,ii)))
      enddo
      enddo
      enddo
      else
      do k=1,n3
      do j=1,n2
      do i=1,n1w
      uc_test(i,j,k)=uc_test(i,j,k)+ccy2(iposit(ii))*ucw(i,j,k,ii)
      uc_test2(i,j,k)=uc_test2(i,j,k)+fact2*ccy2(iposit(ii))*
     &   ucw(i,j,k,ii)
      enddo
      enddo
      enddo
      endif
      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      sum1=0.d0
      sum2=0.d0
      diff=0.d0
      do k=1,n3
      do j=1,n2
      do i=1,n1w
      sum1=sum1+abs(uc_test(i,j,k))**2
      sum2=sum2+abs(uc(i+nnposit,j,k,ist))**2
      diff=diff+abs(uc(i+nnposit,j,k,ist)-uc_test(i,j,k))**2
ccccc THE CHANGE !
      uc_R(i,j,k,num1)=uc(i+nnposit,j,k,ist)-uc_test(i,j,k)      ! uc_test2 does not include the evanescent states
      enddo
      enddo
      enddo
      err=abs(1.d0-sum1/sum2)
      ibad=0
      if(err.gt.err_max) err_max=err
      if(err.gt.0.01.or.err*sum2.gt.0.001) then 
      num_bad=num_bad+1
      ibad=1
      endif
      weight(num1)=sum1/sum2
      aI_tmp(num1)=aIcurrent
      weight(num1+num_stWell)=sum1/sum2*0.99999999999d0   ! just to dist these two states      
600   continue

      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call mpi_allreduce(err_max,err,1,MPI_REAL8,MPI_MAX,
     &   MPI_COMM_WORLD,ierr)

      if(inode.eq.1) then
      write(6,*) "***********************************************"
      write(6,*) "******** iteration", num_iter, "***************"
      write(6,306) num_tot,num_run,err
      endif

306   format("num_tot,num_run,err=",2(i3,1x),2(f6.2,1x))

      if(num_bad.eq.0.or.num_iter.gt.num_iter_evan) goto 3000

2000  continue
3000  continue

      if(num_st.gt.num_st0) num_st=num_st0

      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call mpi_allreduce(ccy2_st,tmp_complex,nwellm*nwellm,
     $     MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
      ccy2_st=tmp_complex
      call mpi_allreduce(weight,tmp_real,nwellm,
     $     MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      weight=tmp_real

      mr=n1w*n2*n3
      mr_n=mr/mnodes
      allocate(uc_tmp(mr_n))
      allocate(uc_tmp2(mr_n))
      call mpi_barrier(MPI_COMM_WORLD,ierr)

      do ii1=1,num_stWell
        do iread=1,mnodes
         do ii=1,mr_n
         jj=ii+(iread-1)*mr_n
         i=(jj-1)/(n2*n3)+1
         j=(jj-1-(i-1)*n2*n3)/n3+1
         k=jj-(i-1)*n2*n3-(j-1)*n3
         uc_tmp(ii)=uc_R(i,j,k,ii1)
         enddo
         call mpi_barrier(MPI_COMM_WORLD,ierr)
         call mpi_allreduce(uc_tmp,uc_tmp2,mr_n,MPI_DOUBLE_COMPLEX,
     $       MPI_SUM,MPI_COMM_WORLD,ierr)
         do ii=1,mr_n
         jj=ii+(iread-1)*mr_n
         i=(jj-1)/(n2*n3)+1
         j=(jj-1-(i-1)*n2*n3)/n3+1
         k=jj-(i-1)*n2*n3-(j-1)*n3
         uc_R(i,j,k,ii1)=uc_tmp2(ii)
         enddo
        enddo
      enddo

      call mpi_barrier(MPI_COMM_WORLD,ierr)
      deallocate(uc_tmp)
      deallocate(uc_tmp2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      sum=0.d0
      do k=1,n3
      do j=1,n2
      do i=1,n1w
      sum=sum+abs(ucw(i,j,k,1))**2
      enddo
      enddo
      enddo

      fact=1.d0/sum

      do num1=num_wr_dis(1),num_wr_dis(2)
      do num2=1,num_stWell
      cc=dcmplx(0.d0,0.d0)
      cc1=dcmplx(0.d0,0.d0)
      do k=1,n3
      do j=1,n2
      do i=1,n1w
      cc=cc+uc_R(i,j,k,num1)*dconjg(uc_R(i,j,k,num2))
      cc1=cc1+dconjg(uc_R(i,j,k,num1))*dconjg(uc_R(i,j,k,num2))
      enddo
      enddo
      enddo
      cc=cc*fact
      cc1=cc1*fact
      cc_R(num1,num2)=cc
      cc_R(num1+num_stWell,num2+num_stWell)=dconjg(cc)
      cc_R(num1+num_stWell,num2)=cc1
      cc_R(num1,num2+num_stWell)=dconjg(cc1)
      enddo
      enddo
      cc_R=dconjg(cc_R)     ! important
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call mpi_allreduce(cc_R,tmp_complex,nwellm*nwellm,
     $     MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
      cc_R=tmp_complex

      deallocate(uc_test)
      deallocate(uc_test2)
      deallocate(uc_R)
ccccccc output is ccy2_st and ucw
      return
      end 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
