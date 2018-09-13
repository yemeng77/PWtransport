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
      integer ikpt_st1w(nstm),ikpt_st2w(nstm),
     &   ikpt_st3w(nstm),
     &   i_st1w(nstm),
     &  i_st2w(nstm),i_st3w(nstm),
     &  num_mx(nstm),imx(10,nstm)
      real*8 x_st1w(nstm),x_st2w(nstm),x_st3w(nstm)
      complex*16 cphase_ucw(nstm)

      real*8 kk1,kk2,kk3,k1k2,k2k3,k2k1,k3k1,k3k2,k1k3

      complex*16, allocatable, dimension (:,:,:) :: uc_test,
     &   uc_test2
      complex*16, allocatable, dimension (:,:,:,:) :: uc_R
      complex*16, allocatable, dimension (:) :: uc_tmp

      real*8 weight(nwellm)
      complex*16 ccy2_st(nwellm,nwellm),cc_R(nwellm,nwellm)
      complex*16 ccm(nwellm,nwellm),ccA(nwellm,nwellm)
      complex*16 ccy(nwellm,10),ccy2(nwellm),ccy3(nwellm)
      integer  ipiv(nwellm),iposit(nwellm),idble(nwellm),ievan(nwellm)

      complex*16 cc,cc1,cc2,cc3,cc_st

      character*7 fileh

      integer num_wr_dis(2)
      integer inode,nnodes

      common /mpi_data/inode,nnodes
 
      dE_min=0.d0
      if(inode.eq.1) write(6,*) "Ew=",Ew

      fileh="wr.new."
      open(77,file=trim(fileh)//"001",form="unformatted")
      rewind(77)
      read(77) n1t,n2t,n3t,mnodes
      close(77)
      if(n1.ne.m1.or.n2.ne.m2.or.n3.ne.m3) then
      if(inode.eq.1) write(6,*) "n1w,n2,n3.ne.n1t,n2t,n3t,stop",
     &   n1w,n2,n3,n1t,n2t,n3t
      call mpi_abort(MPI_COMM_WORLD,1,ierr)
      endif
      mr=n1w*n2*n3
      mr_n=mr/mnodes
      allocate(uc_tmp(mr_n))

      pi=4*datan(1.d0)

      iused_evan=0
      num_st=0
      iposit_next=1
      num_iter=0
3000  continue
      num_iter=num_iter+1
ccccccccccccccccccccccccccccccccccccccc
      if(num_iter.gt.1) goto 1000
ccccccccc  find the running waves
      if(inode.eq.1) then
      num_st_old=0

      do j=1,nstw    ! go through different band-line

      do i=1,numw(j)-1    ! the number of k-points in each band-line

      if((E_linew(i,j).le.Ew.and.E_linew(i+1,j).gt.Ew).or.
     &  (E_linew(i+1,j).lt.Ew.and.E_linew(i,j).ge.Ew)) then


      num_st=num_st+1
      iposit(num_st)=iposit_next
      idble(num_st)=2
      ievan(num_st)=0

      if(i.eq.1) then
       i33=i+2
      elseif(i.eq.(numw(j)-1)) then
       i33=i-1
      else
       if(abs(E_linew(i+2,j)-Ew).gt.abs(E_linew(i-1,j)-Ew)) then
          i33=i-1
       else
          i33=i+2
       endif
      endif

cccccccccccccccccccccccccccccccccccccccccc
cccc note, assuming (it is correct), ikpt_linew(i,j) and i has the same order. They can only be shifted by a number
      if(i33.eq.i+2) then
      i1=i
      i2=i+1
      i3=i+2
      else
      i1=i-1
      i2=i
      i3=i+1
      endif
ccccccccccccccccccccccccccccccccccccccccc
      E1=E_linew(i1,j)
      E2=E_linew(i2,j)
      E3=E_linew(i3,j)
ccccccccccccccccccccccccccccccccccccc
cccc assume E(k)=a1*k^2+a2*k+a3, and k=0 at i2, -1 at i1, 1 at i3 : E(-1)=E1,E(0)=E2,E(1)=E3
ccccc then
      a1=(E3+E1-2*E2)/2
      a2=(E3-E1)/2
      a3=E2
ccccc Then: E(k)=E1*(k^2-k)/2+E2*(1-k^2)+E3*(k^2+k)/2=E1*w1+E2*w2+E3*w3   ! the linear coeff for w1,w2,w3
ccccccccccccccccccccccccccc
cccc Then solve Ew=a1*ak^2+a2*ak+a3 to find ak
cccc ak is the kpoint distance from ikpt_linew(i2,j)
cccccccc
      ak1=(-a2+dsqrt(abs(a2**2-4*a1*(a3-Ew))))/(2*a1)
      ak2=(-a2-dsqrt(abs(a2**2-4*a1*(a3-Ew))))/(2*a1)
      if(abs(ak1).lt.abs(ak2)) then
      ak=ak1
      else
      ak=ak2
      endif

      if(num_st.gt.1.and.ist_linew(i1,j).eq.i_st1w(num_st-1).and.
     &  ikpt_linew(i1,j).eq.ikpt_st1w(num_st-1).and.ist_linew(i2,j).eq.
     &  i_st2w(num_st-1).and.ikpt_linew(i2,j).eq.ikpt_st2w(num_st-1)
     &  .and.ist_linew(i3,j).eq.i_st3w(num_st-1).and.ikpt_linew(i3,j)
     &  .eq.ikpt_st3w(num_st-1)) then
         if(abs(ak1).gt.abs(ak2)) then
         ak=ak1
         else
         ak=ak2
         endif
      endif

      if(abs(ak).gt.1) then
      write(6,*) "ak.gt.1, strange, stop",ak
      stop
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccc
      dE_dk3=(2*a1*ak+a2)/(pi/a11/(nkptw-1))
      dE_dk(num_st)=dE_dk3
      w1=(ak**2-ak)/2
      w2=1.d0-ak**2
      w3=(ak**2+ak)/2

cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ikpt_st1w(num_st)=ikpt_linew(i1,j)
      ikpt_st2w(num_st)=ikpt_linew(i2,j) 
      ikpt_st3w(num_st)=ikpt_linew(i3,j) 

      i_st1w(num_st)=ist_linew(i1,j)
      i_st2w(num_st)=ist_linew(i2,j)
      i_st3w(num_st)=ist_linew(i3,j)

      x_st1w(num_st) = w1
      x_st2w(num_st) = w2
      x_st3w(num_st) = w3
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ak_w(num_st)=ak+ikpt_linew(i2,j)

      nline_w(num_st)=j

      iposit_next=iposit_next+idble(num_st)

      cphase_ucw(num_st)=
     &  cdexp(dcmplx(0.d0,pi*(ak_w(num_st)-1.d0)/(nkptw-1)))

      endif
      enddo
      enddo

      write(6,*) "XXXXXXXXXXXXXXXXXXXXXX"
      write(6,*) "The number of running waves =", num_st
c      do i=1,num_st
c        write(6,*) i_st1w(i),ak_w(i),dE_dk(i)
c      enddo
      write(6,301) (ikpt_st1w(i),i=1,num_st)
      write(6,302) (i_st1w(i),i=1,num_st)
      write(6,*) "XXXXXXXXXXXXXXXXXXXXXX"
301   format("ikpt= ",20(i4,1x)) 
302   format("i_st= ",20(i4,1x))
      endif ! end for inode.eq.1

      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call mpi_bcast(num_st,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(num_st.eq.0) then
      goto 1000      ! go straight to calculate the evanescence states
      else
      num_run=num_st
      goto 1001       ! jump over the following part for evanescence state calculation
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccc
1000  continue        ! calculate the evanescence states

      num_st_old=num_st
      dE_min=1000.d0

      if(inode.eq.1) then
      do j=1,num_evan 
      dE=Ew-E_evan(j)
       if(iused_evan(j).eq.0.and.
     &  dE*E_evanC(j).le.0) then
       if(abs(dE).lt.dE_min) then
       dE_min=abs(dE)
       j_min=j
       endif
       endif
      enddo
      endif

      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call mpi_bcast(dE_min,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

      if(dE_min.gt.dE_evan) goto 5000      !  end
cccccccccccccccccccccccc
      if(inode.eq.1) then
      iused_evan(j_min)=1
      num_st=num_st+1
      iposit(num_st)=iposit_next
      ievan(num_st)=1
      if(iGX_evan(j_min).eq.2) then
      idble(num_st)=2
      else
      idble(num_st)=1
      endif
      iposit_next=iposit_next+idble(num_st)
      
      ikpt_st1w(num_st)=ikpt_evan(j_min)
      if(ikpt_st1w(num_st).lt.nkptw) then
      ikpt_st2w(num_st)=ikpt_st1w(num_st)+1
      else
      ikpt_st2w(num_st)=ikpt_st1w(num_st)-1
      endif
      i_st1w(num_st)=ist_evan(j_min)
      i_st2w(num_st)=ist_evan(j_min)
      x_st1w(num_st)=1.d0
      x_st2w(num_st)=0.d0
      ak_w(num_st)=ikpt_evan(j_min)
      nline_w(num_st)=j_min
      dE_dk(num_st)=0.d0
      ak=ikpt_evan(j_min)
      cphase_ucw(num_st)=cdexp(dcmplx(0.d0,pi*(ak-1.d0)/(nkptw-1)))

      write(6,*) "YYYYYYYYYYYYYYYYYY"
      write(*,*) "evanescent state,ind,ist,iGX,E_evan0,dE_min "
      write(6,501)  num_st-num_run,ist_evan(j_min),iGX_evan(j_min),
     &  E_evan(j_min),dE_min
      write(6,*) "YYYYYYYYYYYYYYYYYY"
501   format("evan_state,ind,ist,iGX,E_evan,dE_min ",
     &  3(i3,1x),2(f9.6,1x))
      endif

cccccccccccccccccccccccccccccccccccccccccccccc
1001  continue       ! end calc. of running wave and evanescence states

ccccccccccccccccccccccccccccccc
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call mpi_bcast(iposit_next,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      num_tot=iposit_next-1

      if(inode.eq.1) then
      write(6,503) num_tot,num_run,num_tot-num_run,num_iter
503   format("num_tot,num_run,num_evan,num_iter",
     &  4(i3,1x))
      endif     

      if(num_iter.eq.1) then
      ccm=dcmplx(0.d0,0.d0)
      ccA=dcmplx(0.d0,0.d0)
      ccy=dcmplx(0.d0,0.d0)
      ccy2=dcmplx(0.d0,0.d0)
      ccy3=dcmplx(0.d0,0.d0)
      ipiv=0
      endif

      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call mpi_bcast(num_st,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(ipiv,num_st,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(iposit,num_st,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(idble,num_st,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(ievan,num_st,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(dE_dk,num_st,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

cccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccc
      if(inode.eq.1) then
      do 1020 ii=num_st_old+1,num_st
      if(num_iter.eq.1) then
      call wave_electrode_3kpts(ikpt_st1w(ii),
     & ikpt_st2w(ii),
     & ikpt_st3w(ii),
     & i_st1w(ii),i_st2w(ii),i_st3w(ii),
     & x_st1w(ii),x_st2w(ii),x_st3w(ii),
     & ucw(1,1,1,ii),n1w,n2,n3,fileh,1,cphase,phase,nintep)
      else
      call wave_electrode(ikpt_st1w(ii),ikpt_st2w(ii),
     & i_st1w(ii),i_st2w(ii),x_st1w(ii),x_st2w(ii),
     & ucw(1,1,1,ii),n1w,n2,n3,fileh,1,cphase,phase,nintep)
      endif

        if(idble(ii).eq.1) then     ! make ucw real, Note, even for X-point evanescent states, we can make it real
        cc=dcmplx(0.d0,0.d0)
        sum=0.d0
        do k=1,n3
        do j=1,n2
        do i=1,n1w
        cc=cc+ucw(i,j,k,ii)**2
        sum=sum+abs(ucw(i,j,k,ii))**2
        enddo
        enddo
        enddo
        cc=1.d0/sqrt(cc/sum)
        sum1=0.d0
        do k=1,n3
        do j=1,n2
        do i=1,n1w
        ucw(i,j,k,ii)=dcmplx(dreal(cc*ucw(i,j,k,ii)),0.d0)
        sum1=sum1+abs(ucw(i,j,k,ii))**2
        enddo
        enddo
        enddo
        fact=dsqrt(sum/sum1)
        do k=1,n3
        do j=1,n2
        do i=1,n1w
        ucw(i,j,k,ii)=fact*ucw(i,j,k,ii)
        enddo
        enddo
        enddo
       
        endif      ! idble(ii).eq.1


1020  continue
ccccccccccccccccccccccccccccccccccc
 
ccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc TEST ORTHOGONALITY
cccccccccc    test
        do ii1=num_st_old+1,num_st
        do ii2=1,ii1-1
        sum2=0.d0
        cc=dcmplx(0.d0,0.d0)
        do k=1,n3
        do j=1,n2
        do i=1,n1w
        sum2=sum2+abs(ucw(i,j,k,ii2))**2
        cc=cc+ucw(i,j,k,ii1)*dconjg(ucw(i,j,k,ii2))
        enddo
        enddo
        enddo
        if(ii2.eq.1) sum000=sum2
        cc=cc/sum2
        orth1=abs(cc)

        if(abs(ak_w(ii1)-ak_w(ii2)).lt.1.E-3) then      ! ak_w, from 1 to 101, for the same ak, make sure they are orth
        do k=1,n3
        do j=1,n2
        do i=1,n1w
        ucw(i,j,k,ii1)=ucw(i,j,k,ii1)-cc*ucw(i,j,k,ii2)
        enddo
        enddo
        enddo
        orth1=0.d0
        endif

        enddo    ! ii2

        sum=0.d0
        do k=1,n3
        do j=1,n2
        do i=1,n1w
        sum=sum+abs(ucw(i,j,k,ii1))**2
        enddo
        enddo
        enddo
        if(ii1.eq.1) then
          sum000=sum
        else
          fact=dsqrt(sum000/sum)
          do k=1,n3
          do j=1,n2
          do i=1,n1w
          ucw(i,j,k,ii1)=ucw(i,j,k,ii1)*fact
          enddo
          enddo
          enddo
        endif
      
        enddo     ! ii1
      endif


      do ist=num_st_old+1,num_st
        do iread=1,mnodes
        if(inode.eq.1) then
         do ii=1,mr_n
         jj=ii+(iread-1)*mr_n
         i=(jj-1)/(n2*n3)+1
         j=(jj-1-(i-1)*n2*n3)/n3+1
         k=jj-(i-1)*n2*n3-(j-1)*n3
         uc_tmp(ii)=ucw(i,j,k,ist)
         enddo
        endif
        call mpi_barrier(MPI_COMM_WORLD,ierr)
        call mpi_bcast(uc_tmp,mr_n,MPI_DOUBLE_COMPLEX,0,
     &       MPI_COMM_WORLD,ierr)
        if(inode.ne.1) then
         do ii=1,mr_n
         jj=ii+(iread-1)*mr_n
         i=(jj-1)/(n2*n3)+1
         j=(jj-1-(i-1)*n2*n3)/n3+1
         k=jj-(i-1)*n2*n3-(j-1)*n3
         ucw(i,j,k,ist)=uc_tmp(ii)
         enddo
        endif
        enddo
      enddo

ccccccccccccccccccccccccccccccccccc
      do ii1=num_st_old+1,num_st
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

c       open(56,file='decomp_systemst.out')
ccccccccccccccccccccccccccccccccccccccccccccccc
      if(num_iter.eq.1) then
      allocate(uc_test(n1w,n2,n3))
      allocate(uc_test2(n1w,n2,n3))
      nst1=num_wr_dis(2)-num_wr_dis(1)+1
      allocate(uc_R(n1w,n2,n3,nst1))
      endif

      num_bad=0

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
      if(err.gt.0.01.or.err*sum2.gt.0.001) then 
      num_bad=num_bad+1
      ibad=1
      endif
      if(ist.eq.1.and.inode.eq.1) then
      write(56,*) "***********************************************"
      write(56,*) "******** iteration", num_iter, "***************"
      endif
      if(inode.eq.1) write(56,306) num_iter,num1,sum1/sum2,
     &          aIcurrent,(abs(ccy2(i)),i=1,num_tot)

      
      weight(num1)=sum1/sum2
      aI_tmp(num1)=aIcurrent
      weight(num1+num_stWell)=sum1/sum2*0.99999999999d0   ! just to dist these two states
303   format("nnposit,ibad,sum,diff,abs_sum", 
     &  i4,2x,i2,f10.7,1x,2(E10.4,1x))
304   format("num_st,sum_fract,diff,sum_absolute,Current", 
     &  i4,2x,f10.7,1x,2(E10.4,1x),2x,E14.8)
305   format(i3,1x,f10.7,1x,E14.8,3x,20(E6.1,1x))
306   format(i3,i3,1x,f10.7,1x,E14.8,3x,20(E6.1,1x))
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
600   continue

      if(num_bad.gt.0.and.num_iter.lt.num_iter_evan.and.
     &   dabs(dE_min).lt.dE_evan) goto 3000
5000  continue

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      sum=0.d0
      do k=1,n3
      do j=1,n2
      do i=1,n1w
      sum=sum+abs(ucw(i,j,k,1))**2
      enddo
      enddo
      enddo

      fact=1.d0/sum

      do num1=1,num_stWell
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
      call mpi_bcast(nline_w,nwellm,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(ak_w,nwellm,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(dE_dk,nwellm,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(cphase_ucw,nwellm,MPI_DOUBLE_COMPLEX,0,
     &    MPI_COMM_WORLD,ierr)

      deallocate(uc_test)
      deallocate(uc_tmp)
ccccccc output is ccy2_st and ucw

      return
      end 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
