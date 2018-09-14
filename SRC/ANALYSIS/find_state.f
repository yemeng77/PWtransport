      subroutine find_state(Ew,num_tot,num_st,num_run,
     & idble,iposit,ievan,ucw,n1w,n2,n3,cphase_ucw,dE_dk,ak_w,nline_w,
     & numw,E_linew,ist_linew,ikpt_linew,nstw,nkptw,nintep,
     & num_mx,imx,cphase,phase,fileh,mnodes,
     & E_evan,E_evanC,ist_evan,ikpt_evan,iGX_evan,num_evan,dE_evan)

      implicit double precision (a-h,o-z)
      include 'mpif.h'

      integer status(MPI_STATUS_SIZE)

      parameter (nm=1200)
      parameter (nevan=4000)
      parameter (nstm=200)

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

      integer  iposit(nstm),idble(nstm),ievan(nstm)

      complex*16 cc
      character*7 fileh

      complex*16, allocatable, dimension (:) :: uc_tmp,uc_tmp2

      integer num_wr_dis(2)
      integer inode,nnodes

      common /mpi_data/inode,nnodes

      if(inode.eq.1) then
      dE_min=0.d0
      write(6,*) "Ew=",Ew
      pi=4*datan(1.d0)
      iused_evan=0
      num_st=0
      iposit_next=1

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
      call mpi_abort(MPI_COMM_WORLD,1,ierr)
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

      num_run=num_st
      write(6,*) "XXXXXXXXXXXXXXXXXXXXXX"
      write(6,*) "The number of running waves =", num_st
      write(6,301) (ikpt_st1w(i),i=1,num_st)
      write(6,302) (i_st1w(i),i=1,num_st)
      write(6,*) "XXXXXXXXXXXXXXXXXXXXXX"
301   format("ikpt= ",20(i4,1x)) 
302   format("i_st= ",20(i4,1x))


3000  continue
      dE_min=1000.d0
      do j=1,num_evan 
      dE=Ew-E_evan(j)
       if(iused_evan(j).eq.0.and.dE*E_evanC(j).le.0.and.
     &   abs(dE).lt.dE_min) then
       dE_min=abs(dE)
       j_min=j
       endif
      enddo

      if(dE_min.gt.dE_evan) goto 5000      !  end
cccccccccccccccccccccccc
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
      goto 3000

5000  continue
cccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccc
      num_tot=iposit_next-1
      endif ! inode.eq.1
cccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccc
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call mpi_bcast(num_tot,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(num_st,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(num_run,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(idble,num_st,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(iposit,num_st,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(ievan,num_st,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(ikpt_st1w,num_st,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(ikpt_st2w,num_st,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(ikpt_st3w,num_st,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(i_st1w,num_st,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(i_st2w,num_st,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(i_st3w,num_st,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(x_st1w,num_st,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(x_st2w,num_st,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(x_st3w,num_st,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(nline_w,num_st,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(ak_w,num_st,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(dE_dk,num_st,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(cphase_ucw,num_st,MPI_DOUBLE_COMPLEX,0,
     &    MPI_COMM_WORLD,ierr)

      ucw=dcmplx(0.d0,0.d0)

      n_tmp1=num_st/nnodes
      n_tmp2=num_st-nnodes*n_tmp1
      if(inode.le.n_tmp2) then
        num_wr_dis(1)=(inode-1)*(n_tmp1+1)+1
        num_wr_dis(2)=inode*(n_tmp1+1)
      else
        num_wr_dis(1)=(inode-1)*n_tmp1+n_tmp2+1
        num_wr_dis(2)=inode*n_tmp1+n_tmp2
      endif

      do 600 ii=num_wr_dis(1),num_wr_dis(2)
      if(ii.le.num_run) then
      call wave_electrode_3kpts(ikpt_st1w(ii),
     & ikpt_st2w(ii),ikpt_st3w(ii),
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

600   continue
ccccccccccccccccccccccccccccccccccc
      mr=n1w*n2*n3
      mr_n=mr/mnodes
      allocate(uc_tmp(mr_n))
      allocate(uc_tmp2(mr_n))
      call mpi_barrier(MPI_COMM_WORLD,ierr)

      do ist=1,num_st
        do iread=1,mnodes
         do ii=1,mr_n
         jj=ii+(iread-1)*mr_n
         i=(jj-1)/(n2*n3)+1
         j=(jj-1-(i-1)*n2*n3)/n3+1
         k=jj-(i-1)*n2*n3-(j-1)*n3
         uc_tmp(ii)=ucw(i,j,k,ist)
         enddo
         call mpi_barrier(MPI_COMM_WORLD,ierr)
         call mpi_allreduce(uc_tmp,uc_tmp2,mr_n,MPI_DOUBLE_COMPLEX,
     $       MPI_SUM,MPI_COMM_WORLD,ierr)
         do ii=1,mr_n
         jj=ii+(iread-1)*mr_n
         i=(jj-1)/(n2*n3)+1
         j=(jj-1-(i-1)*n2*n3)/n3+1
         k=jj-(i-1)*n2*n3-(j-1)*n3
         ucw(i,j,k,ist)=uc_tmp2(ii)
         enddo
        enddo
      enddo
 
ccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc TEST ORTHOGONALITY
cccccccccc    test
      do ii1=1,num_st
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

      deallocate(uc_tmp)
      deallocate(uc_tmp2)
      call mpi_barrier(MPI_COMM_WORLD,ierr)

      return
      end