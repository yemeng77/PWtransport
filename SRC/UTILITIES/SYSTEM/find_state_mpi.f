      program find_state

      implicit double precision (a-h,o-z)
      include "mpif.h"

      integer status(MPI_STATUS_SIZE)

      parameter (nm=400)

      real*8 ALw(3,3),AL(3,3)
cccccc AL is the lattice of the large cell
cccccc ALw is the lattice of the electrode unit cell
      integer ikpt_st1w(nm),ikpt_st2w(nm),ikpt_st3w(nm)
      integer i_st1w(nm),i_st2w(nm),i_st3w(nm)
      real*8 x_st1w(nm),x_st2w(nm),x_st3w(nm)
      real*8 dE_dk(nm),ak_w(nm)
      integer itag(nm)
      complex*16 cphase_ucw(nm)
      integer num_st_dis(2)

      character*7 fileh
      character*20 fwr_all

      integer, allocatable, dimension (:) :: numw
      real*8, allocatable, dimension (:,:) :: E_linew
      integer, allocatable, dimension (:,:) :: ist_linew,ikpt_linew

      integer, allocatable, dimension (:) :: ist_evan,ikpt_evan
      integer, allocatable, dimension (:) :: iGX_evan,iband_evan
      real*8, allocatable, dimension (:) :: E_evan,E_evanC
      integer, allocatable, dimension (:) :: itag_evan

      integer, allocatable, dimension (:) :: iband_max,iband_min
      integer, allocatable, dimension (:) :: inum_max,inum_min
      real*8, allocatable, dimension (:) :: E_max,E_min

      integer, allocatable, dimension (:) :: num1,num2
      integer, allocatable, dimension (:,:) :: ik1,ik2,ikw
      real*8, allocatable, dimension (:,:) :: Ekw

      real*8, allocatable, dimension (:,:,:) :: phase
      complex*16, allocatable, dimension (:,:,:) :: cphase
      complex*16, allocatable, dimension (:,:,:) :: uc
      real*8, allocatable, dimension (:) :: ucR,uCI

      call mpi_init(ierr)
      call mpi_comm_size(MPI_COMM_WORLD,nnodes,ierr)
      call mpi_comm_rank(MPI_COMM_WORLD,inode,ierr)      
      inode=inode+1

      if(inode.eq.1) then
      open(10,file="find_wrl.input")
      rewind(10)
      read(10,*) n1,n2,n3,n1w,nnodes0
      read(10,*) Ew,dV,dE_evan,dk_same
      close(10)
      open(21,file="wr.new.001",form="unformatted")
      rewind(21)
      read(21) n1t,n2t,n3t,nnodesw,ispin_i,ispin_f,iw_i,iw_f
      read(21) ALw
      close(21)
      if(n1w.ne.n1t.or.n2.ne.n2t.or.n3.ne.n3t) then
      write(6,*) "n1,n2,n3 changed in wr.new.***, stop",
     & n1w,n2,n3,n1t,n2t,n3t
      call mpi_abort(MPI_COMM_WORLD,1,ierr)
      endif
      endif

      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call mpi_bcast(nkptw,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(n1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(n2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(n3,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(n1w,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(nnodesw,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(inode.eq.1) then
      AL(:,1)=ALw(:,1)*n1*1.d0/n1w
      AL(:,2:3)=ALw(:,2:3) 
      dE_min=0.d0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Ew1=Ew-(dV/2)
      Ew2=Ew+(dV/2)
      write(6,*) "Ew,Ew1,Ew2",Ew,Ew1,Ew2
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc Real band structure in GX
      open(14,file="E_line_W2K2")
      read(14,*) nbands,nkptw,nintep
      allocate(numw(nbands))
      allocate(E_linew(nkptw,nbands))
      allocate(ist_linew(nkptw,nbands))
      allocate(ikpt_linew(nkptw,nbands))
      do j=1,nbands
      read(14,*) j1,numw(j)
      read(14,*) (E_linew(i,j),i=1,numw(j))
      read(14,*) (ist_linew(i,j),i=1,numw(j))
      read(14,*) (ikpt_linew(i,j),i=1,numw(j))
      enddo
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc End of preprocessing
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call mpi_bcast(nintep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(nkptw,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      allocate(cphase(n1w,n2,n3))
      allocate(phase(n1w,n2,n3))
      pi=4*datan(1.d0)
      do k=1,n3
      do j=1,n2
      do i=1,n1w
      phase(i,j,k)=pi*(i-1.d0)/(nkptw-1.d0)/n1w
      cphase(i,j,k)=cdexp(-dcmplx(0.d0,phase(i,j,k)))
      enddo
      enddo
      enddo

ccccccccccccccccccccccccccccccccccccccc
ccccccccc  find the running waves
      if(inode.eq.1) then
      allocate(num1(nbands))
      allocate(num2(nbands))
      allocate(ik1(nbands,20))
      allocate(ik2(nbands,20))
      ik1=0
      ik2=0
ccc find intersections with RBS for E-dV/2 and E+dV/2
      do j=1,nbands    ! go through different band-line
      ij1=0
      ij2=0
      do i=1,numw(j)-1    ! the number of k-points in each band-line
      if((E_linew(i,j).le.Ew1.and.E_linew(i+1,j).gt.Ew1).or.
     &  (E_linew(i+1,j).lt.Ew1.and.E_linew(i,j).ge.Ew1)) then
         ij1=ij1+1
         ik1(j,ij1)=i
      endif
      if((E_linew(i,j).le.Ew2.and.E_linew(i+1,j).gt.Ew2).or.
     &  (E_linew(i+1,j).lt.Ew2.and.E_linew(i,j).ge.Ew2)) then
         ij2=ij2+1
         ik2(j,ij2)=i
      endif
      enddo !end loop over i(kpts)
      num1(j)=ij1
      num2(j)=ij2
      enddo !end loop over j(bands)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      allocate(ikw(nbands,40))
      allocate(Ekw(nbands,40))
      num_left=0
      num_right=0
      do j=1,nbands    ! go through different band-line
        n1k=num1(j)
        n2k=num2(j)
        num_left=num_left+n1k
        num_right=num_right+n1k
        nj=0
ccccccccc first, identify the same kpt wave from Ew1,Ew2
        if((n1k.ne.0).and.(n2k.ne.0)) then
        do i1=1,n1k
        do i2=1,n2k
          mk12=abs(ik1(j,i1)-ik2(j,i2))
          dmk12=mk12/(nkptw-1.0)
          if(dmk12.lt.dk_same) ik2(j,i2)=-ik2(j,i2) !make a tag
        enddo 
        enddo 
        endif
        do i1=1,n1k
        nj=nj+1
        ikw(j,nj)=ik1(j,i1)
        Ekw(j,nj)=Ew1
        enddo
        do i2=1,n2k
        nj=nj+1
        ikw(j,nj)=ik2(j,i2)! for right 
        Ekw(j,nj)=Ew2
        enddo
        num2(j)=nj
      enddo
      deallocate(ik1)
      deallocate(ik2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      num_st=0
      do j=1,nbands    ! go through different band-line
      do ij=1,num2(j)
      i=ikw(j,ij)
      Ewin=Ekw(j,ij)
      num_st=num_st+1
      if(ij.le.num1(j)) then
        itag(num_st)=0 ! For left electrode
      else if(i.gt.0) then
        itag(num_st)=1 ! For right electrode, and use to generate wl
      else
        itag(num_st)=2 ! For right electrode, and not use to generate wl
        i=-i
      endif
ccccccccc Now, find out the i1,i2,i3 for three point interpolation
ccccccccc Actually, for generating wrl, this might not be necessary, one point should
ccccccccc be good enough
      if(i.eq.1) then
       i33=i+2
      else if(i.eq.(numw(j)-1)) then
       i33=i-1
      else
       if(abs(E_linew(i+2,j)-Ewin).gt.abs(E_linew(i-1,j)-Ewin)) then
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
c      write(6,*) num_st,j,i1,i2,i3,itag(num_st)
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
      ak1=(-a2+dsqrt(abs(a2**2-4*a1*(a3-Ewin))))/(2*a1)
      ak2=(-a2-dsqrt(abs(a2**2-4*a1*(a3-Ewin))))/(2*a1)
      if(abs(ak1).lt.abs(ak2)) then
      ak=ak1
      else
      ak=ak2
      endif
      if(num_st.gt.1.and.((itag(num_st).eq.0.and.itag(num_st-1).eq.0)
     &  .or.((itag(num_st).gt.0.and.itag(num_st-1).gt.0)))) then
      if(ist_linew(i1,j).eq.i_st1w(num_st-1).and.ikpt_linew(i1,j).eq.
     &  ikpt_st1w(num_st-1).and.ist_linew(i2,j).eq.i_st2w(num_st-1).and.
     &  ikpt_linew(i2,j).eq.ikpt_st2w(num_st-1).and.ist_linew(i3,j).eq.
     &  i_st3w(num_st-1).and.ikpt_linew(i3,j).eq.ikpt_st3w(num_st-1)) 
     & then
         if(abs(ak1).gt.abs(ak2)) then
         ak=ak1
         else
         ak=ak2
         endif
      endif
      endif
      if(abs(ak).gt.1) then
        write(6,*) "ak.gt.1, strange, stop",ak
        call mpi_abort(MPI_COMM_WORLD,2,ierr)
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
      x_st1w(num_st)=w1
      x_st2w(num_st)=w2
      x_st3w(num_st)=w3
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ak_w(num_st)=ak+ikpt_linew(i2,j)
      cphase_ucw(num_st)=
     &  cdexp(dcmplx(0.d0,pi*(ak_w(num_st)-1.d0)/(nkptw-1)))
      enddo       ! different state within band
      enddo       !  different bands
ccccccccccccccccccccccccccccccccc
      write(6,*) "XXXXXXXXXXXXXXXXXXXXXX"
      write(6,*) "The number of running waves =", num_st
      write(6,301) (ikpt_st1w(i),i=1,num_st)
      write(6,302) (i_st1w(i),i=1,num_st)
      write(6,*) "XXXXXXXXXXXXXXXXXXXXXX"
301   format("ikpt= ",20(i4,1x)) 
302   format("i_st= ",20(i4,1x))
      write(6,*) 
      deallocate(num1)
      deallocate(num2)
      deallocate(ikw)
      deallocate(Ekw)
      write(6,*) "The number of running waves in left and right
     & elctrodes =",num_left,num_right
      endif !for inode=1

      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call mpi_bcast(num_st,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(ikpt_st1w,nm,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(ikpt_st2w,nm,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(ikpt_st3w,nm,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(i_st1w,nm,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(i_st2w,nm,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(i_st3w,nm,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(x_st1w,nm,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(x_st2w,nm,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(x_st3w,nm,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

      n_tmp1=num_st/nnodes
      n_tmp2=num_st-nnodes*n_tmp1
      if(inode.le.n_tmp2) then
        num_st_dis(1)=(inode-1)*(n_tmp1+1)+1
        num_st_dis(2)=inode*(n_tmp1+1)
      else
        num_st_dis(1)=(inode-1)*n_tmp1+n_tmp2+1
        num_st_dis(2)=inode*n_tmp1+n_tmp2
      endif
      
      allocate(uc(n1w,n2,n3))
      fileh="wr.new."

      do ist=num_st_dis(1),num_st_dis(2)
        call wave_electrode_3kpts_LW(ikpt_st1w(ist),ikpt_st2w(ist),
     &    ikpt_st3w(ist),i_st1w(ist),i_st2w(ist),i_st3w(ist),
     &    x_st1w(ist),x_st2w(ist),x_st3w(ist),uc,n1w,n2,n3,fileh,1,
     &    cphase,phase,nintep)
        ist_100=ist/100
        ist_10=(ist-ist_100*100)/10
        ist_1=ist-ist_100*100-ist_10*10
        fwr_all="wr.run"//char(48+ist_100)//char(48+ist_10)//
     &   char(48+ist_1)
        open(11,file=fwr_all,form="unformatted")
        rewind(11)
        write(11) uc
        close(11)
      enddo

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      num_wl=0
      if(inode.eq.1) then
        open(22,file="wr.left_run",form="unformatted")
        rewind(22)
        write(22) n1w,n2,n3,nnodesw,num_left
        write(22) ALw
        open(23,file="wr.right_run",form="unformatted")
        rewind(23)
        write(23) n1w,n2,n3,nnodesw,num_right
        write(23) ALw
        nr=n1w*n2*n3
        nr_n=nr/nnodesw
        allocate(ucR(nr_n))
        allocate(ucI(nr_n))
        num_tot=0
        do ist=1,num_st
          ist_100=ist/100
          ist_10=(ist-ist_100*100)/10
          ist_1=ist-ist_100*100-ist_10*10
          fwr_all="wr.run"//char(48+ist_100)//char(48+ist_10)//
     &        char(48+ist_1)
          open(11,file=fwr_all,form="unformatted")
          rewind(11)
          read(11) uc
          close(11,status='delete')
          if(itag(ist).eq.0) then
            iout=22
          else
            iout=23
          endif
          write(iout) ikpt_st1w(ist),i_st1w(ist),ak_w(ist),dE_dk(ist)
     &         ,cphase_ucw(ist)
          do iproc=1,nnodesw
          do ii=1,nr_n
            jj=ii+(iproc-1)*nr_n
            i=(jj-1)/(n2*n3)+1
            j=(jj-1-(i-1)*n2*n3)/n3+1
            k=jj-(i-1)*n2*n3-(j-1)*n3
            ucR(ii)=dreal(uc(i,j,k))
            ucI(ii)=dimag(uc(i,j,k))
          enddo
          write(iout) (ucR(i),i=1,nr_n), (ucI(i),i=1,nr_n)
          enddo !end loop over iproc
          if(itag(ist).ne.2) then
            num_wl=num_wl+1
            call write_wl(n1,n2,n3,n1w,nnodesw,uc,AL,num_wl)
          endif
        enddo
      close(22)
      close(23)
      endif

      num_run=num_st
ccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc   calculate the evanescence states
      if(inode.eq.1) then
      num_st=0
      read(14,*) 
      read(14,*) nevan1
      allocate(ikpt_evan(nevan1))
      allocate(ist_evan(nevan1))
      allocate(iGX_evan(nevan1))
      allocate(iband_evan(nevan1))
      allocate(E_evan(nevan1))
      allocate(E_evanC(nevan1))
      do ii=1,nevan1
        read(14,*) ii2,ikpt_evan(ii),ist_evan(ii),
     &  iGX_evan(ii),iband_evan(ii),E_evan(ii),E_evanC(ii)
      enddo
      allocate(itag_evan(nevan1))
      itag_evan=0
      num_left=0
      num_right=0
      do ii=1,nevan1
      dE1=Ew1-E_evan(ii)
      if(dabs(dE1).lt.dE_evan.and.dE1*E_evanC(ii).le.0.d0) then
        num_left=num_left+1
        itag_evan(ii)=itag_evan(ii)+1 ! use for left electorde
      endif
      dE2=Ew2-E_evan(ii)
      if(dabs(dE2).lt.dE_evan.and.dE2*E_evanC(ii).le.0.d0) then
        num_right=num_right+1
        itag_evan(ii)=itag_evan(ii)+2 ! use for right electorde
      endif
      enddo

      do ii=1,nevan1
        if(itag_evan(ii).eq.0) cycle
        num_st=num_st+1
        ikpt_st1w(num_st)=ikpt_evan(ii)
        i_st1w(num_st)=ist_evan(ii)
        itag(num_st)=itag_evan(ii)
        ak=ikpt_evan(ii)
        ak_w(num_st)=ak
        dE_dk(num_st)=0.d0
        cphase_ucw(num_st)=cdexp(dcmplx(0.d0,pi*(ak-1.d0)/(nkptw-1)))
      enddo

      deallocate(ikpt_evan)
      deallocate(ist_evan)
      deallocate(iGX_evan)
      deallocate(iband_evan)
      deallocate(E_evan)
      deallocate(E_evanC)
      deallocate(itag_evan)
cccccccccccccccccccccccccccccccccccccccccccccc
      read(14,*)
      read(14,*) nevan2
      allocate(iband_max(nevan2))
      allocate(iband_min(nevan2))
      allocate(inum_max(nevan2))
      allocate(inum_min(nevan2))
      allocate(E_max(nevan2))
      allocate(E_min(nevan2))
      do ii=1,nevan2
        read(14,*) ii2,inum_max(ii),iband_max(ii),
     &  inum_min(ii),iband_min(ii),E_max(ii),E_min(ii)
      enddo
      close(14)
      allocate(itag_evan(nevan2))
      itag_evan=0
      do ii=1,nevan2
      if(E_max(ii).lt.Ew1.and.E_min(ii).gt.Ew1) then
        num_left=num_left+2
        itag_evan(ii)=itag_evan(ii)+1 ! use for left electorde
      endif
      if(E_max(ii).lt.Ew2.and.E_min(ii).gt.Ew2) then
        num_right=num_right+2
        itag_evan(ii)=itag_evan(ii)+2 ! use for right electorde
      endif
      enddo
      do ii=1,nevan2
        if(itag_evan(ii).eq.0) cycle
        num_st=num_st+2
        ikpt_st1w(num_st-1)=ikpt_linew(inum_max(ii),iband_max(ii))
        i_st1w(num_st-1)=ist_linew(inum_max(ii),iband_max(ii))
        itag(num_st-1)=itag_evan(ii)
        ak=ikpt_st1w(num_st-1)
        ak_w(num_st-1)=ak
        dE_dk(num_st-1)=0.d0
        cphase_ucw(num_st-1)=cdexp(dcmplx(0.d0,pi*(ak-1.d0)/(nkptw-1)))
        ikpt_st1w(num_st)=ikpt_linew(inum_min(ii),iband_min(ii))
        i_st1w(num_st)=ist_linew(inum_min(ii),iband_min(ii))
        itag(num_st)=itag_evan(ii)
        ak=ikpt_st1w(num_st)
        ak_w(num_st)=ak
        dE_dk(num_st)=0.d0
        cphase_ucw(num_st)=cdexp(dcmplx(0.d0,pi*(ak-1.d0)/(nkptw-1)))
      enddo
      deallocate(iband_max)
      deallocate(iband_min)
      deallocate(inum_max)
      deallocate(inum_min)
      deallocate(E_max)
      deallocate(E_min)
      deallocate(itag_evan)
      do ii=1,num_st
        write(6,*) "YYYYYYYYYYYYYYYYYY"
        write(6,501) ii,ikpt_st1w(ii),i_st1w(ii)
      enddo
501   format("evanescent state:ind,ikpt,ist: ",3(i4,1x))
      write(6,*) "The number of evanescent states in left and right
     & elctrodes =",num_left,num_right
      deallocate(numw)
      deallocate(E_linew)
      deallocate(ist_linew)
      deallocate(ikpt_linew)
      endif

cccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccc
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call mpi_bcast(num_st,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(ikpt_st1w,nm,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(i_st1w,nm,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      n_tmp1=num_st/nnodes
      n_tmp2=num_st-nnodes*n_tmp1
      if(inode.le.n_tmp2) then
        num_st_dis(1)=(inode-1)*(n_tmp1+1)+1
        num_st_dis(2)=inode*(n_tmp1+1)
      else
        num_st_dis(1)=(inode-1)*n_tmp1+n_tmp2+1
        num_st_dis(2)=inode*n_tmp1+n_tmp2
      endif

      do ist=num_st_dis(1),num_st_dis(2)
        call wave_electrode_interp(ikpt_st1w(ist),i_st1w(ist),uc,
     &      n1w,n2,n3,fileh,nintep)
        ist_100=ist/100
        ist_10=(ist-ist_100*100)/10
        ist_1=ist-ist_100*100-ist_10*10
        fwr_all="wr.evan"//char(48+ist_100)//char(48+ist_10)//
     &   char(48+ist_1)
        open(11,file=fwr_all,form="unformatted")
        rewind(11)
        write(11) uc
        close(11)
      enddo

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      if(inode.eq.1) then
        open(22,file="wr.left_evan",form="unformatted")
        rewind(22)
        write(22) n1w,n2,n3,nnodesw,num_left
        write(22) ALw
        open(23,file="wr.right_evan",form="unformatted")
        rewind(23)
        write(23) n1w,n2,n3,nnodesw,num_right
        write(23) ALw
        do ist=1,num_st
          ist_100=ist/100
          ist_10=(ist-ist_100*100)/10
          ist_1=ist-ist_100*100-ist_10*10
          fwr_all="wr.evan"//char(48+ist_100)//char(48+ist_10)//
     &        char(48+ist_1)
          open(11,file=fwr_all,form="unformatted")
          rewind(11)
          read(11) uc
          close(11,status='delete')
          num_wl=num_wl+1
          call write_wl(n1,n2,n3,n1w,nnodesw,uc,AL,num_wl)
          if(itag(ist).eq.1.or.itag(ist).eq.3) then
          write(22) ikpt_st1w(ist),i_st1w(ist),ak_w(ist),dE_dk(ist)
     &         ,cphase_ucw(ist)
          endif
          if(itag(ist).eq.2.or.itag(ist).eq.3) then
          write(23) ikpt_st1w(ist),i_st1w(ist),ak_w(ist),dE_dk(ist)
     &         ,cphase_ucw(ist)
          endif
          do iproc=1,nnodesw
          do ii=1,nr_n
            jj=ii+(iproc-1)*nr_n
            i=(jj-1)/(n2*n3)+1
            j=(jj-1-(i-1)*n2*n3)/n3+1
            k=jj-(i-1)*n2*n3-(j-1)*n3
            ucR(ii)=dreal(uc(i,j,k))
            ucI(ii)=dimag(uc(i,j,k))
          enddo
          if(itag(ist).eq.1.or.itag(ist).eq.3) then
          write(22) (ucR(i),i=1,nr_n), (ucI(i),i=1,nr_n)
          endif
          if(itag(ist).eq.2.or.itag(ist).eq.3) then
          write(23) (ucR(i),i=1,nr_n), (ucI(i),i=1,nr_n)
          endif
          enddo !end loop over iproc
        enddo
      close(22)
      close(23)
      deallocate(ucR)
      deallocate(uCI)
      endif

      call mpi_barrier(MPI_COMM_WORLD,ierr)
      deallocate(uc)
      deallocate(phase)
      deallocate(cphase)

ccccccccccccccccccccccccccccccccccc
      call mpi_finalize(ierr)

      stop
      end